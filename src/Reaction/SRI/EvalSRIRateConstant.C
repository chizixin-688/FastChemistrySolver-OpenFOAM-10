/*---------------------------------------------------------------------------*\
  Description
      Computing the forward rate constant of SRI reactions

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"

//=============================================================================//

void FastChemistry::OptReaction::evalSRIRateConstant()const noexcept
{

    unsigned int remainFO = this->n_SRIFO%4;
    __m256d invLog10v = _mm256_set1_pd(0.43429448190325182765112891891661);
    __m256d SRILimiterv = _mm256_set1_pd(FastChemistry::SRILimiter);
    __m256d onev = _mm256_set1_pd(1);
    __m256d logTv = _mm256_set1_pd(this->logT);

    for (unsigned int i = 0; i<this->n_SRIFO-remainFO; i=i+4)
    {
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

        //const double Kinf = this->Kf_[j+this->offset_kinf];
        __m256d Kinfv = _mm256_loadu_pd(&this->Kf_[j+this->offset_kinf]);


        //const double invKinf = 1.0/Kinf;
        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);


        //const double K0 = this->Kf_[j];
        __m256d K0v = _mm256_loadu_pd(&this->Kf_[j]);



        //const double M = this->tmp_M[m];
        __m256d Mv = _mm256_loadu_pd(&this->tmp_M[m]);




        //const double Pr = K0*M*invKinf; 
        __m256d Prv = _mm256_mul_pd(K0v,Mv);
        Prv = _mm256_mul_pd(Prv,invKinfv);


        //const double invLog10 = _mm_cvtsd_f64(_mm256_castpd256_pd128(invLog10v));
        //const double logPr = std::log(max(Pr, small))*invLog10;

        __m256d logPrv = _mm256_max_pd(Prv,SRILimiterv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);


        //const double expbT = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3];
        //const double expTc = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI];
        __m256d expbTv = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3]);
        __m256d expTcv = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI]);



        //const double X = 1/(1 + (logPr*logPr));
        __m256d Xv = _mm256_fmadd_pd(logPrv,logPrv,onev);
        Xv = _mm256_div_pd(onev,Xv);

        //const double a = this->a_[i];
        __m256d av = _mm256_loadu_pd(&this->a_[i]);


        //const double psi = a*expbT + expTc;
        __m256d psiv = _mm256_fmadd_pd(av,expbTv,expTcv);


        //const double logPsi = std::log(psi);
        __m256d logPsiv = vec256_logd(psiv);

        //const double d = this->d_[i];
        //const double e = this->e_[i];
        __m256d dv = _mm256_loadu_pd(&this->d_[i]);
        __m256d ev = _mm256_loadu_pd(&this->e_[i]);


        //const double logT = this->logT;
        //const double F = d*std::pow(psi, X)*std::pow(T, e);
        //const double F = d*std::exp(X*logPsi+e*logT);
        __m256d elogTv = _mm256_mul_pd(ev,logTv);
        __m256d Fv = _mm256_fmadd_pd(Xv,logPsiv,elogTv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);

        //const double invOnePlusPr = 1.0/(1.0+Pr);
        //const double N  = invOnePlusPr*F*K0;

        __m256d Nv = _mm256_add_pd(Prv,onev);
        Nv = _mm256_div_pd(Fv,Nv);
        Nv = _mm256_mul_pd(K0v,Nv);
        __m256d MNv = _mm256_mul_pd(Mv,Nv);

        //this->Kf_[j] = M*N; 
        _mm256_storeu_pd(&this->Kf_[j],MNv);
        //std::cout<<this->Kf_[j+0]<<" "<<this->Kf_[j+1]<<" "<<this->Kf_[j+2]<<" "<<this->Kf_[j+3]<<std::endl;
        //std::cout<<this->Kf_[j+0]<<std::endl;
        //std::cout<<this->Kf_[j+1]<<std::endl;
        //std::cout<<this->Kf_[j+2]<<std::endl;
        //std::cout<<this->Kf_[j+3]<<std::endl;
    }
    if(remainFO==1)
    {
        unsigned int i = this->n_SRIFO-1;
        const unsigned int j = this->SRIFO[i];
        //const unsigned int k = j - this->Ikf[4];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

        const double Kinf = this->Kf_[j+this->offset_kinf];
        const double invKinf = 1.0/Kinf;
        const double K0 = this->Kf_[j];
        const double M = tmp_M[m];

        const double Pr = K0*M*invKinf; 
        const double invLog10 = _mm_cvtsd_f64(_mm256_castpd256_pd128(invLog10v));
        const double logPr = std::log(std::max(Pr, FastChemistry::SRILimiter))*invLog10;

        const double expbT = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3];
        const double expTc = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI];

        const double X = 1/(1 + (logPr*logPr));
        const double a = this->a_[i];
        const double psi = a*expbT + expTc;
        const double logPsi = std::log(psi);

        const double d = this->d_[i];
        const double e = this->e_[i];
        const double logT = this->logT;
        //const double F = d*std::pow(psi, X)*std::pow(T, e);
        const double F = d*std::exp(X*logPsi+e*logT);

        const double invOnePlusPr = 1.0/(1.0+Pr);
        
        const double N  = invOnePlusPr*F*K0;
        this->Kf_[j] = M*N; 
        //std::cout<<this->Kf_[j+0]<<std::endl;
    }
    else if(remainFO==2)
    {
        unsigned int i = this->n_SRIFO-2;
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

        //const double Kinf = this->Kf_[j+this->offset_kinf];
        __m256d Kinfv = _mm256_setzero_pd();
        {
            __m128d tmp = _mm_loadu_pd(&this->Kf_[j+this->offset_kinf]);
            //Kinfv = _mm256_set_m128d(tmp,tmp);
            Kinfv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp), tmp, 1);
        }

        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);

        __m256d K0v = _mm256_setzero_pd();
        {
            __m128d tmp = _mm_loadu_pd(&this->Kf_[j]);
            //K0v = _mm256_set_m128d(tmp,tmp);
            K0v = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp), tmp, 1);
        }

        __m256d Mv = _mm256_setzero_pd();
        {
            __m128d tmp = _mm_loadu_pd(&this->tmp_M[m]);
            //Mv = _mm256_set_m128d(tmp,tmp);
            Mv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp), tmp, 1);
        }

        __m256d Prv = _mm256_mul_pd(K0v,Mv);
        Prv = _mm256_mul_pd(Prv,invKinfv);

        __m256d logPrv = _mm256_max_pd(Prv,SRILimiterv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);

        __m256d expbTv = _mm256_setzero_pd();
        __m256d expTcv = _mm256_setzero_pd();
        {
            __m128d tmp0 = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3]);
            //expbTv = _mm256_set_m128d(tmp0,tmp0);
            expbTv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp0), tmp0, 1);

            __m128d tmp1 = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI]);
            //expTcv = _mm256_set_m128d(tmp1,tmp1);
            expTcv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp1), tmp1, 1);
        }

        __m256d Xv = _mm256_fmadd_pd(logPrv,logPrv,onev);
        Xv = _mm256_div_pd(onev,Xv);

        __m256d av = _mm256_setzero_pd();
        {
            __m128d tmp = _mm_loadu_pd(&this->a_[i]);
            //av = _mm256_set_m128d(tmp,tmp);
            av = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp), tmp, 1);
        }

        __m256d psiv = _mm256_fmadd_pd(av,expbTv,expTcv);

        __m256d logPsiv = vec256_logd(psiv);

        __m256d dv = _mm256_setzero_pd();
        __m256d ev = _mm256_setzero_pd();
        {
            __m128d tmp0 = _mm_loadu_pd(&this->d_[i]);
            //dv = _mm256_set_m128d(tmp0,tmp0);
            dv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp0), tmp0, 1);

            __m128d tmp1 = _mm_loadu_pd(&this->e_[i]);
            //ev = _mm256_set_m128d(tmp1,tmp1);
            ev = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp1), tmp1, 1);
        }

        __m256d elogTv = _mm256_mul_pd(ev,logTv);
        __m256d Fv = _mm256_fmadd_pd(Xv,logPsiv,elogTv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);


        __m256d Nv = _mm256_add_pd(Prv,onev);
        Nv = _mm256_div_pd(Fv,Nv);
        Nv = _mm256_mul_pd(K0v,Nv);
        __m256d MNv = _mm256_mul_pd(Mv,Nv);

        _mm_storeu_pd(&this->Kf_[j],_mm256_castpd256_pd128(MNv));

        //std::cout<<this->Kf_[j+0]<<std::endl;
        //std::cout<<this->Kf_[j+1]<<std::endl;
    }
    else if(remainFO==3)
    {
        unsigned int i = this->n_SRIFO-3;
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

        //const double Kinf = this->Kf_[j+this->offset_kinf];
        const double Kinf0 = this->Kf_[j+this->offset_kinf+0];
        const double Kinf1 = this->Kf_[j+this->offset_kinf+1];
        const double Kinf2 = this->Kf_[j+this->offset_kinf+2];
        __m256d Kinfv = _mm256_setr_pd(Kinf0,Kinf1,Kinf2,Kinf2);


        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);


        const double K00 = this->Kf_[j+0];
        const double K01 = this->Kf_[j+1];
        const double K02 = this->Kf_[j+2];
        __m256d K0v = _mm256_setr_pd(K00,K01,K02,K02);

        const double M0 = this->tmp_M[m+0];
        const double M1 = this->tmp_M[m+1];
        const double M2 = this->tmp_M[m+2];
        __m256d Mv = _mm256_setr_pd(M0,M1,M2,M2);


        __m256d Prv = _mm256_mul_pd(K0v,Mv);
        Prv = _mm256_mul_pd(Prv,invKinfv);

        __m256d logPrv = _mm256_max_pd(Prv,SRILimiterv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);


        const double expbT0 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+0];
        const double expbT1 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+1];
        const double expbT2 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+2];
        __m256d expbTv = _mm256_setr_pd(expbT0,expbT1,expbT2,expbT2);

        const double expTc0 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+0];
        const double expTc1 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+1];
        const double expTc2 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+2];
        __m256d expTcv = _mm256_setr_pd(expTc0,expTc1,expTc2,expTc2);


        __m256d Xv = _mm256_fmadd_pd(logPrv,logPrv,onev);
        Xv = _mm256_div_pd(onev,Xv);

        const double a0 = this->a_[i+0];
        const double a1 = this->a_[i+1];
        const double a2 = this->a_[i+2];
        __m256d av = _mm256_setr_pd(a0,a1,a2,a2);


        __m256d psiv = _mm256_fmadd_pd(av,expbTv,expTcv);

        __m256d logPsiv = vec256_logd(psiv);

        const double d0 = this->d_[i+0];
        const double d1 = this->d_[i+1];
        const double d2 = this->d_[i+2];
        __m256d dv = _mm256_setr_pd(d0,d1,d2,d2);

        const double e0 = this->e_[i+0];
        const double e1 = this->e_[i+1];
        const double e2 = this->e_[i+2];
        __m256d ev = _mm256_setr_pd(e0,e1,e2,e2);


        __m256d elogTv = _mm256_mul_pd(ev,logTv);
        __m256d Fv = _mm256_fmadd_pd(Xv,logPsiv,elogTv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);


        __m256d Nv = _mm256_add_pd(Prv,onev);
        Nv = _mm256_div_pd(Fv,Nv);
        Nv = _mm256_mul_pd(K0v,Nv);
        __m256d MNv = _mm256_mul_pd(Mv,Nv);

        this->Kf_[j+0] = this->get_elem0(MNv);
        this->Kf_[j+1] = this->get_elem1(MNv);
        this->Kf_[j+2] = this->get_elem2(MNv);
        //std::cout<<this->Kf_[j+0]<<std::endl;
        //std::cout<<this->Kf_[j+1]<<std::endl;
        //std::cout<<this->Kf_[j+2]<<std::endl;
    }


    unsigned int remainCA = this->n_SRICA%4;
    for (unsigned int i = 0; i<this->n_SRICA-remainCA; i=i+4)
    {
        const unsigned int j = this->SRICA[i];
        const unsigned int m = j-this->Ikf[5]+this->Itbr[3];

        //const double Kinf = this->Kf_[j+this->offset_kinf];
        __m256d Kinfv = _mm256_loadu_pd(&this->Kf_[j-Ikf[5]+Ikf[8]]);


        //const double invKinf = 1.0/Kinf;
        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);


        //const double K0 = this->Kf_[j];
        __m256d K0v = _mm256_loadu_pd(&this->Kf_[j]);



        //const double M = this->tmp_M[m];
        __m256d Mv = _mm256_loadu_pd(&this->tmp_M[m]);




        //const double Pr = K0*M*invKinf; 
        __m256d Prv = _mm256_mul_pd(K0v,Mv);
        Prv = _mm256_mul_pd(Prv,invKinfv);


        //const double invLog10 = _mm_cvtsd_f64(_mm256_castpd256_pd128(invLog10v));
        //const double logPr = std::log(max(Pr, small))*invLog10;

        __m256d logPrv = _mm256_max_pd(Prv,SRILimiterv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);


        //const double expbT = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3];
        //const double expTc = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI];
        __m256d expbTv = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRIFO]);
        __m256d expTcv = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+this->n_SRIFO]);



        //const double X = 1/(1 + (logPr*logPr));
        __m256d Xv = _mm256_fmadd_pd(logPrv,logPrv,onev);
        Xv = _mm256_div_pd(onev,Xv);

        //const double a = this->a_[i];
        __m256d av = _mm256_loadu_pd(&this->a_[i+this->n_SRIFO]);


        //const double psi = a*expbT + expTc;
        __m256d psiv = _mm256_fmadd_pd(av,expbTv,expTcv);


        //const double logPsi = std::log(psi);
        __m256d logPsiv = vec256_logd(psiv);

        //const double d = this->d_[i];
        //const double e = this->e_[i];
        __m256d dv = _mm256_loadu_pd(&this->d_[i+this->n_SRIFO]);
        __m256d ev = _mm256_loadu_pd(&this->e_[i+this->n_SRIFO]);


        //const double logT = this->logT;
        //const double F = d*std::pow(psi, X)*std::pow(T, e);
        //const double F = d*std::exp(X*logPsi+e*logT);
        __m256d elogTv = _mm256_mul_pd(ev,logTv);
        __m256d Fv = _mm256_fmadd_pd(Xv,logPsiv,elogTv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);

        //const double invOnePlusPr = 1.0/(1.0+Pr);
        //const double N  = invOnePlusPr*F*K0;

        __m256d Nv = _mm256_add_pd(Prv,onev);
        Nv = _mm256_div_pd(Fv,Nv);
        Nv = _mm256_mul_pd(K0v,Nv);

        //this->Kf_[j] = M*N; 
        _mm256_storeu_pd(&this->Kf_[j],Nv);
        //std::cout<<this->Kf_[j+0]<<std::endl;
        //std::cout<<this->Kf_[j+1]<<std::endl;
        //std::cout<<this->Kf_[j+2]<<std::endl;
        //std::cout<<this->Kf_[j+3]<<std::endl;
    }
    if(remainCA==1)
    {
        unsigned int i = this->n_SRICA-1;
        const unsigned int j = this->SRICA[i];
        //const unsigned int k = j - this->Ikf[4];
        const unsigned int m = j-this->Ikf[5]+this->Itbr[3];

        const double Kinf = this->Kf_[j-Ikf[5]+Ikf[8]];
        const double invKinf = 1.0/Kinf;
        const double K0 = this->Kf_[j];
        const double M = tmp_M[m];

        const double Pr = K0*M*invKinf; 
        const double invLog10 = _mm_cvtsd_f64(_mm256_castpd256_pd128(invLog10v));
        const double logPr = std::log(std::max(Pr, FastChemistry::SRILimiter))*invLog10;

        const double expbT = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRIFO];
        const double expTc = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+this->n_SRIFO];

        const double X = 1/(1 + (logPr*logPr));
        const double a = this->a_[i+this->n_SRIFO];
        const double psi = a*expbT + expTc;
        const double logPsi = std::log(psi);

        const double d = this->d_[i+this->n_SRIFO];
        const double e = this->e_[i+this->n_SRIFO];
        const double logT = this->logT;
        //const double F = d*std::pow(psi, X)*std::pow(T, e);
        const double F = d*std::exp(X*logPsi+e*logT);

        const double invOnePlusPr = 1.0/(1.0+Pr);
        
        const double N  = invOnePlusPr*F*K0;
        this->Kf_[j] = N; 
        //std::cout<<this->Kf_[j+0]<<std::endl;
    }
    else if(remainCA==2)
    {
        unsigned int i = this->n_SRICA-2;
        const unsigned int j = this->SRICA[i];
        const unsigned int m = j-this->Ikf[5]+this->Itbr[3];

        //const double Kinf0 = this->Kf_[j+0-Ikf[5]+Ikf[8]];
        //const double Kinf1 = this->Kf_[j+1-Ikf[5]+Ikf[8]];
        __m256d Kinfv = _mm256_setzero_pd();
        {
            __m128d tmp = _mm_loadu_pd(&this->Kf_[j-Ikf[5]+Ikf[8]]);
            //Kinfv = _mm256_set_m128d(tmp,tmp);
            Kinfv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp), tmp, 1);
        }

        //const double invKinf0 = 1.0/Kinf0;
        //const double invKinf1 = 1.0/Kinf1;
        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);

        //const double K00 = this->Kf_[j+0];
        //const double K01 = this->Kf_[j+1];
        __m256d K0v = _mm256_setzero_pd();
        {
            __m128d tmp = _mm_loadu_pd(&this->Kf_[j]);
            //K0v = _mm256_set_m128d(tmp,tmp);
            K0v = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp), tmp, 1);
        }

        //const double M0 = tmp_M[m+0];
        //const double M1 = tmp_M[m+1];
        __m256d Mv = _mm256_setzero_pd();
        {
            __m128d tmp = _mm_loadu_pd(&this->tmp_M[m]);
            //Mv = _mm256_set_m128d(tmp,tmp);
            Mv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp), tmp, 1);
        }

        //const double Pr0 = K00*M0*invKinf0; 
        //const double Pr1 = K01*M1*invKinf1; 
        __m256d Prv = _mm256_mul_pd(K0v,Mv);
        Prv = _mm256_mul_pd(Prv,invKinfv);


        //const double invLog10 = _mm_cvtsd_f64(_mm256_castpd256_pd128(invLog10v));
        //const double logPr0 = std::log(max(Pr0, small))*invLog10;
        //const double logPr1 = std::log(max(Pr1, small))*invLog10;
        __m256d logPrv = _mm256_max_pd(Prv,SRILimiterv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);


        //const double expbT0 = this->tmp_Exp[i+0+this->nSpecies+this->n_Troe*3+this->n_SRIFO];
        //const double expbT1 = this->tmp_Exp[i+1+this->nSpecies+this->n_Troe*3+this->n_SRIFO];
        //const double expTc0 = this->tmp_Exp[i+0+this->nSpecies+this->n_Troe*3+this->n_SRI+this->n_SRIFO];
        //const double expTc1 = this->tmp_Exp[i+1+this->nSpecies+this->n_Troe*3+this->n_SRI+this->n_SRIFO];
        __m256d expbTv = _mm256_setzero_pd();
        __m256d expTcv = _mm256_setzero_pd();
        {
            __m128d tmp0 = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRIFO]);
            //expbTv = _mm256_set_m128d(tmp0,tmp0);
            expbTv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp0), tmp0, 1);

            __m128d tmp1 = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+this->n_SRIFO]);
            //expTcv = _mm256_set_m128d(tmp1,tmp1);
            expTcv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp1), tmp1, 1);
        }

        //const double X0 = 1/(1 + (logPr0*logPr0));
        //const double X1 = 1/(1 + (logPr1*logPr1));
        __m256d Xv = _mm256_fmadd_pd(logPrv,logPrv,onev);
        Xv = _mm256_div_pd(onev,Xv);

        //const double a0 = this->a_[i+0+this->n_SRIFO];
        //const double a1 = this->a_[i+1+this->n_SRIFO];
        __m256d av = _mm256_setzero_pd();
        {
            __m128d tmp = _mm_loadu_pd(&this->a_[i+this->n_SRIFO]);
            //av = _mm256_set_m128d(tmp,tmp);
            av = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp), tmp, 1);
        }

        //const double psi0 = a0*expbT0 + expTc0;
        //const double psi1 = a1*expbT1 + expTc1;
        __m256d psiv = _mm256_fmadd_pd(av,expbTv,expTcv);

        //const double logPsi0 = std::log(psi0);
        //const double logPsi1 = std::log(psi1);
        __m256d logPsiv = vec256_logd(psiv);


        //const double d0 = this->d_[i+0+this->n_SRIFO];
        //const double d1 = this->d_[i+1+this->n_SRIFO];
        //const double e0 = this->e_[i+0+this->n_SRIFO];
        //const double e1 = this->e_[i+1+this->n_SRIFO];
        __m256d dv = _mm256_setzero_pd();
        __m256d ev = _mm256_setzero_pd();
        {
            __m128d tmp0 = _mm_loadu_pd(&this->d_[i+this->n_SRIFO]);
            //dv = _mm256_set_m128d(tmp0,tmp0);
            dv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp0), tmp0, 1);
            __m128d tmp1 = _mm_loadu_pd(&this->e_[i+this->n_SRIFO]);
            //ev = _mm256_set_m128d(tmp1,tmp1);
            ev = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmp1), tmp1, 1);
        }

        //const double F0 = d0*std::exp(X0*logPsi0+e0*logT);
        //const double F1 = d1*std::exp(X1*logPsi1+e1*logT);
        __m256d elogTv = _mm256_mul_pd(ev,logTv);
        __m256d Fv = _mm256_fmadd_pd(Xv,logPsiv,elogTv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);


        
        //const double N0  = F0*K00/(1.0+Pr0);
        //const double N1  = F1*K01/(1.0+Pr1);
        __m256d Nv = _mm256_add_pd(Prv,onev);
        Nv = _mm256_div_pd(Fv,Nv);
        Nv = _mm256_mul_pd(K0v,Nv);

        _mm_storeu_pd(&this->Kf_[j],_mm256_castpd256_pd128(Nv));
        //std::cout<<this->Kf_[j+0]<<std::endl;
        //std::cout<<this->Kf_[j+1]<<std::endl;

    }
    else if(remainCA==3)
    {
        unsigned int i = this->n_SRICA-3;
        const unsigned int j = this->SRICA[i];
        const unsigned int m = j - this->Ikf[5] + this->Itbr[3];

        //const double Kinf = this->Kf_[j+this->offset_kinf];
        const double Kinf0 = this->Kf_[j-Ikf[5]+Ikf[8]+0];
        const double Kinf1 = this->Kf_[j-Ikf[5]+Ikf[8]+1];
        const double Kinf2 = this->Kf_[j-Ikf[5]+Ikf[8]+2];
        __m256d Kinfv = _mm256_setr_pd(Kinf0,Kinf1,Kinf2,Kinf2);

            //const double invKinf0 = 1.0/Kinf0;
            //const double invKinf1 = 1.0/Kinf1;
            //const double invKinf2 = 1.0/Kinf2;
        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);

        const double K00 = this->Kf_[j+0];
        const double K01 = this->Kf_[j+1];
        const double K02 = this->Kf_[j+2];
        __m256d K0v = _mm256_setr_pd(K00,K01,K02,K02);

        const double M0 = this->tmp_M[m+0];
        const double M1 = this->tmp_M[m+1];
        const double M2 = this->tmp_M[m+2];
        __m256d Mv = _mm256_setr_pd(M0,M1,M2,M2);

            //const double Pr0 = K00*M0*invKinf0; 
            //const double Pr1 = K01*M1*invKinf1; 
            //const double Pr2 = K02*M2*invKinf2; 
        __m256d Prv = _mm256_mul_pd(K0v,Mv);
        Prv = _mm256_mul_pd(Prv,invKinfv);


            //const double invLog10 = _mm_cvtsd_f64(_mm256_castpd256_pd128(invLog10v));
            //const double logPr0 = std::log(max(Pr0, small))*invLog10;
            //const double logPr1 = std::log(max(Pr1, small))*invLog10;
            //const double logPr2 = std::log(max(Pr2, small))*invLog10;
        __m256d logPrv = _mm256_max_pd(Prv,SRILimiterv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);

        const double expbT0 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+0+this->n_SRIFO];
        const double expbT1 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+1+this->n_SRIFO];
        const double expbT2 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+2+this->n_SRIFO];
        __m256d expbTv = _mm256_setr_pd(expbT0,expbT1,expbT2,expbT2);

        const double expTc0 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+0+this->n_SRIFO];
        const double expTc1 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+1+this->n_SRIFO];
        const double expTc2 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+2+this->n_SRIFO];
        __m256d expTcv = _mm256_setr_pd(expTc0,expTc1,expTc2,expTc2);

            //const double X0 = 1/(1 + (logPr0*logPr0));
            //const double X1 = 1/(1 + (logPr1*logPr1));
            //const double X2 = 1/(1 + (logPr2*logPr2));
        __m256d Xv = _mm256_fmadd_pd(logPrv,logPrv,onev);
        Xv = _mm256_div_pd(onev,Xv);

        const double a0 = this->a_[i+0+this->n_SRIFO];
        const double a1 = this->a_[i+1+this->n_SRIFO];
        const double a2 = this->a_[i+2+this->n_SRIFO];
        __m256d av = _mm256_setr_pd(a0,a1,a2,a2);

            //const double psi0 = a0*expbT0 + expTc0;
            //const double psi1 = a1*expbT1 + expTc1;
            //const double psi2 = a2*expbT2 + expTc2;
        __m256d psiv = _mm256_fmadd_pd(av,expbTv,expTcv);

            //const double logPsi0 = std::log(psi0);
            //const double logPsi1 = std::log(psi1);
            //const double logPsi2 = std::log(psi2);
        __m256d logPsiv = vec256_logd(psiv);

        const double d0 = this->d_[i+0+this->n_SRIFO];
        const double d1 = this->d_[i+1+this->n_SRIFO];
        const double d2 = this->d_[i+2+this->n_SRIFO];
        __m256d dv = _mm256_setr_pd(d0,d1,d2,d2);

        const double e0 = this->e_[i+0+this->n_SRIFO];
        const double e1 = this->e_[i+1+this->n_SRIFO];
        const double e2 = this->e_[i+2+this->n_SRIFO];
        __m256d ev = _mm256_setr_pd(e0,e1,e2,e2);

            //const double F0 = d0*std::exp(X0*logPsi0+e0*logT);
            //const double F1 = d1*std::exp(X1*logPsi1+e1*logT);
            //const double F2 = d2*std::exp(X2*logPsi2+e2*logT);
        __m256d elogTv = _mm256_mul_pd(ev,logTv);
        __m256d Fv = _mm256_fmadd_pd(Xv,logPsiv,elogTv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);

            //const double N0  = F0*K00/(1.0+Pr0);
            //const double N1  = F1*K01/(1.0+Pr1);
            //const double N2  = F2*K02/(1.0+Pr2);
        __m256d Nv = _mm256_add_pd(Prv,onev);
        Nv = _mm256_div_pd(Fv,Nv);
        Nv = _mm256_mul_pd(K0v,Nv);

        this->Kf_[j+0] = this->get_elem0(Nv);
        this->Kf_[j+1] = this->get_elem1(Nv);
        this->Kf_[j+2] = this->get_elem2(Nv);
        //std::cout<<this->Kf_[j+0]<<std::endl;
        //std::cout<<this->Kf_[j+1]<<std::endl;
        //std::cout<<this->Kf_[j+2]<<std::endl;

    }

}