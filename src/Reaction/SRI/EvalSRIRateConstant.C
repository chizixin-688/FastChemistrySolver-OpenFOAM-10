#include "OptReaction.H"
#include <immintrin.h>  

void OptReaction::evalSRIRateConstant(double T)const noexcept
{

    unsigned int remain = this->n_SRI%4;
    __m256d invLog10v = _mm256_set1_pd(0.43429448190325182765112891891661);
    __m256d SRILimiterv = _mm256_set1_pd(SRILimiter);
    __m256d onev = _mm256_set1_pd(1);
    //__m256d invTv = _mm256_set1_pd(this->invT);
    //__m256d Tv = _mm256_set1_pd(T);
    __m256d logTv = _mm256_set1_pd(this->logT);

    for (unsigned int i = 0; i<this->n_SRI-remain; i=i+4)
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
    }
    if(remain==1)
    {
        unsigned int i = this->n_SRIFO-2;
        const unsigned int j = this->SRIFO[i];
        //const unsigned int k = j - this->Ikf[4];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

        const double Kinf = this->Kf_[j+this->offset_kinf];
        const double invKinf = 1.0/Kinf;
        const double K0 = this->Kf_[j];
        const double M = tmp_M[m];

        const double Pr = K0*M*invKinf; 
        const double invLog10 = _mm_cvtsd_f64(_mm256_castpd256_pd128(invLog10v));
        const double logPr = std::log(max(Pr, small))*invLog10;

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

    }
    else if(remain==2)
    {
        unsigned int i = this->n_SRIFO-2;
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

        //const double Kinf = this->Kf_[j+this->offset_kinf];
        __m256d Kinfv = _mm256_setzero_pd();
        {
            __m128d tmp = _mm_loadu_pd(&this->Kf_[j+this->offset_kinf]);
            Kinfv = _mm256_set_m128d(tmp,tmp);
        }

        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);

        __m256d K0v = _mm256_setzero_pd();
        {
            __m128d tmp = _mm_loadu_pd(&this->Kf_[j]);
            K0v = _mm256_set_m128d(tmp,tmp);
        }

        __m256d Mv = _mm256_setzero_pd();
        {
            __m128d tmp = _mm_loadu_pd(&this->tmp_M[m]);
            Mv = _mm256_set_m128d(tmp,tmp);
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
            expbTv = _mm256_set_m128d(tmp0,tmp0);
            __m128d tmp1 = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI]);
            expTcv = _mm256_set_m128d(tmp1,tmp1);
        }

        __m256d Xv = _mm256_fmadd_pd(logPrv,logPrv,onev);
        Xv = _mm256_div_pd(onev,Xv);

        __m256d av = _mm256_setzero_pd();
        {
            __m128d tmp = _mm_loadu_pd(&this->a_[i]);
            av = _mm256_set_m128d(tmp,tmp);
        }

        __m256d psiv = _mm256_fmadd_pd(av,expbTv,expTcv);

        __m256d logPsiv = vec256_logd(psiv);

        __m256d dv = _mm256_setzero_pd();
        __m256d ev = _mm256_setzero_pd();
        {
            __m128d tmp0 = _mm_loadu_pd(&this->d_[i]);
            dv = _mm256_set_m128d(tmp0,tmp0);
            __m128d tmp1 = _mm_loadu_pd(&this->e_[i]);
            ev = _mm256_set_m128d(tmp1,tmp1);
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
    }
    else if(remain==3)
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
    }

}