/*---------------------------------------------------------------------------*\
  Description
      Computing the forward rate constant and the partial derivatives of SRI
      reactions. Including Kf, dKfdT and dKfdC

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"

//=============================================================================//

void FastChemistry::OptReaction::evalSRIPartialDerivative()const noexcept
{
    __m256d onev = _mm256_set1_pd(1.0);
    __m256d smallv = _mm256_set1_pd(FastChemistry::SRILimiter);
    __m256d invLog10v = _mm256_set1_pd(0.43429448190325182765112891891661);
    __m256d logTv = _mm256_set1_pd(this->logT);
    __m256d invTv = _mm256_set1_pd(this->invT);
    const double invLog10 = 0.43429448190325182765112891891661;


    const unsigned int remainFO = this->n_SRIFO%4;
    for (unsigned int i = 0;i<this->n_SRIFO-remainFO;i=i+4)
    {
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

        //const double Kinf = this->Kf_[j+this->offset_kinf];
        __m256d Kinfv = _mm256_loadu_pd(&this->Kf_[j+this->offset_kinf]);

        //const double invKinf = 1.0/Kinf;
        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);

        //const double K0 = this->Kf_[j];
        __m256d K0v = _mm256_loadu_pd(&this->Kf_[j]);

        //const double M = tmp_M[m];

        //const double Pr = K0*M*invKinf; 
        __m256d Prv = _mm256_mul_pd(K0v,invKinfv);
        {
            __m256d Mv = _mm256_loadu_pd(&this->tmp_M[m]);
            Prv = _mm256_mul_pd(Prv,Mv);
        }

        //const double logPr = std::log(max(Pr, small))*invLog10;
        __m256d logPrv = _mm256_max_pd(Prv, smallv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);

        //const double expbT = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3];
        __m256d expbTv = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3]);


        //const double expTc = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI];
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
        __m256d dv = _mm256_loadu_pd(&this->d_[i]);

        //const double e = this->e_[i];
        __m256d ev = _mm256_loadu_pd(&this->e_[i]);

        //const double logT = this->logT;
            //const double F = d*std::pow(psi, X)*std::pow(T, e);
        //const double F = d*std::exp(X*logPsi+e*logT);
        __m256d Fv = _mm256_mul_pd(ev,logTv);
        Fv = _mm256_fmadd_pd(Xv,logPsiv,Fv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);

        //const double b = this->b_[i];
        __m256d bv = _mm256_loadu_pd(&this->b_[i]);

        //const double invc = this->invc_[i];
        __m256d invcv = _mm256_loadu_pd(&this->invc_[i]);

        //const double invT = this->invT;

        //const double dpsidT = a*b*invT*invT*expbT - invc*expTc;
        __m256d dpsidTv = _mm256_mul_pd(av,bv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        __m256d invcexpTcv = _mm256_mul_pd(invcv,expTcv);
        dpsidTv = _mm256_fmsub_pd(dpsidTv,expbTv,invcexpTcv);

        //const double dlogPrdPr = Pr >= small ? invLog10*1/Pr : 0;
        __m256d cmp = _mm256_cmp_pd(Prv,smallv,_CMP_GE_OQ);
        Prv = _mm256_add_pd(Prv,_mm256_set1_pd(1e-100));
        __m256d dlogPrdPrv = _mm256_div_pd(invLog10v,Prv);
        dlogPrdPrv = _mm256_blendv_pd(_mm256_setzero_pd(),dlogPrdPrv,cmp);

        //const double dK0dT =  this->dKfdT_[j]; 
        __m256d dK0dTv = _mm256_loadu_pd(&this->dKfdT_[j]);
        //const double invOnePlusPr = 1.0/(1.0+Pr);
        __m256d invOnePlusPrv = _mm256_div_pd(onev,_mm256_add_pd(onev,Prv));

        //const double dKinfdT = this->dKfdT_[j+this->offset_kinf];
        __m256d dKinfdTv = _mm256_loadu_pd(&this->dKfdT_[j+this->offset_kinf]);

        __m256d Mv = _mm256_loadu_pd(&this->tmp_M[m]);
        {
            //const double N  = invOnePlusPr*F*K0;
            __m256d Nv = _mm256_mul_pd(invOnePlusPrv,Fv);
            Nv = _mm256_mul_pd(Nv,K0v);
            __m256d MNv = _mm256_mul_pd(Nv,Mv);
            //this->Kf_[j] = M*N; 
            _mm256_storeu_pd(&this->Kf_[j],MNv);

        }
        {
            //const double dKdT   = Pr*dKinfdT;
            __m256d dKdTv = _mm256_mul_pd(Prv,dKinfdTv);

            //const double dFdT = F*(X/psi*dpsidT + e*invT);
            __m256d dFdTv = _mm256_div_pd(Xv,psiv);
            __m256d einvTv = _mm256_mul_pd(ev,invTv);
            dFdTv = _mm256_fmadd_pd(dFdTv,dpsidTv,einvTv);
            dFdTv = _mm256_mul_pd(Fv,dFdTv);

            //const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;
            __m256d PrdKinfdTv = _mm256_mul_pd(Prv,dKinfdTv);
            __m256d dPrdTv = _mm256_fmsub_pd(Mv,dK0dTv,PrdKinfdTv);
            dPrdTv = _mm256_mul_pd(dPrdTv,invKinfv);

            //this->dKfdT_[j] = F*invOnePlusPr*dKdT + F*invOnePlusPr*invOnePlusPr*dPrdT*Kinf + K0*invOnePlusPr*dFdT*M;
            //invOnePlusPr*(F*(dKdT+invOnePlusPr*dPrdT*Kinf) + K0*dFdT*M)
            __m256d dKfdTv = _mm256_mul_pd(invOnePlusPrv,dPrdTv);
            __m256d K0dFdTMv = _mm256_mul_pd(K0v,dFdTv);
            K0dFdTMv = _mm256_mul_pd(K0dFdTMv,Mv);
            dKfdTv = _mm256_fmadd_pd(dKfdTv,Kinfv,dKdTv);
            dKfdTv = _mm256_fmadd_pd(dKfdTv,Fv,K0dFdTMv);
            dKfdTv = _mm256_mul_pd(dKfdTv,invOnePlusPrv);
            _mm256_storeu_pd(&this->dKfdT_[j],dKfdTv);
        
        }
        {
            //const double dXdPr = -(X*X)*2*logPr*dlogPrdPr;
            __m256d dXdPrv = _mm256_mul_pd(Xv,Xv);
            __m256d logPrdlogPrdPrv = _mm256_mul_pd(logPrv,dlogPrdPrv);
            dXdPrv = _mm256_mul_pd(dXdPrv,_mm256_set1_pd(-2));
            dXdPrv = _mm256_mul_pd(dXdPrv,logPrdlogPrdPrv);

            //const double dFdPr = F*logPsi*dXdPr;
            __m256d dFdPrv = _mm256_mul_pd(Fv,logPsiv);
            dFdPrv = _mm256_mul_pd(dFdPrv,dXdPrv);

            //const double N2     = Pr*dFdPr;
            __m256d N2v = _mm256_mul_pd(Prv,dFdPrv);
            //const double N1     = F*invOnePlusPr;
            __m256d N1v = _mm256_mul_pd(Fv,invOnePlusPrv);

            __m256d dKfdCv = _mm256_mul_pd(K0v,invOnePlusPrv);
            __m256d N1N2v = _mm256_add_pd(N1v,N2v);
            dKfdCv = _mm256_mul_pd(dKfdCv,N1N2v);
            _mm256_storeu_pd(&this->dKfdC_[m],dKfdCv);

            //this->dKfdC_[m] =  K0*invOnePlusPr*(N1 + N2); 
        }
    }
    if(remainFO==1)
    {
        const unsigned int i = this->n_SRIFO-1;
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

        const double Kinf = this->Kf_[j+this->offset_kinf];

        const double invKinf = 1.0/Kinf;

        const double K0 = this->Kf_[j];

        const double M = tmp_M[m];

        const double Pr = K0*M*invKinf; 

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

        const double F = d*std::exp(X*logPsi+e*logT);

        const double b = this->b_[i];

        const double invc = this->invc_[i];

        const double invT = this->invT;

        const double dpsidT = a*b*invT*invT*expbT - invc*expTc;

        const double dlogPrdPr = Pr >= FastChemistry::SRILimiter ? invLog10*1/Pr : 0;

        const double dK0dT =  this->dKfdT_[j]; 

        const double invOnePlusPr = 1.0/(1.0+Pr);

        const double dKinfdT = this->dKfdT_[j+this->offset_kinf];
        {
            const double N  = invOnePlusPr*F*K0;
            this->Kf_[j] = M*N; 
        }
        {
            const double dKdT   = Pr*dKinfdT;

            const double dFdT = F*(X/psi*dpsidT + e*invT);

            const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;

            this->dKfdT_[j] = F*invOnePlusPr*dKdT + F*invOnePlusPr*invOnePlusPr*dPrdT*Kinf + K0*invOnePlusPr*dFdT*M;
            //invOnePlusPr*(F*(dKdT+invOnePlusPr*dPrdT*Kinf) + K0*dFdT*M)
        }
        {
            const double dXdPr = -(X*X)*2*logPr*dlogPrdPr;
            const double dFdPr = F*logPsi*dXdPr;

            const double N2     = Pr*dFdPr;

            const double N1     = F*invOnePlusPr;

            this->dKfdC_[m] =  K0*invOnePlusPr*(N1 + N2); 
        }
    }
    else if (remainFO==2)
    {
        const unsigned int i = this->n_SRIFO-2;
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

        //const double Kinf = this->Kf_[j+this->offset_kinf];
        __m256d Kinfv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->Kf_[j+this->offset_kinf]);
            //Kinfv = _mm256_set_m128d(tmpv,tmpv);
            Kinfv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }

        //const double invKinf = 1.0/Kinf;
        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);

        //const double K0 = this->Kf_[j];
        __m256d K0v = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->Kf_[j]);
            //K0v = _mm256_set_m128d(tmpv,tmpv);
            K0v = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        //const double M = tmp_M[m];



        //const double Pr = K0*M*invKinf; 
        __m256d Prv = _mm256_mul_pd(K0v,invKinfv);
        {
            __m256d Mv = _mm256_setzero_pd();
            {
                __m128d tmpv = _mm_loadu_pd(&this->tmp_M[m]);
                //Mv = _mm256_set_m128d(tmpv,tmpv);
                Mv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
            }
            Prv = _mm256_mul_pd(Prv,Mv);
        }



        //const double logPr = std::log(max(Pr, small))*invLog10;
        __m256d logPrv = _mm256_max_pd(Prv, smallv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);

        //const double expbT = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3];
        __m256d expbTv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3]);
            //expbTv = _mm256_set_m128d(tmpv,tmpv);
            expbTv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }

        //const double expTc = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI];
        __m256d expTcv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI]);
            //expTcv = _mm256_set_m128d(tmpv,tmpv);
            expTcv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }

        //const double X = 1/(1 + (logPr*logPr));
        __m256d Xv = _mm256_fmadd_pd(logPrv,logPrv,onev);
        Xv = _mm256_div_pd(onev,Xv);

        //const double a = this->a_[i];
        __m256d av = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->a_[i]);
            //av = _mm256_set_m128d(tmpv,tmpv);
            av = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        //const double psi = a*expbT + expTc;
        __m256d psiv = _mm256_fmadd_pd(av,expbTv,expTcv);
        //const double logPsi = std::log(psi);
        __m256d logPsiv = vec256_logd(psiv);

        //const double d = this->d_[i];
        __m256d dv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->d_[i]);
            //dv = _mm256_set_m128d(tmpv,tmpv);
            dv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        //const double e = this->e_[i];
        __m256d ev = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->e_[i]);
            //ev = _mm256_set_m128d(tmpv,tmpv);
            ev = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        //const double logT = this->logT;
            //const double F = d*std::pow(psi, X)*std::pow(T, e);
        //const double F = d*std::exp(X*logPsi+e*logT);
        __m256d Fv = _mm256_mul_pd(ev,logTv);
        Fv = _mm256_fmadd_pd(Xv,logPsiv,Fv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);

        //const double b = this->b_[i];
        __m256d bv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->b_[i]);
            //bv = _mm256_set_m128d(tmpv,tmpv);
            bv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        //const double invc = this->invc_[i];
        __m256d invcv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->invc_[i]);
            //invcv = _mm256_set_m128d(tmpv,tmpv);
            invcv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        //const double invT = this->invT;

        //const double dpsidT = a*b*invT*invT*expbT - invc*expTc;
        __m256d dpsidTv = _mm256_mul_pd(av,bv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        __m256d invcexpTcv = _mm256_mul_pd(invcv,expTcv);
        dpsidTv = _mm256_fmsub_pd(dpsidTv,expbTv,invcexpTcv);

        //const double dlogPrdPr = Pr >= small ? invLog10*1/Pr : 0;
        __m256d cmp = _mm256_cmp_pd(Prv,smallv,_CMP_GE_OQ);
        Prv = _mm256_add_pd(Prv,_mm256_set1_pd(1e-100));
        __m256d dlogPrdPrv = _mm256_div_pd(invLog10v,Prv);
        dlogPrdPrv = _mm256_blendv_pd(_mm256_setzero_pd(),dlogPrdPrv,cmp);

        //const double dK0dT =  this->dKfdT_[j]; 
        __m256d dK0dTv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->dKfdT_[j]);
            //dK0dTv = _mm256_set_m128d(tmpv,tmpv);
            dK0dTv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }

        //const double invOnePlusPr = 1.0/(1.0+Pr);
        __m256d invOnePlusPrv = _mm256_div_pd(onev,_mm256_add_pd(onev,Prv));

        //const double dKinfdT = this->dKfdT_[j+this->offset_kinf];
        __m256d dKinfdTv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->dKfdT_[j+this->offset_kinf]);
            //dKinfdTv = _mm256_set_m128d(tmpv,tmpv);
            dKinfdTv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        __m256d Mv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->tmp_M[m]);
            //Mv = _mm256_set_m128d(tmpv,tmpv);
            Mv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
            //const double N  = invOnePlusPr*F*K0;
            __m256d Nv = _mm256_mul_pd(invOnePlusPrv,Fv);
            Nv = _mm256_mul_pd(Nv,K0v);
            __m256d MNv = _mm256_mul_pd(Nv,Mv);
            //this->Kf_[j] = M*N; 
            _mm_storeu_pd(&this->Kf_[j],_mm256_castpd256_pd128(MNv));

        }
        {
            //const double dKdT   = Pr*dKinfdT;
            __m256d dKdTv = _mm256_mul_pd(Prv,dKinfdTv);

            //const double dFdT = F*(X/psi*dpsidT + e*invT);
            __m256d dFdTv = _mm256_div_pd(Xv,psiv);
            __m256d einvTv = _mm256_mul_pd(ev,invTv);
            dFdTv = _mm256_fmadd_pd(dFdTv,dpsidTv,einvTv);
            dFdTv = _mm256_mul_pd(Fv,dFdTv);

            //const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;
            __m256d PrdKinfdTv = _mm256_mul_pd(Prv,dKinfdTv);
            __m256d dPrdTv = _mm256_fmsub_pd(Mv,dK0dTv,PrdKinfdTv);
            dPrdTv = _mm256_mul_pd(dPrdTv,invKinfv);

            //this->dKfdT_[j] = F*invOnePlusPr*dKdT + F*invOnePlusPr*invOnePlusPr*dPrdT*Kinf + K0*invOnePlusPr*dFdT*M;
            //invOnePlusPr*(F*(dKdT+invOnePlusPr*dPrdT*Kinf) + K0*dFdT*M)
            __m256d dKfdTv = _mm256_mul_pd(invOnePlusPrv,dPrdTv);
            __m256d K0dFdTMv = _mm256_mul_pd(K0v,dFdTv);
            K0dFdTMv = _mm256_mul_pd(K0dFdTMv,Mv);
            dKfdTv = _mm256_fmadd_pd(dKfdTv,Kinfv,dKdTv);
            dKfdTv = _mm256_fmadd_pd(dKfdTv,Fv,K0dFdTMv);
            dKfdTv = _mm256_mul_pd(dKfdTv,invOnePlusPrv);
            _mm_storeu_pd(&this->dKfdT_[j],_mm256_castpd256_pd128(dKfdTv));
        
        }
        {
            //const double dXdPr = -(X*X)*2*logPr*dlogPrdPr;
            __m256d dXdPrv = _mm256_mul_pd(Xv,Xv);
            __m256d logPrdlogPrdPrv = _mm256_mul_pd(logPrv,dlogPrdPrv);
            dXdPrv = _mm256_mul_pd(dXdPrv,_mm256_set1_pd(-2));
            dXdPrv = _mm256_mul_pd(dXdPrv,logPrdlogPrdPrv);

            //const double dFdPr = F*logPsi*dXdPr;
            __m256d dFdPrv = _mm256_mul_pd(Fv,logPsiv);
            dFdPrv = _mm256_mul_pd(dFdPrv,dXdPrv);

            //const double N2     = Pr*dFdPr;
            __m256d N2v = _mm256_mul_pd(Prv,dFdPrv);
            //const double N1     = F*invOnePlusPr;
            __m256d N1v = _mm256_mul_pd(Fv,invOnePlusPrv);

            __m256d dKfdCv = _mm256_mul_pd(K0v,invOnePlusPrv);
            __m256d N1N2v = _mm256_add_pd(N1v,N2v);
            dKfdCv = _mm256_mul_pd(dKfdCv,N1N2v);
            _mm_storeu_pd(&this->dKfdC_[m],_mm256_castpd256_pd128(dKfdCv));

            //this->dKfdC_[m] =  K0*invOnePlusPr*(N1 + N2); 
        }
    }
    else if(remainFO==3)
    {
        const unsigned int i = this->n_SRIFO-3;
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

        //const double Kinf = this->Kf_[j+this->offset_kinf];
        const double Kinf0 = this->Kf_[j+this->offset_kinf+0];
        const double Kinf1 = this->Kf_[j+this->offset_kinf+1];
        const double Kinf2 = this->Kf_[j+this->offset_kinf+2];
        __m256d Kinfv = _mm256_setr_pd(Kinf0,Kinf1,Kinf2,Kinf2);

        //const double invKinf = 1.0/Kinf;
        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);

        //const double K0 = this->Kf_[j];
        const double K00 = this->Kf_[j+0];
        const double K01 = this->Kf_[j+1];
        const double K02 = this->Kf_[j+2];
        __m256d K0v = _mm256_setr_pd(K00,K01,K02,K02);
 
        //const double M = tmp_M[m];



        //const double Pr = K0*M*invKinf; 
        __m256d Prv = _mm256_mul_pd(K0v,invKinfv);
        {
            const double M0 = this->tmp_M[m+0];
            const double M1 = this->tmp_M[m+1];
            const double M2 = this->tmp_M[m+2];
            __m256d Mv = _mm256_setr_pd(M0,M1,M2,M2);
            Prv = _mm256_mul_pd(Prv,Mv);
        }



        //const double logPr = std::log(max(Pr, small))*invLog10;
        __m256d logPrv = _mm256_max_pd(Prv, smallv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);

        //const double expbT = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3];
        const double expbT0 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+0];
        const double expbT1 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+1];
        const double expbT2 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+2];
        __m256d expbTv = _mm256_setr_pd(expbT0,expbT1,expbT2,expbT2);

        const double expTc0 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+0];
        const double expTc1 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+1];
        const double expTc2 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+2];
        __m256d expTcv = _mm256_setr_pd(expTc0,expTc1,expTc2,expTc2);


        //const double X = 1/(1 + (logPr*logPr));
        __m256d Xv = _mm256_fmadd_pd(logPrv,logPrv,onev);
        Xv = _mm256_div_pd(onev,Xv);

        const double a0 = this->a_[i+0];
        const double a1 = this->a_[i+1];
        const double a2 = this->a_[i+2];
        __m256d av = _mm256_setr_pd(a0,a1,a2,a2);
        //const double psi = a*expbT + expTc;
        __m256d psiv = _mm256_fmadd_pd(av,expbTv,expTcv);
        //const double logPsi = std::log(psi);
        __m256d logPsiv = vec256_logd(psiv);

        const double d0 = this->d_[i+0];
        const double d1 = this->d_[i+1];
        const double d2 = this->d_[i+2];
        __m256d dv = _mm256_setr_pd(d0,d1,d2,d2);


        const double e0 = this->e_[i+0];
        const double e1 = this->e_[i+1];
        const double e2 = this->e_[i+2];

        __m256d ev = _mm256_setr_pd(e0,e1,e2,e2);

        //const double logT = this->logT;
            //const double F = d*std::pow(psi, X)*std::pow(T, e);
        //const double F = d*std::exp(X*logPsi+e*logT);
        __m256d Fv = _mm256_mul_pd(ev,logTv);
        Fv = _mm256_fmadd_pd(Xv,logPsiv,Fv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);

        const double b0 = this->b_[i+0];
        const double b1 = this->b_[i+1];
        const double b2 = this->b_[i+2];
        __m256d bv = _mm256_setr_pd(b0,b1,b2,b2);

        const double invc0 = this->invc_[i+0];
        const double invc1 = this->invc_[i+1];
        const double invc2 = this->invc_[i+2];
        __m256d invcv = _mm256_setr_pd(invc0,invc1,invc2,invc2);

        //const double invT = this->invT;

        //const double dpsidT = a*b*invT*invT*expbT - invc*expTc;
        __m256d dpsidTv = _mm256_mul_pd(av,bv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        __m256d invcexpTcv = _mm256_mul_pd(invcv,expTcv);
        dpsidTv = _mm256_fmsub_pd(dpsidTv,expbTv,invcexpTcv);

        //const double dlogPrdPr = Pr >= small ? invLog10*1/Pr : 0;
        __m256d cmp = _mm256_cmp_pd(Prv,smallv,_CMP_GE_OQ);
        Prv = _mm256_add_pd(Prv,_mm256_set1_pd(1e-100));
        __m256d dlogPrdPrv = _mm256_div_pd(invLog10v,Prv);
        dlogPrdPrv = _mm256_blendv_pd(_mm256_setzero_pd(),dlogPrdPrv,cmp);

        const double dK0dT0 =  this->dKfdT_[j+0];
        const double dK0dT1 =  this->dKfdT_[j+1];
        const double dK0dT2 =  this->dKfdT_[j+2];
        __m256d dK0dTv = _mm256_setr_pd(dK0dT0,dK0dT1,dK0dT2,dK0dT2);


        //const double invOnePlusPr = 1.0/(1.0+Pr);
        __m256d invOnePlusPrv = _mm256_div_pd(onev,_mm256_add_pd(onev,Prv));

        const double dKinfdT0 = this->dKfdT_[j+this->offset_kinf+0];
        const double dKinfdT1 = this->dKfdT_[j+this->offset_kinf+1];
        const double dKinfdT2 = this->dKfdT_[j+this->offset_kinf+2];
        __m256d dKinfdTv = _mm256_setr_pd(dKinfdT0,dKinfdT1,dKinfdT2,dKinfdT2);


        const double M0 = this->tmp_M[m+0];
        const double M1 = this->tmp_M[m+1];
        const double M2 = this->tmp_M[m+2];
        __m256d Mv = _mm256_setr_pd(M0,M1,M2,M2);
        {

            //const double N  = invOnePlusPr*F*K0;
            __m256d Nv = _mm256_mul_pd(invOnePlusPrv,Fv);
            Nv = _mm256_mul_pd(Nv,K0v);
            __m256d MNv = _mm256_mul_pd(Nv,Mv);
            //this->Kf_[j] = M*N; 
            this->Kf_[j+0] = this->get_elem0(MNv);
            this->Kf_[j+1] = this->get_elem1(MNv);
            this->Kf_[j+2] = this->get_elem2(MNv);
        }
        {
            //const double dKdT   = Pr*dKinfdT;
            __m256d dKdTv = _mm256_mul_pd(Prv,dKinfdTv);

            //const double dFdT = F*(X/psi*dpsidT + e*invT);
            __m256d dFdTv = _mm256_div_pd(Xv,psiv);
            __m256d einvTv = _mm256_mul_pd(ev,invTv);
            dFdTv = _mm256_fmadd_pd(dFdTv,dpsidTv,einvTv);
            dFdTv = _mm256_mul_pd(Fv,dFdTv);

            //const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;
            __m256d PrdKinfdTv = _mm256_mul_pd(Prv,dKinfdTv);
            __m256d dPrdTv = _mm256_fmsub_pd(Mv,dK0dTv,PrdKinfdTv);
            dPrdTv = _mm256_mul_pd(dPrdTv,invKinfv);

            //this->dKfdT_[j] = F*invOnePlusPr*dKdT + F*invOnePlusPr*invOnePlusPr*dPrdT*Kinf + K0*invOnePlusPr*dFdT*M;
            //invOnePlusPr*(F*(dKdT+invOnePlusPr*dPrdT*Kinf) + K0*dFdT*M)
            __m256d dKfdTv = _mm256_mul_pd(invOnePlusPrv,dPrdTv);
            __m256d K0dFdTMv = _mm256_mul_pd(K0v,dFdTv);
            K0dFdTMv = _mm256_mul_pd(K0dFdTMv,Mv);
            dKfdTv = _mm256_fmadd_pd(dKfdTv,Kinfv,dKdTv);
            dKfdTv = _mm256_fmadd_pd(dKfdTv,Fv,K0dFdTMv);
            dKfdTv = _mm256_mul_pd(dKfdTv,invOnePlusPrv);
            this->dKfdT_[j+0] = this->get_elem0(dKfdTv);
            this->dKfdT_[j+1] = this->get_elem1(dKfdTv);
            this->dKfdT_[j+2] = this->get_elem2(dKfdTv);
        }
        {
            //const double dXdPr = -(X*X)*2*logPr*dlogPrdPr;
            __m256d dXdPrv = _mm256_mul_pd(Xv,Xv);
            __m256d logPrdlogPrdPrv = _mm256_mul_pd(logPrv,dlogPrdPrv);
            dXdPrv = _mm256_mul_pd(dXdPrv,_mm256_set1_pd(-2));
            dXdPrv = _mm256_mul_pd(dXdPrv,logPrdlogPrdPrv);

            //const double dFdPr = F*logPsi*dXdPr;
            __m256d dFdPrv = _mm256_mul_pd(Fv,logPsiv);
            dFdPrv = _mm256_mul_pd(dFdPrv,dXdPrv);

            //const double N2     = Pr*dFdPr;
            __m256d N2v = _mm256_mul_pd(Prv,dFdPrv);
            //const double N1     = F*invOnePlusPr;
            __m256d N1v = _mm256_mul_pd(Fv,invOnePlusPrv);

            __m256d dKfdCv = _mm256_mul_pd(K0v,invOnePlusPrv);
            __m256d N1N2v = _mm256_add_pd(N1v,N2v);
            dKfdCv = _mm256_mul_pd(dKfdCv,N1N2v);
            this->dKfdC_[m+0] = this->get_elem0(dKfdCv);
            this->dKfdC_[m+1] = this->get_elem1(dKfdCv);
            this->dKfdC_[m+2] = this->get_elem2(dKfdCv);
            //_mm_storeu_pd(&this->dKfdC_[m],_mm256_castpd256_pd128(dKfdCv));

            //this->dKfdC_[m] =  K0*invOnePlusPr*(N1 + N2); 
        }
    }


   const unsigned int remainCA = this->n_SRICA%4;
    for (unsigned int i = 0;i<this->n_SRICA-remainCA;i=i+4)
    {
        const unsigned int j = this->SRICA[i];
        const unsigned int m = j - this->Ikf[5] + this->Itbr[3];

        //const double Kinf = this->Kf_[j+this->offset_kinf];
        __m256d Kinfv = _mm256_loadu_pd(&this->Kf_[j-Ikf[5]+Ikf[8]]);

        //const double invKinf = 1.0/Kinf;
        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);

        //const double K0 = this->Kf_[j];
        __m256d K0v = _mm256_loadu_pd(&this->Kf_[j]);

        //const double M = tmp_M[m];

        //const double Pr = K0*M*invKinf; 
        __m256d Prv = _mm256_mul_pd(K0v,invKinfv);
        {
            __m256d Mv = _mm256_loadu_pd(&this->tmp_M[m]);
            Prv = _mm256_mul_pd(Prv,Mv);
        }

        //const double logPr = std::log(max(Pr, small))*invLog10;
        __m256d logPrv = _mm256_max_pd(Prv, smallv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);

        //const double expbT = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3];
        __m256d expbTv = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRIFO]);


        //const double expTc = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI];
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
        __m256d dv = _mm256_loadu_pd(&this->d_[i+this->n_SRIFO]);

        //const double e = this->e_[i];
        __m256d ev = _mm256_loadu_pd(&this->e_[i+this->n_SRIFO]);

        //const double logT = this->logT;
            //const double F = d*std::pow(psi, X)*std::pow(T, e);
        //const double F = d*std::exp(X*logPsi+e*logT);
        __m256d Fv = _mm256_mul_pd(ev,logTv);
        Fv = _mm256_fmadd_pd(Xv,logPsiv,Fv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);

        //const double b = this->b_[i];
        __m256d bv = _mm256_loadu_pd(&this->b_[i+this->n_SRIFO]);

        //const double invc = this->invc_[i];
        __m256d invcv = _mm256_loadu_pd(&this->invc_[i+this->n_SRIFO]);

        //const double invT = this->invT;

        //const double dpsidT = a*b*invT*invT*expbT - invc*expTc;
        __m256d dpsidTv = _mm256_mul_pd(av,bv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        __m256d invcexpTcv = _mm256_mul_pd(invcv,expTcv);
        dpsidTv = _mm256_fmsub_pd(dpsidTv,expbTv,invcexpTcv);

        //const double dlogPrdPr = Pr >= small ? invLog10*1/Pr : 0;
        __m256d cmp = _mm256_cmp_pd(Prv,smallv,_CMP_GE_OQ);
        Prv = _mm256_add_pd(Prv,_mm256_set1_pd(1e-100));
        __m256d dlogPrdPrv = _mm256_div_pd(invLog10v,Prv);
        dlogPrdPrv = _mm256_blendv_pd(_mm256_setzero_pd(),dlogPrdPrv,cmp);

        //const double dK0dT =  this->dKfdT_[j]; 
        __m256d dK0dTv = _mm256_loadu_pd(&this->dKfdT_[j]);
        //const double invOnePlusPr = 1.0/(1.0+Pr);
        __m256d invOnePlusPrv = _mm256_div_pd(onev,_mm256_add_pd(onev,Prv));

        //const double dKinfdT = this->dKfdT_[j+this->offset_kinf];
        __m256d dKinfdTv = _mm256_loadu_pd(&this->dKfdT_[j-Ikf[5]+Ikf[8]]);

        __m256d Mv = _mm256_loadu_pd(&this->tmp_M[m]);
        {
            //const double N  = invOnePlusPr*F*K0;
            __m256d Nv = _mm256_mul_pd(invOnePlusPrv,Fv);
            Nv = _mm256_mul_pd(Nv,K0v);

            //this->Kf_[j] = N; 
            _mm256_storeu_pd(&this->Kf_[j],Nv);

        }
        {
            //const double dKdT   = dK0dT;
            __m256d dKdTv = dK0dTv;

            //const double dFdT = F*(X/psi*dpsidT + e*invT);
            __m256d dFdTv = _mm256_div_pd(Xv,psiv);
            __m256d einvTv = _mm256_mul_pd(ev,invTv);
            dFdTv = _mm256_fmadd_pd(dFdTv,dpsidTv,einvTv);
            dFdTv = _mm256_mul_pd(Fv,dFdTv);

            //const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;
            __m256d PrdKinfdTv = _mm256_mul_pd(Prv,dKinfdTv);
            __m256d dPrdTv = _mm256_fmsub_pd(Mv,dK0dTv,PrdKinfdTv);
            dPrdTv = _mm256_mul_pd(dPrdTv,invKinfv);

            //this->dKfdT_[j] = F*invOnePlusPr*dKdT 
            //+ F*invOnePlusPr*invOnePlusPr*dPrdT*K0
            //+ K0*invOnePlusPr*dFdT;

            //this->dKfdT_[j] = 1/(1+Pr)*(F*(dKdt+dPrdT*K0/(1+Pr))+dFdT*K0)


            __m256d dKfdTv = _mm256_mul_pd(invOnePlusPrv,dPrdTv);
            __m256d K0dFdTMv = _mm256_mul_pd(K0v,dFdTv);
            dKfdTv = _mm256_fmadd_pd(dKfdTv,K0v,dKdTv);
            dKfdTv = _mm256_fmadd_pd(dKfdTv,Fv,K0dFdTMv);
            dKfdTv = _mm256_mul_pd(dKfdTv,invOnePlusPrv);
            _mm256_storeu_pd(&this->dKfdT_[j],dKfdTv);
        
        }
        {
            //const double dXdPr = -(X*X)*2*logPr*dlogPrdPr;
            __m256d dXdPrv = _mm256_mul_pd(Xv,Xv);
            __m256d logPrdlogPrdPrv = _mm256_mul_pd(logPrv,dlogPrdPrv);
            dXdPrv = _mm256_mul_pd(dXdPrv,_mm256_set1_pd(-2));
            dXdPrv = _mm256_mul_pd(dXdPrv,logPrdlogPrdPrv);

            //const double dFdPr = F*logPsi*dXdPr;
            __m256d dFdPrv = _mm256_mul_pd(Fv,logPsiv);
            dFdPrv = _mm256_mul_pd(dFdPrv,dXdPrv);

            //const double N1     = -F*invOnePlusPr;
            __m256d N1v = _mm256_mul_pd(Fv,invOnePlusPrv);

            //const double N2     = dFdPr;
            __m256d N2v = dFdPrv;

            //this->dKfdC_[m] =  K0*invOnePlusPr*K0*invKinf*(N2 - N1); 
            __m256d dKfdCv = _mm256_mul_pd(K0v,invOnePlusPrv);
            __m256d N2N1v = _mm256_sub_pd(N2v,N1v);
            __m256d K0invKfv = _mm256_mul_pd(K0v,invKinfv);
            dKfdCv = _mm256_mul_pd(dKfdCv,N2N1v);
            dKfdCv = _mm256_mul_pd(dKfdCv,K0invKfv);
            _mm256_storeu_pd(&this->dKfdC_[m],dKfdCv);
        }
    }
    if(remainCA==1)
    {
        const unsigned int i = this->n_SRICA-1;
        const unsigned int j = this->SRICA[i];
        const unsigned int m = j - this->Ikf[5] + this->Itbr[3];

        const double Kinf = this->Kf_[j-Ikf[5]+Ikf[8]];

        const double invKinf = 1.0/Kinf;

        const double K0 = this->Kf_[j];

        const double M = tmp_M[m];

        const double Pr = K0*M*invKinf; 

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

        const double F = d*std::exp(X*logPsi+e*logT);

        const double b = this->b_[i+this->n_SRIFO];

        const double invc = this->invc_[i+this->n_SRIFO];

        const double invT = this->invT;

        const double dpsidT = a*b*invT*invT*expbT - invc*expTc;

        const double dlogPrdPr = Pr >= FastChemistry::SRILimiter ? invLog10*1/Pr : 0;

        const double dK0dT =  this->dKfdT_[j]; 

        const double invOnePlusPr = 1.0/(1.0+Pr);

        const double dKinfdT = this->dKfdT_[j-Ikf[5]+Ikf[8]];
        {
            const double N  = invOnePlusPr*F*K0;
            this->Kf_[j] = N; 
        }
        {
            //const double dKdT   = k<this->n_Fall_Off_Reaction?Pr*dKinfdT:dK0dT;

            const double dFdT = F*(X/psi*dpsidT + e*invT);

            const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;

            this->dKfdT_[j] = F*invOnePlusPr*dK0dT + F*invOnePlusPr*invOnePlusPr*dPrdT*K0 + K0*invOnePlusPr*dFdT;
        }
        {
            const double dXdPr = -(X*X)*2*logPr*dlogPrdPr;
            const double dFdPr = F*logPsi*dXdPr;

            const double N1     = -F*invOnePlusPr;
            const double N2     = dFdPr;

            this->dKfdC_[m] =  K0*invOnePlusPr*K0*invKinf*(N1 + N2); 
        }
    }
    else if (remainCA==2)
    {
        const unsigned int i = this->n_SRICA-2;
        const unsigned int j = this->SRICA[i];
        const unsigned int m = j - this->Ikf[5] + this->Itbr[3];

        //const double Kinf = this->Kf_[j+this->offset_kinf];
        __m256d Kinfv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->Kf_[j-Ikf[5]+Ikf[8]]);
            //Kinfv = _mm256_set_m128d(tmpv,tmpv);
            Kinfv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }

        //const double invKinf = 1.0/Kinf;
        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);

        //const double K0 = this->Kf_[j];
        __m256d K0v = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->Kf_[j]);
            //K0v = _mm256_set_m128d(tmpv,tmpv);
            K0v = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        //const double M = tmp_M[m];



        //const double Pr = K0*M*invKinf; 
        __m256d Prv = _mm256_mul_pd(K0v,invKinfv);
        {
            __m256d Mv = _mm256_setzero_pd();
            {
                __m128d tmpv = _mm_loadu_pd(&this->tmp_M[m]);
                //Mv = _mm256_set_m128d(tmpv,tmpv);
                Mv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
            }
            Prv = _mm256_mul_pd(Prv,Mv);
        }



        //const double logPr = std::log(max(Pr, small))*invLog10;
        __m256d logPrv = _mm256_max_pd(Prv, smallv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);

        //const double expbT = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3];
        __m256d expbTv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRIFO]);
            //expbTv = _mm256_set_m128d(tmpv,tmpv);
            expbTv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }

        //const double expTc = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI];
        __m256d expTcv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+this->n_SRIFO]);
            //expTcv = _mm256_set_m128d(tmpv,tmpv);
            expTcv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }

        //const double X = 1/(1 + (logPr*logPr));
        __m256d Xv = _mm256_fmadd_pd(logPrv,logPrv,onev);
        Xv = _mm256_div_pd(onev,Xv);

        //const double a = this->a_[i];
        __m256d av = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->a_[i+this->n_SRIFO]);
            //av = _mm256_set_m128d(tmpv,tmpv);
            av = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        //const double psi = a*expbT + expTc;
        __m256d psiv = _mm256_fmadd_pd(av,expbTv,expTcv);
        //const double logPsi = std::log(psi);
        __m256d logPsiv = vec256_logd(psiv);

        //const double d = this->d_[i];
        __m256d dv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->d_[i+this->n_SRIFO]);
            //dv = _mm256_set_m128d(tmpv,tmpv);
            dv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        //const double e = this->e_[i];
        __m256d ev = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->e_[i+this->n_SRIFO]);
            //ev = _mm256_set_m128d(tmpv,tmpv);
            ev = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        //const double logT = this->logT;
            //const double F = d*std::pow(psi, X)*std::pow(T, e);
        //const double F = d*std::exp(X*logPsi+e*logT);
        __m256d Fv = _mm256_mul_pd(ev,logTv);
        Fv = _mm256_fmadd_pd(Xv,logPsiv,Fv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);

        //const double b = this->b_[i];
        __m256d bv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->b_[i+this->n_SRIFO]);
            //bv = _mm256_set_m128d(tmpv,tmpv);
            bv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        //const double invc = this->invc_[i];
        __m256d invcv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->invc_[i+this->n_SRIFO]);
            //invcv = _mm256_set_m128d(tmpv,tmpv);
            invcv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        //const double invT = this->invT;

        //const double dpsidT = a*b*invT*invT*expbT - invc*expTc;
        __m256d dpsidTv = _mm256_mul_pd(av,bv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        __m256d invcexpTcv = _mm256_mul_pd(invcv,expTcv);
        dpsidTv = _mm256_fmsub_pd(dpsidTv,expbTv,invcexpTcv);

        //const double dlogPrdPr = Pr >= small ? invLog10*1/Pr : 0;
        __m256d cmp = _mm256_cmp_pd(Prv,smallv,_CMP_GE_OQ);
        Prv = _mm256_add_pd(Prv,_mm256_set1_pd(1e-100));
        __m256d dlogPrdPrv = _mm256_div_pd(invLog10v,Prv);
        dlogPrdPrv = _mm256_blendv_pd(_mm256_setzero_pd(),dlogPrdPrv,cmp);

        //const double dK0dT =  this->dKfdT_[j]; 
        __m256d dK0dTv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->dKfdT_[j]);
            //dK0dTv = _mm256_set_m128d(tmpv,tmpv);
            dK0dTv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }

        //const double invOnePlusPr = 1.0/(1.0+Pr);
        __m256d invOnePlusPrv = _mm256_div_pd(onev,_mm256_add_pd(onev,Prv));

        //const double dKinfdT = this->dKfdT_[j+this->offset_kinf];
        __m256d dKinfdTv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->dKfdT_[j+this->offset_kinf]);
            //dKinfdTv = _mm256_set_m128d(tmpv,tmpv);
            dKinfdTv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);
        }
        __m256d Mv = _mm256_setzero_pd();
        {
            __m128d tmpv = _mm_loadu_pd(&this->tmp_M[m]);
            //Mv = _mm256_set_m128d(tmpv,tmpv);
            Mv = _mm256_insertf128_pd (_mm256_castpd128_pd256 (tmpv), tmpv, 1);

            //const double N  = invOnePlusPr*F*K0;
            __m256d Nv = _mm256_mul_pd(invOnePlusPrv,Fv);
            Nv = _mm256_mul_pd(Nv,K0v);

            //this->Kf_[j] = N; 
            _mm_storeu_pd(&this->Kf_[j],_mm256_castpd256_pd128(Nv));

        }
        {
            //const double dKdT   = Pr*dKinfdT;
            __m256d dKdTv = dK0dTv;

            //const double dFdT = F*(X/psi*dpsidT + e*invT);
            __m256d dFdTv = _mm256_div_pd(Xv,psiv);
            __m256d einvTv = _mm256_mul_pd(ev,invTv);
            dFdTv = _mm256_fmadd_pd(dFdTv,dpsidTv,einvTv);
            dFdTv = _mm256_mul_pd(Fv,dFdTv);

            //const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;
            __m256d PrdKinfdTv = _mm256_mul_pd(Prv,dKinfdTv);
            __m256d dPrdTv = _mm256_fmsub_pd(Mv,dK0dTv,PrdKinfdTv);
            dPrdTv = _mm256_mul_pd(dPrdTv,invKinfv);

            //this->dKfdT_[j] = F*invOnePlusPr*dKdT + F*invOnePlusPr*invOnePlusPr*dPrdT*Kinf + K0*invOnePlusPr*dFdT*M;
            //invOnePlusPr*(F*(dKdT+invOnePlusPr*dPrdT*Kinf) + K0*dFdT*M)
            __m256d dKfdTv = _mm256_mul_pd(invOnePlusPrv,dPrdTv);
            __m256d K0dFdTMv = _mm256_mul_pd(K0v,dFdTv);
            dKfdTv = _mm256_fmadd_pd(dKfdTv,K0v,dKdTv);
            dKfdTv = _mm256_fmadd_pd(dKfdTv,Fv,K0dFdTMv);
            dKfdTv = _mm256_mul_pd(dKfdTv,invOnePlusPrv);
            _mm_storeu_pd(&this->dKfdT_[j],_mm256_castpd256_pd128(dKfdTv));
        
        }
        {
            //const double dXdPr = -(X*X)*2*logPr*dlogPrdPr;
            __m256d dXdPrv = _mm256_mul_pd(Xv,Xv);
            __m256d logPrdlogPrdPrv = _mm256_mul_pd(logPrv,dlogPrdPrv);
            dXdPrv = _mm256_mul_pd(dXdPrv,_mm256_set1_pd(-2));
            dXdPrv = _mm256_mul_pd(dXdPrv,logPrdlogPrdPrv);

            //const double dFdPr = F*logPsi*dXdPr;
            __m256d dFdPrv = _mm256_mul_pd(Fv,logPsiv);
            dFdPrv = _mm256_mul_pd(dFdPrv,dXdPrv);

            //const double N1     = -F*invOnePlusPr;
            __m256d N1v = _mm256_mul_pd(Fv,invOnePlusPrv);
            //const double N2     = dFdPr;
            __m256d N2v = dFdPrv;


            __m256d dKfdCv = _mm256_mul_pd(K0v,invOnePlusPrv);
            __m256d N2N1v = _mm256_sub_pd(N2v,N1v);
            __m256d K0invKinfv = _mm256_mul_pd(K0v,invKinfv);
            dKfdCv = _mm256_mul_pd(dKfdCv,N2N1v);
            dKfdCv = _mm256_mul_pd(dKfdCv,K0invKinfv);

            _mm_storeu_pd(&this->dKfdC_[m],_mm256_castpd256_pd128(dKfdCv));

            //this->dKfdC_[m] =  K0*invOnePlusPr*K0*invKinf*(N2 - N1); 
        }
    }
    else if(remainCA==3)
    {
        const unsigned int i = this->n_SRICA-3;
        const unsigned int j = this->SRICA[i];
        const unsigned int m = j - this->Ikf[5] + this->Itbr[3];

        //const double Kinf = this->Kf_[j+this->offset_kinf];
        const double Kinf0 = this->Kf_[j-Ikf[5]+Ikf[8]+0];
        const double Kinf1 = this->Kf_[j-Ikf[5]+Ikf[8]+1];
        const double Kinf2 = this->Kf_[j-Ikf[5]+Ikf[8]+2];
        __m256d Kinfv = _mm256_setr_pd(Kinf0,Kinf1,Kinf2,Kinf2);

        //const double invKinf = 1.0/Kinf;
        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);

        //const double K0 = this->Kf_[j];
        const double K00 = this->Kf_[j+0];
        const double K01 = this->Kf_[j+1];
        const double K02 = this->Kf_[j+2];
        __m256d K0v = _mm256_setr_pd(K00,K01,K02,K02);
 
        //const double M = tmp_M[m];



        //const double Pr = K0*M*invKinf; 
        __m256d Prv = _mm256_mul_pd(K0v,invKinfv);
        {
            const double M0 = this->tmp_M[m+0];
            const double M1 = this->tmp_M[m+1];
            const double M2 = this->tmp_M[m+2];
            __m256d Mv = _mm256_setr_pd(M0,M1,M2,M2);
            Prv = _mm256_mul_pd(Prv,Mv);
        }



        //const double logPr = std::log(max(Pr, small))*invLog10;
        __m256d logPrv = _mm256_max_pd(Prv, smallv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);

        //const double expbT = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3];
        const double expbT0 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+0+this->n_SRIFO];
        const double expbT1 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+1+this->n_SRIFO];
        const double expbT2 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+2+this->n_SRIFO];
        __m256d expbTv = _mm256_setr_pd(expbT0,expbT1,expbT2,expbT2);

        const double expTc0 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+0+this->n_SRIFO];
        const double expTc1 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+1+this->n_SRIFO];
        const double expTc2 = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI+2+this->n_SRIFO];
        __m256d expTcv = _mm256_setr_pd(expTc0,expTc1,expTc2,expTc2);


        //const double X = 1/(1 + (logPr*logPr));
        __m256d Xv = _mm256_fmadd_pd(logPrv,logPrv,onev);
        Xv = _mm256_div_pd(onev,Xv);

        const double a0 = this->a_[i+0+this->n_SRIFO];
        const double a1 = this->a_[i+1+this->n_SRIFO];
        const double a2 = this->a_[i+2+this->n_SRIFO];
        __m256d av = _mm256_setr_pd(a0,a1,a2,a2);
        //const double psi = a*expbT + expTc;
        __m256d psiv = _mm256_fmadd_pd(av,expbTv,expTcv);
        //const double logPsi = std::log(psi);
        __m256d logPsiv = vec256_logd(psiv);

        const double d0 = this->d_[i+0+this->n_SRIFO];
        const double d1 = this->d_[i+1+this->n_SRIFO];
        const double d2 = this->d_[i+2+this->n_SRIFO];
        __m256d dv = _mm256_setr_pd(d0,d1,d2,d2);


        const double e0 = this->e_[i+0+this->n_SRIFO];
        const double e1 = this->e_[i+1+this->n_SRIFO];
        const double e2 = this->e_[i+2+this->n_SRIFO];

        __m256d ev = _mm256_setr_pd(e0,e1,e2,e2);

        //const double logT = this->logT;
            //const double F = d*std::pow(psi, X)*std::pow(T, e);
        //const double F = d*std::exp(X*logPsi+e*logT);
        __m256d Fv = _mm256_mul_pd(ev,logTv);
        Fv = _mm256_fmadd_pd(Xv,logPsiv,Fv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);

        const double b0 = this->b_[i+0+this->n_SRIFO];
        const double b1 = this->b_[i+1+this->n_SRIFO];
        const double b2 = this->b_[i+2+this->n_SRIFO];
        __m256d bv = _mm256_setr_pd(b0,b1,b2,b2);

        const double invc0 = this->invc_[i+0+this->n_SRIFO];
        const double invc1 = this->invc_[i+1+this->n_SRIFO];
        const double invc2 = this->invc_[i+2+this->n_SRIFO];
        __m256d invcv = _mm256_setr_pd(invc0,invc1,invc2,invc2);

        //const double invT = this->invT;

        //const double dpsidT = a*b*invT*invT*expbT - invc*expTc;
        __m256d dpsidTv = _mm256_mul_pd(av,bv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        __m256d invcexpTcv = _mm256_mul_pd(invcv,expTcv);
        dpsidTv = _mm256_fmsub_pd(dpsidTv,expbTv,invcexpTcv);

        //const double dlogPrdPr = Pr >= small ? invLog10*1/Pr : 0;
        __m256d cmp = _mm256_cmp_pd(Prv,smallv,_CMP_GE_OQ);
        Prv = _mm256_add_pd(Prv,_mm256_set1_pd(1e-100));
        __m256d dlogPrdPrv = _mm256_div_pd(invLog10v,Prv);
        dlogPrdPrv = _mm256_blendv_pd(_mm256_setzero_pd(),dlogPrdPrv,cmp);

        const double dK0dT0 =  this->dKfdT_[j+0];
        const double dK0dT1 =  this->dKfdT_[j+1];
        const double dK0dT2 =  this->dKfdT_[j+2];
        __m256d dK0dTv = _mm256_setr_pd(dK0dT0,dK0dT1,dK0dT2,dK0dT2);


        //const double invOnePlusPr = 1.0/(1.0+Pr);
        __m256d invOnePlusPrv = _mm256_div_pd(onev,_mm256_add_pd(onev,Prv));

        const double dKinfdT0 = this->dKfdT_[j-Ikf[5]+Ikf[8]+0];
        const double dKinfdT1 = this->dKfdT_[j-Ikf[5]+Ikf[8]+1];
        const double dKinfdT2 = this->dKfdT_[j-Ikf[5]+Ikf[8]+2];
        __m256d dKinfdTv = _mm256_setr_pd(dKinfdT0,dKinfdT1,dKinfdT2,dKinfdT2);


        const double M0 = this->tmp_M[m+0];
        const double M1 = this->tmp_M[m+1];
        const double M2 = this->tmp_M[m+2];
        __m256d Mv = _mm256_setr_pd(M0,M1,M2,M2);
        {

            //const double N  = invOnePlusPr*F*K0;
            __m256d Nv = _mm256_mul_pd(invOnePlusPrv,Fv);
            Nv = _mm256_mul_pd(Nv,K0v);
            //this->Kf_[j] = M*N; 
            this->Kf_[j+0] = this->get_elem0(Nv);
            this->Kf_[j+1] = this->get_elem1(Nv);
            this->Kf_[j+2] = this->get_elem2(Nv);
        }
        {
            //const double dKdT   = Pr*dKinfdT;
            __m256d dKdTv = dK0dTv;

            //const double dFdT = F*(X/psi*dpsidT + e*invT);
            __m256d dFdTv = _mm256_div_pd(Xv,psiv);
            __m256d einvTv = _mm256_mul_pd(ev,invTv);
            dFdTv = _mm256_fmadd_pd(dFdTv,dpsidTv,einvTv);
            dFdTv = _mm256_mul_pd(Fv,dFdTv);

            //const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;
            __m256d PrdKinfdTv = _mm256_mul_pd(Prv,dKinfdTv);
            __m256d dPrdTv = _mm256_fmsub_pd(Mv,dK0dTv,PrdKinfdTv);
            dPrdTv = _mm256_mul_pd(dPrdTv,invKinfv);

            //this->dKfdT_[j] = F*invOnePlusPr*dK0dT 
            //+ F*invOnePlusPr*invOnePlusPr*dPrdT*K 
            //+ K0*invOnePlusPr*dFdT*MM;

            __m256d dKfdTv = _mm256_mul_pd(invOnePlusPrv,dPrdTv);
            __m256d K0dFdTv = _mm256_mul_pd(K0v,dFdTv);
            dKfdTv = _mm256_fmadd_pd(dKfdTv,K0v,dKdTv);
            dKfdTv = _mm256_fmadd_pd(dKfdTv,Fv,K0dFdTv);
            dKfdTv = _mm256_mul_pd(dKfdTv,invOnePlusPrv);
            this->dKfdT_[j+0] = this->get_elem0(dKfdTv);
            this->dKfdT_[j+1] = this->get_elem1(dKfdTv);
            this->dKfdT_[j+2] = this->get_elem2(dKfdTv);
        }
        {
            //const double dXdPr = -(X*X)*2*logPr*dlogPrdPr;
            __m256d dXdPrv = _mm256_mul_pd(Xv,Xv);
            __m256d logPrdlogPrdPrv = _mm256_mul_pd(logPrv,dlogPrdPrv);
            dXdPrv = _mm256_mul_pd(dXdPrv,_mm256_set1_pd(-2));
            dXdPrv = _mm256_mul_pd(dXdPrv,logPrdlogPrdPrv);

            //const double dFdPr = F*logPsi*dXdPr;
            __m256d dFdPrv = _mm256_mul_pd(Fv,logPsiv);
            dFdPrv = _mm256_mul_pd(dFdPrv,dXdPrv);


            //const double N1     = -F*invOnePlusPr;
            __m256d N1v = _mm256_mul_pd(Fv,invOnePlusPrv);

            //const double N2     = dFdPr;
            __m256d N2v = dFdPrv;

            __m256d dKfdCv = _mm256_mul_pd(K0v,invOnePlusPrv);
            __m256d N2N1v = _mm256_sub_pd(N2v,N1v);
            __m256d K0invKinfv = _mm256_mul_pd(K0v,invKinfv);
            dKfdCv = _mm256_mul_pd(dKfdCv,N2N1v);
            dKfdCv = _mm256_mul_pd(dKfdCv,K0invKinfv);
            this->dKfdC_[m+0] = this->get_elem0(dKfdCv);
            this->dKfdC_[m+1] = this->get_elem1(dKfdCv);
            this->dKfdC_[m+2] = this->get_elem2(dKfdCv);
            //_mm_storeu_pd(&this->dKfdC_[m],_mm256_castpd256_pd128(dKfdCv));

            //this->dKfdC_[m] =  K0*invOnePlusPr*K0*invKinf*(N2 - N1); 
        }
    }
}