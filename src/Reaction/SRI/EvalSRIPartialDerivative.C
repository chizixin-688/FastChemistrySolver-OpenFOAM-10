#include "OptReaction.H"
#include <immintrin.h>  

void OptReaction::evalSRIPartialDerivative(double T)const noexcept
{
    __m256d onev = _mm256_set1_pd(1.0);
    __m256d smallv = _mm256_set1_pd(TroeLimiter);
    __m256d invLog10v = _mm256_set1_pd(0.43429448190325182765112891891661);
    __m256d logTv = _mm256_set1_pd(this->logT);
    __m256d invTv = _mm256_set1_pd(this->invT);
    const double invLog10 = 0.43429448190325182765112891891661;


    const unsigned int remainFO = this->n_SRIFO%4;

    for (unsigned int i = 0;i<this->n_SRIFO-remainFO;i=i+4)
    {
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

        const double Kinf = this->Kf_[j+this->offset_kinf];
        __m256d Kinfv = _mm256_loadu_pd(&this->Kf_[j+this->offset_kinf]);

        const double invKinf = 1.0/Kinf;
        __m256d invKinfv = _mm256_div_pd(onev,Kinfv);

        const double K0 = this->Kf_[j];
        __m256d K0v = _mm256_loadu_pd(&this->Kf_[j]);

        const double M = tmp_M[m];



        const double Pr = K0*M*invKinf; 
        __m256d Prv = _mm256_mul_pd(K0v,invKinfv);
        {
            __m256d Mv = _mm256_loadu_pd(&this->tmp_M[m]);
            Prv = _mm256_mul_pd(Prv,Mv);
        }



        const double logPr = std::log(max(Pr, small))*invLog10;
        __m256d logPrv = _mm256_max_pd(Prv, smallv);
        logPrv = vec256_logd(logPrv);
        logPrv = _mm256_mul_pd(logPrv,invLog10v);

        const double expbT = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3];
        __m256d expbTv = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3]);


        const double expTc = this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI];
        __m256d expTcv = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*3+this->n_SRI]);

        const double X = 1/(1 + (logPr*logPr));
        __m256d Xv = _mm256_fmadd_pd(logPrv,logPrv,onev);
        Xv = _mm256_div_pd(onev,Xv);

        const double a = this->a_[i];
        __m256d av = _mm256_loadu_pd(&this->a_[i]);

        const double psi = a*expbT + expTc;
        __m256d psiv = _mm256_fmadd_pd(av,expbTv,expTcv);
        const double logPsi = std::log(psi);
        __m256d logPsiv = vec256_logd(psiv);

        const double d = this->d_[i];
        __m256d dv = _mm256_loadu_pd(&this->d_[i]);

        const double e = this->e_[i];
        __m256d ev = _mm256_loadu_pd(&this->e_[i]);

        const double logT = this->logT;
            //const double F = d*std::pow(psi, X)*std::pow(T, e);
        const double F = d*std::exp(X*logPsi+e*logT);
        __m256d Fv = _mm256_mul_pd(ev,logTv);
        Fv = _mm256_fmadd_pd(Xv,logPsiv,Fv);
        Fv = vec256_expd(Fv);
        Fv = _mm256_mul_pd(Fv,dv);

        const double b = this->b_[i];
        __m256d bv = _mm256_loadu_pd(&this->b_[i]);

        const double invc = this->invc_[i];
        __m256d invcv = _mm256_loadu_pd(&this->invc_[i]);

        const double invT = this->invT;

        const double dpsidT = a*b*invT*invT*expbT - invc*expTc;
        __m256d dpsidTv = _mm256_mul_pd(av,bv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        dpsidTv = _mm256_mul_pd(dpsidTv,invTv);
        __m256d invcexpTcv = _mm256_mul_pd(invcv,expTcv);
        dpsidTv = _mm256_fmsub_pd(dpsidTv,expbTv,invcexpTcv);

        const double dlogPrdPr = Pr >= small ? invLog10*1/Pr : 0;
        __m256d cmp = _mm256_cmp_pd(Prv,smallv,_CMP_GE_OQ);
        Prv = _mm256_add_pd(Prv,_mm256_set1_pd(1e-100));
        __m256d dlogPrdPrv = _mm256_div_pd(invLog10v,Prv);
        dlogPrdPrv = _mm256_blendv_pd(_mm256_setzero_pd(),dlogPrdPrv,cmp);

        const double dK0dT =  this->dKfdT_[j]; 
        __m256d dK0dTv = _mm256_loadu_pd(&this->dKfdT_[j]);
        const double invOnePlusPr = 1.0/(1.0+Pr);
        __m256d invOnePlusPrv = _mm256_div_pd(onev,_mm256_add_pd(onev,Prv));

        const double dKinfdT = this->dKfdT_[j+this->offset_kinf];
        __m256d dKinfdTv = _mm256_loadu_pd(&this->dKfdT_[j+this->offset_kinf]);

        __m256d Mv = _mm256_loadu_pd(&this->tmp_M[m]);
        {
            const double N  = invOnePlusPr*F*K0;
            __m256d Nv = _mm256_mul_pd(invOnePlusPrv,Fv);
            Nv = _mm256_mul_pd(Nv,K0v);
            __m256d MNv = _mm256_mul_pd(Nv,Mv);
            this->Kf_[j] = M*N; 
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
        std::cout<<this->Kf_[j+0]<<" "<<this->dKfdT_[j+0]<<" "<<this->dKfdC_[m+0]<<std::endl;
        std::cout<<this->Kf_[j+1]<<" "<<this->dKfdT_[j+1]<<" "<<this->dKfdC_[m+1]<<std::endl;
        std::cout<<this->Kf_[j+2]<<" "<<this->dKfdT_[j+2]<<" "<<this->dKfdC_[m+2]<<std::endl;
        std::cout<<this->Kf_[j+3]<<" "<<this->dKfdT_[j+3]<<" "<<this->dKfdC_[m+3]<<std::endl;
    }
    std::exit(0);

}