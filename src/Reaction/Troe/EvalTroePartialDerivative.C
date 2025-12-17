#include "OptReaction.H"
#include <immintrin.h>  

void OptReaction::evalTroePartialDerivative()const noexcept
{
    unsigned int remainFO = (this->n_TroeFO)%4;
    const double invLog10 = 4.3429448190325200e-01;

    __m256d one = _mm256_set1_pd(1.0);
    __m256d logTen = _mm256_set1_pd(2.3025850929940500);
    double f0 = 0.67;
    double f1 = 0.14;
    double logTens = 2.3025850929940500;
    for (unsigned int i = 0;i<this->n_TroeFO-remainFO;i=i+4)
    {
        const unsigned int j0 = this->TroeFO[i+0];
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];

        __m256d tmpPr = _mm256_setzero_pd();
        __m256d invKinf = _mm256_setzero_pd();
        __m256d tmplogPr = _mm256_setzero_pd();
        {
            __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
            invKinf = _mm256_div_pd(one,Kinf);
            __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);           
            __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
            tmpPr = _mm256_mul_pd(_mm256_mul_pd(M,K0),invKinf);
            __m256d small = _mm256_set1_pd(TroeLimiter);
            tmpPr = _mm256_add_pd(tmpPr,_mm256_set1_pd(1e-100));
            tmplogPr = _mm256_mul_pd(vec256_logd(_mm256_max_pd(small,tmpPr)),_mm256_set1_pd(invLog10));
        }
        const __m256d Pr = tmpPr;
        const __m256d logPr_ = tmplogPr;

        __m256d expTTss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe]);
        __m256d expTTs = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*2]);
        __m256d alpha = _mm256_loadu_pd(&this->alpha_[i]);
        __m256d Fcent = _mm256_setzero_pd();
        {
            __m256d expTTsss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies]);
            Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha),expTTsss);
        }
        __m256d x2 = _mm256_setzero_pd();
        __m256d invx1 = _mm256_setzero_pd();
        __m256d invx3 = _mm256_setzero_pd();
        __m256d x4 = _mm256_setzero_pd();
        __m256d F = _mm256_setzero_pd();
        __m256d invFcent = _mm256_setzero_pd();
        {
            __m256d small = _mm256_set1_pd(TroeLimiter);
            Fcent = _mm256_fmadd_pd(alpha,expTTs,Fcent);
            Fcent = _mm256_add_pd(expTTss,Fcent);
            __m256d tmp0 = _mm256_set1_pd(1e-100);
            tmp0 = _mm256_add_pd(Fcent,tmp0);
            invFcent = _mm256_div_pd(one,tmp0);
            __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),_mm256_set1_pd(invLog10));
            __m256d cc = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(f0),_mm256_set1_pd(0.4));
            __m256d n = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(-1.27),_mm256_set1_pd(0.75));
            __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(cc,logPr_),_mm256_set1_pd(f1),n);
            invx1 = _mm256_div_pd(one,x1);
            x2 = _mm256_mul_pd(_mm256_sub_pd(logPr_,cc),invx1);
            __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
            invx3 = _mm256_div_pd(one,x3);
            x4 = _mm256_mul_pd(logFcent,invx3);
            F = _mm256_mul_pd(x4,logTen);
            F = vec256_expd(F);
        }

        __m256d dFdT = _mm256_setzero_pd();
        __m256d dlogPrdPr = _mm256_mul_pd(Pr,logTen);
        {
            dlogPrdPr = _mm256_div_pd(one,dlogPrdPr);
            __m256d dFcentdT = _mm256_sub_pd(alpha,one);
            __m256d InvTsss = _mm256_loadu_pd(&this->invTsss_[i]);
            dFcentdT = _mm256_mul_pd(dFcentdT,InvTsss);
            __m256d expTTsss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies]);
            dFcentdT = _mm256_mul_pd(dFcentdT,expTTsss);
            __m256d InvTs = _mm256_loadu_pd(&this->invTs_[i]);
            __m256d tmp0 = _mm256_mul_pd(alpha,InvTs);
            tmp0 = _mm256_mul_pd(tmp0,expTTs);
            dFcentdT = _mm256_sub_pd(dFcentdT,tmp0);
            __m256d invT2 = _mm256_set1_pd(invT*invT);
            __m256d small = _mm256_set1_pd(TroeLimiter);
            __m256d Tss = _mm256_loadu_pd(&this->Tss_[i]);
            tmp0 = _mm256_mul_pd(Tss,invT2);
            dFcentdT = _mm256_fmadd_pd(expTTss,tmp0,dFcentdT);
            __m256d cmp2 = _mm256_cmp_pd(Fcent,small,_CMP_GE_OQ);
            __m256d invlogTen = _mm256_set1_pd(invLog10);
            __m256d dlogFcentdT = _mm256_mul_pd(dFcentdT,invFcent);
            dlogFcentdT = _mm256_mul_pd(dlogFcentdT,invlogTen);
            dlogFcentdT = _mm256_blendv_pd(_mm256_setzero_pd(), dlogFcentdT, cmp2);
            __m256d c0 = _mm256_set1_pd(-1.1762);
            dFdT = _mm256_mul_pd(c0,dlogFcentdT);
            dFdT = _mm256_mul_pd(x2,dFdT);
            __m256d c1 = _mm256_set1_pd(-f0);
            __m256d dcdT = _mm256_mul_pd(dlogFcentdT,c1);
            dFdT = _mm256_sub_pd(dcdT,dFdT);
            __m256d two = _mm256_add_pd(one,one);
            dFdT = _mm256_mul_pd(dFdT,invx1);
            dFdT = _mm256_mul_pd(x2,dFdT);
            dFdT = _mm256_mul_pd(dFdT,two);
            dFdT = _mm256_mul_pd(x4,dFdT);
            dFdT = _mm256_sub_pd(dlogFcentdT,dFdT);
            dFdT = _mm256_mul_pd(dFdT,invx3);
            dFdT = _mm256_mul_pd(F,dFdT);
            dFdT = _mm256_mul_pd(logTen,dFdT);
        }

        __m256d dFdPr = _mm256_setzero_pd();
        __m256d invOnePlusPr = _mm256_add_pd(one,Pr);
        invOnePlusPr = _mm256_div_pd(one,invOnePlusPr);
        {
            __m256d small = _mm256_set1_pd(TroeLimiter);
            __m256d cmp_result_Pr = _mm256_cmp_pd(Pr,small,_CMP_GE_OQ);
            dlogPrdPr = _mm256_blendv_pd(_mm256_setzero_pd(), dlogPrdPr, cmp_result_Pr);
            __m256d c2 = _mm256_set1_pd(-f1);
            dFdPr = _mm256_mul_pd(dlogPrdPr,c2);
            dFdPr = _mm256_mul_pd(x2,dFdPr);
            dFdPr = _mm256_sub_pd(dlogPrdPr,dFdPr);
            dFdPr = _mm256_mul_pd(dFdPr,invx1);
            dFdPr = _mm256_mul_pd(x2,dFdPr);
            __m256d two = _mm256_add_pd(one,one);
            dFdPr = _mm256_mul_pd(dFdPr,two);
            dFdPr = _mm256_mul_pd(x4,dFdPr);
            dFdPr = _mm256_mul_pd(invx3,dFdPr);
            dFdPr = -dFdPr;
            dFdPr = _mm256_mul_pd(logTen,dFdPr);
            dFdPr = _mm256_mul_pd(F,dFdPr);
        }

        __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
        __m256d dKinfdT = _mm256_loadu_pd(&this->dKfdT_[j0+this->offset_kinf]);
        __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);           
        __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
        __m256d dK0dT = _mm256_loadu_pd(&this->dKfdT_[j0]);   



        {
            __m256d dKdT = (_mm256_mul_pd(Pr,dKinfdT));
            __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(F,invOnePlusPr),dKdT);
            __m256d dPrdT = _mm256_mul_pd(_mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT)),invKinf);
            __m256d K = (Kinf);
            dKfdT = _mm256_fmadd_pd(K,_mm256_mul_pd(_mm256_mul_pd(F,_mm256_mul_pd(invOnePlusPr,invOnePlusPr)),dPrdT),dKfdT);
            dFdT = _mm256_fmadd_pd(dFdPr,dPrdT,dFdT);
            __m256d MM = (M);
            dKfdT = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),dFdT),MM,dKfdT); 
            _mm256_storeu_pd(&this->dKfdT_[j0],dKfdT);
        }
        {
            __m256d KK = (one);
            __m256d N2 = (_mm256_mul_pd(Pr,dFdPr));
            __m256d N1 = _mm256_mul_pd(F,invOnePlusPr);
            __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),KK),_mm256_add_pd(N1,N2));
            _mm256_storeu_pd(&this->dKfdC_[m0],dKfdC);
        }
        {
            __m256d N = _mm256_mul_pd(_mm256_mul_pd(F,K0),invOnePlusPr);
            __m256d KF = (_mm256_mul_pd(N,M));
            _mm256_storeu_pd(&this->Kf_[j0],KF);
        }
    }
    if(remainFO==1)
    {
        const unsigned int i = this->n_TroeFO-1;
        const unsigned int j = this->TroeFO[i];

        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];
        const double Kinf = this->Kf_[j+this->offset_kinf];
        const double dKinfdT = this->dKfdT_[j+this->offset_kinf];
        const double K0 = this->Kf_[j];
        const double M = this->tmp_M[m];
        const double Pr = K0*M/Kinf; 
        const double invKinf = 1.0/Kinf;
        const double logPr = std::log10(std::max(Pr, small));
        const double expTTsss = this->tmp_Exp[i+this->nSpecies];
        const double expTTss  = this->tmp_Exp[i+this->nSpecies+this->n_Troe];
        const double expTTs   = this->tmp_Exp[i+this->nSpecies+this->n_Troe*2];
        const double Fcent = (1 - this->alpha_[i])*expTTsss + this->alpha_[i]*expTTs + expTTss;

        double x2 = 0;
        double x4 = 0;
        double invx1 = 0;
        double invx3 = 0;
        double F = 0;
        {
            const double logFcent = std::log10(std::max(Fcent, small));
            const double c = -0.4 - 0.67*logFcent;
            const double n = 0.75 - 1.27*logFcent;
            const double x1 = n - 0.14*(logPr + c);
            invx1 = 1.0/x1;
            x2 = (logPr + c)*invx1;
            const double x3 = 1 + (x2*x2);
            invx3 = 1.0/x3;
            x4 = logFcent*invx3;
            F = std::exp(x4*logTens);
        }
        //const double F = std::pow(10, x4);

        const double dFcentdT = - (1 - this->alpha_[i])*this->invTsss_[i]*expTTsss
        - this->alpha_[i]*this->invTs_[i]*expTTs
        + this->Tss_[i]*invT*invT*expTTss;

        double dFdT = 0;
        {
            const double dlogFcentdT = Fcent >= small ? dFcentdT/Fcent*invLog10 : 0;
            const double dcdT = -0.67*dlogFcentdT;
            const double dndT = - 1.27*dlogFcentdT;
            const double dx1dT = dndT - 0.14*dcdT;
            const double dx2dT = (dcdT - x2*dx1dT)*invx1;
            const double dx3dT = 2*x2*dx2dT;
            const double dx4dT = (dlogFcentdT - x4*dx3dT)*invx3;
            dFdT = dx4dT;
        }

        double dFdPr = 0;
        {
            const double dlogPrdPr = Pr >= small ? invLog10/Pr : 0;
            const double dx1dPr = -0.14*dlogPrdPr;
            const double dx2dPr = (dlogPrdPr - x2*dx1dPr)*invx1;
            const double dx3dPr = 2*x2*dx2dPr;
            const double dx4dPr = -x4*dx3dPr*invx3;
            dFdPr = logTens*F*dx4dPr;
        }

        const double dK0dT =  this->dKfdT_[j]; 
        const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;
        dFdT = logTens*F*dFdT + dFdPr*dPrdT;

        const double dKdT   = Pr*dKinfdT;
        const double K      = Kinf      ;
        const double MM     = M         ;
        const double invOneplusPr = 1.0/(1+Pr);
        const double KK     = K0*invKinf;
        const double N1     = -F*invOneplusPr;
        const double N2     = dFdPr;
        const double N  = invOneplusPr*F*K0;
        this->dKfdT_[j] = F*invOneplusPr*dKdT + F*invOneplusPr*invOneplusPr*dPrdT*K + K0*invOneplusPr*dFdT*MM;
        this->dKfdC_[m] =  K0*invOneplusPr*KK*(N1 + N2); 
        this->Kf_[j] = M*N;   
    
    }
    else if(remainFO==2)
    {
        const unsigned int i = this->n_TroeFO-2;
        const unsigned int j0 = this->TroeFO[i+0];
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
        __m128d one128 = _mm256_castpd256_pd128(one);
        __m128d logTen128 = _mm256_castpd256_pd128(logTen);
        __m128d tmpPr = _mm_setzero_pd();
        __m128d invKinf = _mm_setzero_pd();
        __m128d tmplogPr = _mm_setzero_pd();
        {

            __m128d Kinf = _mm_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
            invKinf = _mm_div_pd(one128,Kinf);
            __m128d K0 = _mm_loadu_pd(&this->Kf_[j0]);           
            __m128d M = _mm_loadu_pd(&this->tmp_M[m0]);
            tmpPr = _mm_mul_pd(_mm_mul_pd(M,K0),invKinf);
            __m128d small = _mm_set1_pd(TroeLimiter);
            tmpPr = _mm_add_pd(tmpPr,_mm_set1_pd(1e-100));

            __m256d r = _mm256_set_m128d(_mm_max_pd(small,tmpPr),_mm_max_pd(small,tmpPr));
            r = vec256_logd(r);
            __m128d r1 = _mm256_castpd256_pd128(r);
            tmplogPr = _mm_mul_pd(r1,_mm_set1_pd(invLog10));
        }
        const __m128d Pr = tmpPr;
        const __m128d logPr_ = tmplogPr;

        __m128d expTTss = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe]);
        __m128d expTTs = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*2]);
        __m128d alpha = _mm_loadu_pd(&this->alpha_[i]);
        __m128d Fcent = _mm_setzero_pd();
        {

            __m128d expTTsss = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies]);
            Fcent  = _mm_mul_pd(_mm_sub_pd(one128,alpha),expTTsss);
        }
        __m128d x2 = _mm_setzero_pd();
        __m128d invx1 = _mm_setzero_pd();
        __m128d invx3 = _mm_setzero_pd();
        __m128d x4 = _mm_setzero_pd();
        __m128d F = _mm_setzero_pd();
        __m128d invFcent = _mm_setzero_pd();
        {
            __m128d small = _mm_set1_pd(TroeLimiter);
            Fcent = _mm_fmadd_pd(alpha,expTTs,Fcent);
            Fcent = _mm_add_pd(expTTss,Fcent);
            __m128d tmp0 = _mm_set1_pd(1e-100);
            tmp0 = _mm_add_pd(Fcent,tmp0);

            invFcent = _mm_div_pd(one128,tmp0);
            __m128d r0 = _mm_max_pd(Fcent,small);

            __m256d r1 = _mm256_set_m128d(r0,r0);
            r1 = vec256_logd(r1);
            __m128d r2 = _mm256_castpd256_pd128(r1);
            __m128d logFcent = _mm_mul_pd(r2,_mm_set1_pd(invLog10));
            __m128d cc = _mm_fmadd_pd(logFcent,_mm_set1_pd(f0),_mm_set1_pd(0.4));
            __m128d n = _mm_fmadd_pd(logFcent,_mm_set1_pd(-1.27),_mm_set1_pd(0.75));
            __m128d x1 = _mm_fmadd_pd(_mm_sub_pd(cc,logPr_),_mm_set1_pd(f1),n);
            invx1 = _mm_div_pd(one128,x1);
            x2 = _mm_mul_pd(_mm_sub_pd(logPr_,cc),invx1);
            __m128d x3 = _mm_fmadd_pd(x2,x2,one128);
            invx3 = _mm_div_pd(one128,x3);
            x4 = _mm_mul_pd(logFcent,invx3);
            F = _mm_mul_pd(x4,logTen128);
            __m256d F256 = _mm256_set_m128d(F,F);
            F256 = vec256_expd(F256);
            F = _mm256_castpd256_pd128(F256);
        }

        __m128d dFdT = _mm_setzero_pd();
        __m128d dlogPrdPr = _mm_mul_pd(Pr,logTen128);
        {

            dlogPrdPr = _mm_div_pd(one128,dlogPrdPr);
            __m128d dFcentdT = _mm_sub_pd(alpha,one128);
            __m128d InvTsss = _mm_loadu_pd(&this->invTsss_[i]);
            dFcentdT = _mm_mul_pd(dFcentdT,InvTsss);
            __m128d expTTsss = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies]);
            dFcentdT = _mm_mul_pd(dFcentdT,expTTsss);
            __m128d InvTs = _mm_loadu_pd(&this->invTs_[i]);
            __m128d tmp0 = _mm_mul_pd(alpha,InvTs);
            tmp0 = _mm_mul_pd(tmp0,expTTs);
            dFcentdT = _mm_sub_pd(dFcentdT,tmp0);
            __m128d invT2 = _mm_set1_pd(invT*invT);
            __m128d small = _mm_set1_pd(TroeLimiter);
            __m128d Tss = _mm_loadu_pd(&this->Tss_[i]);
            tmp0 = _mm_mul_pd(Tss,invT2);
            dFcentdT = _mm_fmadd_pd(expTTss,tmp0,dFcentdT);
            __m128d cmp2 = _mm_cmp_pd(Fcent,small,_CMP_GE_OQ);
            __m128d invlogTen = _mm_set1_pd(invLog10);
            __m128d dlogFcentdT = _mm_mul_pd(dFcentdT,invFcent);
            dlogFcentdT = _mm_mul_pd(dlogFcentdT,invlogTen);
            dlogFcentdT = _mm_blendv_pd(_mm_setzero_pd(), dlogFcentdT, cmp2);
            __m128d c0 = _mm_set1_pd(-1.1762);
            dFdT = _mm_mul_pd(c0,dlogFcentdT);
            dFdT = _mm_mul_pd(x2,dFdT);
            __m128d c1 = _mm_set1_pd(-f0);
            __m128d dcdT = _mm_mul_pd(dlogFcentdT,c1);
            dFdT = _mm_sub_pd(dcdT,dFdT);
            __m128d two = _mm_add_pd(one128,one128);
            dFdT = _mm_mul_pd(dFdT,invx1);
            dFdT = _mm_mul_pd(x2,dFdT);
            dFdT = _mm_mul_pd(dFdT,two);
            dFdT = _mm_mul_pd(x4,dFdT);
            dFdT = _mm_sub_pd(dlogFcentdT,dFdT);
            dFdT = _mm_mul_pd(dFdT,invx3);
            dFdT = _mm_mul_pd(F,dFdT);
            dFdT = _mm_mul_pd(logTen128,dFdT);
        }

        __m128d dFdPr = _mm_setzero_pd();
        __m128d invOnePlusPr = _mm_add_pd(one128,Pr);
        invOnePlusPr = _mm_div_pd(one128,invOnePlusPr);
        {
            __m128d small = _mm_set1_pd(TroeLimiter);
            __m128d cmp_result_Pr = _mm_cmp_pd(Pr,small,_CMP_GE_OQ);
            dlogPrdPr = _mm_blendv_pd(_mm_setzero_pd(), dlogPrdPr, cmp_result_Pr);
            __m128d c2 = _mm_set1_pd(-f1);
            dFdPr = _mm_mul_pd(dlogPrdPr,c2);
            dFdPr = _mm_mul_pd(x2,dFdPr);
            dFdPr = _mm_sub_pd(dlogPrdPr,dFdPr);
            dFdPr = _mm_mul_pd(dFdPr,invx1);
            dFdPr = _mm_mul_pd(x2,dFdPr);

            __m128d two = _mm_add_pd(one128,one128);
            dFdPr = _mm_mul_pd(dFdPr,two);
            dFdPr = _mm_mul_pd(x4,dFdPr);
            dFdPr = _mm_mul_pd(invx3,dFdPr);
            dFdPr = -dFdPr;
            dFdPr = _mm_mul_pd(logTen128,dFdPr);
            dFdPr = _mm_mul_pd(F,dFdPr);
        }

        __m128d Kinf = _mm_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
        __m128d dKinfdT = _mm_loadu_pd(&this->dKfdT_[j0+this->offset_kinf]);
        __m128d K0 = _mm_loadu_pd(&this->Kf_[j0]);           
        __m128d M = _mm_loadu_pd(&this->tmp_M[m0]);
        __m128d dK0dT = _mm_loadu_pd(&this->dKfdT_[j0]);   



        {
            __m128d dKdT = (_mm_mul_pd(Pr,dKinfdT));
            __m128d dKfdT = _mm_mul_pd(_mm_mul_pd(F,invOnePlusPr),dKdT);
            __m128d dPrdT = _mm_mul_pd(_mm_fmsub_pd(M,dK0dT,_mm_mul_pd(Pr,dKinfdT)),invKinf);
            __m128d K = (Kinf);
            dKfdT = _mm_fmadd_pd(K,_mm_mul_pd(_mm_mul_pd(F,_mm_mul_pd(invOnePlusPr,invOnePlusPr)),dPrdT),dKfdT);
            dFdT = _mm_fmadd_pd(dFdPr,dPrdT,dFdT);
            __m128d MM = (M);
            dKfdT = _mm_fmadd_pd(_mm_mul_pd(_mm_mul_pd(K0,invOnePlusPr),dFdT),MM,dKfdT); 
            _mm_storeu_pd(&this->dKfdT_[j0],dKfdT);
        }
        {
            __m128d KK = (one128);
            __m128d N2 = (_mm_mul_pd(Pr,dFdPr));
            __m128d N1 = _mm_mul_pd(F,invOnePlusPr);   
            __m128d dKfdC = _mm_mul_pd(_mm_mul_pd(_mm_mul_pd(K0,invOnePlusPr),KK),_mm_add_pd(N1,N2));
            _mm_storeu_pd(&this->dKfdC_[m0],dKfdC);
        }
        {
            __m128d N = _mm_mul_pd(_mm_mul_pd(F,K0),invOnePlusPr);
            __m128d KF = (_mm_mul_pd(N,M));
            _mm_storeu_pd(&this->Kf_[j0],KF);
        }
    }
    else if(remainFO==3)
    {
        const unsigned int i = this->n_TroeFO-3;
        const unsigned int j0 = this->TroeFO[i+0];
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];

        __m256d tmpPr = _mm256_setzero_pd();
        __m256d invKinf = _mm256_setzero_pd();
        __m256d tmplogPr = _mm256_setzero_pd();
        {
            const double Kinf0 = this->Kf_[j0+0+this->offset_kinf];
            const double Kinf1 = this->Kf_[j0+1+this->offset_kinf];
            const double Kinf2 = this->Kf_[j0+2+this->offset_kinf];
            __m256d Kinf = _mm256_setr_pd(Kinf0,Kinf1,Kinf2,1);

            invKinf = _mm256_div_pd(one,Kinf);
            const double kfj0 = this->Kf_[j0+0];
            const double kfj1 = this->Kf_[j0+1];
            const double kfj2 = this->Kf_[j0+2];
            __m256d K0 = _mm256_setr_pd(kfj0,kfj1,kfj2,1);
            const double m0s = this->tmp_M[m0+0];
            const double m1s = this->tmp_M[m0+1];
            const double m2s = this->tmp_M[m0+2];
            __m256d M = _mm256_setr_pd(m0s,m1s,m2s,1);
            tmpPr = _mm256_mul_pd(_mm256_mul_pd(M,K0),invKinf);
            __m256d small = _mm256_set1_pd(TroeLimiter);
            tmpPr = _mm256_add_pd(tmpPr,_mm256_set1_pd(1e-100));
            tmplogPr = _mm256_mul_pd(vec256_logd(_mm256_max_pd(small,tmpPr)),_mm256_set1_pd(invLog10));
        }
        const __m256d Pr = tmpPr;
        const __m256d logPr_ = tmplogPr;

        double tmpi0 = this->tmp_Exp[i+0+this->nSpecies+this->n_Troe];
        double tmpi1 = this->tmp_Exp[i+1+this->nSpecies+this->n_Troe];
        double tmpi2 = this->tmp_Exp[i+2+this->nSpecies+this->n_Troe];
        __m256d expTTss = _mm256_setr_pd(tmpi0,tmpi1,tmpi2,1);
        tmpi0 = this->tmp_Exp[i+0+this->nSpecies+this->n_Troe*2];
        tmpi1 = this->tmp_Exp[i+1+this->nSpecies+this->n_Troe*2];
        tmpi2 = this->tmp_Exp[i+2+this->nSpecies+this->n_Troe*2];
        __m256d expTTs = _mm256_setr_pd(tmpi0,tmpi1,tmpi2,1);

        const double alpha0 = this->alpha_[i+0];
        const double alpha1 = this->alpha_[i+1];
        const double alpha2 = this->alpha_[i+2];
        __m256d alpha = _mm256_setr_pd(alpha0,alpha1,alpha2,1);

        __m256d Fcent = _mm256_setzero_pd();
        {
            const double tmpi0 = this->tmp_Exp[i+0+this->nSpecies];
            const double tmpi1 = this->tmp_Exp[i+1+this->nSpecies];
            const double tmpi2 = this->tmp_Exp[i+2+this->nSpecies];

            __m256d expTTsss = _mm256_setr_pd(tmpi0,tmpi1,tmpi2,1);
            Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha),expTTsss);
        }
        __m256d x2 = _mm256_setzero_pd();
        __m256d invx1 = _mm256_setzero_pd();
        __m256d invx3 = _mm256_setzero_pd();
        __m256d x4 = _mm256_setzero_pd();
        __m256d F = _mm256_setzero_pd();
        __m256d invFcent = _mm256_setzero_pd();
        {
            __m256d small = _mm256_set1_pd(TroeLimiter);
            Fcent = _mm256_fmadd_pd(alpha,expTTs,Fcent);
            Fcent = _mm256_add_pd(expTTss,Fcent);
            __m256d tmp0 = _mm256_set1_pd(1e-100);
            tmp0 = _mm256_add_pd(Fcent,tmp0);
            invFcent = _mm256_div_pd(one,tmp0);
            __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),_mm256_set1_pd(invLog10));
            __m256d cc = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(f0),_mm256_set1_pd(0.4));
            __m256d n = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(-1.27),_mm256_set1_pd(0.75));
            __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(cc,logPr_),_mm256_set1_pd(f1),n);
            invx1 = _mm256_div_pd(one,x1);
            x2 = _mm256_mul_pd(_mm256_sub_pd(logPr_,cc),invx1);
            __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
            invx3 = _mm256_div_pd(one,x3);
            x4 = _mm256_mul_pd(logFcent,invx3);
            F = _mm256_mul_pd(x4,logTen);
            F = vec256_expd(F);
        }

        __m256d dFdT = _mm256_setzero_pd();
        __m256d dlogPrdPr = _mm256_mul_pd(Pr,logTen);
        {
            dlogPrdPr = _mm256_div_pd(one,dlogPrdPr);
            __m256d dFcentdT = _mm256_sub_pd(alpha,one);
            double tmps0 = this->invTsss_[i+0];
            double tmps1 = this->invTsss_[i+1];
            double tmps2 = this->invTsss_[i+2];
            __m256d InvTsss = _mm256_setr_pd(tmps0,tmps1,tmps2,1);
            dFcentdT = _mm256_mul_pd(dFcentdT,InvTsss);

            tmps0 = this->tmp_Exp[i+0+this->nSpecies];
            tmps1 = this->tmp_Exp[i+1+this->nSpecies];
            tmps2 = this->tmp_Exp[i+2+this->nSpecies];
            __m256d expTTsss = _mm256_setr_pd(tmps0,tmps1,tmps2,1);
            dFcentdT = _mm256_mul_pd(dFcentdT,expTTsss);

            tmps0 = this->invTs_[i+0];
            tmps1 = this->invTs_[i+1];
            tmps2 = this->invTs_[i+2];
            __m256d InvTs = _mm256_setr_pd(tmps0,tmps1,tmps2,1);
            __m256d tmp0 = _mm256_mul_pd(alpha,InvTs);
            tmp0 = _mm256_mul_pd(tmp0,expTTs);
            dFcentdT = _mm256_sub_pd(dFcentdT,tmp0);
            __m256d invT2 = _mm256_set1_pd(invT*invT);
            __m256d small = _mm256_set1_pd(TroeLimiter);

            tmps0 = this->Tss_[i+0];
            tmps1 = this->Tss_[i+1];
            tmps2 = this->Tss_[i+2];

            __m256d Tss = _mm256_setr_pd(tmps0,tmps1,tmps2,1);
            tmp0 = _mm256_mul_pd(Tss,invT2);
            dFcentdT = _mm256_fmadd_pd(expTTss,tmp0,dFcentdT);
            __m256d cmp2 = _mm256_cmp_pd(Fcent,small,_CMP_GE_OQ);
            __m256d invlogTen = _mm256_set1_pd(invLog10);
            __m256d dlogFcentdT = _mm256_mul_pd(dFcentdT,invFcent);
            dlogFcentdT = _mm256_mul_pd(dlogFcentdT,invlogTen);
            dlogFcentdT = _mm256_blendv_pd(_mm256_setzero_pd(), dlogFcentdT, cmp2);
            __m256d c0 = _mm256_set1_pd(-1.1762);
            dFdT = _mm256_mul_pd(c0,dlogFcentdT);
            dFdT = _mm256_mul_pd(x2,dFdT);
            __m256d c1 = _mm256_set1_pd(-f0);
            __m256d dcdT = _mm256_mul_pd(dlogFcentdT,c1);
            dFdT = _mm256_sub_pd(dcdT,dFdT);
            __m256d two = _mm256_add_pd(one,one);
            dFdT = _mm256_mul_pd(dFdT,invx1);
            dFdT = _mm256_mul_pd(x2,dFdT);
            dFdT = _mm256_mul_pd(dFdT,two);
            dFdT = _mm256_mul_pd(x4,dFdT);
            dFdT = _mm256_sub_pd(dlogFcentdT,dFdT);
            dFdT = _mm256_mul_pd(dFdT,invx3);
            dFdT = _mm256_mul_pd(F,dFdT);
            dFdT = _mm256_mul_pd(logTen,dFdT);
        }

        __m256d dFdPr = _mm256_setzero_pd();
        __m256d invOnePlusPr = _mm256_add_pd(one,Pr);
        invOnePlusPr = _mm256_div_pd(one,invOnePlusPr);
        {
            __m256d small = _mm256_set1_pd(TroeLimiter);
            __m256d cmp_result_Pr = _mm256_cmp_pd(Pr,small,_CMP_GE_OQ);
            dlogPrdPr = _mm256_blendv_pd(_mm256_setzero_pd(), dlogPrdPr, cmp_result_Pr);
            __m256d c2 = _mm256_set1_pd(-f1);
            dFdPr = _mm256_mul_pd(dlogPrdPr,c2);
            dFdPr = _mm256_mul_pd(x2,dFdPr);
            dFdPr = _mm256_sub_pd(dlogPrdPr,dFdPr);
            dFdPr = _mm256_mul_pd(dFdPr,invx1);
            dFdPr = _mm256_mul_pd(x2,dFdPr);
            __m256d two = _mm256_add_pd(one,one);
            dFdPr = _mm256_mul_pd(dFdPr,two);
            dFdPr = _mm256_mul_pd(x4,dFdPr);
            dFdPr = _mm256_mul_pd(invx3,dFdPr);
            dFdPr = -dFdPr;
            dFdPr = _mm256_mul_pd(logTen,dFdPr);
            dFdPr = _mm256_mul_pd(F,dFdPr);
        }

        const double kinfj0 = this->Kf_[j0+0+this->offset_kinf];
        const double kinfj1 = this->Kf_[j0+1+this->offset_kinf];
        const double kinfj2 = this->Kf_[j0+2+this->offset_kinf];
        __m256d Kinf = _mm256_setr_pd(kinfj0,kinfj1,kinfj2,1);

        const double dKinfdTj0 = this->dKfdT_[j0+0+this->offset_kinf];
        const double dKinfdTj1 = this->dKfdT_[j0+1+this->offset_kinf];
        const double dKinfdTj2 = this->dKfdT_[j0+2+this->offset_kinf];
        __m256d dKinfdT = _mm256_setr_pd(dKinfdTj0,dKinfdTj1,dKinfdTj2,1);

        const double Kfj0 = this->Kf_[j0+0];
        const double Kfj1 = this->Kf_[j0+1];
        const double Kfj2 = this->Kf_[j0+2];
        __m256d K0 = _mm256_setr_pd(Kfj0,Kfj1,Kfj2,1);

        const double m0s = this->tmp_M[m0+0];
        const double m1s = this->tmp_M[m0+1];
        const double m2s = this->tmp_M[m0+2];
        __m256d M = _mm256_setr_pd(m0s,m1s,m2s,1);

        const double dKfdTj0 = this->dKfdT_[j0+0];
        const double dKfdTj1 = this->dKfdT_[j0+1];
        const double dKfdTj2 = this->dKfdT_[j0+2];
        __m256d dK0dT = _mm256_setr_pd(dKfdTj0,dKfdTj1,dKfdTj2,1);




        {
            __m256d dKdT = ( _mm256_mul_pd(Pr,dKinfdT));
            __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(F,invOnePlusPr),dKdT);
            __m256d dPrdT = _mm256_mul_pd(_mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT)),invKinf);
            __m256d K = (Kinf);
            dKfdT = _mm256_fmadd_pd(K,_mm256_mul_pd(_mm256_mul_pd(F,_mm256_mul_pd(invOnePlusPr,invOnePlusPr)),dPrdT),dKfdT);
            dFdT = _mm256_fmadd_pd(dFdPr,dPrdT,dFdT);
            __m256d MM = (M);
            dKfdT = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),dFdT),MM,dKfdT); 
            this->dKfdT_[j0+0] = this->get_elem0(dKfdT);
            this->dKfdT_[j0+1] = this->get_elem1(dKfdT);
            this->dKfdT_[j0+2] = this->get_elem2(dKfdT);

        }
        {
            __m256d KK = (one);
            __m256d N2 = (_mm256_mul_pd(Pr,dFdPr));
            __m256d N1 = _mm256_mul_pd(F,invOnePlusPr);
   
            __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),KK),_mm256_add_pd(N1,N2));
            this->dKfdC_[m0+0] = this->get_elem0(dKfdC);
            this->dKfdC_[m0+1] = this->get_elem1(dKfdC);
            this->dKfdC_[m0+2] = this->get_elem2(dKfdC);

        }
        {
            __m256d N = _mm256_mul_pd(_mm256_mul_pd(F,K0),invOnePlusPr);
            __m256d KF = (_mm256_mul_pd(N,M));
            this->Kf_[j0+0] = this->get_elem0(KF);
            this->Kf_[j0+1] = this->get_elem1(KF);
            this->Kf_[j0+2] = this->get_elem2(KF);
        }
    }














    unsigned int remainCA = (this->n_TroeCA)%4;

    for (unsigned int i = 0;i<this->n_TroeCA-remainCA;i=i+4)
    {
        const unsigned int j0 = this->TroeCA[i+0];
        const unsigned int m0 = j0 - this->Ikf[5] + this->Itbr[3];

        __m256d tmpPr = _mm256_setzero_pd();
        __m256d invKinf = _mm256_setzero_pd();
        __m256d tmplogPr = _mm256_setzero_pd();
        {
            __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
            invKinf = _mm256_div_pd(one,Kinf);
            __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);           
            __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
            tmpPr = _mm256_mul_pd(_mm256_mul_pd(M,K0),invKinf);
            __m256d small = _mm256_set1_pd(TroeLimiter);
            tmpPr = _mm256_add_pd(tmpPr,_mm256_set1_pd(1e-100));
            tmplogPr = _mm256_mul_pd(vec256_logd(_mm256_max_pd(small,tmpPr)),_mm256_set1_pd(invLog10));
        }
        const __m256d Pr = tmpPr;
        const __m256d logPr_ = tmplogPr;

        __m256d expTTss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe+this->n_TroeFO]);
        __m256d expTTs = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*2+this->n_TroeFO]);
        __m256d alpha = _mm256_loadu_pd(&this->alpha_[i+this->n_TroeFO]);
        __m256d Fcent = _mm256_setzero_pd();
        {
            __m256d expTTsss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_TroeFO]);
            Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha),expTTsss);
        }
        __m256d x2 = _mm256_setzero_pd();
        __m256d invx1 = _mm256_setzero_pd();
        __m256d invx3 = _mm256_setzero_pd();
        __m256d x4 = _mm256_setzero_pd();
        __m256d F = _mm256_setzero_pd();
        __m256d invFcent = _mm256_setzero_pd();
        {
            __m256d small = _mm256_set1_pd(TroeLimiter);
            Fcent = _mm256_fmadd_pd(alpha,expTTs,Fcent);
            Fcent = _mm256_add_pd(expTTss,Fcent);
            __m256d tmp0 = _mm256_set1_pd(1e-100);
            tmp0 = _mm256_add_pd(Fcent,tmp0);
            invFcent = _mm256_div_pd(one,tmp0);
            __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),_mm256_set1_pd(invLog10));
            __m256d cc = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(f0),_mm256_set1_pd(0.4));
            __m256d n = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(-1.27),_mm256_set1_pd(0.75));
            __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(cc,logPr_),_mm256_set1_pd(f1),n);
            invx1 = _mm256_div_pd(one,x1);
            x2 = _mm256_mul_pd(_mm256_sub_pd(logPr_,cc),invx1);
            __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
            invx3 = _mm256_div_pd(one,x3);
            x4 = _mm256_mul_pd(logFcent,invx3);
            F = _mm256_mul_pd(x4,logTen);
            F = vec256_expd(F);
        }

        __m256d dFdT = _mm256_setzero_pd();
        __m256d dlogPrdPr = _mm256_mul_pd(Pr,logTen);
        {
            dlogPrdPr = _mm256_div_pd(one,dlogPrdPr);
            __m256d dFcentdT = _mm256_sub_pd(alpha,one);
            __m256d InvTsss = _mm256_loadu_pd(&this->invTsss_[i+this->n_TroeFO]);
            dFcentdT = _mm256_mul_pd(dFcentdT,InvTsss);
            __m256d expTTsss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_TroeFO]);
            dFcentdT = _mm256_mul_pd(dFcentdT,expTTsss);
            __m256d InvTs = _mm256_loadu_pd(&this->invTs_[i+this->n_TroeFO]);
            __m256d tmp0 = _mm256_mul_pd(alpha,InvTs);
            tmp0 = _mm256_mul_pd(tmp0,expTTs);
            dFcentdT = _mm256_sub_pd(dFcentdT,tmp0);
            __m256d invT2 = _mm256_set1_pd(invT*invT);
            __m256d small = _mm256_set1_pd(TroeLimiter);
            __m256d Tss = _mm256_loadu_pd(&this->Tss_[i+this->n_TroeFO]);
            tmp0 = _mm256_mul_pd(Tss,invT2);
            dFcentdT = _mm256_fmadd_pd(expTTss,tmp0,dFcentdT);
            __m256d cmp2 = _mm256_cmp_pd(Fcent,small,_CMP_GE_OQ);
            __m256d invlogTen = _mm256_set1_pd(invLog10);
            __m256d dlogFcentdT = _mm256_mul_pd(dFcentdT,invFcent);
            dlogFcentdT = _mm256_mul_pd(dlogFcentdT,invlogTen);
            dlogFcentdT = _mm256_blendv_pd(_mm256_setzero_pd(), dlogFcentdT, cmp2);
            __m256d c0 = _mm256_set1_pd(-1.1762);
            dFdT = _mm256_mul_pd(c0,dlogFcentdT);
            dFdT = _mm256_mul_pd(x2,dFdT);
            __m256d c1 = _mm256_set1_pd(-f0);
            __m256d dcdT = _mm256_mul_pd(dlogFcentdT,c1);
            dFdT = _mm256_sub_pd(dcdT,dFdT);
            __m256d two = _mm256_add_pd(one,one);
            dFdT = _mm256_mul_pd(dFdT,invx1);
            dFdT = _mm256_mul_pd(x2,dFdT);
            dFdT = _mm256_mul_pd(dFdT,two);
            dFdT = _mm256_mul_pd(x4,dFdT);
            dFdT = _mm256_sub_pd(dlogFcentdT,dFdT);
            dFdT = _mm256_mul_pd(dFdT,invx3);
            dFdT = _mm256_mul_pd(F,dFdT);
            dFdT = _mm256_mul_pd(logTen,dFdT);
        }

        __m256d dFdPr = _mm256_setzero_pd();
        __m256d invOnePlusPr = _mm256_add_pd(one,Pr);
        invOnePlusPr = _mm256_div_pd(one,invOnePlusPr);
        {
            __m256d small = _mm256_set1_pd(TroeLimiter);
            __m256d cmp_result_Pr = _mm256_cmp_pd(Pr,small,_CMP_GE_OQ);
            dlogPrdPr = _mm256_blendv_pd(_mm256_setzero_pd(), dlogPrdPr, cmp_result_Pr);
            __m256d c2 = _mm256_set1_pd(-f1);
            dFdPr = _mm256_mul_pd(dlogPrdPr,c2);
            dFdPr = _mm256_mul_pd(x2,dFdPr);
            dFdPr = _mm256_sub_pd(dlogPrdPr,dFdPr);
            dFdPr = _mm256_mul_pd(dFdPr,invx1);
            dFdPr = _mm256_mul_pd(x2,dFdPr);
            __m256d two = _mm256_add_pd(one,one);
            dFdPr = _mm256_mul_pd(dFdPr,two);
            dFdPr = _mm256_mul_pd(x4,dFdPr);
            dFdPr = _mm256_mul_pd(invx3,dFdPr);
            dFdPr = -dFdPr;
            dFdPr = _mm256_mul_pd(logTen,dFdPr);
            dFdPr = _mm256_mul_pd(F,dFdPr);
        }

        //__m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
        __m256d dKinfdT = _mm256_loadu_pd(&this->dKfdT_[j0+this->offset_kinf]);
        __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);           
        __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
        __m256d dK0dT = _mm256_loadu_pd(&this->dKfdT_[j0]);   




        {
            __m256d dKdT = (dK0dT);
            __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(F,invOnePlusPr),dKdT);
            __m256d dPrdT = _mm256_mul_pd(_mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT)),invKinf);
            __m256d K = (K0);
            dKfdT = _mm256_fmadd_pd(K,_mm256_mul_pd(_mm256_mul_pd(F,_mm256_mul_pd(invOnePlusPr,invOnePlusPr)),dPrdT),dKfdT);
            dFdT = _mm256_fmadd_pd(dFdPr,dPrdT,dFdT);
            __m256d MM = (one);
            dKfdT = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),dFdT),MM,dKfdT); 
            _mm256_storeu_pd(&this->dKfdT_[j0],dKfdT);
        }
        {
            __m256d KK = (_mm256_mul_pd(K0,invKinf));
            __m256d N2 = (dFdPr);
            __m256d N1 = _mm256_mul_pd(F,invOnePlusPr);
            N1 = (_mm256_sub_pd(_mm256_setzero_pd(),N1));    
            __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),KK),_mm256_add_pd(N1,N2));
            _mm256_storeu_pd(&this->dKfdC_[m0],dKfdC);
        }
        {
            __m256d N = _mm256_mul_pd(_mm256_mul_pd(F,K0),invOnePlusPr);
            __m256d KF = (N);
            _mm256_storeu_pd(&this->Kf_[j0],KF);
        }
    }
    if(remainCA==1)
    {
        const unsigned int i = this->n_TroeCA-1;
        const unsigned int j = this->TroeCA[i];

        const unsigned int m = j - this->Ikf[5] + this->Itbr[3];
        const double Kinf = this->Kf_[j+this->offset_kinf];
        const double dKinfdT = this->dKfdT_[j+this->offset_kinf];
        const double K0 = this->Kf_[j];
        const double M = this->tmp_M[m];
        const double Pr = K0*M/Kinf; 
        const double invKinf = 1.0/Kinf;
        const double logPr = std::log10(std::max(Pr, small));
        const double expTTsss = this->tmp_Exp[i+this->nSpecies+this->n_TroeFO];
        const double expTTss  = this->tmp_Exp[i+this->nSpecies+this->n_TroeCA+this->n_TroeFO];
        const double expTTs   = this->tmp_Exp[i+this->nSpecies+this->n_TroeCA*2+this->n_TroeFO];
        const double Fcent = (1 - this->alpha_[i+this->n_TroeFO])*expTTsss + 
            this->alpha_[i+this->n_TroeFO]*expTTs + expTTss;

        double x2 = 0;
        double x4 = 0;
        double invx1 = 0;
        double invx3 = 0;
        double F = 0;
        {
            const double logFcent = std::log10(std::max(Fcent, small));
            const double c = -0.4 - 0.67*logFcent;
            const double n = 0.75 - 1.27*logFcent;
            const double x1 = n - 0.14*(logPr + c);
            invx1 = 1.0/x1;
            x2 = (logPr + c)*invx1;
            const double x3 = 1 + (x2*x2);
            invx3 = 1.0/x3;
            x4 = logFcent*invx3;
            F = std::exp(x4*logTens);
        }
        //const double F = std::pow(10, x4);

        const double dFcentdT = - (1 - this->alpha_[i+this->n_TroeFO])*this->invTsss_[i+this->n_TroeFO]*expTTsss
        - this->alpha_[i+this->n_TroeFO]*this->invTs_[i+this->n_TroeFO]*expTTs
        + this->Tss_[i+this->n_TroeFO]*invT*invT*expTTss;

        double dFdT = 0;
        {
            const double dlogFcentdT = Fcent >= small ? dFcentdT/Fcent*invLog10 : 0;
            const double dcdT = -0.67*dlogFcentdT;
            const double dndT = - 1.27*dlogFcentdT;
            const double dx1dT = dndT - 0.14*dcdT;
            const double dx2dT = (dcdT - x2*dx1dT)*invx1;
            const double dx3dT = 2*x2*dx2dT;
            const double dx4dT = (dlogFcentdT - x4*dx3dT)*invx3;
            dFdT = dx4dT;
        }

        double dFdPr = 0;
        {
            const double dlogPrdPr = Pr >= small ? invLog10/Pr : 0;
            const double dx1dPr = -0.14*dlogPrdPr;
            const double dx2dPr = (dlogPrdPr - x2*dx1dPr)*invx1;
            const double dx3dPr = 2*x2*dx2dPr;
            const double dx4dPr = -x4*dx3dPr*invx3;
            dFdPr = logTens*F*dx4dPr;
        }

        const double dK0dT =  this->dKfdT_[j]; 
        const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;
        dFdT = logTens*F*dFdT + dFdPr*dPrdT;

        const double dKdT   = dK0dT;
        const double K      = K0;
        const double MM     = 1;
        const double invOneplusPr = 1.0/(1+Pr);
        const double KK     = K0*invKinf;
        const double N1     = -F*invOneplusPr;
        const double N2     = dFdPr;
        const double N  = invOneplusPr*F*K0;
        this->dKfdT_[j] = F*invOneplusPr*dKdT + F*invOneplusPr*invOneplusPr*dPrdT*K + K0*invOneplusPr*dFdT*MM;
        this->dKfdC_[m] =  K0*invOneplusPr*KK*(N1 + N2); 
        this->Kf_[j] = N;
    
    }
    else if(remainCA==2)
    {
        const unsigned int i = this->n_TroeCA-2;
        const unsigned int j0 = this->TroeCA[i];
        const unsigned int m0 = j0 - this->Ikf[5] + this->Itbr[3];
        __m128d one128 = _mm256_castpd256_pd128(one);
        __m128d logTen128 = _mm256_castpd256_pd128(logTen);
        __m128d tmpPr = _mm_setzero_pd();
        __m128d invKinf = _mm_setzero_pd();
        __m128d tmplogPr = _mm_setzero_pd();
        {

            __m128d Kinf = _mm_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
            invKinf = _mm_div_pd(one128,Kinf);
            __m128d K0 = _mm_loadu_pd(&this->Kf_[j0]);           
            __m128d M = _mm_loadu_pd(&this->tmp_M[m0]);
            tmpPr = _mm_mul_pd(_mm_mul_pd(M,K0),invKinf);
            __m128d small = _mm_set1_pd(TroeLimiter);
            tmpPr = _mm_add_pd(tmpPr,_mm_set1_pd(1e-100));

            __m256d r = _mm256_set_m128d(_mm_max_pd(small,tmpPr),_mm_max_pd(small,tmpPr));
            r = vec256_logd(r);
            __m128d r1 = _mm256_castpd256_pd128(r);
            tmplogPr = _mm_mul_pd(r1,_mm_set1_pd(invLog10));
        }
        const __m128d Pr = tmpPr;
        const __m128d logPr_ = tmplogPr;

        __m128d expTTss = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_TroeCA+this->n_TroeFO]);
        __m128d expTTs = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_TroeCA*2+this->n_TroeFO]);
        __m128d alpha = _mm_loadu_pd(&this->alpha_[i+this->n_TroeFO]);
        __m128d Fcent = _mm_setzero_pd();
        {

            __m128d expTTsss = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_TroeFO]);
            Fcent  = _mm_mul_pd(_mm_sub_pd(one128,alpha),expTTsss);
        }
        __m128d x2 = _mm_setzero_pd();
        __m128d invx1 = _mm_setzero_pd();
        __m128d invx3 = _mm_setzero_pd();
        __m128d x4 = _mm_setzero_pd();
        __m128d F = _mm_setzero_pd();
        __m128d invFcent = _mm_setzero_pd();
        {
            __m128d small = _mm_set1_pd(TroeLimiter);
            Fcent = _mm_fmadd_pd(alpha,expTTs,Fcent);
            Fcent = _mm_add_pd(expTTss,Fcent);
            __m128d tmp0 = _mm_set1_pd(1e-100);
            tmp0 = _mm_add_pd(Fcent,tmp0);

            invFcent = _mm_div_pd(one128,tmp0);
            __m128d r0 = _mm_max_pd(Fcent,small);

            __m256d r1 = _mm256_set_m128d(r0,r0);
            r1 = vec256_logd(r1);
            __m128d r2 = _mm256_castpd256_pd128(r1);
            __m128d logFcent = _mm_mul_pd(r2,_mm_set1_pd(invLog10));
            __m128d cc = _mm_fmadd_pd(logFcent,_mm_set1_pd(f0),_mm_set1_pd(0.4));
            __m128d n = _mm_fmadd_pd(logFcent,_mm_set1_pd(-1.27),_mm_set1_pd(0.75));
            __m128d x1 = _mm_fmadd_pd(_mm_sub_pd(cc,logPr_),_mm_set1_pd(f1),n);
            invx1 = _mm_div_pd(one128,x1);
            x2 = _mm_mul_pd(_mm_sub_pd(logPr_,cc),invx1);
            __m128d x3 = _mm_fmadd_pd(x2,x2,one128);
            invx3 = _mm_div_pd(one128,x3);
            x4 = _mm_mul_pd(logFcent,invx3);
            F = _mm_mul_pd(x4,logTen128);
            __m256d F256 = _mm256_set_m128d(F,F);
            F256 = vec256_expd(F256);
            F = _mm256_castpd256_pd128(F256);
        }

        __m128d dFdT = _mm_setzero_pd();
        __m128d dlogPrdPr = _mm_mul_pd(Pr,logTen128);
        {

            dlogPrdPr = _mm_div_pd(one128,dlogPrdPr);
            __m128d dFcentdT = _mm_sub_pd(alpha,one128);
            __m128d InvTsss = _mm_loadu_pd(&this->invTsss_[i+this->n_TroeFO]);
            dFcentdT = _mm_mul_pd(dFcentdT,InvTsss);
            __m128d expTTsss = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_TroeFO]);
            dFcentdT = _mm_mul_pd(dFcentdT,expTTsss);
            __m128d InvTs = _mm_loadu_pd(&this->invTs_[i+this->n_TroeFO]);
            __m128d tmp0 = _mm_mul_pd(alpha,InvTs);
            tmp0 = _mm_mul_pd(tmp0,expTTs);
            dFcentdT = _mm_sub_pd(dFcentdT,tmp0);
            __m128d invT2 = _mm_set1_pd(invT*invT);
            __m128d small = _mm_set1_pd(TroeLimiter);
            __m128d Tss = _mm_loadu_pd(&this->Tss_[i+this->n_TroeFO]);
            tmp0 = _mm_mul_pd(Tss,invT2);
            dFcentdT = _mm_fmadd_pd(expTTss,tmp0,dFcentdT);
            __m128d cmp2 = _mm_cmp_pd(Fcent,small,_CMP_GE_OQ);
            __m128d invlogTen = _mm_set1_pd(invLog10);
            __m128d dlogFcentdT = _mm_mul_pd(dFcentdT,invFcent);
            dlogFcentdT = _mm_mul_pd(dlogFcentdT,invlogTen);
            dlogFcentdT = _mm_blendv_pd(_mm_setzero_pd(), dlogFcentdT, cmp2);
            __m128d c0 = _mm_set1_pd(-1.1762);
            dFdT = _mm_mul_pd(c0,dlogFcentdT);
            dFdT = _mm_mul_pd(x2,dFdT);
            __m128d c1 = _mm_set1_pd(-f0);
            __m128d dcdT = _mm_mul_pd(dlogFcentdT,c1);
            dFdT = _mm_sub_pd(dcdT,dFdT);
            __m128d two = _mm_add_pd(one128,one128);
            dFdT = _mm_mul_pd(dFdT,invx1);
            dFdT = _mm_mul_pd(x2,dFdT);
            dFdT = _mm_mul_pd(dFdT,two);
            dFdT = _mm_mul_pd(x4,dFdT);
            dFdT = _mm_sub_pd(dlogFcentdT,dFdT);
            dFdT = _mm_mul_pd(dFdT,invx3);
            dFdT = _mm_mul_pd(F,dFdT);
            dFdT = _mm_mul_pd(logTen128,dFdT);
        }

        __m128d dFdPr = _mm_setzero_pd();
        __m128d invOnePlusPr = _mm_add_pd(one128,Pr);
        invOnePlusPr = _mm_div_pd(one128,invOnePlusPr);
        {
            __m128d small = _mm_set1_pd(TroeLimiter);
            __m128d cmp_result_Pr = _mm_cmp_pd(Pr,small,_CMP_GE_OQ);
            dlogPrdPr = _mm_blendv_pd(_mm_setzero_pd(), dlogPrdPr, cmp_result_Pr);
            __m128d c2 = _mm_set1_pd(-f1);
            dFdPr = _mm_mul_pd(dlogPrdPr,c2);
            dFdPr = _mm_mul_pd(x2,dFdPr);
            dFdPr = _mm_sub_pd(dlogPrdPr,dFdPr);
            dFdPr = _mm_mul_pd(dFdPr,invx1);
            dFdPr = _mm_mul_pd(x2,dFdPr);

            __m128d two = _mm_add_pd(one128,one128);
            dFdPr = _mm_mul_pd(dFdPr,two);
            dFdPr = _mm_mul_pd(x4,dFdPr);
            dFdPr = _mm_mul_pd(invx3,dFdPr);
            dFdPr = -dFdPr;
            dFdPr = _mm_mul_pd(logTen128,dFdPr);
            dFdPr = _mm_mul_pd(F,dFdPr);
        }

        //__m128d Kinf = _mm_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
        __m128d dKinfdT = _mm_loadu_pd(&this->dKfdT_[j0+this->offset_kinf]);
        __m128d K0 = _mm_loadu_pd(&this->Kf_[j0]);           
        __m128d M = _mm_loadu_pd(&this->tmp_M[m0]);
        __m128d dK0dT = _mm_loadu_pd(&this->dKfdT_[j0]);   
        //const unsigned int k0 = j0 - this->Ikf[4];
        //__m128d k = _mm_setr_pd(k0,k0+1);


        {
            __m128d dKdT = (dK0dT);
            __m128d dKfdT = _mm_mul_pd(_mm_mul_pd(F,invOnePlusPr),dKdT);
            __m128d dPrdT = _mm_mul_pd(_mm_fmsub_pd(M,dK0dT,_mm_mul_pd(Pr,dKinfdT)),invKinf);
            __m128d K = (K0);
            dKfdT = _mm_fmadd_pd(K,_mm_mul_pd(_mm_mul_pd(F,_mm_mul_pd(invOnePlusPr,invOnePlusPr)),dPrdT),dKfdT);
            dFdT = _mm_fmadd_pd(dFdPr,dPrdT,dFdT);
            __m128d MM = (one128);
            dKfdT = _mm_fmadd_pd(_mm_mul_pd(_mm_mul_pd(K0,invOnePlusPr),dFdT),MM,dKfdT); 
            _mm_storeu_pd(&this->dKfdT_[j0],dKfdT);
        }
        {
            __m128d KK = (_mm_mul_pd(K0,invKinf));
            __m128d N2 = (dFdPr);
            __m128d N1 = _mm_mul_pd(F,invOnePlusPr);
            N1 = (_mm_sub_pd(_mm_setzero_pd(),N1));    
            __m128d dKfdC = _mm_mul_pd(_mm_mul_pd(_mm_mul_pd(K0,invOnePlusPr),KK),_mm_add_pd(N1,N2));
            _mm_storeu_pd(&this->dKfdC_[m0],dKfdC);
        }
        {
            __m128d N = _mm_mul_pd(_mm_mul_pd(F,K0),invOnePlusPr);
            __m128d KF = (N);
            _mm_storeu_pd(&this->Kf_[j0],KF);
        }
    }
    else if(remainCA==3)
    {
        const unsigned int i = this->n_TroeCA-3;
        const unsigned int j0 = this->TroeCA[i];
        const unsigned int m0 = j0 - this->Ikf[5] + this->Itbr[3];

        __m256d tmpPr = _mm256_setzero_pd();
        __m256d invKinf = _mm256_setzero_pd();
        __m256d tmplogPr = _mm256_setzero_pd();
        {
            const double Kinf0 = this->Kf_[j0+0+this->offset_kinf];
            const double Kinf1 = this->Kf_[j0+1+this->offset_kinf];
            const double Kinf2 = this->Kf_[j0+2+this->offset_kinf];
            __m256d Kinf = _mm256_setr_pd(Kinf0,Kinf1,Kinf2,1);

            invKinf = _mm256_div_pd(one,Kinf);
            const double kfj0 = this->Kf_[j0+0];
            const double kfj1 = this->Kf_[j0+1];
            const double kfj2 = this->Kf_[j0+2];
            __m256d K0 = _mm256_setr_pd(kfj0,kfj1,kfj2,1);
            const double m0s = this->tmp_M[m0+0];
            const double m1s = this->tmp_M[m0+1];
            const double m2s = this->tmp_M[m0+2];
            __m256d M = _mm256_setr_pd(m0s,m1s,m2s,1);
            tmpPr = _mm256_mul_pd(_mm256_mul_pd(M,K0),invKinf);
            __m256d small = _mm256_set1_pd(TroeLimiter);
            tmpPr = _mm256_add_pd(tmpPr,_mm256_set1_pd(1e-100));
            tmplogPr = _mm256_mul_pd(vec256_logd(_mm256_max_pd(small,tmpPr)),_mm256_set1_pd(invLog10));
        }
        const __m256d Pr = tmpPr;
        const __m256d logPr_ = tmplogPr;

        double tmpi0 = this->tmp_Exp[i+0+this->nSpecies+this->n_TroeCA+this->n_TroeFO];
        double tmpi1 = this->tmp_Exp[i+1+this->nSpecies+this->n_TroeCA+this->n_TroeFO];
        double tmpi2 = this->tmp_Exp[i+2+this->nSpecies+this->n_TroeCA+this->n_TroeFO];
        __m256d expTTss = _mm256_setr_pd(tmpi0,tmpi1,tmpi2,1);
        tmpi0 = this->tmp_Exp[i+0+this->nSpecies+this->n_TroeCA*2+this->n_TroeFO];
        tmpi1 = this->tmp_Exp[i+1+this->nSpecies+this->n_TroeCA*2+this->n_TroeFO];
        tmpi2 = this->tmp_Exp[i+2+this->nSpecies+this->n_TroeCA*2+this->n_TroeFO];
        __m256d expTTs = _mm256_setr_pd(tmpi0,tmpi1,tmpi2,1);

        const double alpha0 = this->alpha_[i+0+this->n_TroeFO];
        const double alpha1 = this->alpha_[i+1+this->n_TroeFO];
        const double alpha2 = this->alpha_[i+2+this->n_TroeFO];
        __m256d alpha = _mm256_setr_pd(alpha0,alpha1,alpha2,1);

        __m256d Fcent = _mm256_setzero_pd();
        {
            const double tmpi0 = this->tmp_Exp[i+0+this->nSpecies+this->n_TroeFO];
            const double tmpi1 = this->tmp_Exp[i+1+this->nSpecies+this->n_TroeFO];
            const double tmpi2 = this->tmp_Exp[i+2+this->nSpecies+this->n_TroeFO];

            __m256d expTTsss = _mm256_setr_pd(tmpi0,tmpi1,tmpi2,1);
            Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha),expTTsss);
        }
        __m256d x2 = _mm256_setzero_pd();
        __m256d invx1 = _mm256_setzero_pd();
        __m256d invx3 = _mm256_setzero_pd();
        __m256d x4 = _mm256_setzero_pd();
        __m256d F = _mm256_setzero_pd();
        __m256d invFcent = _mm256_setzero_pd();
        {
            __m256d small = _mm256_set1_pd(TroeLimiter);
            Fcent = _mm256_fmadd_pd(alpha,expTTs,Fcent);
            Fcent = _mm256_add_pd(expTTss,Fcent);
            __m256d tmp0 = _mm256_set1_pd(1e-100);
            tmp0 = _mm256_add_pd(Fcent,tmp0);
            invFcent = _mm256_div_pd(one,tmp0);
            __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),_mm256_set1_pd(invLog10));
            __m256d cc = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(f0),_mm256_set1_pd(0.4));
            __m256d n = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(-1.27),_mm256_set1_pd(0.75));
            __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(cc,logPr_),_mm256_set1_pd(f1),n);
            invx1 = _mm256_div_pd(one,x1);
            x2 = _mm256_mul_pd(_mm256_sub_pd(logPr_,cc),invx1);
            __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
            invx3 = _mm256_div_pd(one,x3);
            x4 = _mm256_mul_pd(logFcent,invx3);
            F = _mm256_mul_pd(x4,logTen);
            F = vec256_expd(F);
        }

        __m256d dFdT = _mm256_setzero_pd();
        __m256d dlogPrdPr = _mm256_mul_pd(Pr,logTen);
        {
            dlogPrdPr = _mm256_div_pd(one,dlogPrdPr);
            __m256d dFcentdT = _mm256_sub_pd(alpha,one);
            double tmps0 = this->invTsss_[i+0+this->n_TroeFO];
            double tmps1 = this->invTsss_[i+1+this->n_TroeFO];
            double tmps2 = this->invTsss_[i+2+this->n_TroeFO];
            __m256d InvTsss = _mm256_setr_pd(tmps0,tmps1,tmps2,1);
            dFcentdT = _mm256_mul_pd(dFcentdT,InvTsss);

            tmps0 = this->tmp_Exp[i+0+this->nSpecies+this->n_TroeFO];
            tmps1 = this->tmp_Exp[i+1+this->nSpecies+this->n_TroeFO];
            tmps2 = this->tmp_Exp[i+2+this->nSpecies+this->n_TroeFO];
            __m256d expTTsss = _mm256_setr_pd(tmps0,tmps1,tmps2,1);
            dFcentdT = _mm256_mul_pd(dFcentdT,expTTsss);

            tmps0 = this->invTs_[i+0+this->n_TroeFO];
            tmps1 = this->invTs_[i+1+this->n_TroeFO];
            tmps2 = this->invTs_[i+2+this->n_TroeFO];
            __m256d InvTs = _mm256_setr_pd(tmps0,tmps1,tmps2,1);
            __m256d tmp0 = _mm256_mul_pd(alpha,InvTs);
            tmp0 = _mm256_mul_pd(tmp0,expTTs);
            dFcentdT = _mm256_sub_pd(dFcentdT,tmp0);
            __m256d invT2 = _mm256_set1_pd(invT*invT);
            __m256d small = _mm256_set1_pd(TroeLimiter);

            tmps0 = this->Tss_[i+0+this->n_TroeFO];
            tmps1 = this->Tss_[i+1+this->n_TroeFO];
            tmps2 = this->Tss_[i+2+this->n_TroeFO];

            __m256d Tss = _mm256_setr_pd(tmps0,tmps1,tmps2,1);
            tmp0 = _mm256_mul_pd(Tss,invT2);
            dFcentdT = _mm256_fmadd_pd(expTTss,tmp0,dFcentdT);
            __m256d cmp2 = _mm256_cmp_pd(Fcent,small,_CMP_GE_OQ);
            __m256d invlogTen = _mm256_set1_pd(invLog10);
            __m256d dlogFcentdT = _mm256_mul_pd(dFcentdT,invFcent);
            dlogFcentdT = _mm256_mul_pd(dlogFcentdT,invlogTen);
            dlogFcentdT = _mm256_blendv_pd(_mm256_setzero_pd(), dlogFcentdT, cmp2);
            __m256d c0 = _mm256_set1_pd(-1.1762);
            dFdT = _mm256_mul_pd(c0,dlogFcentdT);
            dFdT = _mm256_mul_pd(x2,dFdT);
            __m256d c1 = _mm256_set1_pd(-f0);
            __m256d dcdT = _mm256_mul_pd(dlogFcentdT,c1);
            dFdT = _mm256_sub_pd(dcdT,dFdT);
            __m256d two = _mm256_add_pd(one,one);
            dFdT = _mm256_mul_pd(dFdT,invx1);
            dFdT = _mm256_mul_pd(x2,dFdT);
            dFdT = _mm256_mul_pd(dFdT,two);
            dFdT = _mm256_mul_pd(x4,dFdT);
            dFdT = _mm256_sub_pd(dlogFcentdT,dFdT);
            dFdT = _mm256_mul_pd(dFdT,invx3);
            dFdT = _mm256_mul_pd(F,dFdT);
            dFdT = _mm256_mul_pd(logTen,dFdT);
        }

        __m256d dFdPr = _mm256_setzero_pd();
        __m256d invOnePlusPr = _mm256_add_pd(one,Pr);
        invOnePlusPr = _mm256_div_pd(one,invOnePlusPr);
        {
            __m256d small = _mm256_set1_pd(TroeLimiter);
            __m256d cmp_result_Pr = _mm256_cmp_pd(Pr,small,_CMP_GE_OQ);
            dlogPrdPr = _mm256_blendv_pd(_mm256_setzero_pd(), dlogPrdPr, cmp_result_Pr);
            __m256d c2 = _mm256_set1_pd(-f1);
            dFdPr = _mm256_mul_pd(dlogPrdPr,c2);
            dFdPr = _mm256_mul_pd(x2,dFdPr);
            dFdPr = _mm256_sub_pd(dlogPrdPr,dFdPr);
            dFdPr = _mm256_mul_pd(dFdPr,invx1);
            dFdPr = _mm256_mul_pd(x2,dFdPr);
            __m256d two = _mm256_add_pd(one,one);
            dFdPr = _mm256_mul_pd(dFdPr,two);
            dFdPr = _mm256_mul_pd(x4,dFdPr);
            dFdPr = _mm256_mul_pd(invx3,dFdPr);
            dFdPr = -dFdPr;
            dFdPr = _mm256_mul_pd(logTen,dFdPr);
            dFdPr = _mm256_mul_pd(F,dFdPr);
        }

        const double kinfj0 = this->Kf_[j0+0+this->offset_kinf];
        const double kinfj1 = this->Kf_[j0+1+this->offset_kinf];
        const double kinfj2 = this->Kf_[j0+2+this->offset_kinf];
        __m256d Kinf = _mm256_setr_pd(kinfj0,kinfj1,kinfj2,1);

        const double dKinfdTj0 = this->dKfdT_[j0+0+this->offset_kinf];
        const double dKinfdTj1 = this->dKfdT_[j0+1+this->offset_kinf];
        const double dKinfdTj2 = this->dKfdT_[j0+2+this->offset_kinf];
        __m256d dKinfdT = _mm256_setr_pd(dKinfdTj0,dKinfdTj1,dKinfdTj2,1);

        const double Kfj0 = this->Kf_[j0+0];
        const double Kfj1 = this->Kf_[j0+1];
        const double Kfj2 = this->Kf_[j0+2];
        __m256d K0 = _mm256_setr_pd(Kfj0,Kfj1,Kfj2,1);

        const double m0s = this->tmp_M[m0+0];
        const double m1s = this->tmp_M[m0+1];
        const double m2s = this->tmp_M[m0+2];
        __m256d M = _mm256_setr_pd(m0s,m1s,m2s,1);

        const double dKfdTj0 = this->dKfdT_[j0+0];
        const double dKfdTj1 = this->dKfdT_[j0+1];
        const double dKfdTj2 = this->dKfdT_[j0+2];
        __m256d dK0dT = _mm256_setr_pd(dKfdTj0,dKfdTj1,dKfdTj2,1);

        //const unsigned int k0 = j0 - this->Ikf[4];
        //__m256d k = _mm256_setr_pd(k0,k0+1,k0+2,k0+3);


        {
            __m256d dKdT = (dK0dT);
            __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(F,invOnePlusPr),dKdT);
            __m256d dPrdT = _mm256_mul_pd(_mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT)),invKinf);
            __m256d K = (K0);
            dKfdT = _mm256_fmadd_pd(K,_mm256_mul_pd(_mm256_mul_pd(F,_mm256_mul_pd(invOnePlusPr,invOnePlusPr)),dPrdT),dKfdT);
            dFdT = _mm256_fmadd_pd(dFdPr,dPrdT,dFdT);
            __m256d MM = (one);
            dKfdT = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),dFdT),MM,dKfdT); 
            this->dKfdT_[j0+0] = this->get_elem0(dKfdT);
            this->dKfdT_[j0+1] = this->get_elem1(dKfdT);
            this->dKfdT_[j0+2] = this->get_elem2(dKfdT);

        }
        {
            __m256d KK = (_mm256_mul_pd(K0,invKinf));
            __m256d N2 = (dFdPr);
            __m256d N1 = _mm256_mul_pd(F,invOnePlusPr);
            N1 = (_mm256_sub_pd(_mm256_setzero_pd(),N1));    
            __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),KK),_mm256_add_pd(N1,N2));
            this->dKfdC_[m0+0] = this->get_elem0(dKfdC);
            this->dKfdC_[m0+1] = this->get_elem1(dKfdC);
            this->dKfdC_[m0+2] = this->get_elem2(dKfdC);

        }
        {
            __m256d N = _mm256_mul_pd(_mm256_mul_pd(F,K0),invOnePlusPr);
            __m256d KF = (N);
            this->Kf_[j0+0] = this->get_elem0(KF);
            this->Kf_[j0+1] = this->get_elem1(KF);
            this->Kf_[j0+2] = this->get_elem2(KF);
        }
    }
}