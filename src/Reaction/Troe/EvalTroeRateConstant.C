/*---------------------------------------------------------------------------*\
  Description
      Computing the forward rate constant of Troe reactions

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"

//=============================================================================//

void FastChemistry::OptReaction::evalTroeRateConstant()const noexcept
{
    unsigned int remainFO = this->n_TroeFO%4;


    __m256d one = _mm256_set1_pd(1.0);
    __m256d invLog10v = _mm256_set1_pd(0.43429448190325182765112891891661);
    __m256d n0v = _mm256_set1_pd(0.67);
    __m256d n1 = _mm256_set1_pd(0.4);
    __m256d n2 = _mm256_set1_pd(-1.27);
    __m256d n3 = _mm256_set1_pd(0.75);
    __m256d n4 = _mm256_set1_pd(0.14);
    __m256d n5 = _mm256_set1_pd(10);
    __m256d small = _mm256_set1_pd(FastChemistry::TroeLimiter);
        for(unsigned int i = 0; i < this->n_TroeFO-remainFO;i=i+4)
        {
            const unsigned int j = this->TroeFO[i+0];
            const unsigned int m0 = j - this->Ikf[4] + Itbr[2];
            __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j+this->offset_kinf]);
            __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
            __m256d K0 = _mm256_loadu_pd(&this->Kf_[j]);    
            __m256d Pr_ = _mm256_div_pd(_mm256_mul_pd(K0,M),Kinf);
            Pr_ = _mm256_max_pd(small,Pr_);
            __m256d logPr_ = _mm256_mul_pd(vec256_logd(Pr_),invLog10v);
            __m256d alpha = _mm256_loadu_pd(&this->alpha_[i]);
            __m256d invOnePlusPr = _mm256_div_pd(one,_mm256_add_pd(one,Pr_));
            __m256d expTTsss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies]);
            __m256d expTTss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe]);
            __m256d expTTs = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*2]);

            __m256d Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha), expTTsss);
            Fcent = _mm256_fmadd_pd(alpha, expTTs,Fcent);
            Fcent = _mm256_add_pd(expTTss,Fcent);

            __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),invLog10v);
            __m256d cc = _mm256_fmadd_pd(logFcent,n0v,n1);
            __m256d n = _mm256_fmadd_pd(logFcent,n2,n3);
            __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(cc,logPr_),n4,n);
            __m256d x2 = _mm256_div_pd(_mm256_sub_pd(logPr_,cc),x1);
            __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
            __m256d x4 = _mm256_div_pd(logFcent,x3);
            __m256d F_ = vec256_powd(n5,x4);
            __m256d N = _mm256_mul_pd(_mm256_mul_pd(K0,F_),invOnePlusPr);

            __m256d Kf = _mm256_mul_pd(M,N);
            _mm256_storeu_pd(&this->Kf_[j],Kf);
        }        
        if(remainFO==1)       
        {
            const unsigned int i = this->n_TroeFO-1;
            const unsigned int j = this->TroeFO[i+0];
            const unsigned int m0 = j - this->Ikf[4] + this->Itbr[2];
            const double Kinf0 = this->Kf_[j+this->offset_kinf];
            double M0 = this->tmp_M[m0];   
            const double K00 = this->Kf_[j];
            const double Pr0 = K00*M0/Kinf0;  
            const double logPr0 = std::log10(std::max(Pr0, FastChemistry::TroeLimiter));
            const double expTTsss0 = this->tmp_Exp[i + 0 + this->nSpecies];
            const double expTTss0  = this->tmp_Exp[i + 0 + this->nSpecies + this->n_Troe];
            const double expTTs0   = this->tmp_Exp[i + 0 + this->nSpecies + this->n_Troe*2];
            const double Fcent0 =(1 - this->alpha_[i+0])*expTTsss0 + this->alpha_[i+0]*expTTs0 + expTTss0;
            const double logFcent0 = std::log10(std::max(Fcent0, FastChemistry::TroeLimiter));
            const double c0 = -0.4 - 0.67*logFcent0;
            const double n0 = 0.75 - 1.27*logFcent0;
            const double x1_0 = n0 - 0.14*(logPr0 + c0);
            const double x2_0 = (logPr0 + c0)/x1_0;
            const double x3_0 = 1 + (x2_0*x2_0);
            const double x4_0 = logFcent0/x3_0;
            const double F0 = std::pow(10, x4_0);
            const double N0  = 1/(1+Pr0)*F0*K00;
            this->Kf_[j] = M0*N0;   
        }
        else if(remainFO==2)  
        {
            const unsigned int i = this->n_TroeFO-2;
            const unsigned int j = this->TroeFO[i+0];  
            const unsigned int m0 = j - this->Ikf[4] + Itbr[2];
            __m128d one128 = _mm256_castpd256_pd128(one);
            __m128d Kinf = _mm_loadu_pd(&this->Kf_[j+this->offset_kinf]);
            __m128d M = _mm_loadu_pd(&this->tmp_M[m0]);
            __m128d K0 = _mm_loadu_pd(&this->Kf_[j]);    
            __m256d Pr_ = _mm256_insertf128_pd (_mm256_setzero_pd (), _mm_div_pd(_mm_mul_pd(K0,M),Kinf), 0);

            Pr_ = _mm256_max_pd(small,Pr_);

            __m128d logPr_ = _mm256_castpd256_pd128(_mm256_mul_pd(vec256_logd((Pr_)),(invLog10v)));
            __m128d alpha = _mm_loadu_pd(&this->alpha_[i]);
            __m128d invOnePlusPr = _mm_div_pd(one128,_mm_add_pd(one128,_mm256_castpd256_pd128(Pr_)));
            __m128d expTTsss = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies]);
            __m128d expTTss = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe]);
            __m128d expTTs = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*2]);

            __m128d Fcent  = _mm_mul_pd(_mm_sub_pd(one128,alpha), expTTsss);
            Fcent = _mm_fmadd_pd(alpha, expTTs,Fcent);
            Fcent = _mm_add_pd(expTTss,Fcent);

            __m128d logFcent = _mm_mul_pd
            (
                _mm256_castpd256_pd128(vec256_logd
                (
                    _mm256_max_pd
                    (
                        _mm256_insertf128_pd (_mm256_setzero_pd (), Fcent, 0),
                        small
                    )
                )),
                _mm256_castpd256_pd128(invLog10v)
            );
            __m128d cc = _mm_fmadd_pd(logFcent,_mm256_castpd256_pd128(n0v),_mm256_castpd256_pd128(n1));
            __m128d n = _mm_fmadd_pd(logFcent,_mm256_castpd256_pd128(n2),_mm256_castpd256_pd128(n3));
            __m128d x1 = _mm_fmadd_pd(_mm_sub_pd(cc,logPr_),_mm256_castpd256_pd128(n4),n);
            __m128d x2 = _mm_div_pd(_mm_sub_pd(logPr_,cc),x1);
            __m128d x3 = _mm_fmadd_pd(x2,x2,one128);
            __m128d x4 = _mm_div_pd(logFcent,x3);

            __m256d F_ = vec256_powd
                        (
                            n5,
                            _mm256_insertf128_pd (_mm256_setzero_pd (), x4, 0)
                        );
            
            __m128d N = _mm_mul_pd(_mm_mul_pd(K0,_mm256_castpd256_pd128(F_)),invOnePlusPr);

            __m128d Kf = _mm_mul_pd(M,N);
            _mm_storeu_pd(&this->Kf_[j],Kf);

        }
        else if(remainFO==3)  
        {
            const unsigned int i = this->n_TroeFO-3;

            const unsigned int j = this->TroeFO[i+0];

            const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

            const double Kinf0 = this->Kf_[j+0+this->offset_kinf];
            const double Kinf1 = this->Kf_[j+1+this->offset_kinf];
            const double Kinf2 = this->Kf_[j+2+this->offset_kinf];
            __m256d Kinf = _mm256_setr_pd(Kinf0,Kinf1,Kinf2,1);

            double M0 = this->tmp_M[m+0];   
            double M1 = this->tmp_M[m+1]; 
            double M2 = this->tmp_M[m+2]; 
            __m256d M = _mm256_setr_pd(M0,M1,M2,1);

            const double K00 = this->Kf_[j+0];
            const double K01 = this->Kf_[j+1];
            const double K02 = this->Kf_[j+2];
            __m256d K0 = _mm256_setr_pd(K00,K01,K02,1);

            __m256d Pr_ = _mm256_div_pd(_mm256_mul_pd(K0,M),Kinf);
            Pr_ = _mm256_max_pd(small,Pr_);
            __m256d logPr_ = _mm256_mul_pd(vec256_logd(Pr_),invLog10v);

            const double alpha0 = this->alpha_[i+0];
            const double alpha1 = this->alpha_[i+1];
            const double alpha2 = this->alpha_[i+2];
            __m256d alpha = _mm256_setr_pd(alpha0,alpha1,alpha2,1);

            __m256d invOnePlusPr = _mm256_div_pd(one,_mm256_add_pd(one,Pr_));

            const double expTTsss0 = this->tmp_Exp[i+0+this->nSpecies];
            const double expTTsss1 = this->tmp_Exp[i+1+this->nSpecies];
            const double expTTsss2 = this->tmp_Exp[i+2+this->nSpecies];
            __m256d expTTsss = _mm256_setr_pd(expTTsss0,expTTsss1,expTTsss2,1);

            const double expTTss0  = this->tmp_Exp[i+0+this->nSpecies+this->n_Troe];
            const double expTTss1  = this->tmp_Exp[i+1+this->nSpecies+this->n_Troe];
            const double expTTss2  = this->tmp_Exp[i+2+this->nSpecies+this->n_Troe];
            __m256d expTTss = _mm256_setr_pd(expTTss0,expTTss1,expTTss2,1);

            const double expTTs0   = this->tmp_Exp[i+0+this->nSpecies+this->n_Troe*2];
            const double expTTs1   = this->tmp_Exp[i+1+this->nSpecies+this->n_Troe*2];
            const double expTTs2   = this->tmp_Exp[i+2+this->nSpecies+this->n_Troe*2];
            __m256d expTTs = _mm256_setr_pd(expTTs0,expTTs1,expTTs2,1);

            __m256d Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha), expTTsss);
            Fcent = _mm256_fmadd_pd(alpha, expTTs,Fcent);
            Fcent = _mm256_add_pd(expTTss,Fcent);

            __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),invLog10v);
            __m256d cc = _mm256_fmadd_pd(logFcent,n0v,n1);
            __m256d n = _mm256_fmadd_pd(logFcent,n2,n3);
            __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(cc,logPr_),n4,n);
            __m256d x2 = _mm256_div_pd(_mm256_sub_pd(logPr_,cc),x1);
            __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
            __m256d x4 = _mm256_div_pd(logFcent,x3);
            __m256d F_ = vec256_powd(n5,x4);
            __m256d N = _mm256_mul_pd(_mm256_mul_pd(K0,F_),invOnePlusPr);

            __m256d Kf = _mm256_mul_pd(M,N);
            this->Kf_[j+0] = get0(Kf);   
            this->Kf_[j+1] = get1(Kf);
            this->Kf_[j+2] = get2(Kf);   
        }


    unsigned int remainCA = this->n_TroeCA%4;
        for(unsigned int i = 0; i < this->n_TroeCA-remainCA;i=i+4)
        {
            const unsigned int j = this->TroeCA[i+0];
            const unsigned int m0 = j - this->Ikf[5] + Itbr[3];
            __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j-Ikf[5]+Ikf[8]]);
            __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
            __m256d K0 = _mm256_loadu_pd(&this->Kf_[j]);
            __m256d Pr_ = _mm256_div_pd(_mm256_mul_pd(K0,M),Kinf);
            Pr_ = _mm256_max_pd(small,Pr_);
            __m256d logPr_ = _mm256_mul_pd(vec256_logd(Pr_),invLog10v);
            __m256d alpha = _mm256_loadu_pd(&this->alpha_[i+this->n_TroeFO]);
            __m256d invOnePlusPr = _mm256_div_pd(one,_mm256_add_pd(one,Pr_));
            __m256d expTTsss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_TroeFO]);
            __m256d expTTss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe+this->n_TroeFO]);
            __m256d expTTs = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*2+this->n_TroeFO]);

            __m256d Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha), expTTsss);
            Fcent = _mm256_fmadd_pd(alpha, expTTs,Fcent);
            Fcent = _mm256_add_pd(expTTss,Fcent);

            __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),invLog10v);
            __m256d cc = _mm256_fmadd_pd(logFcent,n0v,n1);
            __m256d n = _mm256_fmadd_pd(logFcent,n2,n3);
            __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(cc,logPr_),n4,n);
            __m256d x2 = _mm256_div_pd(_mm256_sub_pd(logPr_,cc),x1);
            __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
            __m256d x4 = _mm256_div_pd(logFcent,x3);
            __m256d F_ = vec256_powd(n5,x4);
            __m256d N = _mm256_mul_pd(_mm256_mul_pd(K0,F_),invOnePlusPr);
            __m256d Kf = (N);
            _mm256_storeu_pd(&this->Kf_[j],Kf);
        }        
        if(remainCA==1)       
        {
            const unsigned int i = this->n_TroeCA-1;
            const unsigned int j = this->TroeCA[i+0];
            const unsigned int m0 = j - this->Ikf[5] + this->Itbr[3];
            const double Kinf0 = this->Kf_[j-Ikf[5]+Ikf[8]];
            double M0 = this->tmp_M[m0];   
            const double K00 = this->Kf_[j];
            const double Pr0 = K00*M0/Kinf0;  
            const double logPr0 = std::log10(std::max(Pr0, FastChemistry::TroeLimiter));
            const double expTTsss0 = this->tmp_Exp[i+this->nSpecies+this->n_TroeFO];
            const double expTTss0  = this->tmp_Exp[i+this->nSpecies+this->n_Troe+this->n_TroeFO];
            const double expTTs0   = this->tmp_Exp[i+this->nSpecies+this->n_Troe*2+this->n_TroeFO];
            const double Fcent0 =(1 - this->alpha_[i+this->n_TroeFO])*expTTsss0 + this->alpha_[i+this->n_TroeFO]*expTTs0 + expTTss0;
            const double logFcent0 = std::log10(std::max(Fcent0, FastChemistry::TroeLimiter));
            const double c0 = -0.4 - 0.67*logFcent0;
            const double n0 = 0.75 - 1.27*logFcent0;
            const double x1_0 = n0 - 0.14*(logPr0 + c0);
            const double x2_0 = (logPr0 + c0)/x1_0;
            const double x3_0 = 1 + (x2_0*x2_0);
            const double x4_0 = logFcent0/x3_0;
            const double F0 = std::pow(10, x4_0);
            const double N0  = 1/(1+Pr0)*F0*K00;
            this->Kf_[j] = N0;   
        }
        else if(remainCA==2)  
        {
            const unsigned int i = this->n_TroeCA-2;
            const unsigned int j = this->TroeCA[i+0];
            const unsigned int m0 = j - this->Ikf[5] + this->Itbr[3];
            __m128d one128 = _mm256_castpd256_pd128(one);
            __m128d Kinf = _mm_loadu_pd(&this->Kf_[j-Ikf[5]+Ikf[8]]);
            __m128d M = _mm_loadu_pd(&this->tmp_M[m0]);
            __m128d K0 = _mm_loadu_pd(&this->Kf_[j]);    
            __m256d Pr_ = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_div_pd(_mm_mul_pd(K0,M),Kinf), 0);
            Pr_ = _mm256_max_pd(small,Pr_);

            __m128d logPr_ = _mm256_castpd256_pd128(_mm256_mul_pd(vec256_logd((Pr_)),(invLog10v)));
            __m128d alpha = _mm_loadu_pd(&this->alpha_[i+this->n_TroeFO]);
            __m128d invOnePlusPr = _mm_div_pd(one128,_mm_add_pd(one128,_mm256_castpd256_pd128(Pr_)));
            __m128d expTTsss = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_TroeFO]);
            __m128d expTTss = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe+this->n_TroeFO]);
            __m128d expTTs = _mm_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*2+this->n_TroeFO]);

            __m128d Fcent  = _mm_mul_pd(_mm_sub_pd(one128,alpha), expTTsss);
            Fcent = _mm_fmadd_pd(alpha, expTTs,Fcent);
            Fcent = _mm_add_pd(expTTss,Fcent);

            __m128d logFcent = _mm_mul_pd
            (
                _mm256_castpd256_pd128(vec256_logd
                (
                    _mm256_max_pd
                    (
                        _mm256_insertf128_pd (_mm256_setzero_pd (), Fcent, 0),
                        small
                    )
                )),
                _mm256_castpd256_pd128(invLog10v)
            );
            __m128d cc = _mm_fmadd_pd(logFcent,_mm256_castpd256_pd128(n0v),_mm256_castpd256_pd128(n1));
            __m128d n = _mm_fmadd_pd(logFcent,_mm256_castpd256_pd128(n2),_mm256_castpd256_pd128(n3));
            __m128d x1 = _mm_fmadd_pd(_mm_sub_pd(cc,logPr_),_mm256_castpd256_pd128(n4),n);
            __m128d x2 = _mm_div_pd(_mm_sub_pd(logPr_,cc),x1);
            __m128d x3 = _mm_fmadd_pd(x2,x2,one128);
            __m128d x4 = _mm_div_pd(logFcent,x3);

            __m256d F_ = vec256_powd
                        (
                            n5,
                            _mm256_insertf128_pd (_mm256_setzero_pd (), x4, 0)
                        );
            
            __m128d N = _mm_mul_pd(_mm_mul_pd(K0,_mm256_castpd256_pd128(F_)),invOnePlusPr);
            _mm_storeu_pd(&this->Kf_[j],N);
        }
        else if(remainCA==3)  
        {
            const unsigned int i = this->n_TroeCA-3;
            const unsigned int j = this->TroeCA[i+0];
            const unsigned int m = j - this->Ikf[5] + Itbr[3];

            const double Kinf0 = this->Kf_[j+0-Ikf[5]+Ikf[8]];
            const double Kinf1 = this->Kf_[j+1-Ikf[5]+Ikf[8]];
            const double Kinf2 = this->Kf_[j+2-Ikf[5]+Ikf[8]];
            __m256d Kinf = _mm256_setr_pd(Kinf0,Kinf1,Kinf2,1);

            double M0 = this->tmp_M[m+0];   
            double M1 = this->tmp_M[m+1]; 
            double M2 = this->tmp_M[m+2]; 
            __m256d M = _mm256_setr_pd(M0,M1,M2,1);

            const double K00 = this->Kf_[j+0];
            const double K01 = this->Kf_[j+1];
            const double K02 = this->Kf_[j+2];
            __m256d K0 = _mm256_setr_pd(K00,K01,K02,1);


            __m256d Pr_ = _mm256_div_pd(_mm256_mul_pd(K0,M),Kinf);
            Pr_ = _mm256_max_pd(small,Pr_);
            __m256d logPr_ = _mm256_mul_pd(vec256_logd(Pr_),invLog10v);


            double alpha0 = this->alpha_[i+0+this->n_TroeFO];
            double alpha1 = this->alpha_[i+1+this->n_TroeFO];
            double alpha2 = this->alpha_[i+2+this->n_TroeFO];
            __m256d alpha = _mm256_setr_pd(alpha0,alpha1,alpha2,1);

            __m256d invOnePlusPr = _mm256_div_pd(one,_mm256_add_pd(one,Pr_));

            const double expTTsss0 = this->tmp_Exp[i+0+this->nSpecies+this->n_TroeFO];
            const double expTTsss1 = this->tmp_Exp[i+1+this->nSpecies+this->n_TroeFO];
            const double expTTsss2 = this->tmp_Exp[i+2+this->nSpecies+this->n_TroeFO];
            __m256d expTTsss = _mm256_setr_pd(expTTsss0,expTTsss1,expTTsss2,1);

            const double expTTss0  = this->tmp_Exp[i+0+this->nSpecies+this->n_Troe+this->n_TroeFO];
            const double expTTss1  = this->tmp_Exp[i+1+this->nSpecies+this->n_Troe+this->n_TroeFO];
            const double expTTss2  = this->tmp_Exp[i+2+this->nSpecies+this->n_Troe+this->n_TroeFO];
            __m256d expTTss = _mm256_setr_pd(expTTss0,expTTss1,expTTss2,1);


            const double expTTs0 = this->tmp_Exp[i+0+this->nSpecies+this->n_Troe*2+this->n_TroeFO];
            const double expTTs1 = this->tmp_Exp[i+1+this->nSpecies+this->n_Troe*2+this->n_TroeFO];
            const double expTTs2 = this->tmp_Exp[i+2+this->nSpecies+this->n_Troe*2+this->n_TroeFO];
            __m256d expTTs = _mm256_setr_pd(expTTs0,expTTs1,expTTs2,1);

            __m256d Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha), expTTsss);
            Fcent = _mm256_fmadd_pd(alpha, expTTs,Fcent);
            Fcent = _mm256_add_pd(expTTss,Fcent);

            __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),invLog10v);
            __m256d cc = _mm256_fmadd_pd(logFcent,n0v,n1);
            __m256d n = _mm256_fmadd_pd(logFcent,n2,n3);
            __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(cc,logPr_),n4,n);
            __m256d x2 = _mm256_div_pd(_mm256_sub_pd(logPr_,cc),x1);
            __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
            __m256d x4 = _mm256_div_pd(logFcent,x3);
            __m256d F_ = vec256_powd(n5,x4);
            __m256d N = _mm256_mul_pd(_mm256_mul_pd(K0,F_),invOnePlusPr);
            __m256d Kf = (N);
            this->Kf_[j+0] = get0(Kf);
            this->Kf_[j+1] = get1(Kf);
            this->Kf_[j+2] = get2(Kf);
        }

}