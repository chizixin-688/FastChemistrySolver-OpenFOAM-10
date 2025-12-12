void OptReaction::evalTroeRateConstant()const noexcept
{
    unsigned int remain = this->Troe.size()%4;


    __m256d one = _mm256_set1_pd(1.0);
    __m256d invLog10v = _mm256_set1_pd(0.43429448190325182765112891891661);
    __m256d n0 = _mm256_set1_pd(0.67);
    __m256d n1 = _mm256_set1_pd(0.4);
    __m256d n2 = _mm256_set1_pd(-1.27);
    __m256d n3 = _mm256_set1_pd(0.75);
    __m256d n4 = _mm256_set1_pd(0.14);
    __m256d n5 = _mm256_set1_pd(10);
    __m256d cond = _mm256_set1_pd(double(this->n_Fall_Off_Reaction));
    __m256d small = _mm256_set1_pd(TroeLimiter);
        for(unsigned int i = 0; i < this->Troe.size()-remain;i=i+4)
        {
            const unsigned int j0 = this->Troe[i+0];
            const unsigned int m0 = j0 - this->Ikf[4] + Itbr[2];
            __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
            __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
            __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);    
            __m256d Pr_ = _mm256_div_pd(_mm256_mul_pd(K0,M),Kinf);
            Pr_ = _mm256_max_pd(small,Pr_);
            __m256d logPr_ = _mm256_mul_pd(vec256_logd(Pr_),invLog10v);
            __m256d alpha = _mm256_loadu_pd(&this->alpha_[i]);
            __m256d invOnePlusPr = _mm256_div_pd(one,_mm256_add_pd(one,Pr_));
            __m256d expTTsss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies]);
            __m256d expTTss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->Troe.size()]);
            __m256d expTTs = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->Troe.size()*2]);

            __m256d Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha), expTTsss);
            Fcent = _mm256_fmadd_pd(alpha, expTTs,Fcent);
            Fcent = _mm256_add_pd(expTTss,Fcent);

            __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),invLog10v);
            __m256d cc = _mm256_fmadd_pd(logFcent,n0,n1);
            __m256d n = _mm256_fmadd_pd(logFcent,n2,n3);
            __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(cc,logPr_),n4,n);
            __m256d x2 = _mm256_div_pd(_mm256_sub_pd(logPr_,cc),x1);
            __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
            __m256d x4 = _mm256_div_pd(logFcent,x3);
            __m256d F_ = vec256_powd(n5,x4);
            __m256d N = _mm256_mul_pd(_mm256_mul_pd(K0,F_),invOnePlusPr);

            __m256d j0v = _mm256_set1_pd(double(j0- this->Ikf[4]));
            //__m256d k = _mm256_setr_pd(k0,k1,k2,k3);
            __m256d inc = _mm256_setr_pd(0.0,1.0,2.0,3.0);
            __m256d k = _mm256_add_pd(j0v,inc);
            __m256d cmp = _mm256_cmp_pd(k,cond,_CMP_LT_OQ);
            __m256d Kf = _mm256_blendv_pd(N,_mm256_mul_pd(M,N),cmp);
            _mm256_storeu_pd(&this->Kf_[j0],Kf);
        }        
        if(remain==1)       
        {
            const size_t i = this->Troe.size()-1;
            const unsigned int j0 = this->Troe[i+0];
            const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
            const unsigned int k0 = j0 - this->Ikf[4];
            const double Kinf0 = this->Kf_[j0+this->offset_kinf];
            double M0 = this->tmp_M[m0];   
            const double K00 = this->Kf_[j0];
            const double Pr0 = K00*M0/Kinf0;  
            const double logPr0 = std::log10(std::max(Pr0, TroeLimiter));
            const double expTTsss0 = this->tmp_Exp[i + 0 + this->nSpecies];
            const double expTTss0  = this->tmp_Exp[i + 0 + this->nSpecies + this->Troe.size()];
            const double expTTs0   = this->tmp_Exp[i + 0 + this->nSpecies + this->Troe.size()*2];
            const double Fcent0 =(1 - this->alpha_[i + 0])*expTTsss0 + this->alpha_[i + 0]*expTTs0 + expTTss0;
            const double logFcent0 = std::log10(std::max(Fcent0, TroeLimiter));
            const double c0 = -0.4 - 0.67*logFcent0;
            const double n0 = 0.75 - 1.27*logFcent0;
            const double x1_0 = n0 - 0.14*(logPr0 + c0);
            const double x2_0 = (logPr0 + c0)/x1_0;
            const double x3_0 = 1 + (x2_0*x2_0);
            const double x4_0 = logFcent0/x3_0;
            const double F0 = std::pow(10, x4_0);
            const double N0  = 1/(1+Pr0)*F0*K00;
            this->Kf_[j0] = k0<this->n_Fall_Off_Reaction ? M0*N0 : N0;   
        }
        else if(remain==2)  
        {
            const size_t i = this->Troe.size()-2;
            
            const unsigned int j0 = this->Troe[i+0];
            const unsigned int j1 = this->Troe[i+1];

            const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
            const unsigned int m1 = j1 - this->Ikf[4] + this->Itbr[2];

            const unsigned int k0 = j0 - this->Ikf[4];
            const unsigned int k1 = j1 - this->Ikf[4];

            const double Kinf0 = this->Kf_[j0+this->offset_kinf];
            const double Kinf1 = this->Kf_[j1+this->offset_kinf];

            double M0 = this->tmp_M[m0];   
            double M1 = this->tmp_M[m1]; 

            const double K00 = this->Kf_[j0];
            const double K01 = this->Kf_[j1];

            const double Pr0 = K00*M0/Kinf0;  
            const double Pr1 = K01*M1/Kinf1;  

            const double logPr0 = std::log10(std::max(Pr0, TroeLimiter));
            const double logPr1 = std::log10(std::max(Pr1, TroeLimiter));

            const double expTTsss0 = this->tmp_Exp[i + 0 + this->nSpecies];
            const double expTTsss1 = this->tmp_Exp[i + 1 + this->nSpecies];

            const double expTTss0  = this->tmp_Exp[i + 0 + this->nSpecies + this->Troe.size()];
            const double expTTss1  = this->tmp_Exp[i + 1 + this->nSpecies + this->Troe.size()];

            const double expTTs0   = this->tmp_Exp[i + 0 + this->nSpecies + this->Troe.size()*2];
            const double expTTs1   = this->tmp_Exp[i + 1 + this->nSpecies + this->Troe.size()*2];

            const double Fcent0 =(1 - this->alpha_[i + 0])*expTTsss0 + this->alpha_[i + 0]*expTTs0 + expTTss0;
            const double Fcent1 =(1 - this->alpha_[i + 1])*expTTsss1 + this->alpha_[i + 1]*expTTs1 + expTTss1;

            const double logFcent0 = std::log10(std::max(Fcent0, TroeLimiter));
            const double logFcent1 = std::log10(std::max(Fcent1, TroeLimiter));

            const double c0 = -0.4 - 0.67*logFcent0;
            const double c1 = -0.4 - 0.67*logFcent1;

            const double n0 = 0.75 - 1.27*logFcent0;
            const double n1 = 0.75 - 1.27*logFcent1;

            const double x1_0 = n0 - 0.14*(logPr0 + c0);
            const double x1_1 = n1 - 0.14*(logPr1 + c1);

            const double x2_0 = (logPr0 + c0)/x1_0;
            const double x2_1 = (logPr1 + c1)/x1_1;

            const double x3_0 = 1 + (x2_0*x2_0);
            const double x3_1 = 1 + (x2_1*x2_1);

            const double x4_0 = logFcent0/x3_0;
            const double x4_1 = logFcent1/x3_1;

            const double F0 = std::pow(10, x4_0);
            const double F1 = std::pow(10, x4_1);

            const double N0  = 1/(1+Pr0)*F0*K00;
            const double N1  = 1/(1+Pr1)*F1*K01;

            this->Kf_[j0] = k0<this->n_Fall_Off_Reaction ? M0*N0 : N0;
            this->Kf_[j1] = k1<this->n_Fall_Off_Reaction ? M1*N1 : N1;
        }
        else if(remain==3)  {this->Troe_F_3();}
}