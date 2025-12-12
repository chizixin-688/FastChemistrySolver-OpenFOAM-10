void OptReaction::evalLindemannRateConstant()const noexcept
{
    __m256d one = _mm256_set1_pd(1.0);
    __m256d condition = _mm256_set1_pd(double(this->n_Fall_Off_Reaction));
    __m256d inc = _mm256_setr_pd(0.0,1.0,2.0,3.0);
    unsigned int remain = this->Lindemann.size()%4;
    for (unsigned int i = 0;i<this->Lindemann.size()-remain;i=i+4)
    {
        const unsigned int j0 = this->Lindemann[i+0]+0;
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
        __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
        __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
        __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);
        __m256d Pr = _mm256_div_pd(_mm256_mul_pd(K0,M),Kinf);
        __m256d N = _mm256_div_pd(K0,_mm256_add_pd(Pr,one));
        __m256d k0v = _mm256_set1_pd(j0-this->Ikf[4]);
        __m256d k = _mm256_add_pd(k0v,inc);
        __m256d cmp = _mm256_cmp_pd(k,condition,_CMP_LT_OQ);
        __m256d Kf = _mm256_blendv_pd(N,_mm256_mul_pd(M,N),cmp);
        _mm256_storeu_pd(&this->Kf_[j0],Kf);
    }
    if(remain==1)
    {
        size_t i = this->Lindemann.size()-1;
        const unsigned int j = this->Lindemann[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];
        const unsigned int k = j - this->Ikf[4];
        const double Kinf = this->Kf_[j+this->offset_kinf];
        double M = this->tmp_M[m];     
        const double K0 = this->Kf_[j];
        const double Pr = K0*M/Kinf;   
        const double N          = 1/(1+Pr)*K0;
        this->Kf_[j] = k<this->n_Fall_Off_Reaction ? M*N : N;            
    }
    else if(remain==2)
    {
        size_t i = this->Lindemann.size()-2;
        const unsigned int j0 = this->Lindemann[i+0]+0;
        const unsigned int j1 = this->Lindemann[i+0]+1;
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
        const unsigned int k0 = j0 - this->Ikf[4];
        const unsigned int k1 = j1 - this->Ikf[4];
        __m128d Kinf = _mm_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
        __m128d M = _mm_loadu_pd(&this->tmp_M[m0]);
        __m128d K0 = _mm_loadu_pd(&this->Kf_[j0]);
        __m128d Pr = _mm_div_pd(_mm_mul_pd(K0,M),Kinf);
        __m128d one128 = _mm256_castpd256_pd128(one);
        __m128d N = _mm_div_pd(K0,_mm_add_pd(Pr,one128));
        __m128d k = _mm_setr_pd(k0,k1);
        __m128d cmp = _mm_cmp_pd(k,_mm_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
        __m128d Kf = _mm_blendv_pd(N,_mm_mul_pd(M,N),cmp);
        _mm_storeu_pd(&this->Kf_[j0],Kf);
    }
    else if(remain==3)
    {
        size_t i = this->Lindemann.size()-3;
        const unsigned int j0 = this->Lindemann[i+0];
        const unsigned int j1 = this->Lindemann[i+1];
        const unsigned int j2 = this->Lindemann[i+2];
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
        const unsigned int m1 = j1 - this->Ikf[4] + this->Itbr[2];
        const unsigned int m2 = j2 - this->Ikf[4] + this->Itbr[2];
        const unsigned int k0 = j0 - this->Ikf[4];
        const unsigned int k1 = j1 - this->Ikf[4];
        const unsigned int k2 = j2 - this->Ikf[4];
        const double Kinf0 = this->Kf_[j0+this->offset_kinf];
        const double Kinf1 = this->Kf_[j1+this->offset_kinf];
        const double Kinf2 = this->Kf_[j2+this->offset_kinf];
        __m256d Kinf = _mm256_setr_pd(Kinf0,Kinf1,Kinf2,1);

        double M0 = this->tmp_M[m0];
        double M1 = this->tmp_M[m1];
        double M2 = this->tmp_M[m2];
        __m256d M = _mm256_setr_pd(M0,M1,M2,0);

        const double K00 = this->Kf_[j0];
        const double K01 = this->Kf_[j1];
        const double K02 = this->Kf_[j2];
        __m256d K0 = _mm256_setr_pd(K00,K01,K02,0);    
        
        __m256d Pr = _mm256_div_pd(_mm256_mul_pd(K0,M),Kinf);
        __m256d N = _mm256_div_pd(K0,_mm256_add_pd(Pr,one));
        __m256d k = _mm256_setr_pd(k0,k1,k2,0);
        __m256d cmp = _mm256_cmp_pd(k,_mm256_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
        __m256d Kf = _mm256_blendv_pd(N,_mm256_mul_pd(M,N),cmp);
        this->Kf_[j0] = get_elem0(Kf);
        this->Kf_[j1] = get_elem1(Kf);
        this->Kf_[j2] = get_elem2(Kf);
    }
}