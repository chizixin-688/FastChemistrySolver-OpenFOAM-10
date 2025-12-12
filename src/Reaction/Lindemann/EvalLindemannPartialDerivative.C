
void OptReaction::evalLindemannPartialDerivative()const noexcept
{
    size_t remain_Lindemann = (Lindemann.size())%4;
    __m256d one = _mm256_set1_pd(1.0);    
    __m256d condition = _mm256_set1_pd(this->n_Fall_Off_Reaction);
    __m256d inc = _mm256_setr_pd(0,1,2,3);
    for (size_t i = 0;i<Lindemann.size()-remain_Lindemann;i=i+4)
    {
        const unsigned int j0 = this->Lindemann[i+0]+0;
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
        __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0 + this->offset_kinf]);    
        __m256d invKinf = _mm256_div_pd(one,Kinf);
        __m256d dKinfdT = _mm256_loadu_pd(&this->dKfdT_[j0 + this->offset_kinf]);
        __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);        
        __m256d dK0dT = _mm256_loadu_pd(&this->dKfdT_[j0]);     
        __m256d M = _mm256_loadu_pd(&tmp_M[m0]);


        __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(K0,M),invKinf);
        __m256d dPrdT = _mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT));
        dPrdT = _mm256_mul_pd(dPrdT,invKinf);

        __m256d k = _mm256_set1_pd(j0-this->Ikf[4]);
        k = _mm256_add_pd(k,inc);
        //__m256d k = _mm256_setr_pd(k0,k1,k2,k3);

        __m256d cmp = _mm256_cmp_pd(k,condition,_CMP_LT_OQ);
        __m256d dKdT = _mm256_blendv_pd(dK0dT,_mm256_mul_pd(Pr,dKinfdT),cmp);
        __m256d K = _mm256_blendv_pd(K0,Kinf,cmp);
        __m256d KK = _mm256_blendv_pd(_mm256_mul_pd(K0,invKinf),one,cmp);
        __m256d tmp = _mm256_div_pd(one,_mm256_add_pd(one,Pr));
        __m256d N1 = _mm256_blendv_pd(-tmp,tmp,cmp);
        __m256d N = _mm256_mul_pd(tmp,K0);
        __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(tmp,tmp),dPrdT),K);
        dKfdT = _mm256_fmadd_pd(tmp,dKdT,dKfdT);
        _mm256_storeu_pd(&this->dKfdT_[j0],dKfdT);
        __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,tmp),KK),N1);
        _mm256_storeu_pd(&this->dKfdC_[m0],dKfdC);
        __m256d KF = _mm256_blendv_pd(N,_mm256_mul_pd(M,N),cmp);
        _mm256_storeu_pd(&this->Kf_[j0],KF);
    }
    if(remain_Lindemann==1)
    {
        size_t i = this->Lindemann.size()-1;
        const unsigned int j0 = this->Lindemann[i+0];
        const unsigned int k0 = j0 - this->Ikf[4];
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
        const double Kinf0 = this->Kf_[j0 + this->offset_kinf];
        const double dKinfdT0 = this->dKfdT_[j0 + this->offset_kinf];
        const double K00 = this->Kf_[j0];
        double M0 = tmp_M[m0];
        const double invKinf0 = 1.0/Kinf0;
        const double Pr0 = K00*M0*invKinf0; 
        const double dK0dT0 =  this->dKfdT_[j0];
        const double dPrdT0 = (M0*dK0dT0-Pr0*dKinfdT0)*invKinf0;
        const double invOnePlusPr0 = 1/(1+Pr0);
        const double dKdT0   = j0 - this->Ikf[4]<this->n_Fall_Off_Reaction?Pr0*dKinfdT0:dK0dT0;
        const double K0      = j0 - this->Ikf[4]<this->n_Fall_Off_Reaction?Kinf0      :K00;
        const double KK0     = j0 - this->Ikf[4]<this->n_Fall_Off_Reaction?1         :K00*invKinf0;
        const double N10     = j0 - this->Ikf[4]<this->n_Fall_Off_Reaction?invOnePlusPr0  :-invOnePlusPr0;
        const double N0  = invOnePlusPr0*1*K00;
        this->dKfdT_[j0] = invOnePlusPr0*dKdT0 + invOnePlusPr0*invOnePlusPr0*dPrdT0*K0;
        this->dKfdC_[m0] =  K00*KK0*(N10)*invOnePlusPr0; 
        this->Kf_[j0] = k0<this->n_Fall_Off_Reaction ? M0*N0 : N0;    
    }
    else if (remain_Lindemann==2)
    {
        size_t i = this->Lindemann.size()-2;
        const unsigned int j0 = this->Lindemann[i+0]+0;
        const unsigned int j1 = this->Lindemann[i+0]+1;
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
        const unsigned int m1 = j1 - this->Ikf[4] + this->Itbr[2];            
        const unsigned int k0 = j0 - this->Ikf[4];
        const unsigned int k1 = j1 - this->Ikf[4];
        __m256d Kinf = _mm256_setr_pd(Kf_[j0 + this->offset_kinf],Kf_[j1 + this->offset_kinf],1,1);    
        __m256d invKinf = _mm256_div_pd(_mm256_set1_pd(1.0),Kinf);
        __m256d dKinfdT = _mm256_setr_pd(dKfdT_[j0 + this->offset_kinf],dKfdT_[j1 + this->offset_kinf],1,1);
        __m256d K0 = _mm256_setr_pd(Kf_[j0],Kf_[j1],1,1);        
        __m256d dK0dT = _mm256_setr_pd(dKfdT_[j0],dKfdT_[j1],1,1);     
        __m256d M = _mm256_setr_pd(tmp_M[m0],tmp_M[m1],1,1);
        __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(K0,M),invKinf);
        __m256d dPrdT = _mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT));
        dPrdT = _mm256_mul_pd(dPrdT,invKinf);
        __m256d k = _mm256_setr_pd(k0,k1,1,1);
        __m256d cmp = _mm256_cmp_pd(k,_mm256_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
        __m256d dKdT = _mm256_blendv_pd(dK0dT,_mm256_mul_pd(Pr,dKinfdT),cmp);
        __m256d K = _mm256_blendv_pd(K0,Kinf,cmp);
        __m256d KK = _mm256_blendv_pd(_mm256_mul_pd(K0,invKinf),one,cmp);
        __m256d tmp = _mm256_div_pd(one,_mm256_add_pd(one,Pr));
        __m256d N1 = _mm256_blendv_pd(-tmp,tmp,cmp);
        __m256d N = _mm256_mul_pd(tmp,K0);
        __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(tmp,tmp),dPrdT),K);
        dKfdT = _mm256_fmadd_pd(tmp,dKdT,dKfdT);
        dKfdT_[j0] = get_elem0(dKfdT);
        dKfdT_[j1] = get_elem1(dKfdT);
        __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,tmp),KK),N1);
        dKfdC_[m0] = get_elem0(dKfdC);
        dKfdC_[m1] = get_elem1(dKfdC);
        __m256d KF = _mm256_blendv_pd(N,_mm256_mul_pd(M,N),cmp);   
        Kf_[j0] = get_elem0(KF);
        Kf_[j1] = get_elem1(KF);
    }
    else if (remain_Lindemann==3)
    {
        size_t i = this->Lindemann.size()-3;
        const unsigned int j0 = this->Lindemann[i+0]+0;
        const unsigned int j1 = this->Lindemann[i+0]+1;
        const unsigned int j2 = this->Lindemann[i+0]+2;
        const unsigned int k0 = j0 - this->Ikf[4];
        const unsigned int k1 = j1 - this->Ikf[4];
        const unsigned int k2 = j2 - this->Ikf[4];
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
        const unsigned int m1 = j1 - this->Ikf[4] + this->Itbr[2];
        const unsigned int m2 = j2 - this->Ikf[4] + this->Itbr[2];
        __m256d Kinf = _mm256_setr_pd(Kf_[j0+this->offset_kinf],Kf_[j1+this->offset_kinf],Kf_[j2+this->offset_kinf],1);    
     
        __m256d invKinf = _mm256_div_pd(_mm256_set1_pd(1.0),Kinf);
        __m256d dKinfdT = _mm256_setr_pd(dKfdT_[j0+this->offset_kinf],dKfdT_[j1+this->offset_kinf],dKfdT_[j2+this->offset_kinf],1);
        __m256d K0 = _mm256_setr_pd(Kf_[j0],Kf_[j1],Kf_[j2],1);        
        __m256d dK0dT = _mm256_setr_pd(dKfdT_[j0],dKfdT_[j1],dKfdT_[j2],1);     
        __m256d M = _mm256_setr_pd(tmp_M[m0],tmp_M[m1],tmp_M[m2],1);
        __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(K0,M),invKinf);
        __m256d dPrdT = _mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT));
        dPrdT = _mm256_mul_pd(dPrdT,invKinf);
        __m256d k = _mm256_setr_pd(k0,k1,k2,1);
        __m256d cmp = _mm256_cmp_pd(k,_mm256_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
        __m256d dKdT = _mm256_blendv_pd(dK0dT,_mm256_mul_pd(Pr,dKinfdT),cmp);
        __m256d K = _mm256_blendv_pd(K0,Kinf,cmp);
        __m256d KK = _mm256_blendv_pd(_mm256_mul_pd(K0,invKinf),one,cmp);
        __m256d tmp = _mm256_div_pd(one,_mm256_add_pd(one,Pr));
        __m256d N1 = _mm256_blendv_pd(-tmp,tmp,cmp);
        __m256d N = _mm256_mul_pd(tmp,K0);
        __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(tmp,tmp),dPrdT),K);
        dKfdT = _mm256_fmadd_pd(tmp,dKdT,dKfdT);
        this->dKfdT_[j0] = get_elem0(dKfdT);
        this->dKfdT_[j1] = get_elem1(dKfdT);
        this->dKfdT_[j2] = get_elem2(dKfdT);
        __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,tmp),KK),N1);
        this->dKfdC_[m0] = get_elem0(dKfdC);
        this->dKfdC_[m1] = get_elem1(dKfdC);
        this->dKfdC_[m2] = get_elem2(dKfdC);
        __m256d KF = _mm256_blendv_pd(N,_mm256_mul_pd(M,N),cmp);
        this->Kf_[j0] = get_elem0(KF);
        this->Kf_[j1] = get_elem1(KF);
        this->Kf_[j2] = get_elem2(KF);
    }   
}