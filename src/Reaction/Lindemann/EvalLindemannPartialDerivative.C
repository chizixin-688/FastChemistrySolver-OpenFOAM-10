#include "OptReaction.H"
void FastChemistry::OptReaction::evalLindemannPartialDerivative()const noexcept
{
    const unsigned int remainFO = (this->n_LindemannFO)%4;
    __m256d one = _mm256_set1_pd(1.0);    

    for (unsigned int i = 0;i<this->n_LindemannFO-remainFO;i=i+4)
    {
        const unsigned int j0 = this->LindemannFO[i+0]+0;
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


        __m256d dKdT = (_mm256_mul_pd(Pr,dKinfdT));
        __m256d K = (Kinf);
        __m256d KK = (one);
        __m256d tmp = _mm256_div_pd(one,_mm256_add_pd(one,Pr));
        __m256d N1 = (tmp);
        __m256d N = _mm256_mul_pd(tmp,K0);
        __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(tmp,tmp),dPrdT),K);
        dKfdT = _mm256_fmadd_pd(tmp,dKdT,dKfdT);
        _mm256_storeu_pd(&this->dKfdT_[j0],dKfdT);
        __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,tmp),KK),N1);
        _mm256_storeu_pd(&this->dKfdC_[m0],dKfdC);
        __m256d KF = _mm256_mul_pd(M,N);
        _mm256_storeu_pd(&this->Kf_[j0],KF);
    }
    if(remainFO==1)
    {
        const unsigned int i = this->n_LindemannFO-1;
        const unsigned int j0 = this->LindemannFO[i+0];
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
        const double dKdT0   = Pr0*dKinfdT0;
        const double K0      = Kinf0;
        const double KK0     = 1;
        const double N10     = invOnePlusPr0;
        const double N0  = invOnePlusPr0*1*K00;
        this->dKfdT_[j0] = invOnePlusPr0*dKdT0 + invOnePlusPr0*invOnePlusPr0*dPrdT0*K0;
        this->dKfdC_[m0] =  K00*KK0*(N10)*invOnePlusPr0; 
        this->Kf_[j0] = M0*N0;    
    }
    else if (remainFO==2)
    {
        const unsigned int i = this->n_LindemannFO-2;
        const unsigned int j0 = this->LindemannFO[i+0]+0;
        const unsigned int j1 = this->LindemannFO[i+0]+1;
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
        const unsigned int m1 = j1 - this->Ikf[4] + this->Itbr[2];            
        //const unsigned int k0 = j0 - this->Ikf[4];
        //const unsigned int k1 = j1 - this->Ikf[4];
        __m256d Kinf = _mm256_setr_pd(Kf_[j0 + this->offset_kinf],Kf_[j1 + this->offset_kinf],1,1);    
        __m256d invKinf = _mm256_div_pd(_mm256_set1_pd(1.0),Kinf);
        __m256d dKinfdT = _mm256_setr_pd(dKfdT_[j0 + this->offset_kinf],dKfdT_[j1 + this->offset_kinf],1,1);
        __m256d K0 = _mm256_setr_pd(Kf_[j0],Kf_[j1],1,1);        
        __m256d dK0dT = _mm256_setr_pd(dKfdT_[j0],dKfdT_[j1],1,1);     
        __m256d M = _mm256_setr_pd(tmp_M[m0],tmp_M[m1],1,1);
        __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(K0,M),invKinf);
        __m256d dPrdT = _mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT));
        dPrdT = _mm256_mul_pd(dPrdT,invKinf);


        __m256d dKdT = _mm256_mul_pd(Pr,dKinfdT);
        __m256d K = (Kinf);
        __m256d KK = (one);
        __m256d tmp = _mm256_div_pd(one,_mm256_add_pd(one,Pr));
        __m256d N1 = (tmp);
        __m256d N = _mm256_mul_pd(tmp,K0);
        __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(tmp,tmp),dPrdT),K);
        dKfdT = _mm256_fmadd_pd(tmp,dKdT,dKfdT);
        dKfdT_[j0] = get_elem0(dKfdT);
        dKfdT_[j1] = get_elem1(dKfdT);
        __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,tmp),KK),N1);
        dKfdC_[m0] = get_elem0(dKfdC);
        dKfdC_[m1] = get_elem1(dKfdC);
        __m256d KF = (_mm256_mul_pd(M,N));   
        Kf_[j0] = get_elem0(KF);
        Kf_[j1] = get_elem1(KF);
    }
    else if (remainFO==3)
    {
        const unsigned int i = this->n_LindemannFO-3;
        const unsigned int j0 = this->LindemannFO[i+0]+0;
        const unsigned int j1 = this->LindemannFO[i+0]+1;
        const unsigned int j2 = this->LindemannFO[i+0]+2;
        //const unsigned int k0 = j0 - this->Ikf[4];
        //const unsigned int k1 = j1 - this->Ikf[4];
        //const unsigned int k2 = j2 - this->Ikf[4];
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


        __m256d dKdT = (_mm256_mul_pd(Pr,dKinfdT));
        __m256d K = (Kinf);
        __m256d KK = (one);
        __m256d tmp = _mm256_div_pd(one,_mm256_add_pd(one,Pr));
        __m256d N1 = (tmp);
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
        __m256d KF = (_mm256_mul_pd(M,N));
        this->Kf_[j0] = get_elem0(KF);
        this->Kf_[j1] = get_elem1(KF);
        this->Kf_[j2] = get_elem2(KF);
    }

    const unsigned int remainCA = (this->n_LindemannCA)%4;
    for (unsigned int i = 0;i<this->n_LindemannCA-remainCA;i=i+4)
    {
        const unsigned int j0 = this->LindemannCA[i+0]+0;
        const unsigned int m0 = j0 - this->Ikf[5] + this->Itbr[3];
        __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0 + this->offset_kinf]);    
        __m256d invKinf = _mm256_div_pd(one,Kinf);
        __m256d dKinfdT = _mm256_loadu_pd(&this->dKfdT_[j0 + this->offset_kinf]);
        __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);        
        __m256d dK0dT = _mm256_loadu_pd(&this->dKfdT_[j0]);
        __m256d M = _mm256_loadu_pd(&tmp_M[m0]);


        __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(K0,M),invKinf);
        __m256d dPrdT = _mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT));
        dPrdT = _mm256_mul_pd(dPrdT,invKinf);




        __m256d dKdT = (dK0dT);
        __m256d K = (K0);
        __m256d KK = (_mm256_mul_pd(K0,invKinf));
        __m256d tmp = _mm256_div_pd(one,_mm256_add_pd(one,Pr));
        __m256d N1 = (-tmp);
        __m256d N = _mm256_mul_pd(tmp,K0);
        __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(tmp,tmp),dPrdT),K);
        dKfdT = _mm256_fmadd_pd(tmp,dKdT,dKfdT);
        _mm256_storeu_pd(&this->dKfdT_[j0],dKfdT);
        __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,tmp),KK),N1);
        _mm256_storeu_pd(&this->dKfdC_[m0],dKfdC);
        __m256d KF = (N);
        _mm256_storeu_pd(&this->Kf_[j0],KF);
    }
    if(remainCA==1)
    {
        const unsigned int i = this->n_LindemannCA-1;
        const unsigned int j0 = this->LindemannCA[i+0];
        const unsigned int m0 = j0 - this->Ikf[5] + this->Itbr[3];
        const double Kinf0 = this->Kf_[j0- Ikf[5] + Ikf[8]];
        const double dKinfdT0 = this->dKfdT_[j0- Ikf[5] + Ikf[8]];
        const double K00 = this->Kf_[j0];
        double M0 = tmp_M[m0];
        const double invKinf0 = 1.0/Kinf0;
        const double Pr0 = K00*M0*invKinf0; 
        const double dK0dT0 =  this->dKfdT_[j0];
        const double dPrdT0 = (M0*dK0dT0-Pr0*dKinfdT0)*invKinf0;
        const double invOnePlusPr0 = 1/(1+Pr0);
        const double dKdT0   = dK0dT0;
        const double K0      = K00;
        const double KK0     = K00*invKinf0;
        const double N10     = -invOnePlusPr0;
        const double N0  = invOnePlusPr0*1*K00;
        // - k0*k0*dMdc/kInf/sqr(1 + Pr)

        const double r0 = invOnePlusPr0*dKdT0 - invOnePlusPr0*invOnePlusPr0*dPrdT0*K0;
        const double r1 = K00*KK0*(N10)*invOnePlusPr0; 
        this->dKfdT_[j0] = r0;
        this->dKfdC_[m0] =  r1; 
        this->Kf_[j0] = N0;
    }
    else if (remainCA==2)
    {
        const unsigned int i = this->n_LindemannCA-2;
        const unsigned int j0 = this->LindemannCA[i+0]+0;
        const unsigned int j1 = this->LindemannCA[i+0]+1;
        const unsigned int m0 = j0 - this->Ikf[5] + this->Itbr[3];
        const unsigned int m1 = j1 - this->Ikf[5] + this->Itbr[3];            


        __m256d Kinf = _mm256_setr_pd(Kf_[j0 + this->offset_kinf],Kf_[j1 + this->offset_kinf],1,1);    
        __m256d invKinf = _mm256_div_pd(_mm256_set1_pd(1.0),Kinf);
        __m256d dKinfdT = _mm256_setr_pd(dKfdT_[j0 + this->offset_kinf],dKfdT_[j1 + this->offset_kinf],1,1);
        __m256d K0 = _mm256_setr_pd(Kf_[j0],Kf_[j1],1,1);        
        __m256d dK0dT = _mm256_setr_pd(dKfdT_[j0],dKfdT_[j1],1,1);     
        __m256d M = _mm256_setr_pd(tmp_M[m0],tmp_M[m1],1,1);
        __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(K0,M),invKinf);
        __m256d dPrdT = _mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT));
        dPrdT = _mm256_mul_pd(dPrdT,invKinf);


        __m256d dKdT = (dK0dT);
        __m256d K = (K0);
        __m256d KK = (_mm256_mul_pd(K0,invKinf));
        __m256d tmp = _mm256_div_pd(one,_mm256_add_pd(one,Pr));
        __m256d N1 = (-tmp);
        __m256d N = _mm256_mul_pd(tmp,K0);
        __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(tmp,tmp),dPrdT),K);
        dKfdT = _mm256_fmadd_pd(tmp,dKdT,dKfdT);
        dKfdT_[j0] = get_elem0(dKfdT);
        dKfdT_[j1] = get_elem1(dKfdT);
        __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,tmp),KK),N1);
        dKfdC_[m0] = get_elem0(dKfdC);
        dKfdC_[m1] = get_elem1(dKfdC);
        __m256d KF = (N);   
        Kf_[j0] = get_elem0(KF);
        Kf_[j1] = get_elem1(KF);
    }
    else if (remainCA==3)
    {
        const unsigned int i = this->n_LindemannCA-3;
        const unsigned int j0 = this->LindemannCA[i+0]+0;
        const unsigned int j1 = this->LindemannCA[i+0]+1;
        const unsigned int j2 = this->LindemannCA[i+0]+2;

        const unsigned int m0 = j0 - this->Ikf[5] + this->Itbr[3];
        const unsigned int m1 = j1 - this->Ikf[5] + this->Itbr[3];
        const unsigned int m2 = j2 - this->Ikf[5] + this->Itbr[3];
        __m256d Kinf = _mm256_setr_pd(Kf_[j0+this->offset_kinf],Kf_[j1+this->offset_kinf],Kf_[j2+this->offset_kinf],1);    
     
        __m256d invKinf = _mm256_div_pd(_mm256_set1_pd(1.0),Kinf);
        __m256d dKinfdT = _mm256_setr_pd(dKfdT_[j0+this->offset_kinf],dKfdT_[j1+this->offset_kinf],dKfdT_[j2+this->offset_kinf],1);
        __m256d K0 = _mm256_setr_pd(Kf_[j0],Kf_[j1],Kf_[j2],1);        
        __m256d dK0dT = _mm256_setr_pd(dKfdT_[j0],dKfdT_[j1],dKfdT_[j2],1);     
        __m256d M = _mm256_setr_pd(tmp_M[m0],tmp_M[m1],tmp_M[m2],1);
        __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(K0,M),invKinf);
        __m256d dPrdT = _mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT));
        dPrdT = _mm256_mul_pd(dPrdT,invKinf);


        __m256d dKdT = (dK0dT);
        __m256d K = (K0);
        __m256d KK = (_mm256_mul_pd(K0,invKinf));
        __m256d tmp = _mm256_div_pd(one,_mm256_add_pd(one,Pr));
        __m256d N1 = (-tmp);
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
        __m256d KF = (N);
        this->Kf_[j0] = get_elem0(KF);
        this->Kf_[j1] = get_elem1(KF);
        this->Kf_[j2] = get_elem2(KF);
    }
}