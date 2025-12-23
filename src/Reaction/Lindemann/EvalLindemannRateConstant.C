#include "OptReaction.H"
void FastChemistry::OptReaction::evalLindemannRateConstant()const noexcept
{
    __m256d one = _mm256_set1_pd(1.0);


    unsigned int remainFO = this->n_LindemannFO%4;
    for (unsigned int i = 0;i<this->n_LindemannFO-remainFO;i=i+4)
    {
        const unsigned int j0 = this->LindemannFO[i+0]+0;
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
        __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
        __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
        __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);
        __m256d Pr = _mm256_div_pd(_mm256_mul_pd(K0,M),Kinf);
        __m256d N = _mm256_div_pd(K0,_mm256_add_pd(Pr,one));
        __m256d Kf = _mm256_mul_pd(M,N);
        _mm256_storeu_pd(&this->Kf_[j0],Kf);
    }
    if(remainFO==1)
    {
        const unsigned int i = this->n_LindemannFO-1;
        const unsigned int j = this->LindemannFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];
        //const unsigned int k = j - this->Ikf[4];
        const double Kinf = this->Kf_[j+this->offset_kinf];
        double M = this->tmp_M[m];     
        const double K0 = this->Kf_[j];
        const double Pr = K0*M/Kinf;
        const double N  = 1/(1+Pr)*K0;
        this->Kf_[j] = M*N;            
    }
    else if(remainFO==2)
    {
        const unsigned int i = this->n_LindemannFO-2;
        const unsigned int j0 = this->LindemannFO[i+0]+0;
        //const unsigned int j1 = this->LindemannFO[i+0]+1;
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
        //const unsigned int k0 = j0 - this->Ikf[4];
        //const unsigned int k1 = j1 - this->Ikf[4];
        __m128d Kinf = _mm_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
        __m128d M = _mm_loadu_pd(&this->tmp_M[m0]);
        __m128d K0 = _mm_loadu_pd(&this->Kf_[j0]);
        __m128d Pr = _mm_div_pd(_mm_mul_pd(K0,M),Kinf);
        __m128d one128 = _mm256_castpd256_pd128(one);
        __m128d N = _mm_div_pd(K0,_mm_add_pd(Pr,one128));
        __m128d Kf = _mm_mul_pd(M,N);
        _mm_storeu_pd(&this->Kf_[j0],Kf);
    }
    else if(remainFO==3)
    {
        const unsigned int i = this->n_LindemannFO-3;
        const unsigned int j0 = this->LindemannFO[i+0];
        const unsigned int j1 = this->LindemannFO[i+1];
        const unsigned int j2 = this->LindemannFO[i+2];
        const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
        const unsigned int m1 = j1 - this->Ikf[4] + this->Itbr[2];
        const unsigned int m2 = j2 - this->Ikf[4] + this->Itbr[2];
        //const unsigned int k0 = j0 - this->Ikf[4];
        //const unsigned int k1 = j1 - this->Ikf[4];
        //const unsigned int k2 = j2 - this->Ikf[4];
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

        __m256d Kf = _mm256_mul_pd(M,N);
        this->Kf_[j0] = get_elem0(Kf);
        this->Kf_[j1] = get_elem1(Kf);
        this->Kf_[j2] = get_elem2(Kf);
    }

    unsigned int remainCA = this->n_LindemannCA%4;
    for (unsigned int i = 0;i<this->n_LindemannCA-remainCA;i=i+4)
    {
        const unsigned int j0 = this->LindemannCA[i+0]+0;
        const unsigned int m0 = j0 - this->Ikf[5] + this->Itbr[3];
        __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
        __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
        __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);
        __m256d Pr = _mm256_div_pd(_mm256_mul_pd(K0,M),Kinf);
        __m256d N = _mm256_div_pd(K0,_mm256_add_pd(Pr,one));
        __m256d Kf = (N);
        _mm256_storeu_pd(&this->Kf_[j0],Kf);
    }
    if(remainCA==1)
    {
        const unsigned int i = this->n_LindemannCA-1;
        const unsigned int j = this->LindemannCA[i];
        const unsigned int m = j - this->Ikf[5] + this->Itbr[3];

        const double Kinf = this->Kf_[j- Ikf[5] + Ikf[8]];
        double M = this->tmp_M[m];     
        const double K0 = this->Kf_[j];
        const double Pr = K0*M/Kinf;
        const double N  = 1/(1+Pr)*K0;
        this->Kf_[j] = N;            
    }
    else if(remainCA==2)
    {
        const unsigned int i = this->n_LindemannCA-2;
        const unsigned int j0 = this->LindemannCA[i+0]+0;
        //const unsigned int j1 = this->LindemannCA[i+0]+1;
        const unsigned int m0 = j0 - this->Ikf[5] + this->Itbr[3];

        __m128d Kinf = _mm_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
        __m128d M = _mm_loadu_pd(&this->tmp_M[m0]);
        __m128d K0 = _mm_loadu_pd(&this->Kf_[j0]);
        __m128d Pr = _mm_div_pd(_mm_mul_pd(K0,M),Kinf);
        __m128d one128 = _mm256_castpd256_pd128(one);
        __m128d N = _mm_div_pd(K0,_mm_add_pd(Pr,one128));
        __m128d Kf = (N);
        _mm_storeu_pd(&this->Kf_[j0],Kf);
        
    }
    else if(remainCA==3)
    {
        const unsigned int i = this->n_LindemannCA-3;
        const unsigned int j0 = this->LindemannCA[i+0];
        const unsigned int j1 = this->LindemannCA[i+1];
        const unsigned int j2 = this->LindemannCA[i+2];
        const unsigned int m0 = j0 - this->Ikf[5] + this->Itbr[3];
        const unsigned int m1 = j1 - this->Ikf[5] + this->Itbr[3];
        const unsigned int m2 = j2 - this->Ikf[5] + this->Itbr[3];

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

        __m256d Kf = (N);
        this->Kf_[j0] = get_elem0(Kf);
        this->Kf_[j1] = get_elem1(Kf);
        this->Kf_[j2] = get_elem2(Kf);
    }
}