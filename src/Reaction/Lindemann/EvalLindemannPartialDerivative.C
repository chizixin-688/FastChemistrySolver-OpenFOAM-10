/*---------------------------------------------------------------------------*\
  Description
      Computing the forward rate constant and the partial derivatives of 
      Lindemann reactions. Including Kf and dKfdT

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"

//=============================================================================//

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

        __m256d tmp = _mm256_div_pd(one,_mm256_add_pd(one,Pr));
        __m256d N1 = (tmp);
        __m256d N = _mm256_mul_pd(tmp,K0);
        __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(tmp,tmp),dPrdT),K);
        dKfdT = _mm256_fmadd_pd(tmp,dKdT,dKfdT);
        _mm256_storeu_pd(&this->dKfdT_[j0],dKfdT);
        __m256d dKfdC = _mm256_mul_pd((_mm256_mul_pd(K0,tmp)),N1);
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
        const unsigned int j = this->LindemannFO[i+0]+0;
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];        

        //__m256d Kinf = _mm256_setr_pd(Kf_[j + this->offset_kinf],Kf_[j1 + this->offset_kinf],1,1);  
        __m128d Kinf = _mm_loadu_pd(&Kf_[j+this->offset_kinf]);  

        __m128d invKinf = _mm_div_pd(_mm256_castpd256_pd128(one),Kinf);

        //__m128d dKinfdT = _mm256_setr_pd(dKfdT_[j + this->offset_kinf],dKfdT_[j1 + this->offset_kinf],1,1);
        __m128d dKinfdT = _mm_loadu_pd(&dKfdT_[j+this->offset_kinf]);

        __m128d K0 = _mm_loadu_pd(&Kf_[j]); 

        __m128d dK0dT = _mm_loadu_pd(&dKfdT_[j]);     
        __m128d M = _mm_loadu_pd(&tmp_M[m]);
        __m128d Pr = _mm_mul_pd(_mm_mul_pd(K0,M),invKinf);
        __m128d dPrdT = _mm_fmsub_pd(M,dK0dT,_mm_mul_pd(Pr,dKinfdT));
        dPrdT = _mm_mul_pd(dPrdT,invKinf);


        __m128d dKdT = _mm_mul_pd(Pr,dKinfdT);
        __m128d K = (Kinf);
        __m128d tmp = _mm_div_pd(_mm256_castpd256_pd128(one),_mm_add_pd(_mm256_castpd256_pd128(one),Pr));
        __m128d N1 = (tmp);
        __m128d N = _mm_mul_pd(tmp,K0);
        __m128d dKfdT = _mm_mul_pd(_mm_mul_pd(_mm_mul_pd(tmp,tmp),dPrdT),K);
        dKfdT = _mm_fmadd_pd(tmp,dKdT,dKfdT);

        //dKfdT_[j] = get_elem0(dKfdT);
        //dKfdT_[j1] = get_elem1(dKfdT);
        _mm_storeu_pd(&dKfdT_[j],dKfdT);

        __m128d dKfdC = _mm_mul_pd((_mm_mul_pd(K0,tmp)),N1);
        //dKfdC_[m] = get_elem0(dKfdC);
        //dKfdC_[m1] = get_elem1(dKfdC);
        _mm_storeu_pd(&dKfdC_[m],dKfdC);

        __m128d KF = (_mm_mul_pd(M,N));   
        //Kf_[j] = get_elem0(KF);
        //Kf_[j1] = get_elem1(KF);
        _mm_storeu_pd(&Kf_[j],KF);
    }
    else if (remainFO==3)
    {
        const unsigned int i = this->n_LindemannFO-3;
        const unsigned int j = this->LindemannFO[i+0]+0;

        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];

        const double Kinf0 = this->Kf_[j+0+this->offset_kinf];
        const double Kinf1 = this->Kf_[j+1+this->offset_kinf];
        const double Kinf2 = this->Kf_[j+2+this->offset_kinf];
        __m256d Kinfv = _mm256_setr_pd(Kinf0,Kinf1,Kinf2,1);    
     
        __m256d invKinfv = _mm256_div_pd(one,Kinfv);


        const double dKinfdT0 = this->dKfdT_[j+0+this->offset_kinf];
        const double dKinfdT1 = this->dKfdT_[j+1+this->offset_kinf];
        const double dKinfdT2 = this->dKfdT_[j+2+this->offset_kinf];
        __m256d dKinfdTv = _mm256_setr_pd(dKinfdT0,dKinfdT1,dKinfdT2,1);
        
        const double k00 = this->Kf_[j+0];
        const double k01 = this->Kf_[j+1];
        const double k02 = this->Kf_[j+2];
        __m256d K0v = _mm256_setr_pd(k00,k01,k02,1);

        
        const double dK0dT0 = this->dKfdT_[j+0];
        const double dK0dT1 = this->dKfdT_[j+1];
        const double dK0dT2 = this->dKfdT_[j+2];
        __m256d dK0dTv = _mm256_setr_pd(dK0dT0,dK0dT1,dK0dT2,1);     


        const double M0 = this->tmp_M[m+0];
        const double M1 = this->tmp_M[m+1];
        const double M2 = this->tmp_M[m+2];
        __m256d Mv = _mm256_setr_pd(M0,M1,M2,1);


        __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(K0v,Mv),invKinfv);
        __m256d dPrdT = _mm256_fmsub_pd(Mv,dK0dTv,_mm256_mul_pd(Pr,dKinfdTv));
        dPrdT = _mm256_mul_pd(dPrdT,invKinfv);


        __m256d dKdT = (_mm256_mul_pd(Pr,dKinfdTv));
        __m256d K = (Kinfv);
        __m256d tmp = _mm256_div_pd(one,_mm256_add_pd(one,Pr));
        __m256d N1 = (tmp);
        __m256d N = _mm256_mul_pd(tmp,K0v);
        __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(tmp,tmp),dPrdT),K);
        dKfdT = _mm256_fmadd_pd(tmp,dKdT,dKfdT);
        this->dKfdT_[j+0] = get_elem0(dKfdT);
        this->dKfdT_[j+1] = get_elem1(dKfdT);
        this->dKfdT_[j+2] = get_elem2(dKfdT);
        __m256d dKfdC = _mm256_mul_pd((_mm256_mul_pd(K0v,tmp)),N1);
        this->dKfdC_[m+00] = get_elem0(dKfdC);
        this->dKfdC_[m+1] = get_elem1(dKfdC);
        this->dKfdC_[m+2] = get_elem2(dKfdC);
        __m256d KF = (_mm256_mul_pd(Mv,N));
        this->Kf_[j+0] = get_elem0(KF);
        this->Kf_[j+1] = get_elem1(KF);
        this->Kf_[j+2] = get_elem2(KF);
    }

    const unsigned int remainCA = (this->n_LindemannCA)%4;
    for (unsigned int i = 0;i<this->n_LindemannCA-remainCA;i=i+4)
    {
        const unsigned int j0 = this->LindemannCA[i+0]+0;
        const unsigned int m0 = j0 - this->Ikf[5] + this->Itbr[3];
        __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0-this->Ikf[5]+this->Ikf[8]]);    
        __m256d invKinf = _mm256_div_pd(one,Kinf);
        __m256d dKinfdT = _mm256_loadu_pd(&this->dKfdT_[j0-this->Ikf[5]+this->Ikf[8]]);
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
        const double invOnePlusPr0 = 1.0/(1.0+Pr0);
        const double dKdT0   = dK0dT0;
        const double K0      = K00;
        const double KK0     = K00*invKinf0;
        const double N10     = -invOnePlusPr0;
        const double N0  = invOnePlusPr0*1*K00;
        // - k0*k0*dMdc/kInf/sqr(1 + Pr)

        const double r0 = invOnePlusPr0*dKdT0 - invOnePlusPr0*invOnePlusPr0*dPrdT0*K0;
        const double r1 = K00*KK0*(N10)*invOnePlusPr0; 
        this->dKfdT_[j0] = r0;
        this->dKfdC_[m0] = r1; 
        this->Kf_[j0] = N0;
    }
    else if (remainCA==2)
    {
        const unsigned int i = this->n_LindemannCA-2;
        const unsigned int j = this->LindemannCA[i+0]+0;
        const unsigned int m = j - this->Ikf[5] + this->Itbr[3];

        //__m256d Kinf = _mm256_setr_pd(Kf_[j0-this->Ikf[5]+this->Ikf[8]],Kf_[j1-this->Ikf[5]+this->Ikf[8]],1,1);  
        __m128d Kinfv = _mm_loadu_pd(&Kf_[j-this->Ikf[5]+this->Ikf[8]]);

        //__m256d invKinf = _mm256_div_pd(_mm256_set1_pd(1.0),Kinf);
        __m128d invKinfv = _mm_div_pd(_mm256_castpd256_pd128(one),Kinfv);

        //__m256d dKinfdT = _mm256_setr_pd(dKfdT_[j0-this->Ikf[5]+this->Ikf[8]],dKfdT_[j1-this->Ikf[5]+this->Ikf[8]],1,1);
        __m128d dKinfdTv = _mm_loadu_pd(&dKfdT_[j-this->Ikf[5]+this->Ikf[8]]);
        
        //__m256d K0 = _mm256_setr_pd(Kf_[j0],Kf_[j1],1,1);        
        __m128d K0v = _mm_loadu_pd(&Kf_[j]);

        //__m256d dK0dT = _mm256_setr_pd(dKfdT_[j0],dKfdT_[j1],1,1);     
        __m128d dK0dTv = _mm_loadu_pd(&dKfdT_[j]);

        //__m256d M = _mm256_setr_pd(tmp_M[m0],tmp_M[m1],1,1);
        __m128d Mv = _mm_loadu_pd(&tmp_M[m]);

        __m128d Prv = _mm_mul_pd(_mm_mul_pd(K0v,Mv),invKinfv);

        __m128d dPrdTv = _mm_fmsub_pd(Mv,dK0dTv,_mm_mul_pd(Prv,dKinfdTv));
        dPrdTv = _mm_mul_pd(dPrdTv,invKinfv);


        __m128d dKdTv = (dK0dTv);
        __m128d Kv = (K0v);
        __m128d KKv = (_mm_mul_pd(K0v,invKinfv));
        __m128d tmpv = _mm_div_pd(_mm256_castpd256_pd128(one),_mm_add_pd(_mm256_castpd256_pd128(one),Prv));
        __m128d N1v = (-tmpv);
        __m128d Nv = _mm_mul_pd(tmpv,K0v);
        __m128d dKfdTv = _mm_mul_pd(_mm_mul_pd(_mm_mul_pd(tmpv,tmpv),dPrdTv),Kv);
        dKfdTv = _mm_fmadd_pd(tmpv,dKdTv,dKfdTv);

        //dKfdT_[j0] = get_elem0(dKfdT);
        //dKfdT_[j1] = get_elem1(dKfdT);
        _mm_storeu_pd(&dKfdT_[j],dKfdTv);

        __m128d dKfdCv = _mm_mul_pd(_mm_mul_pd(_mm_mul_pd(K0v,tmpv),KKv),N1v);
        //dKfdC_[m0] = get_elem0(dKfdC);
        //dKfdC_[m1] = get_elem1(dKfdC);
        _mm_storeu_pd(&dKfdC_[m],dKfdCv);

        __m128d kfv = (Nv);   
        //Kf_[j0] = get_elem0(KF);
        //Kf_[j1] = get_elem1(KF);
        _mm_storeu_pd(&Kf_[j],kfv);

    }
    else if (remainCA==3)
    {
        const unsigned int i = this->n_LindemannCA-3;
        const unsigned int j = this->LindemannCA[i+0]+0;

        const unsigned int m = j+0-this->Ikf[5] + this->Itbr[3];

        const double Kinf0 = this->Kf_[j+0-this->Ikf[5]+this->Ikf[8]];
        const double Kinf1 = this->Kf_[j+1-this->Ikf[5]+this->Ikf[8]];
        const double Kinf2 = this->Kf_[j+2-this->Ikf[5]+this->Ikf[8]];
        __m256d Kinfv = _mm256_setr_pd(Kinf0,Kinf1,Kinf2,1);    
     
        __m256d invKinfv = _mm256_div_pd(one,Kinfv);

        const double dKinfdT0 = this->dKfdT_[j+0-this->Ikf[5]+this->Ikf[8]];
        const double dKinfdT1 = this->dKfdT_[j+1-this->Ikf[5]+this->Ikf[8]];
        const double dKinfdT2 = this->dKfdT_[j+2-this->Ikf[5]+this->Ikf[8]];
        __m256d dKinfdTv = _mm256_setr_pd(dKinfdT0,dKinfdT1,dKinfdT2,1);

        const double K0 = this->Kf_[j+0];
        const double K1 = this->Kf_[j+1];
        const double K2 = this->Kf_[j+2];
        __m256d K0v = _mm256_setr_pd(K0,K1,K2,1);        

        const double dK0dT0 = this->dKfdT_[j+0];
        const double dK0dT1 = this->dKfdT_[j+1];
        const double dK0dT2 = this->dKfdT_[j+2];
        __m256d dK0dTv = _mm256_setr_pd(dK0dT0,dK0dT1,dK0dT2,1);  
        
        const double M0 = this->tmp_M[m+0];
        const double M1 = this->tmp_M[m+1];
        const double M2 = this->tmp_M[m+2];

        __m256d M = _mm256_setr_pd(M0,M1,M2,1);

        __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(K0v,M),invKinfv);
        __m256d dPrdT = _mm256_fmsub_pd(M,dK0dTv,_mm256_mul_pd(Pr,dKinfdTv));
        dPrdT = _mm256_mul_pd(dPrdT,invKinfv);


        __m256d dKdT = (dK0dTv);
        __m256d K = (K0v);
        __m256d KK = (_mm256_mul_pd(K0v,invKinfv));
        __m256d tmp = _mm256_div_pd(one,_mm256_add_pd(one,Pr));
        __m256d N1 = (-tmp);
        __m256d N = _mm256_mul_pd(tmp,K0v);
        __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(tmp,tmp),dPrdT),K);
        dKfdT = _mm256_fmadd_pd(tmp,dKdT,dKfdT);
        this->dKfdT_[j+0] = get_elem0(dKfdT);
        this->dKfdT_[j+1] = get_elem1(dKfdT);
        this->dKfdT_[j+2] = get_elem2(dKfdT);
        __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0v,tmp),KK),N1);
        this->dKfdC_[m+0] = get_elem0(dKfdC);
        this->dKfdC_[m+1] = get_elem1(dKfdC);
        this->dKfdC_[m+2] = get_elem2(dKfdC);
        __m256d KF = (N);
        this->Kf_[j+0] = get_elem0(KF);
        this->Kf_[j+1] = get_elem1(KF);
        this->Kf_[j+2] = get_elem2(KF);
    }
}