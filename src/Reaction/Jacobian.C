#include "OptReaction.H"
#include <immintrin.h>  


void 
OptReaction::ddNdtByVdcTp
(
    double p,
    double Temperature,
    double* __restrict__ Phi,
    double* __restrict__ c,
    double* __restrict__ dNdtByV,
    double* __restrict__ dBdT,
    double* __restrict__ dCpdT,
    double* __restrict__ Cp,
    double* __restrict__ Ha,
    double* __restrict__ rhoMByRhoi,
    double* __restrict__ WiByrhoM,
    double* __restrict__ ddNdtByVdcT
) const noexcept
{
    Temperature = Temperature<TlowMin?TlowMin:Temperature;
    Temperature = Temperature>ThighMax?ThighMax:Temperature;
    this->logT = std::log(Temperature);
    this->T = Temperature;
    this->invT = 1/Temperature;

    this->setPtrCoeffs(Temperature);

    {
        this->JacobianThermo
        (
            p,
            Temperature,
            Phi,
            c,
            &this->tmp_Exp[0],
            dBdT,
            dCpdT,
            Cp,
            Ha,
            rhoMByRhoi,
            WiByrhoM
        );
    }

    this->update_Pow_pByRT_SumVki2(Temperature);

    for(size_t i = 0; i <this->Troe.size();i++)
    {
        size_t j0 = i + this->nSpecies;
        size_t j1 = i + this->nSpecies + this->Troe.size();
        size_t j2 = i + this->nSpecies + this->Troe.size()*2;         
        this->tmp_Exp[j0] = -Temperature*this->invTsss_[i];            
        this->tmp_Exp[j1] = -this->Tss_[i]*invT; 
        this->tmp_Exp[j2] = -Temperature*this->invTs_[i];            
    }
    
    for(size_t i = 0; i <this->SRI.size();i++)
    {
        size_t j0 = i + this->nSpecies + this->Troe.size()*3;
        size_t j1 = i + this->nSpecies + this->Troe.size()*3 + this->SRI.size();
        this->tmp_Exp[j0] = -this->b_[i]*invT;
        this->tmp_Exp[j1] = -Temperature*this->invc_[i];            
    }   



    {
        size_t remain = this->tmp_ExpSize%4;
        for(size_t i = 0; i < this->tmp_ExpSize-remain;i=i+4)
        {
            __m256d tmp = _mm256_loadu_pd(&this->tmp_Exp[i]);
            tmp = vec256_expd(tmp);
            _mm256_storeu_pd(&this->tmp_Exp[i],tmp);
        }
        if(remain==1)
        {
            size_t i = this->tmp_ExpSize-1;
            this->tmp_Exp[i] = std::exp(this->tmp_Exp[i]);
        }
        else if(remain==2)
        {
            size_t i0 = this->tmp_ExpSize-2;
            size_t i1 = this->tmp_ExpSize-1;
            __m256d tmp = _mm256_setr_pd(tmp_Exp[i0],tmp_Exp[i1],0,0);
            tmp = vec256_expd(tmp);
            this->tmp_Exp[i0] = get_elem0(tmp);
            this->tmp_Exp[i1] = get_elem1(tmp);
        }
        else if(remain==3)
        {
            size_t i0 = this->tmp_ExpSize-3;
            size_t i1 = this->tmp_ExpSize-2;
            size_t i2 = this->tmp_ExpSize-1;
            __m256d tmp = _mm256_setr_pd(tmp_Exp[i0],tmp_Exp[i1],tmp_Exp[i2],0);
            tmp = vec256_expd(tmp);
            this->tmp_Exp[i0] = get_elem0(tmp);
            this->tmp_Exp[i1] = get_elem1(tmp);
            this->tmp_Exp[i2] = get_elem2(tmp);
        }
    }

    if(this->n_PlogReaction>0)
    {
        this->logP = std::log(p);
        for(unsigned int i = 0; i< this->n_PlogReaction; i ++)
        {
            const size_t length = this->Prange[i].size();
            if(p<=this->Prange[i][0])
            {
                double A0 = this->APlog[i][0];
                double beta0 = this->betaPlog[i][0];
                double Ta0 = this->TaPlog[i][0];

                A[i+this->Ikf[6]] = A0;
                A[i+this->Ikf[11]] = A0;
                beta[i+this->Ikf[6]] = beta0;
                beta[i+this->Ikf[11]] = beta0;
                Ta[i+this->Ikf[6]] = Ta0;
                Ta[i+this->Ikf[11]] = Ta0;
                this->Pindex[i] = 0;
            }
            else if(p>=this->Prange[i][length-1])
            {
                double A1 = this->APlog[i][length-1];
                double beta1 = this->betaPlog[i][length-1];
                double Ta1 = this->TaPlog[i][length-1];
                
                A[i+this->Ikf[6]] = A1;
                A[i+this->Ikf[11]] = A1;
                beta[i+this->Ikf[6]] = beta1;
                beta[i+this->Ikf[11]] = beta1;
                Ta[i+this->Ikf[6]] = Ta1;
                Ta[i+this->Ikf[11]] = Ta1;
                this->Pindex[i] = static_cast<unsigned int>(length-1);
            }
            else
            {
                unsigned int index = 0;
                for(unsigned int j = 0; j < length-1;j++)
                {
                    if(this->Prange[i][j]<=p && p<this->Prange[i][j+1])
                    {
                        index = j;
                        break;
                    }
                }
                this->Pindex[i] = index;
                double A0 = this->APlog[i][index+0];
                double A1 = this->APlog[i][index+1];
                double beta0 = this->betaPlog[i][index+0];
                double beta1 = this->betaPlog[i][index+1];
                double Ta0 = this->TaPlog[i][index+0];
                double Ta1 = this->TaPlog[i][index+1];
                A[i+this->Ikf[6]] = A0;
                A[i+this->Ikf[11]] = A1;
                beta[i+this->Ikf[6]] = beta0;
                beta[i+this->Ikf[11]] = beta1;
                Ta[i+this->Ikf[6]] = Ta0;
                Ta[i+this->Ikf[11]] = Ta1;
            }
        }
    }

    {
        __m256d LogTv = _mm256_set1_pd(logT);
        __m256d InvTv = _mm256_set1_pd(invT);        
        unsigned int remain = (this->Ikf[12]-this->n_Temperature_Independent_Reaction)%4;
        unsigned int times = (this->Ikf[12]-this->n_Temperature_Independent_Reaction)/4;
        for(unsigned int z = 0; z <times;z=z+1)
        {
            unsigned int i = z*4 + this->n_Temperature_Independent_Reaction;
            __m256d betav = _mm256_loadu_pd(&this->beta[i]);
            __m256d Tav = _mm256_loadu_pd(&this->Ta[i]);

            __m256d Kfv = _mm256_mul_pd(Tav,-InvTv);
            Kfv = _mm256_fmadd_pd(betav,LogTv,Kfv);
            __m256d A_ = _mm256_loadu_pd(&this->A[i]);
            Kfv = vec256_expd(Kfv);
            Kfv = _mm256_mul_pd(A_,Kfv);
            _mm256_storeu_pd(&this->Kf_[i],Kfv);
            __m256d dKfdT = _mm256_mul_pd(_mm256_fmadd_pd(Tav,InvTv,betav),InvTv);
            dKfdT = _mm256_mul_pd(dKfdT,Kfv); 
            _mm256_storeu_pd(&this->dKfdT_[i+0],dKfdT);           
        }
        if(remain==1)
        {
            unsigned int i = this->Ikf[12]-1;
            this->Kf_[i] = this->A[i]*std::exp(this->beta[i+0]*logT-this->Ta[i+0]*invT);   
            this->dKfdT_[i+0] = this->Kf_[i+0]*(this->beta[i+0]+this->Ta[i+0]*invT)*invT;  
        }
        else if(remain==2)
        {
            unsigned int i0 = this->Ikf[12]-2;
            unsigned int i1 = this->Ikf[12]-1;    
            __m256d betav = _mm256_setr_pd(beta[i0],beta[i1],0,0);
            __m256d Av = _mm256_setr_pd(A[i0],A[i1],0,0);
            __m256d Tav = _mm256_setr_pd(Ta[i0],Ta[i1],0,0);
            __m256d tmp = _mm256_fmsub_pd(betav,LogTv,_mm256_mul_pd(Tav,InvTv));
            tmp = vec256_expd(tmp);
            __m256d Kfv = _mm256_mul_pd(Av,tmp);
            this->Kf_[i0] = get_elem0(Kfv);
            this->Kf_[i1] = get_elem1(Kfv);
            tmp = _mm256_fmadd_pd(Tav,InvTv,betav);
            __m256d dKfdTv = _mm256_mul_pd(Kfv,_mm256_mul_pd(tmp,InvTv));
            this->dKfdT_[i0] = get_elem0(dKfdTv);
            this->dKfdT_[i1] = get_elem1(dKfdTv);
        }
        else if(remain==3)
        {
            unsigned int i0 = this->Ikf[12]-3;
            unsigned int i1 = this->Ikf[12]-2;
            unsigned int i2 = this->Ikf[12]-1;    
            __m256d betav = _mm256_setr_pd(beta[i0],beta[i1],beta[i2],0);
            __m256d Av = _mm256_setr_pd(A[i0],A[i1],A[i2],0);
            __m256d Tav = _mm256_setr_pd(Ta[i0],Ta[i1],Ta[i2],0);
            __m256d tmp = _mm256_fmsub_pd(betav,LogTv,_mm256_mul_pd(Tav,InvTv));
            tmp = vec256_expd(tmp);
            __m256d Kfv = _mm256_mul_pd(Av,tmp);
            this->Kf_[i0] = get_elem0(Kfv);
            this->Kf_[i1] = get_elem1(Kfv);
            this->Kf_[i2] = get_elem2(Kfv);
            tmp = _mm256_fmadd_pd(Tav,InvTv,betav);
            __m256d dKfdTv = _mm256_mul_pd(Kfv,_mm256_mul_pd(tmp,InvTv));
            this->dKfdT_[i0] = get_elem0(dKfdTv);  
            this->dKfdT_[i1] = get_elem1(dKfdTv); 
            this->dKfdT_[i2] = get_elem2(dKfdTv); 
        }
    }



    if(this->n_PlogReaction>0)
    {
        for(unsigned int i = 0; i< this->n_PlogReaction; i ++)
        {
            const size_t length = this->Prange[i].size();
            if(this->Pindex[i] == 0 || this->Pindex[i] == length-1)
            {
                continue;
            }
            else
            {
                unsigned int index = this->Pindex[i];
                double weight = (this->logP - this->logPi[i][index])*this->rDeltaP_[i][index];
                double Kf0 = this->Kf_[i+this->Ikf[6]];
                double Kf1 = this->Kf_[i+this->Ikf[11]];
                double Kf = Kf0*std::pow(Kf1/Kf0,weight);
                this->Kf_[i+this->Ikf[6]] = Kf;

                double beta0 = this->beta[i+this->Ikf[6]];
                double beta1 = this->beta[i+this->Ikf[11]];
                double Ta0 = this->Ta[i+this->Ikf[6]];
                double Ta1 = this->Ta[i+this->Ikf[11]];
                double invt = this->invT;

                double dKfdT = Kf*invt*(beta0 + Ta0*invt + (beta1-beta0+(Ta1-Ta0)*invt)*weight);
                this->dKfdT_[i+this->Ikf[6]] = dKfdT;
            }
        }
    }
    {
        
        unsigned int Tremain = (this->Itbr[5])%4;
        for(unsigned int i = 0; i < this->Itbr[5]-Tremain; i=i+4)
        {
            __m256d arrM_0 = _mm256_setzero_pd();
            __m256d arrM_1 = _mm256_setzero_pd();
            __m256d arrM_2 = _mm256_setzero_pd();
            __m256d arrM_3 = _mm256_setzero_pd();
            double* __restrict__ TBF1DRowi0 = &ThirdBodyFactor1D[(i+0)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi1 = &ThirdBodyFactor1D[(i+1)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi2 = &ThirdBodyFactor1D[(i+2)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi3 = &ThirdBodyFactor1D[(i+3)*this->AlignSpecies];
            for(unsigned int j  = 0;j<this->AlignSpecies;j=j+4)
            {
                __m256d Factor0 = _mm256_loadu_pd(&TBF1DRowi0[j+0]);
                __m256d Factor1 = _mm256_loadu_pd(&TBF1DRowi1[j+0]);
                __m256d Factor2 = _mm256_loadu_pd(&TBF1DRowi2[j+0]);
                __m256d Factor3 = _mm256_loadu_pd(&TBF1DRowi3[j+0]);
                __m256d C_ = _mm256_loadu_pd(&c[j+0]);
                
                arrM_0 = _mm256_fmadd_pd(Factor0,C_,arrM_0);
                arrM_1 = _mm256_fmadd_pd(Factor1,C_,arrM_1);
                arrM_2 = _mm256_fmadd_pd(Factor2,C_,arrM_2);
                arrM_3 = _mm256_fmadd_pd(Factor3,C_,arrM_3);
            }
            __m256d s0h = _mm256_hadd_pd(arrM_0, arrM_1); 
            __m256d s1h = _mm256_hadd_pd(arrM_2, arrM_3); 
            s0h = _mm256_permute4x64_pd(s0h, 0b11011000);
            s1h = _mm256_permute4x64_pd(s1h, 0b11011000);
            __m256d sum_all = _mm256_hadd_pd(s0h, s1h); 
            sum_all = _mm256_permute4x64_pd(sum_all, 0b11011000);

           _mm256_storeu_pd(&this->tmp_M[i+0],sum_all);

        }
        if(Tremain==3)
        {
            unsigned int i =(this->Itbr[5]) -3;
            double* __restrict__ TBF1DRowi0 = &ThirdBodyFactor1D[(i+0)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi1 = &ThirdBodyFactor1D[(i+1)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi2 = &ThirdBodyFactor1D[(i+2)*this->AlignSpecies];
            double M0 = 0;
            double M1 = 0;           
            double M2 = 0; 
            __m256d arrM_0 = _mm256_setzero_pd();
            __m256d arrM_1 = _mm256_setzero_pd();
            __m256d arrM_2 = _mm256_setzero_pd();
            for(unsigned int j  = 0;j<this->AlignSpecies;j=j+4)
            {
                __m256d Factor0 = _mm256_loadu_pd(&TBF1DRowi0[j+0]);
                __m256d Factor1 = _mm256_loadu_pd(&TBF1DRowi1[j+0]);
                __m256d Factor2 = _mm256_loadu_pd(&TBF1DRowi2[j+0]);
                __m256d C_ = _mm256_loadu_pd(&c[j+0]);
                arrM_0 = _mm256_fmadd_pd(Factor0,C_,arrM_0);
                arrM_1 = _mm256_fmadd_pd(Factor1,C_,arrM_1);
                arrM_2 = _mm256_fmadd_pd(Factor2,C_,arrM_2);
            }

            M0 = M0 + hsum4(arrM_0);
            M1 = M1 + hsum4(arrM_1);
            M2 = M2 + hsum4(arrM_2);

            this->tmp_M[i+0] = M0;
            this->tmp_M[i+1] = M1;
            this->tmp_M[i+2] = M2;
        }
        else if(Tremain==2)
        {
            unsigned int i =(this->Itbr[5]) -2;
            double* __restrict__ TBF1DRowi0 = &ThirdBodyFactor1D[(i+0)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi1 = &ThirdBodyFactor1D[(i+1)*this->AlignSpecies];
            double M0 = 0;
            double M1 = 0;           
            __m256d arrM_0 = _mm256_setzero_pd();
            __m256d arrM_1 = _mm256_setzero_pd();
            for(unsigned int j  = 0;j<this->AlignSpecies;j=j+4)
            {
                __m256d Factor0 = _mm256_loadu_pd(&TBF1DRowi0[j+0]);
                __m256d Factor1 = _mm256_loadu_pd(&TBF1DRowi1[j+0]);
                __m256d C_ = _mm256_loadu_pd(&c[j+0]);
                arrM_0 = _mm256_fmadd_pd(Factor0,C_,arrM_0);
                arrM_1 = _mm256_fmadd_pd(Factor1,C_,arrM_1);
            }

            M0 = M0 + hsum4(arrM_0);
            M1 = M1 + hsum4(arrM_1);

            this->tmp_M[i+0] = M0;
            this->tmp_M[i+1] = M1;
        }
        else if(Tremain==1)
        {
            unsigned int i =(this->Itbr[5]) -1;
            double* __restrict__ TBF1DRowi0 = &ThirdBodyFactor1D[(i+0)*this->AlignSpecies];            
            double M0 = 0;
            __m256d arrM_0 = _mm256_setzero_pd();
            for(unsigned int j  = 0;j<this->AlignSpecies;j=j+4)
            {
                __m256d Factor0 = _mm256_loadu_pd(&TBF1DRowi0[j+0]);
                __m256d C_ = _mm256_loadu_pd(&c[j+0]);
                arrM_0 = _mm256_fmadd_pd(Factor0,C_,arrM_0);
            }
            M0 = M0 + hsum4(arrM_0);
            this->tmp_M[i+0] = M0;
        }
    }

    for(unsigned int i = 0; i < this->n_ThirdBodyReaction; i++)
    {
        double M = this->tmp_M[i+this->Itbr[1]];
        const unsigned int j = i + this->Ikf[3];
        this->dKfdC_[i+this->Itbr[1]] = this->Kf_[j];   
        this->Kf_[j] = this->Kf_[j]*M;
        this->dKfdT_[j] = this->dKfdT_[j]*M;
    }
    


    for(unsigned int i = 0; i < this->n_NonEquilibriumThirdBodyReaction; i++)
    {
        double Mfwd = this->tmp_M[i];
        double Mrev = this->tmp_M[this->Itbr[4]+i];
        this->dKfdC_[i] = this->Kf_[this->Ikf[2]+i];
        this->dKfdC_[this->Itbr[4]+i] = this->Kf_[this->Ikf[10]+i];
        this->Kf_[this->Ikf[2]+i] = this->Kf_[this->Ikf[2]+i]*Mfwd;
        this->dKfdT_[this->Ikf[2]+i] = this->dKfdT_[this->Ikf[2]+i]*Mfwd;
        this->Kf_[this->Ikf[10]+i] = this->Kf_[this->Ikf[10]+i]*Mrev;
        this->dKfdT_[this->Ikf[10]+i] = this->dKfdT_[this->Ikf[10]+i]*Mrev;
    } 


    {
        size_t remain_Lindemann = (Lindemann.size())%4;
        for (size_t i = 0;i<Lindemann.size()-remain_Lindemann;i=i+4)
        {
            const unsigned int j0 = this->Lindemann[i+0]+0;
            const unsigned int j1 = this->Lindemann[i+0]+1;
            const unsigned int j2 = this->Lindemann[i+0]+2;
            const unsigned int j3 = this->Lindemann[i+0]+3;
            const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
            const unsigned int k0 = j0 - this->Ikf[4];
            const unsigned int k1 = j1 - this->Ikf[4];
            const unsigned int k2 = j2 - this->Ikf[4];
            const unsigned int k3 = j3 - this->Ikf[4];
            __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0 + this->offset_kinf]);    
            __m256d one = _mm256_set1_pd(1.0);       
            __m256d invKinf = _mm256_div_pd(_mm256_set1_pd(1.0),Kinf);
            __m256d dKinfdT = _mm256_loadu_pd(&this->dKfdT_[j0 + this->offset_kinf]);
            __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);        
            __m256d dK0dT = _mm256_loadu_pd(&this->dKfdT_[j0]);     
            __m256d M = _mm256_loadu_pd(&tmp_M[m0]);
            __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(K0,M),invKinf);
            __m256d dPrdT = _mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT));
            dPrdT = _mm256_mul_pd(dPrdT,invKinf);
            __m256d k = _mm256_setr_pd(k0,k1,k2,k3);
            __m256d cmp = _mm256_cmp_pd(k,_mm256_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
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
            __m256d one = _mm256_set1_pd(1.0);       
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
            __m256d one = _mm256_set1_pd(1.0);       
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


    {
        size_t remain_Troe = (Troe.size())%4;
        for (size_t i = 0;i<Troe.size()-remain_Troe;i=i+4)
        {
            const unsigned int j0 = this->Troe[i+0];
            const unsigned int j1 = this->Troe[i+1];
            const unsigned int j2 = this->Troe[i+2];
            const unsigned int j3 = this->Troe[i+3];
            const unsigned int k0 = j0 - this->Ikf[4];
            const unsigned int k1 = j1 - this->Ikf[4];
            const unsigned int k2 = j2 - this->Ikf[4];
            const unsigned int k3 = j3 - this->Ikf[4];
            const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
            __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
            __m256d invKinf = _mm256_div_pd(_mm256_set1_pd(1.0),Kinf);
            __m256d dKinfdT = _mm256_loadu_pd(&this->dKfdT_[j0+this->offset_kinf]);
            __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);           
            __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
            __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(M,K0),invKinf);
            __m256d small = _mm256_set1_pd(2.2e-16);
            __m256d cmp_result_Pr = _mm256_cmp_pd(Pr,small,_CMP_GE_OQ);
            Pr = _mm256_add_pd(Pr,_mm256_set1_pd(1e-100));
            const double invLog10 = 1.0/std::log(10);
            __m256d logPr_ = _mm256_mul_pd(vec256_logd(_mm256_max_pd(small,Pr)),_mm256_set1_pd(invLog10));
            __m256d InvTsss = _mm256_loadu_pd(&this->invTsss_[i]);
            __m256d InvTs = _mm256_loadu_pd(&this->invTs_[i]);
            __m256d Tss = _mm256_loadu_pd(&this->Tss_[i]);
            __m256d expTTsss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies]);
            __m256d expTTss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->Troe.size()]);
            __m256d expTTs = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->Troe.size()*2]);
            __m256d one = _mm256_set1_pd(1.0);
            __m256d alpha = _mm256_loadu_pd(&this->alpha_[i]);
            __m256d Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha),expTTsss);
            Fcent = _mm256_fmadd_pd(alpha,expTTs,Fcent);
            Fcent = _mm256_add_pd(expTTss,Fcent);
            __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),_mm256_set1_pd(invLog10));
            __m256d cc = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(0.67),_mm256_set1_pd(0.4));
            __m256d n = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(-1.27),_mm256_set1_pd(0.75));
            __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(cc,logPr_),_mm256_set1_pd(0.14),n);
            __m256d invx1 = _mm256_div_pd(one,x1);
            __m256d x2 = _mm256_mul_pd(_mm256_sub_pd(logPr_,cc),invx1);
            __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
            __m256d invx3 = _mm256_div_pd(one,x3);
            __m256d x4 = _mm256_mul_pd(logFcent,invx3);
            __m256d  F = vec256_powd(_mm256_set1_pd(10),x4);
            __m256d logTen = _mm256_set1_pd(std::log(10));
            __m256d dFcentdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_sub_pd(alpha,one),InvTsss),expTTsss);
            dFcentdT = _mm256_sub_pd(dFcentdT,_mm256_mul_pd(_mm256_mul_pd(alpha,InvTs),expTTs));
            __m256d invT2 = _mm256_set1_pd(invT*invT);
            dFcentdT = _mm256_fmadd_pd(expTTss,_mm256_mul_pd(Tss,invT2),dFcentdT);
            __m256d cmp2 = _mm256_cmp_pd(Fcent,small,_CMP_GE_OQ);
            __m256d dlogFcentdT = _mm256_div_pd(_mm256_div_pd(dFcentdT,_mm256_max_pd(Fcent,small)),logTen);
            dlogFcentdT = _mm256_blendv_pd(_mm256_setzero_pd(), dlogFcentdT, cmp2);
            __m256d dcdT = _mm256_mul_pd(dlogFcentdT,_mm256_set1_pd(-0.67));
            __m256d dndT = _mm256_mul_pd(dlogFcentdT,_mm256_set1_pd(-1.27));
            __m256d dx1dT = _mm256_fmadd_pd(dcdT,_mm256_set1_pd(-0.14),dndT);
            __m256d dx2dT = _mm256_mul_pd(_mm256_sub_pd(dcdT,_mm256_mul_pd(x2,dx1dT)),invx1);
            __m256d dx3dT = _mm256_mul_pd(_mm256_mul_pd(x2,dx2dT),_mm256_set1_pd(2.0));
            __m256d dx4dT = _mm256_mul_pd(_mm256_sub_pd(dlogFcentdT,_mm256_mul_pd(x4,dx3dT)),invx3);
            __m256d dFdT = _mm256_mul_pd(logTen,_mm256_mul_pd(F,dx4dT));
            __m256d dlogPrdPr = _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(Pr,logTen));
            dlogPrdPr = _mm256_blendv_pd(_mm256_setzero_pd(), dlogPrdPr, cmp_result_Pr);
            __m256d dx1dPr = _mm256_mul_pd(dlogPrdPr,_mm256_set1_pd(-0.14));
            __m256d dx2dPr = _mm256_mul_pd(_mm256_sub_pd(dlogPrdPr,_mm256_mul_pd(x2,dx1dPr)),invx1);
            __m256d dx3dPr = _mm256_mul_pd(_mm256_mul_pd(x2,dx2dPr),_mm256_set1_pd(2.0));
            __m256d dx4dPr = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(-1.0),x4),dx3dPr),invx3);
            __m256d dFdPr  = _mm256_mul_pd(_mm256_mul_pd(logTen,F),dx4dPr);    
            __m256d dK0dT = _mm256_loadu_pd(&this->dKfdT_[j0]);            
            __m256d dPrdT = _mm256_mul_pd(_mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT)),invKinf);
            dFdT = _mm256_fmadd_pd(dFdPr,dPrdT,dFdT);
            __m256d k = _mm256_setr_pd(k0,k1,k2,k3);
            __m256d cmp = _mm256_cmp_pd(k,_mm256_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
            __m256d dKdT = _mm256_blendv_pd(dK0dT, _mm256_mul_pd(Pr,dKinfdT), cmp);
            __m256d K = _mm256_blendv_pd(K0, Kinf, cmp);
            __m256d MM = _mm256_blendv_pd(_mm256_set1_pd(1), M, cmp);
            __m256d KK = _mm256_blendv_pd(_mm256_mul_pd(K0,invKinf), _mm256_set1_pd(1), cmp);
            __m256d invOnePlusPr = _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),Pr));
            __m256d N1 = _mm256_mul_pd(F,invOnePlusPr);
            N1 = _mm256_blendv_pd(_mm256_sub_pd(_mm256_setzero_pd(),N1),N1,cmp);            
            __m256d N2 = _mm256_blendv_pd(dFdPr,_mm256_mul_pd(Pr,dFdPr),cmp);
            __m256d N = _mm256_mul_pd(_mm256_mul_pd(F,K0),invOnePlusPr);
            __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(F,invOnePlusPr),dKdT);
            dKfdT = _mm256_fmadd_pd(K,_mm256_mul_pd(_mm256_mul_pd(F,_mm256_mul_pd(invOnePlusPr,invOnePlusPr)),dPrdT),dKfdT);
            dKfdT = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),dFdT),MM,dKfdT); 
            _mm256_storeu_pd(&this->dKfdT_[j0],dKfdT);           
            __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),KK),_mm256_add_pd(N1,N2));
            _mm256_storeu_pd(&this->dKfdC_[m0],dKfdC);     
            __m256d KF = _mm256_blendv_pd(N,_mm256_mul_pd(N,M),cmp);
            _mm256_storeu_pd(&this->Kf_[j0],KF);
        }
        if(remain_Troe==1)
        {
            this->Troe_Jac_1();
        }
        else if(remain_Troe==2)
        {

            this->Troe_Jac_2();
        }
        else if(remain_Troe==3)
        {
            this->Troe_Jac_3();
        }
    }        


    {
        for (unsigned int i = 0;i<this->SRI.size();i++)
        {
            const unsigned int j = this->SRI[i];
            const unsigned int k = j - this->Ikf[4];
            const unsigned int m = j - this->Ikf[4] + this->Itbr[2];
            const double Kinf = this->Kf_[j+this->offset_kinf];
            const double invKinf = 1.0/Kinf;
            const double K0 = this->Kf_[j];
            const double dKinfdT = this->dKfdT_[j+this->offset_kinf];
            double F ;
            double dFdT;
            double dFdPr;
            double M = tmp_M[m];
            const double Pr = K0*M*invKinf; 
            this->SRI_F_dFdT_dFdPr(Temperature,Pr,i,F,dFdT,dFdPr);
            const double dK0dT =  this->dKfdT_[j]; 
            const double invOnePlusPr = 1.0/(1.0+Pr);
            const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;
            const double dKdT   = k<this->n_Fall_Off_Reaction?Pr*dKinfdT:dK0dT;
            const double K      = k<this->n_Fall_Off_Reaction?Kinf      :K0;
            const double MM     = k<this->n_Fall_Off_Reaction?M         :1;
            const double KK     = k<this->n_Fall_Off_Reaction?1         :K0*invKinf;
            const double N1     = k<this->n_Fall_Off_Reaction?F*invOnePlusPr  :-F*invOnePlusPr;
            const double N2     = k<this->n_Fall_Off_Reaction?Pr*dFdPr  :dFdPr;
            const double N  = invOnePlusPr*F*K0;
            this->dKfdT_[j] = F*invOnePlusPr*dKdT 
            + F*invOnePlusPr*invOnePlusPr*dPrdT*K 
            + K0*invOnePlusPr*dFdT*MM;
            this->dKfdC_[m] =  K0*invOnePlusPr*KK*(N1 + N2); 
            this->Kf_[j] = k<this->n_Fall_Off_Reaction ? M*N : N;   
        }
    }

    for(unsigned int z = 0; z < this->Ikf[7];z++)
    {
        if(this->isGlobal[z]==1)
        {
            this->JFGNI(z,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);            
            continue;
        }        
        const unsigned int i = z ;
        const auto j = lhsOffset[i+1]-lhsOffset[i];
        const auto k = rhsOffset[i+1]-rhsOffset[i];
        if(j==2)
        {
            if(k==2)        {this->JF22(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else if(k==1)   {this->JF21(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else if(k==3)   {this->JF23(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
        }
        else if(j==1)
        {
            if(k==2)        {this->JF12(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else if(k==1)   {this->JF11(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else if(k==3)   {this->JF13(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
        }
        else if(j==3)
        {
            if(k==2)        {this->JF32(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else if(k==1)   {this->JF31(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else if(k==3)   {this->JF33(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
        }
        if(j>3 || k>3){this->JFGI(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
    }
}
void 
OptReaction::ddNdtByVdcTp
(
    double p,
    double Temperature,
    double* __restrict__ Phi,
    double* __restrict__ c,
    double* __restrict__ dNdtByV,
    double* __restrict__ dBdT,
    double* __restrict__ ddNdtByVdcT
) const noexcept
{
    //Temperature = Temperature<TlowMin?TlowMin:Temperature;
    //Temperature = Temperature>ThighMax?ThighMax:Temperature;
    //this->logT = std::log(Temperature);
    this->T = Temperature;
    //this->invT = 1/Temperature;


    this->update_Pow_pByRT_SumVki(Temperature);
    this->update_Pow_pByRT_SumVki2(Temperature);

    for(unsigned int i=0; i<this->nSpecies; i++)
    {
        this->negGstdByRT[i] = this->tmp_Exp[i];
    }
    for(size_t i = 0; i <this->Troe.size();i++)
    {
        size_t j0 = i + this->nSpecies;
        size_t j1 = i + this->nSpecies + this->Troe.size();
        size_t j2 = i + this->nSpecies + this->Troe.size()*2;         
        this->tmp_Exp[j0] = -Temperature*this->invTsss_[i];            
        this->tmp_Exp[j1] = -this->Tss_[i]*invT; 
        this->tmp_Exp[j2] = -Temperature*this->invTs_[i];            
    }
    
    for(size_t i = 0; i <this->SRI.size();i++)
    {
        size_t j0 = i + this->nSpecies + this->Troe.size()*3;
        size_t j1 = i + this->nSpecies + this->Troe.size()*3 + this->SRI.size();
        this->tmp_Exp[j0] = -this->b_[i]*invT;
        this->tmp_Exp[j1] = -Temperature*this->invc_[i];            
    }   



    {
        size_t remain = this->tmp_ExpSize%4;
        for(size_t i = 0; i < this->tmp_ExpSize-remain;i=i+4)
        {
            __m256d tmp = _mm256_loadu_pd(&this->tmp_Exp[i]);
            tmp = vec256_expd(tmp);
            _mm256_storeu_pd(&this->tmp_Exp[i],tmp);
        }
        if(remain==1)
        {
            size_t i = this->tmp_ExpSize-1;
            this->tmp_Exp[i] = std::exp(this->tmp_Exp[i]);
        }
        else if(remain==2)
        {
            size_t i0 = this->tmp_ExpSize-2;
            size_t i1 = this->tmp_ExpSize-1;
            __m256d tmp = _mm256_setr_pd(tmp_Exp[i0],tmp_Exp[i1],0,0);
            tmp = vec256_expd(tmp);
            this->tmp_Exp[i0] = get_elem0(tmp);
            this->tmp_Exp[i1] = get_elem1(tmp);
        }
        else if(remain==3)
        {
            size_t i0 = this->tmp_ExpSize-3;
            size_t i1 = this->tmp_ExpSize-2;
            size_t i2 = this->tmp_ExpSize-1;
            __m256d tmp = _mm256_setr_pd(tmp_Exp[i0],tmp_Exp[i1],tmp_Exp[i2],0);
            tmp = vec256_expd(tmp);
            this->tmp_Exp[i0] = get_elem0(tmp);
            this->tmp_Exp[i1] = get_elem1(tmp);
            this->tmp_Exp[i2] = get_elem2(tmp);
        }
    }

    if(this->n_PlogReaction>0)
    {
        this->logP = std::log(p);
        for(unsigned int i = 0; i< this->n_PlogReaction; i ++)
        {
            const size_t length = this->Prange[i].size();
            if(p<=this->Prange[i][0])
            {
                double A0 = this->APlog[i][0];
                double beta0 = this->betaPlog[i][0];
                double Ta0 = this->TaPlog[i][0];

                A[i+this->Ikf[6]] = A0;
                A[i+this->Ikf[11]] = A0;
                beta[i+this->Ikf[6]] = beta0;
                beta[i+this->Ikf[11]] = beta0;
                Ta[i+this->Ikf[6]] = Ta0;
                Ta[i+this->Ikf[11]] = Ta0;
                this->Pindex[i] = 0;
            }
            else if(p>=this->Prange[i][length-1])
            {
                double A1 = this->APlog[i][length-1];
                double beta1 = this->betaPlog[i][length-1];
                double Ta1 = this->TaPlog[i][length-1];
                
                A[i+this->Ikf[6]] = A1;
                A[i+this->Ikf[11]] = A1;
                beta[i+this->Ikf[6]] = beta1;
                beta[i+this->Ikf[11]] = beta1;
                Ta[i+this->Ikf[6]] = Ta1;
                Ta[i+this->Ikf[11]] = Ta1;
                this->Pindex[i] = static_cast<unsigned int>(length-1);
            }
            else
            {
                unsigned int index = 0;
                for(unsigned int j = 0; j < length-1;j++)
                {
                    if(this->Prange[i][j]<=p && p<this->Prange[i][j+1])
                    {
                        index = j;
                        break;
                    }
                }
                this->Pindex[i] = index;
                double A0 = this->APlog[i][index+0];
                double A1 = this->APlog[i][index+1];
                double beta0 = this->betaPlog[i][index+0];
                double beta1 = this->betaPlog[i][index+1];
                double Ta0 = this->TaPlog[i][index+0];
                double Ta1 = this->TaPlog[i][index+1];
                A[i+this->Ikf[6]] = A0;
                A[i+this->Ikf[11]] = A1;
                beta[i+this->Ikf[6]] = beta0;
                beta[i+this->Ikf[11]] = beta1;
                Ta[i+this->Ikf[6]] = Ta0;
                Ta[i+this->Ikf[11]] = Ta1;
            }
        }
    }

    {
        __m256d LogTv = _mm256_set1_pd(logT);
        __m256d InvTv = _mm256_set1_pd(invT);        
        unsigned int remain = (this->Ikf[12]-this->n_Temperature_Independent_Reaction)%4;
        unsigned int times = (this->Ikf[12]-this->n_Temperature_Independent_Reaction)/4;
        for(unsigned int z = 0; z <times;z=z+1)
        {
            unsigned int i = z*4 + this->n_Temperature_Independent_Reaction;
            __m256d betav = _mm256_loadu_pd(&this->beta[i]);
            __m256d Tav = _mm256_loadu_pd(&this->Ta[i]);

            __m256d Kfv = _mm256_mul_pd(Tav,-InvTv);
            Kfv = _mm256_fmadd_pd(betav,LogTv,Kfv);
            __m256d A_ = _mm256_loadu_pd(&this->A[i]);
            Kfv = vec256_expd(Kfv);
            Kfv = _mm256_mul_pd(A_,Kfv);
            _mm256_storeu_pd(&this->Kf_[i],Kfv);
            __m256d dKfdT = _mm256_mul_pd(_mm256_fmadd_pd(Tav,InvTv,betav),InvTv);
            dKfdT = _mm256_mul_pd(dKfdT,Kfv); 
            _mm256_storeu_pd(&this->dKfdT_[i+0],dKfdT);           
        }
        if(remain==1)
        {
            unsigned int i = this->Ikf[12]-1;
            this->Kf_[i] = this->A[i]*std::exp(this->beta[i+0]*logT-this->Ta[i+0]*invT);   
            this->dKfdT_[i+0] = this->Kf_[i+0]*(this->beta[i+0]+this->Ta[i+0]*invT)*invT;  
        }
        else if(remain==2)
        {
            unsigned int i0 = this->Ikf[12]-2;
            unsigned int i1 = this->Ikf[12]-1;    
            __m256d betav = _mm256_setr_pd(beta[i0],beta[i1],0,0);
            __m256d Av = _mm256_setr_pd(A[i0],A[i1],0,0);
            __m256d Tav = _mm256_setr_pd(Ta[i0],Ta[i1],0,0);
            __m256d tmp = _mm256_fmsub_pd(betav,LogTv,_mm256_mul_pd(Tav,InvTv));
            tmp = vec256_expd(tmp);
            __m256d Kfv = _mm256_mul_pd(Av,tmp);
            this->Kf_[i0] = get_elem0(Kfv);
            this->Kf_[i1] = get_elem1(Kfv);
            tmp = _mm256_fmadd_pd(Tav,InvTv,betav);
            __m256d dKfdTv = _mm256_mul_pd(Kfv,_mm256_mul_pd(tmp,InvTv));
            this->dKfdT_[i0] = get_elem0(dKfdTv);
            this->dKfdT_[i1] = get_elem1(dKfdTv);
        }
        else if(remain==3)
        {
            unsigned int i0 = this->Ikf[12]-3;
            unsigned int i1 = this->Ikf[12]-2;
            unsigned int i2 = this->Ikf[12]-1;    
            __m256d betav = _mm256_setr_pd(beta[i0],beta[i1],beta[i2],0);
            __m256d Av = _mm256_setr_pd(A[i0],A[i1],A[i2],0);
            __m256d Tav = _mm256_setr_pd(Ta[i0],Ta[i1],Ta[i2],0);
            __m256d tmp = _mm256_fmsub_pd(betav,LogTv,_mm256_mul_pd(Tav,InvTv));
            tmp = vec256_expd(tmp);
            __m256d Kfv = _mm256_mul_pd(Av,tmp);
            this->Kf_[i0] = get_elem0(Kfv);
            this->Kf_[i1] = get_elem1(Kfv);
            this->Kf_[i2] = get_elem2(Kfv);
            tmp = _mm256_fmadd_pd(Tav,InvTv,betav);
            __m256d dKfdTv = _mm256_mul_pd(Kfv,_mm256_mul_pd(tmp,InvTv));
            this->dKfdT_[i0] = get_elem0(dKfdTv);  
            this->dKfdT_[i1] = get_elem1(dKfdTv); 
            this->dKfdT_[i2] = get_elem2(dKfdTv); 
        }
    }



    if(this->n_PlogReaction>0)
    {
        for(unsigned int i = 0; i< this->n_PlogReaction; i ++)
        {
            const size_t length = this->Prange[i].size();
            if(this->Pindex[i] == 0 || this->Pindex[i] == length-1)
            {
                continue;
            }
            else
            {
                unsigned int index = this->Pindex[i];
                double weight = (this->logP - this->logPi[i][index])*this->rDeltaP_[i][index];
                double Kf0 = this->Kf_[i+this->Ikf[6]];
                double Kf1 = this->Kf_[i+this->Ikf[11]];
                double Kf = Kf0*std::pow(Kf1/Kf0,weight);
                this->Kf_[i+this->Ikf[6]] = Kf;

                double beta0 = this->beta[i+this->Ikf[6]];
                double beta1 = this->beta[i+this->Ikf[11]];
                double Ta0 = this->Ta[i+this->Ikf[6]];
                double Ta1 = this->Ta[i+this->Ikf[11]];
                double invt = this->invT;

                double dKfdT = Kf*invt*(beta0 + Ta0*invt + (beta1-beta0+(Ta1-Ta0)*invt)*weight);
                this->dKfdT_[i+this->Ikf[6]] = dKfdT;
            }
        }
    }
    {
        
        unsigned int Tremain = (this->Itbr[5])%4;
        for(unsigned int i = 0; i < this->Itbr[5]-Tremain; i=i+4)
        {
            __m256d arrM_0 = _mm256_setzero_pd();
            __m256d arrM_1 = _mm256_setzero_pd();
            __m256d arrM_2 = _mm256_setzero_pd();
            __m256d arrM_3 = _mm256_setzero_pd();
            double* __restrict__ TBF1DRowi0 = &ThirdBodyFactor1D[(i+0)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi1 = &ThirdBodyFactor1D[(i+1)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi2 = &ThirdBodyFactor1D[(i+2)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi3 = &ThirdBodyFactor1D[(i+3)*this->AlignSpecies];
            for(unsigned int j  = 0;j<this->AlignSpecies;j=j+4)
            {
                __m256d Factor0 = _mm256_loadu_pd(&TBF1DRowi0[j+0]);
                __m256d Factor1 = _mm256_loadu_pd(&TBF1DRowi1[j+0]);
                __m256d Factor2 = _mm256_loadu_pd(&TBF1DRowi2[j+0]);
                __m256d Factor3 = _mm256_loadu_pd(&TBF1DRowi3[j+0]);
                __m256d C_ = _mm256_loadu_pd(&c[j+0]);
                
                arrM_0 = _mm256_fmadd_pd(Factor0,C_,arrM_0);
                arrM_1 = _mm256_fmadd_pd(Factor1,C_,arrM_1);
                arrM_2 = _mm256_fmadd_pd(Factor2,C_,arrM_2);
                arrM_3 = _mm256_fmadd_pd(Factor3,C_,arrM_3);
            }
            __m256d s0h = _mm256_hadd_pd(arrM_0, arrM_1); 
            __m256d s1h = _mm256_hadd_pd(arrM_2, arrM_3); 
            s0h = _mm256_permute4x64_pd(s0h, 0b11011000);
            s1h = _mm256_permute4x64_pd(s1h, 0b11011000);
            __m256d sum_all = _mm256_hadd_pd(s0h, s1h); 
            sum_all = _mm256_permute4x64_pd(sum_all, 0b11011000);

           _mm256_storeu_pd(&this->tmp_M[i+0],sum_all);

        }
        if(Tremain==3)
        {
            unsigned int i =(this->Itbr[5]) -3;
            double* __restrict__ TBF1DRowi0 = &ThirdBodyFactor1D[(i+0)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi1 = &ThirdBodyFactor1D[(i+1)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi2 = &ThirdBodyFactor1D[(i+2)*this->AlignSpecies];
            double M0 = 0;
            double M1 = 0;           
            double M2 = 0; 
            __m256d arrM_0 = _mm256_setzero_pd();
            __m256d arrM_1 = _mm256_setzero_pd();
            __m256d arrM_2 = _mm256_setzero_pd();
            for(unsigned int j  = 0;j<this->AlignSpecies;j=j+4)
            {
                __m256d Factor0 = _mm256_loadu_pd(&TBF1DRowi0[j+0]);
                __m256d Factor1 = _mm256_loadu_pd(&TBF1DRowi1[j+0]);
                __m256d Factor2 = _mm256_loadu_pd(&TBF1DRowi2[j+0]);
                __m256d C_ = _mm256_loadu_pd(&c[j+0]);
                arrM_0 = _mm256_fmadd_pd(Factor0,C_,arrM_0);
                arrM_1 = _mm256_fmadd_pd(Factor1,C_,arrM_1);
                arrM_2 = _mm256_fmadd_pd(Factor2,C_,arrM_2);
            }

            M0 = M0 + hsum4(arrM_0);
            M1 = M1 + hsum4(arrM_1);
            M2 = M2 + hsum4(arrM_2);

            this->tmp_M[i+0] = M0;
            this->tmp_M[i+1] = M1;
            this->tmp_M[i+2] = M2;
        }
        else if(Tremain==2)
        {
            unsigned int i =(this->Itbr[5]) -2;
            double* __restrict__ TBF1DRowi0 = &ThirdBodyFactor1D[(i+0)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi1 = &ThirdBodyFactor1D[(i+1)*this->AlignSpecies];
            double M0 = 0;
            double M1 = 0;           
            __m256d arrM_0 = _mm256_setzero_pd();
            __m256d arrM_1 = _mm256_setzero_pd();
            for(unsigned int j  = 0;j<this->AlignSpecies;j=j+4)
            {
                __m256d Factor0 = _mm256_loadu_pd(&TBF1DRowi0[j+0]);
                __m256d Factor1 = _mm256_loadu_pd(&TBF1DRowi1[j+0]);
                __m256d C_ = _mm256_loadu_pd(&c[j+0]);
                arrM_0 = _mm256_fmadd_pd(Factor0,C_,arrM_0);
                arrM_1 = _mm256_fmadd_pd(Factor1,C_,arrM_1);
            }

            M0 = M0 + hsum4(arrM_0);
            M1 = M1 + hsum4(arrM_1);

            this->tmp_M[i+0] = M0;
            this->tmp_M[i+1] = M1;
        }
        else if(Tremain==1)
        {
            unsigned int i =(this->Itbr[5]) -1;
            double* __restrict__ TBF1DRowi0 = &ThirdBodyFactor1D[(i+0)*this->AlignSpecies];            
            double M0 = 0;
            __m256d arrM_0 = _mm256_setzero_pd();
            for(unsigned int j  = 0;j<this->AlignSpecies;j=j+4)
            {
                __m256d Factor0 = _mm256_loadu_pd(&TBF1DRowi0[j+0]);
                __m256d C_ = _mm256_loadu_pd(&c[j+0]);
                arrM_0 = _mm256_fmadd_pd(Factor0,C_,arrM_0);
            }
            M0 = M0 + hsum4(arrM_0);
            this->tmp_M[i+0] = M0;
        }
    }

    for(unsigned int i = 0; i < this->n_ThirdBodyReaction; i++)
    {
        double M = this->tmp_M[i+this->Itbr[1]];
        const unsigned int j = i + this->Ikf[3];
        this->dKfdC_[i+this->Itbr[1]] = this->Kf_[j];   
        this->Kf_[j] = this->Kf_[j]*M;
        this->dKfdT_[j] = this->dKfdT_[j]*M;
    }
    


    for(unsigned int i = 0; i < this->n_NonEquilibriumThirdBodyReaction; i++)
    {
        double Mfwd = this->tmp_M[i];
        double Mrev = this->tmp_M[this->Itbr[4]+i];
        this->dKfdC_[i] = this->Kf_[this->Ikf[2]+i];
        this->dKfdC_[this->Itbr[4]+i] = this->Kf_[this->Ikf[10]+i];
        this->Kf_[this->Ikf[2]+i] = this->Kf_[this->Ikf[2]+i]*Mfwd;
        this->dKfdT_[this->Ikf[2]+i] = this->dKfdT_[this->Ikf[2]+i]*Mfwd;
        this->Kf_[this->Ikf[10]+i] = this->Kf_[this->Ikf[10]+i]*Mrev;
        this->dKfdT_[this->Ikf[10]+i] = this->dKfdT_[this->Ikf[10]+i]*Mrev;
    } 


    {
        size_t remain_Lindemann = (Lindemann.size())%4;
        for (size_t i = 0;i<Lindemann.size()-remain_Lindemann;i=i+4)
        {
            const unsigned int j0 = this->Lindemann[i+0]+0;
            const unsigned int j1 = this->Lindemann[i+0]+1;
            const unsigned int j2 = this->Lindemann[i+0]+2;
            const unsigned int j3 = this->Lindemann[i+0]+3;
            const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
            const unsigned int k0 = j0 - this->Ikf[4];
            const unsigned int k1 = j1 - this->Ikf[4];
            const unsigned int k2 = j2 - this->Ikf[4];
            const unsigned int k3 = j3 - this->Ikf[4];
            __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0 + this->offset_kinf]);    
            __m256d one = _mm256_set1_pd(1.0);       
            __m256d invKinf = _mm256_div_pd(_mm256_set1_pd(1.0),Kinf);
            __m256d dKinfdT = _mm256_loadu_pd(&this->dKfdT_[j0 + this->offset_kinf]);
            __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);        
            __m256d dK0dT = _mm256_loadu_pd(&this->dKfdT_[j0]);     
            __m256d M = _mm256_loadu_pd(&tmp_M[m0]);
            __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(K0,M),invKinf);
            __m256d dPrdT = _mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT));
            dPrdT = _mm256_mul_pd(dPrdT,invKinf);
            __m256d k = _mm256_setr_pd(k0,k1,k2,k3);
            __m256d cmp = _mm256_cmp_pd(k,_mm256_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
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
            __m256d one = _mm256_set1_pd(1.0);       
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
            __m256d one = _mm256_set1_pd(1.0);       
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


    {
        size_t remain_Troe = (Troe.size())%4;
        for (size_t i = 0;i<Troe.size()-remain_Troe;i=i+4)
        {
            const unsigned int j0 = this->Troe[i+0];
            const unsigned int j1 = this->Troe[i+1];
            const unsigned int j2 = this->Troe[i+2];
            const unsigned int j3 = this->Troe[i+3];
            const unsigned int k0 = j0 - this->Ikf[4];
            const unsigned int k1 = j1 - this->Ikf[4];
            const unsigned int k2 = j2 - this->Ikf[4];
            const unsigned int k3 = j3 - this->Ikf[4];
            const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
            __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
            __m256d invKinf = _mm256_div_pd(_mm256_set1_pd(1.0),Kinf);
            __m256d dKinfdT = _mm256_loadu_pd(&this->dKfdT_[j0+this->offset_kinf]);
            __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);           
            __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
            __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(M,K0),invKinf);
            __m256d small = _mm256_set1_pd(2.2e-16);
            __m256d cmp_result_Pr = _mm256_cmp_pd(Pr,small,_CMP_GE_OQ);
            Pr = _mm256_add_pd(Pr,_mm256_set1_pd(1e-100));
            const double invLog10 = 1.0/std::log(10);
            __m256d logPr_ = _mm256_mul_pd(vec256_logd(_mm256_max_pd(small,Pr)),_mm256_set1_pd(invLog10));
            __m256d InvTsss = _mm256_loadu_pd(&this->invTsss_[i]);
            __m256d InvTs = _mm256_loadu_pd(&this->invTs_[i]);
            __m256d Tss = _mm256_loadu_pd(&this->Tss_[i]);
            __m256d expTTsss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies]);
            __m256d expTTss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->Troe.size()]);
            __m256d expTTs = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->Troe.size()*2]);
            __m256d one = _mm256_set1_pd(1.0);
            __m256d alpha = _mm256_loadu_pd(&this->alpha_[i]);
            __m256d Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha),expTTsss);
            Fcent = _mm256_fmadd_pd(alpha,expTTs,Fcent);
            Fcent = _mm256_add_pd(expTTss,Fcent);
            __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),_mm256_set1_pd(invLog10));
            __m256d cc = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(0.67),_mm256_set1_pd(0.4));
            __m256d n = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(-1.27),_mm256_set1_pd(0.75));
            __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(cc,logPr_),_mm256_set1_pd(0.14),n);
            __m256d invx1 = _mm256_div_pd(one,x1);
            __m256d x2 = _mm256_mul_pd(_mm256_sub_pd(logPr_,cc),invx1);
            __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
            __m256d invx3 = _mm256_div_pd(one,x3);
            __m256d x4 = _mm256_mul_pd(logFcent,invx3);
            __m256d  F = vec256_powd(_mm256_set1_pd(10),x4);
            __m256d logTen = _mm256_set1_pd(std::log(10));
            __m256d dFcentdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_sub_pd(alpha,one),InvTsss),expTTsss);
            dFcentdT = _mm256_sub_pd(dFcentdT,_mm256_mul_pd(_mm256_mul_pd(alpha,InvTs),expTTs));
            __m256d invT2 = _mm256_set1_pd(invT*invT);
            dFcentdT = _mm256_fmadd_pd(expTTss,_mm256_mul_pd(Tss,invT2),dFcentdT);
            __m256d cmp2 = _mm256_cmp_pd(Fcent,small,_CMP_GE_OQ);
            __m256d dlogFcentdT = _mm256_div_pd(_mm256_div_pd(dFcentdT,_mm256_max_pd(Fcent,small)),logTen);
            dlogFcentdT = _mm256_blendv_pd(_mm256_setzero_pd(), dlogFcentdT, cmp2);
            __m256d dcdT = _mm256_mul_pd(dlogFcentdT,_mm256_set1_pd(-0.67));
            __m256d dndT = _mm256_mul_pd(dlogFcentdT,_mm256_set1_pd(-1.27));
            __m256d dx1dT = _mm256_fmadd_pd(dcdT,_mm256_set1_pd(-0.14),dndT);
            __m256d dx2dT = _mm256_mul_pd(_mm256_sub_pd(dcdT,_mm256_mul_pd(x2,dx1dT)),invx1);
            __m256d dx3dT = _mm256_mul_pd(_mm256_mul_pd(x2,dx2dT),_mm256_set1_pd(2.0));
            __m256d dx4dT = _mm256_mul_pd(_mm256_sub_pd(dlogFcentdT,_mm256_mul_pd(x4,dx3dT)),invx3);
            __m256d dFdT = _mm256_mul_pd(logTen,_mm256_mul_pd(F,dx4dT));
            __m256d dlogPrdPr = _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(Pr,logTen));
            dlogPrdPr = _mm256_blendv_pd(_mm256_setzero_pd(), dlogPrdPr, cmp_result_Pr);
            __m256d dx1dPr = _mm256_mul_pd(dlogPrdPr,_mm256_set1_pd(-0.14));
            __m256d dx2dPr = _mm256_mul_pd(_mm256_sub_pd(dlogPrdPr,_mm256_mul_pd(x2,dx1dPr)),invx1);
            __m256d dx3dPr = _mm256_mul_pd(_mm256_mul_pd(x2,dx2dPr),_mm256_set1_pd(2.0));
            __m256d dx4dPr = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(-1.0),x4),dx3dPr),invx3);
            __m256d dFdPr  = _mm256_mul_pd(_mm256_mul_pd(logTen,F),dx4dPr);    
            __m256d dK0dT = _mm256_loadu_pd(&this->dKfdT_[j0]);            
            __m256d dPrdT = _mm256_mul_pd(_mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT)),invKinf);
            dFdT = _mm256_fmadd_pd(dFdPr,dPrdT,dFdT);
            __m256d k = _mm256_setr_pd(k0,k1,k2,k3);
            __m256d cmp = _mm256_cmp_pd(k,_mm256_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
            __m256d dKdT = _mm256_blendv_pd(dK0dT, _mm256_mul_pd(Pr,dKinfdT), cmp);
            __m256d K = _mm256_blendv_pd(K0, Kinf, cmp);
            __m256d MM = _mm256_blendv_pd(_mm256_set1_pd(1), M, cmp);
            __m256d KK = _mm256_blendv_pd(_mm256_mul_pd(K0,invKinf), _mm256_set1_pd(1), cmp);
            __m256d invOnePlusPr = _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),Pr));
            __m256d N1 = _mm256_mul_pd(F,invOnePlusPr);
            N1 = _mm256_blendv_pd(_mm256_sub_pd(_mm256_setzero_pd(),N1),N1,cmp);            
            __m256d N2 = _mm256_blendv_pd(dFdPr,_mm256_mul_pd(Pr,dFdPr),cmp);
            __m256d N = _mm256_mul_pd(_mm256_mul_pd(F,K0),invOnePlusPr);
            __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(F,invOnePlusPr),dKdT);
            dKfdT = _mm256_fmadd_pd(K,_mm256_mul_pd(_mm256_mul_pd(F,_mm256_mul_pd(invOnePlusPr,invOnePlusPr)),dPrdT),dKfdT);
            dKfdT = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),dFdT),MM,dKfdT); 
            _mm256_storeu_pd(&this->dKfdT_[j0],dKfdT);           
            __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),KK),_mm256_add_pd(N1,N2));
            _mm256_storeu_pd(&this->dKfdC_[m0],dKfdC);     
            __m256d KF = _mm256_blendv_pd(N,_mm256_mul_pd(N,M),cmp);
            _mm256_storeu_pd(&this->Kf_[j0],KF);
        }
        if(remain_Troe==1)
        {
            this->Troe_Jac_1();
        }
        else if(remain_Troe==2)
        {

            this->Troe_Jac_2();
        }
        else if(remain_Troe==3)
        {
            this->Troe_Jac_3();
        }
    }        


    {
        for (unsigned int i = 0;i<this->SRI.size();i++)
        {
            const unsigned int j = this->SRI[i];
            const unsigned int k = j - this->Ikf[4];
            const unsigned int m = j - this->Ikf[4] + this->Itbr[2];
            const double Kinf = this->Kf_[j+this->offset_kinf];
            const double invKinf = 1.0/Kinf;
            const double K0 = this->Kf_[j];
            const double dKinfdT = this->dKfdT_[j+this->offset_kinf];
            double F ;
            double dFdT;
            double dFdPr;
            double M = tmp_M[m];
            const double Pr = K0*M*invKinf; 
            this->SRI_F_dFdT_dFdPr(Temperature,Pr,i,F,dFdT,dFdPr);
            const double dK0dT =  this->dKfdT_[j]; 
            const double invOnePlusPr = 1.0/(1.0+Pr);
            const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;
            const double dKdT   = k<this->n_Fall_Off_Reaction?Pr*dKinfdT:dK0dT;
            const double K      = k<this->n_Fall_Off_Reaction?Kinf      :K0;
            const double MM     = k<this->n_Fall_Off_Reaction?M         :1;
            const double KK     = k<this->n_Fall_Off_Reaction?1         :K0*invKinf;
            const double N1     = k<this->n_Fall_Off_Reaction?F*invOnePlusPr  :-F*invOnePlusPr;
            const double N2     = k<this->n_Fall_Off_Reaction?Pr*dFdPr  :dFdPr;
            const double N  = invOnePlusPr*F*K0;
            this->dKfdT_[j] = F*invOnePlusPr*dKdT 
            + F*invOnePlusPr*invOnePlusPr*dPrdT*K 
            + K0*invOnePlusPr*dFdT*MM;
            this->dKfdC_[m] =  K0*invOnePlusPr*KK*(N1 + N2); 
            this->Kf_[j] = k<this->n_Fall_Off_Reaction ? M*N : N;   
        }
    }

    for(unsigned int z = 0; z < this->Ikf[7];z++)
    {
        if(this->isGlobal[z]==1)
        {
            this->JFGNI(z,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);            
            continue;
        }        
        const unsigned int i = z ;
        const auto j = lhsOffset[i+1]-lhsOffset[i];
        const auto k = rhsOffset[i+1]-rhsOffset[i];
        if(j==2)
        {
            if(k==2)        {this->JF22(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else if(k==1)   {this->JF21(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else if(k==3)   {this->JF23(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
        }
        else if(j==1)
        {
            if(k==2)        {this->JF12(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else if(k==1)   {this->JF11(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else if(k==3)   {this->JF13(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
        }
        else if(j==3)
        {
            if(k==2)        {this->JF32(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else if(k==1)   {this->JF31(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else if(k==3)   {this->JF33(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
        }
        if(j>3 || k>3){this->JFGI(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
    }
}

void 
OptReaction::ddYdtdY_Vec1_0
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ rhoMByRhoi,
    const double* __restrict__ WiByrhoM,  
    const double* __restrict__ dPhidt,  
    const double* __restrict__ Phi, 
    double* __restrict__ Jac
) const noexcept
{
    double* __restrict__ invWPtr = &invW[0];
    __m256d rhoMv = _mm256_set1_pd(rhoM);
    for(unsigned int j = 0; j < this->nSpecies; j = j + 4)
    {
        __m256d rhoMvj_ = _mm256_loadu_pd(&rhoMByRhoi[j+0]);
        for(unsigned int i = 0; i < this->nSpecies; i=i+4)
        {
            __m256d rhoMByWiYTv = _mm256_mul_pd(_mm256_loadu_pd(&invWPtr[i+0]), _mm256_loadu_pd(&Phi[i+0]));
            rhoMByWiYTv = _mm256_mul_pd(-rhoMv,rhoMByWiYTv);
            _mm256_storeu_pd(&buffer[i*4+0+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0x00),rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+4+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0x55),rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+8+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0xAA),rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+12+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0xFF),rhoMvj_));
        }

        buffer[j*4+0] += rhoM*invWPtr[j+0];
        buffer[j*4+5] += rhoM*invWPtr[j+1];
        buffer[j*4+10] += rhoM*invWPtr[j+2];
        buffer[j*4+15] += rhoM*invWPtr[j+3];

        for(unsigned int i=0; i<this->nSpecies; i=i+4)
        {
            const double Wi0ByrhoM_ = WiByrhoM[i+0];
            const double Wi1ByrhoM_ = WiByrhoM[i+1];
            const double Wi2ByrhoM_ = WiByrhoM[i+2];
            const double Wi3ByrhoM_ = WiByrhoM[i+3];
            const double dYi0dt = dPhidt[i+0]*Wi0ByrhoM_;
            const double dYi1dt = dPhidt[i+1]*Wi1ByrhoM_;
            const double dYi2dt = dPhidt[i+2]*Wi2ByrhoM_;
            const double dYi3dt = dPhidt[i+3]*Wi3ByrhoM_;
            __m256d ddNi0dtByVdYj = _mm256_setzero_pd();
            __m256d ddNi1dtByVdYj = _mm256_setzero_pd();
            __m256d ddNi2dtByVdYj = _mm256_setzero_pd();
            __m256d ddNi3dtByVdYj = _mm256_setzero_pd();
            const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*alignN];
            const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*alignN];
            const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*alignN];
            const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i+3)*alignN];
            for (unsigned int k=0; k<this->nSpecies; k=k+4)
            {
                const double ddNi0dtByVdck0 = JcRowi0[k+0];
                const double ddNi1dtByVdck0 = JcRowi1[k+0];
                const double ddNi2dtByVdck0 = JcRowi2[k+0];
                const double ddNi3dtByVdck0 = JcRowi3[k+0];
                const double ddNi0dtByVdck1 = JcRowi0[k+1];
                const double ddNi1dtByVdck1 = JcRowi1[k+1];
                const double ddNi2dtByVdck1 = JcRowi2[k+1];
                const double ddNi3dtByVdck1 = JcRowi3[k+1];
                const double ddNi0dtByVdck2 = JcRowi0[k+2];
                const double ddNi1dtByVdck2 = JcRowi1[k+2];
                const double ddNi2dtByVdck2 = JcRowi2[k+2];
                const double ddNi3dtByVdck2 = JcRowi3[k+2];
                const double ddNi0dtByVdck3 = JcRowi0[k+3];
                const double ddNi1dtByVdck3 = JcRowi1[k+3];
                const double ddNi2dtByVdck3 = JcRowi2[k+3];
                const double ddNi3dtByVdck3 = JcRowi3[k+3];

                __m256d dCk0dYj = _mm256_loadu_pd(&buffer[k*4]);
                __m256d dCk1dYj = _mm256_loadu_pd(&buffer[k*4+4]);
                __m256d dCk2dYj = _mm256_loadu_pd(&buffer[k*4+8]);
                __m256d dCk3dYj = _mm256_loadu_pd(&buffer[k*4+12]);

                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck0),dCk0dYj,ddNi0dtByVdYj);
                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck1),dCk1dYj,ddNi0dtByVdYj);
                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck2),dCk2dYj,ddNi0dtByVdYj);
                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck3),dCk3dYj,ddNi0dtByVdYj);

                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck0),dCk0dYj,ddNi1dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck1),dCk1dYj,ddNi1dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck2),dCk2dYj,ddNi1dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck3),dCk3dYj,ddNi1dtByVdYj);


                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck0),dCk0dYj,ddNi2dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck1),dCk1dYj,ddNi2dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck2),dCk2dYj,ddNi2dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck3),dCk3dYj,ddNi2dtByVdYj);

                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck0),dCk0dYj,ddNi3dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck1),dCk1dYj,ddNi3dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck2),dCk2dYj,ddNi3dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck3),dCk3dYj,ddNi3dtByVdYj);

            }
             __m256d WiByrhoM_0 = _mm256_set1_pd(Wi0ByrhoM_);
            __m256d dYidtv0 = _mm256_set1_pd(dYi0dt);
            __m256d ddYi0dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv0);
            ddYi0dtdYj = _mm256_fmadd_pd(WiByrhoM_0,ddNi0dtByVdYj,ddYi0dtdYj);
            __m256d WiByrhoM_1 = _mm256_set1_pd(Wi1ByrhoM_);
            __m256d dYidtv1 = _mm256_set1_pd(dYi1dt);
            __m256d ddYi1dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv1);
            ddYi1dtdYj = _mm256_fmadd_pd(WiByrhoM_1,ddNi1dtByVdYj,ddYi1dtdYj);
            __m256d WiByrhoM_2 = _mm256_set1_pd(Wi2ByrhoM_);
            __m256d dYidtv2 = _mm256_set1_pd(dYi2dt);
            __m256d ddYi2dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv2);
            ddYi2dtdYj = _mm256_fmadd_pd(WiByrhoM_2,ddNi2dtByVdYj,ddYi2dtdYj);
            __m256d WiByrhoM_3 = _mm256_set1_pd(Wi3ByrhoM_);
            __m256d dYidtv3 = _mm256_set1_pd(dYi3dt);
            __m256d ddYi3dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv3);
            ddYi3dtdYj = _mm256_fmadd_pd(WiByrhoM_3,ddNi3dtByVdYj,ddYi3dtdYj);
            _mm256_storeu_pd(&Jac[(i+0)*(alignN)+ j+0],ddYi0dtdYj);
            _mm256_storeu_pd(&Jac[(i+1)*(alignN)+ j+0],ddYi1dtdYj);
            _mm256_storeu_pd(&Jac[(i+2)*(alignN)+ j+0],ddYi2dtdYj);
            _mm256_storeu_pd(&Jac[(i+3)*(alignN)+ j+0],ddYi3dtdYj);
        }     
    }     
}


void 
OptReaction::ddYdtdY_Vec1_1
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ rhoMByRhoi,
    const double* __restrict__ WiByrhoM,  
    const double* __restrict__ dPhidt,  
    const double* __restrict__ Phi, 
    double* __restrict__ Jac
) const noexcept
{
    double* __restrict__ invWPtr = &invW[0];
    __m256d rhoMv = _mm256_set1_pd(rhoM);
    for(unsigned int j = 0; j < this->nSpecies - 1; j = j + 4)
    {
        __m256d rhoMvj_ = _mm256_loadu_pd(&rhoMByRhoi[j+0]);

        for(unsigned int i = 0; i < this->nSpecies-1; i=i+4)
        {
            __m256d rhoMByWiYTv = _mm256_mul_pd(_mm256_loadu_pd(&invWPtr[i+0]), _mm256_loadu_pd(&Phi[i+0]));
            rhoMByWiYTv = _mm256_mul_pd(-rhoMv,rhoMByWiYTv);
            _mm256_storeu_pd(&buffer[i*4+0+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0x00),rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+4+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0x55),rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+8+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0xAA),rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+12+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0xFF),rhoMvj_));           
        }

        unsigned int ii = this->nSpecies-1;
        const double rhoMByWiYTPi0 = -rhoM*invWPtr[ii+0]*Phi[ii+0];
        __m256d result1 = _mm256_mul_pd(_mm256_set1_pd(rhoMByWiYTPi0),rhoMvj_);
        _mm256_storeu_pd(&buffer[ii*4],result1);
        buffer[j*4+0] += rhoM*invWPtr[j+0];
        buffer[j*4+5] += rhoM*invWPtr[j+1];
        buffer[j*4+10] += rhoM*invWPtr[j+2];
        buffer[j*4+15] += rhoM*invWPtr[j+3];
  
        for(unsigned int i=0; i<this->nSpecies - 1; i=i+4)
        {
       
            const double Wi0ByrhoM_ = WiByrhoM[i+0];
            const double Wi1ByrhoM_ = WiByrhoM[i+1];
            const double Wi2ByrhoM_ = WiByrhoM[i+2];
            const double Wi3ByrhoM_ = WiByrhoM[i+3];

            const double dYi0dt = dPhidt[i+0]*Wi0ByrhoM_;
            const double dYi1dt = dPhidt[i+1]*Wi1ByrhoM_;
            const double dYi2dt = dPhidt[i+2]*Wi2ByrhoM_;
            const double dYi3dt = dPhidt[i+3]*Wi3ByrhoM_;

            __m256d ddNi0dtByVdYj = _mm256_setzero_pd();
            __m256d ddNi1dtByVdYj = _mm256_setzero_pd();
            __m256d ddNi2dtByVdYj = _mm256_setzero_pd();
            __m256d ddNi3dtByVdYj = _mm256_setzero_pd();

            const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*(alignN)];
            const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*(alignN)];
            const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*(alignN)];
            const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i+3)*(alignN)];
            for (unsigned int k=0; k<this->nSpecies-1; k=k+4)
            {
                const double ddNi0dtByVdck0 = JcRowi0[k+0];
                const double ddNi0dtByVdck1 = JcRowi0[k+1];
                const double ddNi0dtByVdck2 = JcRowi0[k+2];
                const double ddNi0dtByVdck3 = JcRowi0[k+3];

                const double ddNi1dtByVdck0 = JcRowi1[k+0];
                const double ddNi1dtByVdck1 = JcRowi1[k+1];
                const double ddNi1dtByVdck2 = JcRowi1[k+2];
                const double ddNi1dtByVdck3 = JcRowi1[k+3];

                const double ddNi2dtByVdck0 = JcRowi2[k+0];
                const double ddNi2dtByVdck1 = JcRowi2[k+1];
                const double ddNi2dtByVdck2 = JcRowi2[k+2];
                const double ddNi2dtByVdck3 = JcRowi2[k+3];

                const double ddNi3dtByVdck0 = JcRowi3[k+0];
                const double ddNi3dtByVdck1 = JcRowi3[k+1];
                const double ddNi3dtByVdck2 = JcRowi3[k+2];
                const double ddNi3dtByVdck3 = JcRowi3[k+3];

                __m256d dCk0dYj = _mm256_loadu_pd(&buffer[k*4]);
                __m256d dCk1dYj = _mm256_loadu_pd(&buffer[k*4+4]);
                __m256d dCk2dYj = _mm256_loadu_pd(&buffer[k*4+8]);
                __m256d dCk3dYj = _mm256_loadu_pd(&buffer[k*4+12]);

                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck0),dCk0dYj,ddNi0dtByVdYj);
                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck1),dCk1dYj,ddNi0dtByVdYj);
                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck2),dCk2dYj,ddNi0dtByVdYj);
                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck3),dCk3dYj,ddNi0dtByVdYj);

                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck0),dCk0dYj,ddNi1dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck1),dCk1dYj,ddNi1dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck2),dCk2dYj,ddNi1dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck3),dCk3dYj,ddNi1dtByVdYj);

                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck0),dCk0dYj,ddNi2dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck1),dCk1dYj,ddNi2dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck2),dCk2dYj,ddNi2dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck3),dCk3dYj,ddNi2dtByVdYj);

                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck0),dCk0dYj,ddNi3dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck1),dCk1dYj,ddNi3dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck2),dCk2dYj,ddNi3dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck3),dCk3dYj,ddNi3dtByVdYj);
            }
            unsigned int k = this->nSpecies-1;
            const double ddNi0dtByVdck0 = JcRowi0[k+0];
            const double ddNi1dtByVdck0 = JcRowi1[k+0];
            const double ddNi2dtByVdck0 = JcRowi2[k+0];
            const double ddNi3dtByVdck0 = JcRowi3[k+0];
            __m256d dCk0dYj = _mm256_loadu_pd(&buffer[k*4]);
            ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck0),dCk0dYj,ddNi0dtByVdYj);
            ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck0),dCk0dYj,ddNi1dtByVdYj);
            ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck0),dCk0dYj,ddNi2dtByVdYj);
            ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck0),dCk0dYj,ddNi3dtByVdYj);


            __m256d WiByrhoM_0 = _mm256_set1_pd(Wi0ByrhoM_);
            __m256d WiByrhoM_1 = _mm256_set1_pd(Wi1ByrhoM_);
            __m256d WiByrhoM_2 = _mm256_set1_pd(Wi2ByrhoM_);
            __m256d WiByrhoM_3 = _mm256_set1_pd(Wi3ByrhoM_);

            __m256d dYidtv0 = _mm256_set1_pd(dYi0dt);
            __m256d dYidtv1 = _mm256_set1_pd(dYi1dt);
            __m256d dYidtv2 = _mm256_set1_pd(dYi2dt);
            __m256d dYidtv3 = _mm256_set1_pd(dYi3dt);

            __m256d ddYi0dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv0);
            ddYi0dtdYj = _mm256_fmadd_pd(WiByrhoM_0,ddNi0dtByVdYj,ddYi0dtdYj);

            __m256d ddYi1dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv1);
            ddYi1dtdYj = _mm256_fmadd_pd(WiByrhoM_1,ddNi1dtByVdYj,ddYi1dtdYj);
            
            __m256d ddYi2dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv2);
            ddYi2dtdYj = _mm256_fmadd_pd(WiByrhoM_2,ddNi2dtByVdYj,ddYi2dtdYj);
            
            __m256d ddYi3dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv3);
            ddYi3dtdYj = _mm256_fmadd_pd(WiByrhoM_3,ddNi3dtByVdYj,ddYi3dtdYj);  


            _mm256_storeu_pd(&Jac[(i+0)*(alignN)+ j+0],ddYi0dtdYj);
            _mm256_storeu_pd(&Jac[(i+1)*(alignN)+ j+0],ddYi1dtdYj);
            _mm256_storeu_pd(&Jac[(i+2)*(alignN)+ j+0],ddYi2dtdYj);
            _mm256_storeu_pd(&Jac[(i+3)*(alignN)+ j+0],ddYi3dtdYj);
        }     

        const double* __restrict__ JcRowii = &ddNdtByVdcT[ii*alignN];
        const double WiByrhoM_ = WiByrhoM[ii];
        double dYidt = dPhidt[ii];
        dYidt *= WiByrhoM_;
        __m256d ddNidtByVdYjv = _mm256_setzero_pd();
        for (unsigned int k=0; k<this->nSpecies-1; k=k+4)
        {
            __m256d Arrk0 = _mm256_loadu_pd(&buffer[k*4]);
            ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(JcRowii[k+0]),Arrk0,ddNidtByVdYjv);
            __m256d Arrk1 = _mm256_loadu_pd(&buffer[k*4+4]);
            ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(JcRowii[k+1]),Arrk1,ddNidtByVdYjv);
            __m256d Arrk2 = _mm256_loadu_pd(&buffer[k*4+8]);
            ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(JcRowii[k+2]),Arrk2,ddNidtByVdYjv);
            __m256d Arrk3 = _mm256_loadu_pd(&buffer[k*4+12]);
            ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(JcRowii[k+3]),Arrk3,ddNidtByVdYjv);  
        }
        unsigned int k0 = nSpecies-1;
        __m256d Arrk3 = _mm256_loadu_pd(&buffer[k0*4]);
        ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(JcRowii[k0]),Arrk3,ddNidtByVdYjv);          
        __m256d result = _mm256_mul_pd(rhoMvj_,_mm256_set1_pd(dYidt));
        result = _mm256_fmadd_pd(ddNidtByVdYjv,_mm256_set1_pd(WiByrhoM_),result);
        _mm256_storeu_pd(&Jac[ii*(alignN) + j+0],result);
    }
    unsigned int j = this->nSpecies-1;
    const double rhoMvj_ = rhoMByRhoi[j];    
    __m256d rhoMvjv = _mm256_set1_pd(-rhoMvj_);


    for(unsigned int i=0; i<this->nSpecies -1; i=i+4)
    {
        __m256d rhoMByWiv = _mm256_mul_pd(rhoMv,_mm256_loadu_pd(&invWPtr[i]));
        __m256d YTpv = _mm256_loadu_pd(&Phi[i+0]);
        __m256d result = _mm256_mul_pd(rhoMByWiv,_mm256_mul_pd(YTpv,rhoMvjv));
        _mm256_storeu_pd(&buffer[0*4+0+i+0],result);
    }
    buffer[j] = rhoM*invWPtr[j]*(1-rhoMByRhoi[j]*Phi[j]);

    for(unsigned int i=0; i<this->nSpecies-1; i=i+4)
    {
        __m256d WiByrhoMv = _mm256_loadu_pd(&WiByrhoM[i+0]);
        __m256d dYidtv = _mm256_mul_pd(_mm256_loadu_pd(&dPhidt[i+0]),WiByrhoMv);
        __m256d sum0 = _mm256_setzero_pd();
        __m256d sum1 = _mm256_setzero_pd();
        __m256d sum2 = _mm256_setzero_pd();
        __m256d sum3 = _mm256_setzero_pd();
        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*alignN];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*alignN];
        const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*alignN];
        const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i+3)*alignN];

        for (unsigned int k=0; k<this->nSpecies-1; k=k+4)
        {
            __m256d a0 = _mm256_loadu_pd(&JcRowi0[k]);
            __m256d a1 = _mm256_loadu_pd(&JcRowi1[k]);
            __m256d a2 = _mm256_loadu_pd(&JcRowi2[k]);
            __m256d a3 = _mm256_loadu_pd(&JcRowi3[k]);

            __m256d b0 = _mm256_loadu_pd(&buffer[k]);
            sum0 = _mm256_fmadd_pd(a0, b0, sum0); 
            sum1 = _mm256_fmadd_pd(a1, b0, sum1); 
            sum2 = _mm256_fmadd_pd(a2, b0, sum2); 
            sum3 = _mm256_fmadd_pd(a3, b0, sum3); 
        }
        __m256d ddNidtByVdYjv = hsum4x4(sum0,sum1,sum2,sum3);
        unsigned int k = this->nSpecies-1;
        ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_setr_pd(JcRowi0[k],JcRowi1[k],JcRowi2[k],JcRowi3[k]),_mm256_set1_pd(buffer[k]),ddNidtByVdYjv);
        __m256d result = _mm256_mul_pd(WiByrhoMv,ddNidtByVdYjv);
        result = _mm256_fmadd_pd(_mm256_set1_pd(rhoMvj_),dYidtv,result);

        Jac[(i+0)*(alignN) + j] = get_elem0(result);
        Jac[(i+1)*(alignN) + j] = get_elem1(result);
        Jac[(i+2)*(alignN) + j] = get_elem2(result);
        Jac[(i+3)*(alignN) + j] = get_elem3(result);
    }
    {
        unsigned int i = this->nSpecies-1;
        const double WiByrhoM_0 = WiByrhoM[i+0];
        double dYi0dt = WiByrhoM_0*dPhidt[i+0];
        double ddNi0dtByVdYj = 0;
        __m256d sum = _mm256_setzero_pd();
        const double* __restrict__ JcRowi = &ddNdtByVdcT[i*alignN];

        for (unsigned int k=0; k<this->nSpecies-1; k=k+4)
        {
            __m256d a = _mm256_loadu_pd(&JcRowi[k]);
            __m256d b = _mm256_loadu_pd(&buffer[k]);
            sum = _mm256_fmadd_pd(a, b, sum); 
        }
        __m256d tmp = _mm256_permute2f128_pd(sum, sum, 0x01);
        __m256d sum1 = _mm256_add_pd(sum, tmp);
        __m128d lo = _mm256_castpd256_pd128(sum1);            
        __m128d hi = _mm_unpackhi_pd(lo, lo);                 
        __m128d sum2 = _mm_add_pd(lo, hi);
        double dot =  _mm_cvtsd_f64(sum2);
        ddNi0dtByVdYj += dot;
        unsigned int k = nSpecies-1;
        ddNi0dtByVdYj += JcRowi[k]*buffer[k];
        Jac[(i+0)*(alignN) + j] = WiByrhoM_0*ddNi0dtByVdYj + rhoMvj_*dYi0dt;
    }
}

void 
OptReaction::ddYdtdY_Vec1_2
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ rhoMByRhoi,
    const double* __restrict__ WiByrhoM,  
    const double* __restrict__ dPhidt,  
    const double* __restrict__ Phi, 
    double* __restrict__ Jac
) const noexcept
{
    double* __restrict__ invWPtr = &invW[0];
    __m256d rhoMv = _mm256_set1_pd(rhoM);
    unsigned int remain = 2;
    for(unsigned int j = 0; j < this->nSpecies - remain; j = j + 4)
    {
        const double rhoMvj0_ = rhoMByRhoi[j+0];
        const double rhoMvj1_ = rhoMByRhoi[j+1]; 
        const double rhoMvj2_ = rhoMByRhoi[j+2]; 
        const double rhoMvj3_ = rhoMByRhoi[j+3]; 
        __m256d rhoMvj_ = _mm256_loadu_pd(&rhoMByRhoi[j+0]);
        for(unsigned int i = 0; i < this->nSpecies-remain; i=i+4)
        {
            __m256d rhoMByWiYTv = _mm256_mul_pd(_mm256_loadu_pd(&invWPtr[i+0]), _mm256_loadu_pd(&Phi[i+0]));
            rhoMByWiYTv = _mm256_mul_pd(-rhoMv,rhoMByWiYTv);
            _mm256_storeu_pd(&buffer[i*4+0+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0x00),rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+4+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0x55),rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+8+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0xAA),rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+12+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0xFF),rhoMvj_));    
        }

        {
            unsigned int i = this->nSpecies-2;
            const double rhoMByWiYTPi0 = -rhoM*invWPtr[i+0]*Phi[i+0];
            const double rhoMByWiYTPi1 = -rhoM*invWPtr[i+1]*Phi[i+1];
            buffer[i*4+0+0] = rhoMByWiYTPi0*rhoMvj0_;
            buffer[i*4+0+1] = rhoMByWiYTPi0*rhoMvj1_;
            buffer[i*4+0+2] = rhoMByWiYTPi0*rhoMvj2_;
            buffer[i*4+0+3] = rhoMByWiYTPi0*rhoMvj3_;   
            buffer[i*4+4+0] = rhoMByWiYTPi1*rhoMvj0_;
            buffer[i*4+4+1] = rhoMByWiYTPi1*rhoMvj1_;
            buffer[i*4+4+2] = rhoMByWiYTPi1*rhoMvj2_;
            buffer[i*4+4+3] = rhoMByWiYTPi1*rhoMvj3_;     
        }
     
        buffer[j*4+0] += rhoM*invWPtr[j+0];
        buffer[j*4+5] += rhoM*invWPtr[j+1];
        buffer[j*4+10] += rhoM*invWPtr[j+2];
        buffer[j*4+15] += rhoM*invWPtr[j+3];


        for(unsigned int i=0; i<this->nSpecies - remain; i=i+4)
        {
            double Wi0ByrhoM_ = WiByrhoM[i+0];
            double Wi1ByrhoM_ = WiByrhoM[i+1];
            double Wi2ByrhoM_ = WiByrhoM[i+2];
            double Wi3ByrhoM_ = WiByrhoM[i+3];

            double dYi0dt = dPhidt[i+0]*Wi0ByrhoM_;
            double dYi1dt = dPhidt[i+1]*Wi1ByrhoM_;
            double dYi2dt = dPhidt[i+2]*Wi2ByrhoM_;
            double dYi3dt = dPhidt[i+3]*Wi3ByrhoM_;

            __m256d ddNi0dtByVdYj = _mm256_setzero_pd();
            __m256d ddNi1dtByVdYj = _mm256_setzero_pd();
            __m256d ddNi2dtByVdYj = _mm256_setzero_pd();
            __m256d ddNi3dtByVdYj = _mm256_setzero_pd();
            const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*(alignN)];
            const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*(alignN)];
            const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*(alignN)];
            const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i+3)*(alignN)];
            for (unsigned int k=0; k<this->nSpecies-remain; k=k+4)
            {

                const double ddNi0dtByVdck0 = JcRowi0[k+0];
                const double ddNi1dtByVdck0 = JcRowi1[k+0];
                const double ddNi2dtByVdck0 = JcRowi2[k+0];
                const double ddNi3dtByVdck0 = JcRowi3[k+0];

                const double ddNi0dtByVdck1 = JcRowi0[k+1];
                const double ddNi1dtByVdck1 = JcRowi1[k+1];
                const double ddNi2dtByVdck1 = JcRowi2[k+1];
                const double ddNi3dtByVdck1 = JcRowi3[k+1];

                const double ddNi0dtByVdck2 = JcRowi0[k+2];
                const double ddNi1dtByVdck2 = JcRowi1[k+2];
                const double ddNi2dtByVdck2 = JcRowi2[k+2];
                const double ddNi3dtByVdck2 = JcRowi3[k+2];

                const double ddNi0dtByVdck3 = JcRowi0[k+3];
                const double ddNi1dtByVdck3 = JcRowi1[k+3];
                const double ddNi2dtByVdck3 = JcRowi2[k+3];
                const double ddNi3dtByVdck3 = JcRowi3[k+3];

                __m256d dCk0dYj = _mm256_loadu_pd(&buffer[k*4+0]);
                __m256d dCk1dYj = _mm256_loadu_pd(&buffer[k*4+4]);
                __m256d dCk2dYj = _mm256_loadu_pd(&buffer[k*4+8]);
                __m256d dCk3dYj = _mm256_loadu_pd(&buffer[k*4+12]);

                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck0),dCk0dYj,ddNi0dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck0),dCk0dYj,ddNi1dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck0),dCk0dYj,ddNi2dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck0),dCk0dYj,ddNi3dtByVdYj);

                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck1),dCk1dYj,ddNi0dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck1),dCk1dYj,ddNi1dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck1),dCk1dYj,ddNi2dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck1),dCk1dYj,ddNi3dtByVdYj);

                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck2),dCk2dYj,ddNi0dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck2),dCk2dYj,ddNi1dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck2),dCk2dYj,ddNi2dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck2),dCk2dYj,ddNi3dtByVdYj);     
                
                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck3),dCk3dYj,ddNi0dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck3),dCk3dYj,ddNi1dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck3),dCk3dYj,ddNi2dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck3),dCk3dYj,ddNi3dtByVdYj); 
            }
            
            unsigned int h = this->nSpecies-2;
            const double ddNi0dtByVdck0 = JcRowi0[h+0];
            const double ddNi1dtByVdck0 = JcRowi1[h+0];
            const double ddNi2dtByVdck0 = JcRowi2[h+0];
            const double ddNi3dtByVdck0 = JcRowi3[h+0];

            const double ddNi0dtByVdck1 = JcRowi0[h+1];
            const double ddNi1dtByVdck1 = JcRowi1[h+1];
            const double ddNi2dtByVdck1 = JcRowi2[h+1];
            const double ddNi3dtByVdck1 = JcRowi3[h+1];

            __m256d dCk0dYj = _mm256_loadu_pd(&buffer[h*4+0]);
            __m256d dCk1dYj = _mm256_loadu_pd(&buffer[h*4+4]);
            ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck0),dCk0dYj,ddNi0dtByVdYj);
            ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck0),dCk0dYj,ddNi1dtByVdYj);
            ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck0),dCk0dYj,ddNi2dtByVdYj);
            ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck0),dCk0dYj,ddNi3dtByVdYj);
            ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck1),dCk1dYj,ddNi0dtByVdYj);
            ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck1),dCk1dYj,ddNi1dtByVdYj);
            ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck1),dCk1dYj,ddNi2dtByVdYj);
            ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck1),dCk1dYj,ddNi3dtByVdYj); 
            __m256d WiByrhoM_0 = _mm256_set1_pd(Wi0ByrhoM_);
            __m256d dYidtv0 = _mm256_set1_pd(dYi0dt);
            __m256d ddYi0dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv0);
            ddYi0dtdYj = _mm256_fmadd_pd(WiByrhoM_0,ddNi0dtByVdYj,ddYi0dtdYj);
            __m256d WiByrhoM_1 = _mm256_set1_pd(Wi1ByrhoM_);
            __m256d dYidtv1 = _mm256_set1_pd(dYi1dt);
            __m256d ddYi1dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv1);
            ddYi1dtdYj = _mm256_fmadd_pd(WiByrhoM_1,ddNi1dtByVdYj,ddYi1dtdYj);
            __m256d WiByrhoM_2 = _mm256_set1_pd(Wi2ByrhoM_);
            __m256d dYidtv2 = _mm256_set1_pd(dYi2dt);
            __m256d ddYi2dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv2);
            ddYi2dtdYj = _mm256_fmadd_pd(WiByrhoM_2,ddNi2dtByVdYj,ddYi2dtdYj);
            __m256d WiByrhoM_3 = _mm256_set1_pd(Wi3ByrhoM_);
            __m256d dYidtv3 = _mm256_set1_pd(dYi3dt);
            __m256d ddYi3dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv3);
            ddYi3dtdYj = _mm256_fmadd_pd(WiByrhoM_3,ddNi3dtByVdYj,ddYi3dtdYj);
            _mm256_storeu_pd(&Jac[(i+0)*(alignN)+ j+0],ddYi0dtdYj);
            _mm256_storeu_pd(&Jac[(i+1)*(alignN)+ j+0],ddYi1dtdYj);
            _mm256_storeu_pd(&Jac[(i+2)*(alignN)+ j+0],ddYi2dtdYj);
            _mm256_storeu_pd(&Jac[(i+3)*(alignN)+ j+0],ddYi3dtdYj);
            {
                unsigned int i0 = this->nSpecies-2;
                unsigned int i1 = this->nSpecies-1;
                Wi0ByrhoM_ = WiByrhoM[i0];
                Wi1ByrhoM_ = WiByrhoM[i1];
                dYi0dt = dPhidt[i0]*WiByrhoM[i0];
                dYi1dt = dPhidt[i1]*WiByrhoM[i1];
                __m256d ddNi0dtByVdYjv = _mm256_setzero_pd();
                __m256d ddNi1dtByVdYjv = _mm256_setzero_pd();
                const double* __restrict__ JcRowi00 = &ddNdtByVdcT[(i0)*(alignN)];
                const double* __restrict__ JcRowi11 = &ddNdtByVdcT[(i1)*(alignN)];                 
                for (unsigned int k=0; k<this->nSpecies-2; k=k+4)
                {
                    const double ddNi0dtByVdck0a = JcRowi00[k+0];
                    const double ddNi1dtByVdck0b = JcRowi11[k+0];
                    __m256d Arrk00 = _mm256_loadu_pd(&buffer[k*4+0]);
                    ddNi0dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck0a),Arrk00,ddNi0dtByVdYjv);
                    ddNi1dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck0b),Arrk00,ddNi1dtByVdYjv);
                    const double ddNi0dtByVdck1a = JcRowi00[k+1];
                    const double ddNi1dtByVdck1b = JcRowi11[k+1];
                    __m256d Arrk10 = _mm256_loadu_pd(&buffer[k*4+4]);
                    ddNi0dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck1a),Arrk10,ddNi0dtByVdYjv);
                    ddNi1dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck1b),Arrk10,ddNi1dtByVdYjv);
                    const double ddNi0dtByVdck2a = JcRowi00[k+2];
                    const double ddNi1dtByVdck2b = JcRowi11[k+2];
                    __m256d Arrk20 = _mm256_loadu_pd(&buffer[k*4+8]);
                    ddNi0dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck2a),Arrk20,ddNi0dtByVdYjv);
                    ddNi1dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck2b),Arrk20,ddNi1dtByVdYjv);
                    const double ddNi0dtByVdck3a = JcRowi00[k+3];
                    const double ddNi1dtByVdck3b = JcRowi11[k+3];
                    __m256d Arrk30 = _mm256_loadu_pd(&buffer[k*4+12]);
                    ddNi0dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck3a),Arrk30,ddNi0dtByVdYjv);
                    ddNi1dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck3b),Arrk30,ddNi1dtByVdYjv);  
                }
                {
                    unsigned int k0 = this->nSpecies-2;
                    unsigned int k1 = this->nSpecies-1;
                    const double ddNi0dtByVdck0a = JcRowi00[k0];
                    const double ddNi1dtByVdck0b = JcRowi11[k0];
                    __m256d Arrk00 = _mm256_loadu_pd(&buffer[k0*4+0]);
                    ddNi0dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck0a),Arrk00,ddNi0dtByVdYjv);
                    ddNi1dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck0b),Arrk00,ddNi1dtByVdYjv); 

                    const double ddNi0dtByVdck1a = JcRowi00[k1];
                    const double ddNi1dtByVdck1b = JcRowi11[k1];
                    __m256d Arrk10 = _mm256_loadu_pd(&buffer[k1*4+0]);
                    ddNi0dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck1a),Arrk10,ddNi0dtByVdYjv);
                    ddNi1dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck1b),Arrk10,ddNi1dtByVdYjv);
                }
                __m256d r0 = _mm256_mul_pd(ddNi0dtByVdYjv,_mm256_set1_pd(Wi0ByrhoM_));
                r0 = _mm256_fmadd_pd(rhoMvj_,_mm256_set1_pd(dYi0dt),r0);
                _mm256_storeu_pd(&Jac[i0*(alignN) + j+0],r0);
                __m256d r1 = _mm256_mul_pd(ddNi1dtByVdYjv,_mm256_set1_pd(Wi1ByrhoM_));
                r1 = _mm256_fmadd_pd(rhoMvj_,_mm256_set1_pd(dYi1dt),r1);
                _mm256_storeu_pd(&Jac[i1*(alignN) + j+0],r1);
            }
        }
    }
    
    unsigned int j0 = this->nSpecies-2;
    unsigned int j1 = this->nSpecies-1;

    __m128d rhoMvjv = _mm_loadu_pd(&rhoMByRhoi[j0]);
    for(unsigned int i=0; i<this->nSpecies-2; i=i+4)
    {
        const double rhoMByWi0 = rhoM*invWPtr[i+0];
        const double rhoMByWi1 = rhoM*invWPtr[i+1];
        const double rhoMByWi2 = rhoM*invWPtr[i+2];
        const double rhoMByWi3 = rhoM*invWPtr[i+3];

        __m128d Arr0 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-rhoMByWi0),rhoMvjv),_mm_set1_pd(Phi[i+0]));
        _mm_storeu_pd(&buffer[i*4+0],Arr0);
        __m128d Arr1 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-rhoMByWi1),rhoMvjv),_mm_set1_pd(Phi[i+1]));
        _mm_storeu_pd(&buffer[i*4+4],Arr1);
        __m128d Arr2 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-rhoMByWi2),rhoMvjv),_mm_set1_pd(Phi[i+2]));
        _mm_storeu_pd(&buffer[i*4+8],Arr2); 
        __m128d Arr3 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-rhoMByWi3),rhoMvjv),_mm_set1_pd(Phi[i+3]));
        _mm_storeu_pd(&buffer[i*4+12],Arr3);                
    }
    {
        unsigned int i = this->nSpecies-2;
        const double rhoMByWi0 = rhoM*invWPtr[i+0];
        const double rhoMByWi1 = rhoM*invWPtr[i+1];
        __m128d Arr0 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-rhoMByWi0),rhoMvjv),_mm_set1_pd(Phi[i+0]));
        _mm_storeu_pd(&buffer[i*4+0],Arr0);  
        __m128d Arr1 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-rhoMByWi1),rhoMvjv),_mm_set1_pd(Phi[i+1]));
        _mm_storeu_pd(&buffer[i*4+4],Arr1);             
    }
    buffer[j0*4+0] = buffer[j0*4+0] + rhoM*invWPtr[j0];
    buffer[j1*4+1] = buffer[j1*4+1] + rhoM*invWPtr[j1];

    for(unsigned int i=0; i<this->nSpecies-2; i=i+4)
    {
        const double WiByrhoM_0 = WiByrhoM[i+0];
        const double WiByrhoM_1 = WiByrhoM[i+1];
        const double WiByrhoM_2 = WiByrhoM[i+2];
        const double WiByrhoM_3 = WiByrhoM[i+3];

        double dYidt0 = dPhidt[i+0]*WiByrhoM_0;
        double dYidt1 = dPhidt[i+1]*WiByrhoM_1;
        double dYidt2 = dPhidt[i+2]*WiByrhoM_2;
        double dYidt3 = dPhidt[i+3]*WiByrhoM_3;

        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*(alignN)];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*(alignN)];
        const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*(alignN)];
        const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i+3)*(alignN)];

        __m128d ddNi0dtByVdYjv = _mm_setzero_pd();
        __m128d ddNi1dtByVdYjv = _mm_setzero_pd();
        __m128d ddNi2dtByVdYjv = _mm_setzero_pd();
        __m128d ddNi3dtByVdYjv = _mm_setzero_pd();
        for (unsigned int k=0; k<this->nSpecies-2; k=k+4)
        {
            const double ddNi0dtByVdck0 = JcRowi0[k+0];
            const double ddNi0dtByVdck1 = JcRowi0[k+1];
            const double ddNi0dtByVdck2 = JcRowi0[k+2];
            const double ddNi0dtByVdck3 = JcRowi0[k+3];

            const double ddNi1dtByVdck0 = JcRowi1[k+0];
            const double ddNi1dtByVdck1 = JcRowi1[k+1];
            const double ddNi1dtByVdck2 = JcRowi1[k+2];
            const double ddNi1dtByVdck3 = JcRowi1[k+3];

            const double ddNi2dtByVdck0 = JcRowi2[k+0];
            const double ddNi2dtByVdck1 = JcRowi2[k+1];
            const double ddNi2dtByVdck2 = JcRowi2[k+2];
            const double ddNi2dtByVdck3 = JcRowi2[k+3];

            const double ddNi3dtByVdck0 = JcRowi3[k+0];
            const double ddNi3dtByVdck1 = JcRowi3[k+1];
            const double ddNi3dtByVdck2 = JcRowi3[k+2];
            const double ddNi3dtByVdck3 = JcRowi3[k+3];


            __m128d dCk0dYj = _mm_loadu_pd(&buffer[k*4+0]);
            __m128d dCk1dYj = _mm_loadu_pd(&buffer[k*4+4]);
            __m128d dCk2dYj = _mm_loadu_pd(&buffer[k*4+8]);
            __m128d dCk3dYj = _mm_loadu_pd(&buffer[k*4+12]);

            ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck0),dCk0dYj,ddNi0dtByVdYjv);
            ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck1),dCk1dYj,ddNi0dtByVdYjv);
            ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck2),dCk2dYj,ddNi0dtByVdYjv);
            ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck3),dCk3dYj,ddNi0dtByVdYjv);

            ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck0),dCk0dYj,ddNi1dtByVdYjv);
            ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck1),dCk1dYj,ddNi1dtByVdYjv);
            ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck2),dCk2dYj,ddNi1dtByVdYjv);
            ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck3),dCk3dYj,ddNi1dtByVdYjv);   
            
            ddNi2dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi2dtByVdck0),dCk0dYj,ddNi2dtByVdYjv);
            ddNi2dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi2dtByVdck1),dCk1dYj,ddNi2dtByVdYjv);
            ddNi2dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi2dtByVdck2),dCk2dYj,ddNi2dtByVdYjv);
            ddNi2dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi2dtByVdck3),dCk3dYj,ddNi2dtByVdYjv);

            ddNi3dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi3dtByVdck0),dCk0dYj,ddNi3dtByVdYjv);
            ddNi3dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi3dtByVdck1),dCk1dYj,ddNi3dtByVdYjv);
            ddNi3dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi3dtByVdck2),dCk2dYj,ddNi3dtByVdYjv);
            ddNi3dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi3dtByVdck3),dCk3dYj,ddNi3dtByVdYjv);            
        }
        {
            unsigned int k = this->nSpecies-2;
            const double ddNi0dtByVdck0 = JcRowi0[k+0];
            const double ddNi0dtByVdck1 = JcRowi0[k+1];

            const double ddNi1dtByVdck0 = JcRowi1[k+0];
            const double ddNi1dtByVdck1 = JcRowi1[k+1];  

            const double ddNi2dtByVdck0 = JcRowi2[k+0];
            const double ddNi2dtByVdck1 = JcRowi2[k+1];  

            const double ddNi3dtByVdck0 = JcRowi3[k+0];
            const double ddNi3dtByVdck1 = JcRowi3[k+1];
            __m128d dCk0dYj0 = _mm_loadu_pd(&buffer[k*4+0]);
            __m128d dCk1dYj0 = _mm_loadu_pd(&buffer[k*4+4]);
            ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck0),dCk0dYj0,ddNi0dtByVdYjv);
            ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck1),dCk1dYj0,ddNi0dtByVdYjv);

            ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck0),dCk0dYj0,ddNi1dtByVdYjv);
            ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck1),dCk1dYj0,ddNi1dtByVdYjv);
            
            ddNi2dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi2dtByVdck0),dCk0dYj0,ddNi2dtByVdYjv);
            ddNi2dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi2dtByVdck1),dCk1dYj0,ddNi2dtByVdYjv);
            
            ddNi3dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi3dtByVdck0),dCk0dYj0,ddNi3dtByVdYjv);
            ddNi3dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi3dtByVdck1),dCk1dYj0,ddNi3dtByVdYjv);
        }
        __m128d Wi0ByrhoMv = _mm_set1_pd(WiByrhoM_0);
        __m128d Wi1ByrhoMv = _mm_set1_pd(WiByrhoM_1);
        __m128d Wi2ByrhoMv = _mm_set1_pd(WiByrhoM_2);
        __m128d Wi3ByrhoMv = _mm_set1_pd(WiByrhoM_3);

        __m128d rhoMByrhoi = _mm_loadu_pd(&rhoMByRhoi[j0]);

        __m128d dYi0dtv = _mm_set1_pd(dYidt0);
        __m128d dYi1dtv = _mm_set1_pd(dYidt1);
        __m128d dYi2dtv = _mm_set1_pd(dYidt2);
        __m128d dYi3dtv = _mm_set1_pd(dYidt3);

        __m128d result0 = _mm_mul_pd(Wi0ByrhoMv,ddNi0dtByVdYjv);
        result0 = _mm_fmadd_pd(rhoMByrhoi,dYi0dtv,result0);
        _mm_storeu_pd(&Jac[(i+0)*(alignN) + j0],result0);

        __m128d result1 = _mm_mul_pd(Wi1ByrhoMv,ddNi1dtByVdYjv);
        result1 = _mm_fmadd_pd(rhoMByrhoi,dYi1dtv,result1);
        _mm_storeu_pd(&Jac[(i+1)*(alignN) + j0],result1);

        __m128d result2 = _mm_mul_pd(Wi2ByrhoMv,ddNi2dtByVdYjv);
        result2 = _mm_fmadd_pd(rhoMByrhoi,dYi2dtv,result2);
        _mm_storeu_pd(&Jac[(i+2)*(alignN) + j0],result2);

        __m128d result3 = _mm_mul_pd(Wi3ByrhoMv,ddNi3dtByVdYjv);
        result3 = _mm_fmadd_pd(rhoMByrhoi,dYi3dtv,result3);
        _mm_storeu_pd(&Jac[(i+3)*(alignN) + j0],result3);               
    }

    {
        unsigned int i = nSpecies - 2;

        const double WiByrhoM_0 = WiByrhoM[i+0];
        const double WiByrhoM_1 = WiByrhoM[i+1];

        double dYidt0 = dPhidt[i+0]*WiByrhoM_0;
        double dYidt1 = dPhidt[i+1]*WiByrhoM_1;

        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*(alignN)];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*(alignN)];
        __m128d ddNi0dtByVdYjv = _mm_setzero_pd();
        __m128d ddNi1dtByVdYjv = _mm_setzero_pd();
        for (unsigned int k=0; k<this->nSpecies-2; k=k+4)
        {
            const double ddNi0dtByVdck0 = JcRowi0[k+0];
            const double ddNi0dtByVdck1 = JcRowi0[k+1];
            const double ddNi0dtByVdck2 = JcRowi0[k+2];
            const double ddNi0dtByVdck3 = JcRowi0[k+3];
            const double ddNi1dtByVdck0 = JcRowi1[k+0];
            const double ddNi1dtByVdck1 = JcRowi1[k+1];
            const double ddNi1dtByVdck2 = JcRowi1[k+2];
            const double ddNi1dtByVdck3 = JcRowi1[k+3]; 

            __m128d dCk0dYj = _mm_loadu_pd(&buffer[k*4+0]);
            __m128d dCk1dYj = _mm_loadu_pd(&buffer[k*4+4]);
            __m128d dCk2dYj = _mm_loadu_pd(&buffer[k*4+8]);
            __m128d dCk3dYj = _mm_loadu_pd(&buffer[k*4+12]);

            ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck0),dCk0dYj,ddNi0dtByVdYjv);
            ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck1),dCk1dYj,ddNi0dtByVdYjv);
            ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck2),dCk2dYj,ddNi0dtByVdYjv);
            ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck3),dCk3dYj,ddNi0dtByVdYjv);

            ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck0),dCk0dYj,ddNi1dtByVdYjv);
            ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck1),dCk1dYj,ddNi1dtByVdYjv);
            ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck2),dCk2dYj,ddNi1dtByVdYjv);
            ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck3),dCk3dYj,ddNi1dtByVdYjv);            
        }
        {
            unsigned int k = this->nSpecies-2;
            const double ddNi0dtByVdck0 = JcRowi0[k+0];
            const double ddNi0dtByVdck1 = JcRowi0[k+1];

            const double ddNi1dtByVdck0 = JcRowi1[k+0];
            const double ddNi1dtByVdck1 = JcRowi1[k+1];
            __m128d dCk0dYj0 = _mm_loadu_pd(&buffer[k*4+0]);
            __m128d dCk1dYj0 = _mm_loadu_pd(&buffer[k*4+4]);
            ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck0),dCk0dYj0,ddNi0dtByVdYjv);
            ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck1),dCk1dYj0,ddNi0dtByVdYjv);

            ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck0),dCk0dYj0,ddNi1dtByVdYjv);
            ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck1),dCk1dYj0,ddNi1dtByVdYjv);            
        }
        __m128d Wi0ByrhoMv = _mm_set1_pd(WiByrhoM_0);
        __m128d Wi1ByrhoMv = _mm_set1_pd(WiByrhoM_1);
        __m128d rhiMvjv = _mm_loadu_pd(&rhoMByRhoi[j0]);
        __m128d dYi0dtv = _mm_set1_pd(dYidt0);
        __m128d dYi1dtv = _mm_set1_pd(dYidt1);
        __m128d result0 = _mm_mul_pd(Wi0ByrhoMv,ddNi0dtByVdYjv);
        __m128d result1 = _mm_mul_pd(Wi1ByrhoMv,ddNi1dtByVdYjv);
        result0 = _mm_fmadd_pd(rhiMvjv,dYi0dtv,result0);
        result1 = _mm_fmadd_pd(rhiMvjv,dYi1dtv,result1);
        _mm_storeu_pd(&Jac[(i+0)*(alignN) + j0],result0);
        _mm_storeu_pd(&Jac[(i+1)*(alignN) + j0],result1);
    }  
}

void 
OptReaction::ddYdtdY_Vec1_3
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ rhoMByRhoi,
    const double* __restrict__ WiByrhoM,  
    const double* __restrict__ dPhidt,  
    const double* __restrict__ Phi, 
    double* __restrict__ Jac
) const noexcept
{
    double* __restrict__ invWPtr = &invW[0];
    __m256d rhoMv = _mm256_set1_pd(rhoM);
    for(unsigned int j = 0; j < this->nSpecies - 3; j = j + 4)
    {
        __m256d rhoMvj_ = _mm256_loadu_pd(&rhoMByRhoi[j+0]);
        for(unsigned int i = 0; i < this->nSpecies-3; i=i+4)
        {
            __m256d rhoMByWiYTv = _mm256_mul_pd(_mm256_loadu_pd(&invWPtr[i+0]), _mm256_loadu_pd(&Phi[i+0]));
            rhoMByWiYTv = _mm256_mul_pd(-rhoMv,rhoMByWiYTv);
            _mm256_storeu_pd(&buffer[i*4+0+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0x00),rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+4+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0x55),rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+8+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0xAA),rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+12+0],_mm256_mul_pd(_mm256_permute4x64_pd(rhoMByWiYTv, 0xFF),rhoMvj_));              
        }
        {
            unsigned int i = this->nSpecies-3;
            const double rhoMByWiYTPi0 = -rhoM*invWPtr[i+0]*Phi[i+0];
            const double rhoMByWiYTPi1 = -rhoM*invWPtr[i+1]*Phi[i+1];
            const double rhoMByWiYTPi2 = -rhoM*invWPtr[i+2]*Phi[i+2];

            __m256d rhoMByWiYTPi0v = _mm256_set1_pd(rhoMByWiYTPi0);
            __m256d rhoMByWiYTPi1v = _mm256_set1_pd(rhoMByWiYTPi1);
            __m256d rhoMByWiYTPi2v = _mm256_set1_pd(rhoMByWiYTPi2);
            _mm256_storeu_pd(&buffer[i*4],_mm256_mul_pd(rhoMByWiYTPi0v,rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+4],_mm256_mul_pd(rhoMByWiYTPi1v,rhoMvj_));
            _mm256_storeu_pd(&buffer[i*4+8],_mm256_mul_pd(rhoMByWiYTPi2v,rhoMvj_));                        
        }
        buffer[j*4] += rhoM*invWPtr[j+0];
        buffer[j*4+5] += rhoM*invWPtr[j+1];
        buffer[j*4+10] += rhoM*invWPtr[j+2];
        buffer[j*4+15] += rhoM*invWPtr[j+3];  

        for(unsigned int i=0; i<this->nSpecies - 3; i=i+4)
        {
            const double Wi0ByrhoM_ = WiByrhoM[i+0];
            const double Wi1ByrhoM_ = WiByrhoM[i+1];
            const double Wi2ByrhoM_ = WiByrhoM[i+2];
            const double Wi3ByrhoM_ = WiByrhoM[i+3];

            const double dYi0dt = dPhidt[i+0]*Wi0ByrhoM_;
            const double dYi1dt = dPhidt[i+1]*Wi1ByrhoM_;
            const double dYi2dt = dPhidt[i+2]*Wi2ByrhoM_;
            const double dYi3dt = dPhidt[i+3]*Wi3ByrhoM_;

            __m256d ddNi0dtByVdYj = _mm256_setzero_pd();
            __m256d ddNi1dtByVdYj = _mm256_setzero_pd();
            __m256d ddNi2dtByVdYj = _mm256_setzero_pd();
            __m256d ddNi3dtByVdYj = _mm256_setzero_pd();
            const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*(alignN)];
            const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*(alignN)];
            const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*(alignN)];
            const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i+3)*(alignN)];
            for (unsigned int k=0; k<this->nSpecies-3; k=k+4)
            {

                const double ddNi0dtByVdck0 = JcRowi0[k+0];
                const double ddNi1dtByVdck0 = JcRowi1[k+0];
                const double ddNi2dtByVdck0 = JcRowi2[k+0];
                const double ddNi3dtByVdck0 = JcRowi3[k+0];

                const double ddNi0dtByVdck1 = JcRowi0[k+1];
                const double ddNi1dtByVdck1 = JcRowi1[k+1];
                const double ddNi2dtByVdck1 = JcRowi2[k+1];
                const double ddNi3dtByVdck1 = JcRowi3[k+1];

                const double ddNi0dtByVdck2 = JcRowi0[k+2];
                const double ddNi1dtByVdck2 = JcRowi1[k+2];
                const double ddNi2dtByVdck2 = JcRowi2[k+2];
                const double ddNi3dtByVdck2 = JcRowi3[k+2];

                const double ddNi0dtByVdck3 = JcRowi0[k+3];
                const double ddNi1dtByVdck3 = JcRowi1[k+3];
                const double ddNi2dtByVdck3 = JcRowi2[k+3];
                const double ddNi3dtByVdck3 = JcRowi3[k+3];

                __m256d dCk0dYj = _mm256_loadu_pd(&buffer[k*4]);
                __m256d dCk1dYj = _mm256_loadu_pd(&buffer[k*4+4]);
                __m256d dCk2dYj = _mm256_loadu_pd(&buffer[k*4+8]);
                __m256d dCk3dYj = _mm256_loadu_pd(&buffer[k*4+12]);

                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck0),dCk0dYj,ddNi0dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck0),dCk0dYj,ddNi1dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck0),dCk0dYj,ddNi2dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck0),dCk0dYj,ddNi3dtByVdYj);

                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck1),dCk1dYj,ddNi0dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck1),dCk1dYj,ddNi1dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck1),dCk1dYj,ddNi2dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck1),dCk1dYj,ddNi3dtByVdYj);

                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck2),dCk2dYj,ddNi0dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck2),dCk2dYj,ddNi1dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck2),dCk2dYj,ddNi2dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck2),dCk2dYj,ddNi3dtByVdYj);     
                
                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck3),dCk3dYj,ddNi0dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck3),dCk3dYj,ddNi1dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck3),dCk3dYj,ddNi2dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck3),dCk3dYj,ddNi3dtByVdYj); 
            }

            {
                unsigned int k = this->nSpecies-3;
                const double ddNi0dtByVdck0 = JcRowi0[k+0];
                const double ddNi1dtByVdck0 = JcRowi1[k+0];
                const double ddNi2dtByVdck0 = JcRowi2[k+0];
                const double ddNi3dtByVdck0 = JcRowi3[k+0];

                const double ddNi0dtByVdck1 = JcRowi0[k+1];
                const double ddNi1dtByVdck1 = JcRowi1[k+1];
                const double ddNi2dtByVdck1 = JcRowi2[k+1];
                const double ddNi3dtByVdck1 = JcRowi3[k+1];

                const double ddNi0dtByVdck2 = JcRowi0[k+2];
                const double ddNi1dtByVdck2 = JcRowi1[k+2];
                const double ddNi2dtByVdck2 = JcRowi2[k+2];
                const double ddNi3dtByVdck2 = JcRowi3[k+2];

                __m256d dCk0dYj = _mm256_loadu_pd(&buffer[k*4+0]);
                __m256d dCk1dYj = _mm256_loadu_pd(&buffer[k*4+4]);
                __m256d dCk2dYj = _mm256_loadu_pd(&buffer[k*4+8]);

                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck0),dCk0dYj,ddNi0dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck0),dCk0dYj,ddNi1dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck0),dCk0dYj,ddNi2dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck0),dCk0dYj,ddNi3dtByVdYj);

                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck1),dCk1dYj,ddNi0dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck1),dCk1dYj,ddNi1dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck1),dCk1dYj,ddNi2dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck1),dCk1dYj,ddNi3dtByVdYj);

                ddNi0dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck2),dCk2dYj,ddNi0dtByVdYj);
                ddNi1dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck2),dCk2dYj,ddNi1dtByVdYj);
                ddNi2dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck2),dCk2dYj,ddNi2dtByVdYj);
                ddNi3dtByVdYj = _mm256_fmadd_pd(_mm256_set1_pd(ddNi3dtByVdck2),dCk2dYj,ddNi3dtByVdYj);     
            }
            __m256d WiByrhoM_0 = _mm256_set1_pd(Wi0ByrhoM_);
            __m256d dYidtv0 = _mm256_set1_pd(dYi0dt);
            __m256d ddYi0dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv0);
            ddYi0dtdYj = _mm256_fmadd_pd(WiByrhoM_0,ddNi0dtByVdYj,ddYi0dtdYj);
            __m256d WiByrhoM_1 = _mm256_set1_pd(Wi1ByrhoM_);
            __m256d dYidtv1 = _mm256_set1_pd(dYi1dt);
            __m256d ddYi1dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv1);
            ddYi1dtdYj = _mm256_fmadd_pd(WiByrhoM_1,ddNi1dtByVdYj,ddYi1dtdYj);
            __m256d WiByrhoM_2 = _mm256_set1_pd(Wi2ByrhoM_);
            __m256d dYidtv2 = _mm256_set1_pd(dYi2dt);
            __m256d ddYi2dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv2);
            ddYi2dtdYj = _mm256_fmadd_pd(WiByrhoM_2,ddNi2dtByVdYj,ddYi2dtdYj);
            __m256d WiByrhoM_3 = _mm256_set1_pd(Wi3ByrhoM_);
            __m256d dYidtv3 = _mm256_set1_pd(dYi3dt);
            __m256d ddYi3dtdYj = _mm256_mul_pd(rhoMvj_,dYidtv3);
            ddYi3dtdYj = _mm256_fmadd_pd(WiByrhoM_3,ddNi3dtByVdYj,ddYi3dtdYj);
            _mm256_storeu_pd(&Jac[(i+0)*(alignN)+ j+0],ddYi0dtdYj);
            _mm256_storeu_pd(&Jac[(i+1)*(alignN)+ j+0],ddYi1dtdYj);
            _mm256_storeu_pd(&Jac[(i+2)*(alignN)+ j+0],ddYi2dtdYj);
            _mm256_storeu_pd(&Jac[(i+3)*(alignN)+ j+0],ddYi3dtdYj);
        }     
        
        {
            unsigned int i0 = this->nSpecies-3;            
            unsigned int i1 = this->nSpecies-2;
            unsigned int i2 = this->nSpecies-1;
            const double Wi0ByrhoM_ = WiByrhoM[i0];
            const double Wi1ByrhoM_ = WiByrhoM[i1];
            const double Wi2ByrhoM_ = WiByrhoM[i2];
            double dYi0dt = dPhidt[i0]*Wi0ByrhoM_;
            double dYi1dt = dPhidt[i1]*Wi1ByrhoM_;
            double dYi2dt = dPhidt[i2]*Wi2ByrhoM_;

            __m256d ddNi0dtByVdYjv = _mm256_setzero_pd();
            __m256d ddNi1dtByVdYjv = _mm256_setzero_pd();
            __m256d ddNi2dtByVdYjv = _mm256_setzero_pd();
            const double* __restrict__ JcRowi0 = &ddNdtByVdcT[i0*(alignN)];
            const double* __restrict__ JcRowi1 = &ddNdtByVdcT[i1*(alignN)];
            const double* __restrict__ JcRowi2 = &ddNdtByVdcT[i2*(alignN)];
            for (unsigned int k=0; k<this->nSpecies; k++)
            {
                const double ddNi0dtByVdck = JcRowi0[k];
                const double ddNi1dtByVdck = JcRowi1[k];
                const double ddNi2dtByVdck = JcRowi2[k];
                __m256d Arrkv = _mm256_loadu_pd(&buffer[k*4]);
                ddNi0dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi0dtByVdck),Arrkv,ddNi0dtByVdYjv);
                ddNi1dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi1dtByVdck),Arrkv,ddNi1dtByVdYjv);
                ddNi2dtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNi2dtByVdck),Arrkv,ddNi2dtByVdYjv);
            }
            __m256d result0 = _mm256_mul_pd(rhoMvj_,_mm256_set1_pd(dYi0dt));
            result0 = _mm256_fmadd_pd(_mm256_set1_pd(Wi0ByrhoM_),ddNi0dtByVdYjv,result0);
            _mm256_storeu_pd(&Jac[i0*(alignN) + j+0],result0);
            __m256d result1 = _mm256_mul_pd(rhoMvj_,_mm256_set1_pd(dYi1dt));
            result1 = _mm256_fmadd_pd(_mm256_set1_pd(Wi1ByrhoM_),ddNi1dtByVdYjv,result1);
            _mm256_storeu_pd(&Jac[i1*(alignN) + j+0],result1);    
            __m256d result2 = _mm256_mul_pd(rhoMvj_,_mm256_set1_pd(dYi2dt));
            result2 = _mm256_fmadd_pd(_mm256_set1_pd(Wi2ByrhoM_),ddNi2dtByVdYjv,result2);
            _mm256_storeu_pd(&Jac[i2*(alignN) + j+0],result2);             
        }        
    }

    {
        unsigned int j0 = this->nSpecies-3;
        unsigned int j1 = this->nSpecies-2;
        unsigned int j2 = this->nSpecies-1;

        for(unsigned int i=0; i<this->nSpecies; i++)
        {
            const double rhoMByWi = rhoM*invWPtr[i];
            const double rhoMvj0_ = rhoMByRhoi[j0];
            const double rhoMvj1_ = rhoMByRhoi[j1];
            const double rhoMvj2_ = rhoMByRhoi[j2];
            buffer[i*4+0] = rhoMByWi*(0 - rhoMvj0_*Phi[i]);
            buffer[i*4+1] = rhoMByWi*(0 - rhoMvj1_*Phi[i]);
            buffer[i*4+2] = rhoMByWi*(0 - rhoMvj2_*Phi[i]);
            buffer[i*4+3] = 0;
        }
        buffer[j0*4+0] = buffer[j0*4+0] + rhoM*invWPtr[j0];
        buffer[j1*4+1] = buffer[j1*4+1] + rhoM*invWPtr[j1];
        buffer[j2*4+2] = buffer[j2*4+2] + rhoM*invWPtr[j2];
        for(unsigned int i=0; i<this->nSpecies; i++)
        {
            const double WiByrhoM_ = WiByrhoM[i];
            const double rhoMvj0_ = rhoMByRhoi[j0];
            const double rhoMvj1_ = rhoMByRhoi[j1];
            const double rhoMvj2_ = rhoMByRhoi[j2];

            double dYidt = dPhidt[i]*WiByrhoM[i];

            __m256d ddNidtByVdYjv = _mm256_setzero_pd();
            const double* __restrict__ JcRowi = &ddNdtByVdcT[i*(alignN)];
            for (unsigned int k=0; k<this->nSpecies-3; k=k+4)
            {
                const double ddNidtByVdck0 = JcRowi[k+0];
                ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNidtByVdck0),_mm256_loadu_pd(&buffer[k*4]),ddNidtByVdYjv);
                const double ddNidtByVdck1 = JcRowi[k+1];
                ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNidtByVdck1),_mm256_loadu_pd(&buffer[k*4+4]),ddNidtByVdYjv);
                const double ddNidtByVdck2 = JcRowi[k+2];
                ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNidtByVdck2),_mm256_loadu_pd(&buffer[k*4+8]),ddNidtByVdYjv);
                const double ddNidtByVdck3 = JcRowi[k+3];
                ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNidtByVdck3),_mm256_loadu_pd(&buffer[k*4+12]),ddNidtByVdYjv);
            }
            {
                unsigned int k = this->nSpecies-3;
                const double ddNidtByVdck0 = JcRowi[k+0];
                ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNidtByVdck0),_mm256_loadu_pd(&buffer[k*4]),ddNidtByVdYjv);
                const double ddNidtByVdck1 = JcRowi[k+1];
                ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNidtByVdck1),_mm256_loadu_pd(&buffer[k*4+4]),ddNidtByVdYjv);
                const double ddNidtByVdck2 = JcRowi[k+2];   
                ddNidtByVdYjv = _mm256_fmadd_pd(_mm256_set1_pd(ddNidtByVdck2),_mm256_loadu_pd(&buffer[k*4+8]),ddNidtByVdYjv);       
            }

            __m256d rhoMvjv = _mm256_setr_pd(rhoMvj0_,rhoMvj1_,rhoMvj2_,0);
            __m256d result = _mm256_mul_pd(rhoMvjv,_mm256_set1_pd(dYidt));
            result = _mm256_fmadd_pd(ddNidtByVdYjv,_mm256_set1_pd(WiByrhoM_),result);
            __m128d lo = _mm256_castpd256_pd128(result);         
            __m128d hi = _mm256_extractf128_pd(result, 1);       
            double r0 = _mm_cvtsd_f64(lo);
            double r1 = _mm_cvtsd_f64(_mm_unpackhi_pd(lo, lo));
            double r2 = _mm_cvtsd_f64(hi);
            Jac[i*(alignN) + j0] = r0;
            Jac[i*(alignN) + j1] = r1;
            Jac[i*(alignN) + j2] = r2;
        }
    }        
}



void 
OptReaction::ddYdtdTP_Vec_3
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ WiByrhoM,
    const double* __restrict__ c,  
    double* __restrict__ dPhidt, 
    double* __restrict__ J
) const noexcept
{

    double alphavM = this->invT;
    __m256d alphavMv = _mm256_set1_pd(this->invT);
    for (unsigned int i=0; i<this->nSpecies-3; i=i+4)
    {
        unsigned int i0 = i;
        unsigned int i1 = i+1;
        unsigned int i2 = i+2;
        unsigned int i3 = i+3;

        __m256d Wi0ByrhoM_v = _mm256_loadu_pd(&WiByrhoM[i0]);
        __m256d dYTpdtv = _mm256_loadu_pd(&dPhidt[i0]);   
        dYTpdtv = _mm256_mul_pd(dYTpdtv,Wi0ByrhoM_v);   
        _mm256_storeu_pd(&dPhidt[i0],dYTpdtv);
        __m256d dYi0dtv = _mm256_loadu_pd(&dPhidt[i0]);

        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i0)*(alignN)];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i1)*(alignN)];
        const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i2)*(alignN)];
        const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i3)*(alignN)];

        double ddNi0dtByVdT = JcRowi0[nSpecies];
        double ddNi1dtByVdT = JcRowi1[nSpecies];
        double ddNi2dtByVdT = JcRowi2[nSpecies];
        double ddNi3dtByVdT = JcRowi3[nSpecies];
        __m256d sum0 = _mm256_setzero_pd();
        __m256d sum1 = _mm256_setzero_pd();
        __m256d sum2 = _mm256_setzero_pd();
        __m256d sum3 = _mm256_setzero_pd();
        for (unsigned int j=0; j<this->nSpecies-3; j=j+4)
        {
            __m256d Cv = _mm256_loadu_pd(&c[j+0]);
            sum0 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi0[j]),Cv),alphavMv,sum0);
            sum1 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi1[j]),Cv),alphavMv,sum1);
            sum2 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi2[j]),Cv),alphavMv,sum2);
            sum3 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi3[j]),Cv),alphavMv,sum3);
        }
        ddNi0dtByVdT = ddNi0dtByVdT - (hsum4(sum0));
        ddNi1dtByVdT = ddNi1dtByVdT - (hsum4(sum1));
        ddNi2dtByVdT = ddNi2dtByVdT - (hsum4(sum2));
        ddNi3dtByVdT = ddNi3dtByVdT - (hsum4(sum3));
        {
            unsigned int j0 = nSpecies-3;
            unsigned int j1 = nSpecies-2;
            unsigned int j2 = nSpecies-1;
            ddNi0dtByVdT -= JcRowi0[j0]*c[j0]*alphavM;
            ddNi1dtByVdT -= JcRowi1[j0]*c[j0]*alphavM;
            ddNi2dtByVdT -= JcRowi2[j0]*c[j0]*alphavM;
            ddNi3dtByVdT -= JcRowi3[j0]*c[j0]*alphavM;

            ddNi0dtByVdT -= JcRowi0[j1]*c[j1]*alphavM;
            ddNi1dtByVdT -= JcRowi1[j1]*c[j1]*alphavM;
            ddNi2dtByVdT -= JcRowi2[j1]*c[j1]*alphavM;
            ddNi3dtByVdT -= JcRowi3[j1]*c[j1]*alphavM;

            ddNi0dtByVdT -= JcRowi0[j2]*c[j2]*alphavM;
            ddNi1dtByVdT -= JcRowi1[j2]*c[j2]*alphavM;
            ddNi2dtByVdT -= JcRowi2[j2]*c[j2]*alphavM;
            ddNi3dtByVdT -= JcRowi3[j2]*c[j2]*alphavM;
        }
        __m256d ddNi0dtByVdTv = _mm256_setr_pd(ddNi0dtByVdT,ddNi1dtByVdT,ddNi2dtByVdT,ddNi3dtByVdT);
        __m256d result = _mm256_fmadd_pd(Wi0ByrhoM_v,ddNi0dtByVdTv,_mm256_mul_pd(alphavMv,dYi0dtv));
        J[i0*(alignN) + nSpecies] = get_elem0(result);
        J[i1*(alignN) + nSpecies] = get_elem1(result);
        J[i2*(alignN) + nSpecies] = get_elem2(result);
        J[i3*(alignN) + nSpecies] = get_elem3(result);
      
    }
    {
        unsigned int i0 = this->nSpecies-3;
        unsigned int i1 = this->nSpecies-2;
        unsigned int i2 = this->nSpecies-1;
        const double Wi0ByrhoM_ = WiByrhoM[i0];
        const double Wi1ByrhoM_ = WiByrhoM[i1];
        const double Wi2ByrhoM_ = WiByrhoM[i2];

        dPhidt[i0] = dPhidt[i0]*Wi0ByrhoM_;
        dPhidt[i1] = dPhidt[i1]*Wi1ByrhoM_;
        dPhidt[i2] = dPhidt[i2]*Wi2ByrhoM_;

        double dYi0dt = dPhidt[i0];
        double dYi1dt = dPhidt[i1];
        double dYi2dt = dPhidt[i2];

        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i0)*(alignN)];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i1)*(alignN)];
        const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i2)*(alignN)];

        double ddNi0dtByVdT = JcRowi0[nSpecies];
        double ddNi1dtByVdT = JcRowi1[nSpecies];
        double ddNi2dtByVdT = JcRowi2[nSpecies];
        __m256d ddNi0dtByVdTv = _mm256_setzero_pd();
        __m256d ddNi1dtByVdTv = _mm256_setzero_pd();
        __m256d ddNi2dtByVdTv = _mm256_setzero_pd();
        for (unsigned int j=0; j<this->nSpecies-3; j=j+4)
        {
            __m256d Cv = _mm256_loadu_pd(&c[j+0]);
            __m256d JcRowi0v = _mm256_loadu_pd(&JcRowi0[j+0]);
            ddNi0dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(JcRowi0v,Cv),alphavMv,ddNi0dtByVdTv);
            __m256d JcRowi1v = _mm256_loadu_pd(&JcRowi1[j+0]);
            ddNi1dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(JcRowi1v,Cv),alphavMv,ddNi1dtByVdTv);
            __m256d JcRowi2v = _mm256_loadu_pd(&JcRowi2[j+0]);
            ddNi2dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(JcRowi2v,Cv),alphavMv,ddNi2dtByVdTv);                            
        }
        {
            unsigned int j = nSpecies-3;
            ddNi0dtByVdT -= JcRowi0[j+0]*c[j+0]*alphavM;
            ddNi0dtByVdT -= JcRowi0[j+1]*c[j+1]*alphavM;
            ddNi0dtByVdT -= JcRowi0[j+2]*c[j+2]*alphavM;
            ddNi1dtByVdT -= JcRowi1[j+0]*c[j+0]*alphavM;
            ddNi1dtByVdT -= JcRowi1[j+1]*c[j+1]*alphavM;
            ddNi1dtByVdT -= JcRowi1[j+2]*c[j+2]*alphavM;
            ddNi2dtByVdT -= JcRowi2[j+0]*c[j+0]*alphavM;
            ddNi2dtByVdT -= JcRowi2[j+1]*c[j+1]*alphavM; 
            ddNi2dtByVdT -= JcRowi2[j+2]*c[j+2]*alphavM;
        }
        ddNi0dtByVdT -= hsum4(ddNi0dtByVdTv);
        ddNi1dtByVdT -= hsum4(ddNi1dtByVdTv);
        ddNi2dtByVdT -= hsum4(ddNi2dtByVdTv);
        J[i0*(alignN) + nSpecies] = Wi0ByrhoM_*ddNi0dtByVdT + alphavM*dYi0dt;
        J[i1*(alignN) + nSpecies] = Wi1ByrhoM_*ddNi1dtByVdT + alphavM*dYi1dt;
        J[i2*(alignN) + nSpecies] = Wi2ByrhoM_*ddNi2dtByVdT + alphavM*dYi2dt;
    }
}



void 
OptReaction::ddYdtdTP_Vec_2
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ WiByrhoM,
    const double* __restrict__ c,  
    double* __restrict__ dPhidt, 
    double* __restrict__ Jac
) const noexcept
{

    double alphavM = this->invT;
    __m256d alphavMv = _mm256_set1_pd(this->invT);
    for (unsigned int i=0; i<this->nSpecies-2; i=i+4)
    {
        unsigned int i0 = i;
        unsigned int i1 = i+1;
        unsigned int i2 = i+2;
        unsigned int i3 = i+3;

        __m256d Wi0ByrhoM_v = _mm256_loadu_pd(&WiByrhoM[i0]);
        __m256d dYTpdtv = _mm256_loadu_pd(&dPhidt[i0]);   
        dYTpdtv = _mm256_mul_pd(dYTpdtv,Wi0ByrhoM_v);   
        _mm256_storeu_pd(&dPhidt[i0],dYTpdtv);
        __m256d dYi0dtv = _mm256_loadu_pd(&dPhidt[i0]);

        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i0)*(alignN)];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i1)*(alignN)];
        const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i2)*(alignN)];
        const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i3)*(alignN)];

        double ddNi0dtByVdT = JcRowi0[nSpecies];
        double ddNi1dtByVdT = JcRowi1[nSpecies];
        double ddNi2dtByVdT = JcRowi2[nSpecies];
        double ddNi3dtByVdT = JcRowi3[nSpecies];
        __m256d sum0 = _mm256_setzero_pd();
        __m256d sum1 = _mm256_setzero_pd();
        __m256d sum2 = _mm256_setzero_pd();
        __m256d sum3 = _mm256_setzero_pd();
        for (unsigned int j=0; j<this->nSpecies-2; j=j+4)
        {
            __m256d Cv = _mm256_loadu_pd(&c[j+0]);
            sum0 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi0[j]),Cv),alphavMv,sum0);
            sum1 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi1[j]),Cv),alphavMv,sum1);
            sum2 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi2[j]),Cv),alphavMv,sum2);
            sum3 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi3[j]),Cv),alphavMv,sum3);
        }
        ddNi0dtByVdT = ddNi0dtByVdT - (hsum4(sum0));
        ddNi1dtByVdT = ddNi1dtByVdT - (hsum4(sum1));
        ddNi2dtByVdT = ddNi2dtByVdT - (hsum4(sum2));
        ddNi3dtByVdT = ddNi3dtByVdT - (hsum4(sum3));
        {
            unsigned int j0 = nSpecies-2;
            unsigned int j1 = nSpecies-1;
            ddNi0dtByVdT -= JcRowi0[j0]*c[j0]*alphavM;
            ddNi1dtByVdT -= JcRowi1[j0]*c[j0]*alphavM;
            ddNi2dtByVdT -= JcRowi2[j0]*c[j0]*alphavM;
            ddNi3dtByVdT -= JcRowi3[j0]*c[j0]*alphavM;
            ddNi0dtByVdT -= JcRowi0[j1]*c[j1]*alphavM;
            ddNi1dtByVdT -= JcRowi1[j1]*c[j1]*alphavM;
            ddNi2dtByVdT -= JcRowi2[j1]*c[j1]*alphavM;
            ddNi3dtByVdT -= JcRowi3[j1]*c[j1]*alphavM;
        }
        __m256d ddNi0dtByVdTv = _mm256_setr_pd(ddNi0dtByVdT,ddNi1dtByVdT,ddNi2dtByVdT,ddNi3dtByVdT);
        __m256d result = _mm256_fmadd_pd(Wi0ByrhoM_v,ddNi0dtByVdTv,_mm256_mul_pd(alphavMv,dYi0dtv));
        Jac[i0*(alignN) + nSpecies] = get_elem0(result);
        Jac[i1*(alignN) + nSpecies] = get_elem1(result);
        Jac[i2*(alignN) + nSpecies] = get_elem2(result);
        Jac[i3*(alignN) + nSpecies] = get_elem3(result);
    }

    {
        unsigned int i0 = this->nSpecies-2;
        unsigned int i1 = this->nSpecies-1;


        const double Wi0ByrhoM_ = WiByrhoM[i0];
        const double Wi1ByrhoM_ = WiByrhoM[i1];

        dPhidt[i0] = dPhidt[i0]*Wi0ByrhoM_;
        dPhidt[i1] = dPhidt[i1]*Wi1ByrhoM_;

        double dYi0dt = dPhidt[i0];
        double dYi1dt = dPhidt[i1];

        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i0)*(alignN)];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i1)*(alignN)];

        double ddNi0dtByVdT = JcRowi0[nSpecies];
        double ddNi1dtByVdT = JcRowi1[nSpecies];
        __m256d ddNi0dtByVdTv = _mm256_setzero_pd();
        __m256d ddNi1dtByVdTv = _mm256_setzero_pd();
        for (unsigned int j=0; j<this->nSpecies-2; j=j+4)
        {
            __m256d Cv = _mm256_loadu_pd(&c[j+0]);
            __m256d JcRowi0v = _mm256_loadu_pd(&JcRowi0[j+0]);
            ddNi0dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(JcRowi0v,Cv),alphavMv,ddNi0dtByVdTv);
            __m256d JcRowi1v = _mm256_loadu_pd(&JcRowi1[j+0]);
            ddNi1dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(JcRowi1v,Cv),alphavMv,ddNi1dtByVdTv);
        }
        {
            unsigned int j = nSpecies-2;
            ddNi0dtByVdT -= ddNdtByVdcT[i0*(alignN) + j+0]*c[j+0]*alphavM;
            ddNi0dtByVdT -= ddNdtByVdcT[i0*(alignN) + j+1]*c[j+1]*alphavM;
            ddNi1dtByVdT -= ddNdtByVdcT[i1*(alignN) + j+0]*c[j+0]*alphavM;
            ddNi1dtByVdT -= ddNdtByVdcT[i1*(alignN) + j+1]*c[j+1]*alphavM;
        }
        ddNi0dtByVdT -= hsum4(ddNi0dtByVdTv);
        ddNi1dtByVdT -= hsum4(ddNi1dtByVdTv);

        Jac[i0*(alignN) + nSpecies] = Wi0ByrhoM_*ddNi0dtByVdT + alphavM*dYi0dt;
        Jac[i1*(alignN) + nSpecies] = Wi1ByrhoM_*ddNi1dtByVdT + alphavM*dYi1dt;
    }
}


void 
OptReaction::ddYdtdTP_Vec_1
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ WiByrhoM,
    const double* __restrict__ c,  
    double* __restrict__ dPhidt, 
    double* __restrict__ Jac
) const noexcept
{

    double alphavM = this->invT;
    __m256d alphavMv = _mm256_set1_pd(this->invT);
    for (unsigned int i=0; i<this->nSpecies-1; i=i+4)
    {
        unsigned int i0 = i;
        unsigned int i1 = i+1;
        unsigned int i2 = i+2;
        unsigned int i3 = i+3;

        __m256d Wi0ByrhoM_v = _mm256_loadu_pd(&WiByrhoM[i0]);
        __m256d dYTpdtv = _mm256_loadu_pd(&dPhidt[i0]);   
        dYTpdtv = _mm256_mul_pd(dYTpdtv,Wi0ByrhoM_v);   
        _mm256_storeu_pd(&dPhidt[i0],dYTpdtv);
        __m256d dYi0dtv = _mm256_loadu_pd(&dPhidt[i0]);

        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i0)*(alignN)];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i1)*(alignN)];
        const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i2)*(alignN)];
        const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i3)*(alignN)];

        double ddNi0dtByVdT = JcRowi0[nSpecies];
        double ddNi1dtByVdT = JcRowi1[nSpecies];
        double ddNi2dtByVdT = JcRowi2[nSpecies];
        double ddNi3dtByVdT = JcRowi3[nSpecies];

        __m256d sum0 = _mm256_setzero_pd();
        __m256d sum1 = _mm256_setzero_pd();
        __m256d sum2 = _mm256_setzero_pd();
        __m256d sum3 = _mm256_setzero_pd();

        for (unsigned int j=0; j<this->nSpecies-1; j=j+4)
        {
            __m256d Cv = _mm256_loadu_pd(&c[j+0]);
            sum0 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi0[j]),Cv),alphavMv,sum0);
            sum1 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi1[j]),Cv),alphavMv,sum1);
            sum2 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi2[j]),Cv),alphavMv,sum2);
            sum3 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi3[j]),Cv),alphavMv,sum3);
        }

        __m256d sum_all = hsum4x4(sum0,sum1,sum2,sum3);

        __m256d ddNi0dtByVdTv = _mm256_setr_pd(ddNi0dtByVdT,ddNi1dtByVdT,ddNi2dtByVdT,ddNi3dtByVdT);
        ddNi0dtByVdTv = _mm256_sub_pd(ddNi0dtByVdTv,sum_all);
        {

            unsigned int j0 = nSpecies-1;
            __m256d JcRowiv = _mm256_setr_pd(-JcRowi0[j0],-JcRowi1[j0],-JcRowi2[j0],-JcRowi3[j0]);
            __m256d cv = _mm256_set1_pd(c[j0]);
            
            ddNi0dtByVdTv = _mm256_fmadd_pd(JcRowiv,_mm256_mul_pd(cv,alphavMv),ddNi0dtByVdTv);

        }
        __m256d result = _mm256_fmadd_pd(Wi0ByrhoM_v,ddNi0dtByVdTv,_mm256_mul_pd(alphavMv,dYi0dtv));
        Jac[i0*(alignN) + nSpecies] = get_elem0(result);
        Jac[i1*(alignN) + nSpecies] = get_elem1(result);
        Jac[i2*(alignN) + nSpecies] = get_elem2(result);
        Jac[i3*(alignN) + nSpecies] = get_elem3(result);         
    }
    {
        unsigned int i0 = this->nSpecies-1;
        const double Wi0ByrhoM_ = WiByrhoM[i0];
        dPhidt[i0] = dPhidt[i0]*Wi0ByrhoM_;
        double dYi0dt = dPhidt[i0];
        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[i0*(alignN)];
        double ddNi0dtByVdT = JcRowi0[nSpecies];
        __m256d ddNi0dtByVdTv = _mm256_setzero_pd();
        for (unsigned int j=0; j<this->nSpecies-1; j=j+4)
        {
            __m256d JcRowi0v = _mm256_loadu_pd(&JcRowi0[j+0]);
            __m256d Cv = _mm256_loadu_pd(&c[j+0]);
            ddNi0dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(JcRowi0v,Cv),alphavMv,ddNi0dtByVdTv);
        }
        {
            unsigned int j = this->nSpecies-1;
            ddNi0dtByVdT = ddNi0dtByVdT - JcRowi0[j]*c[j]*alphavM - hsum4(ddNi0dtByVdTv);
        }
        
        Jac[i0*(alignN) + nSpecies] = Wi0ByrhoM_*ddNi0dtByVdT + alphavM*dYi0dt;
    }
}


void 
OptReaction::ddYdtdTP_Vec_0
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ WiByrhoM,
    const double* __restrict__ c,  
    double* __restrict__ dPhidt, 
    double* __restrict__ Jac
) const noexcept
{

    __m256d alphavMv = _mm256_set1_pd(this->invT);
    for (unsigned int i=0; i<this->nSpecies; i=i+4)
    {
        unsigned int i0 = i;
        unsigned int i1 = i+1;
        unsigned int i2 = i+2;
        unsigned int i3 = i+3;

        __m256d Wi0ByrhoM_v = _mm256_loadu_pd(&WiByrhoM[i0]);
        __m256d dYTpdtv = _mm256_loadu_pd(&dPhidt[i0]);   
        dYTpdtv = _mm256_mul_pd(dYTpdtv,Wi0ByrhoM_v);   
        _mm256_storeu_pd(&dPhidt[i0],dYTpdtv);
        __m256d dYi0dtv = _mm256_loadu_pd(&dPhidt[i0]);

        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i0)*(alignN)];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i1)*(alignN)];
        const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i2)*(alignN)];
        const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i3)*(alignN)];

        double ddNi0dtByVdT = JcRowi0[nSpecies];
        double ddNi1dtByVdT = JcRowi1[nSpecies];
        double ddNi2dtByVdT = JcRowi2[nSpecies];
        double ddNi3dtByVdT = JcRowi3[nSpecies];
        __m256d sum0 = _mm256_setzero_pd();
        __m256d sum1 = _mm256_setzero_pd();
        __m256d sum2 = _mm256_setzero_pd();
        __m256d sum3 = _mm256_setzero_pd();
        for (unsigned int j=0; j<this->nSpecies; j=j+4)
        {
            __m256d Cv = _mm256_loadu_pd(&c[j+0]);
            sum0 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi0[j]),Cv),alphavMv,sum0);
            sum1 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi1[j]),Cv),alphavMv,sum1);
            sum2 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi2[j]),Cv),alphavMv,sum2);
            sum3 =_mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(&JcRowi3[j]),Cv),alphavMv,sum3);
        }

        ddNi0dtByVdT = ddNi0dtByVdT - (hsum4(sum0));
        ddNi1dtByVdT = ddNi1dtByVdT - (hsum4(sum1));
        ddNi2dtByVdT = ddNi2dtByVdT - (hsum4(sum2));
        ddNi3dtByVdT = ddNi3dtByVdT - (hsum4(sum3));

        __m256d ddNi0dtByVdTv = _mm256_setr_pd(ddNi0dtByVdT,ddNi1dtByVdT,ddNi2dtByVdT,ddNi3dtByVdT);
        __m256d result = _mm256_fmadd_pd(Wi0ByrhoM_v,ddNi0dtByVdTv,_mm256_mul_pd(alphavMv,dYi0dtv));

        Jac[i0*(alignN) + nSpecies] = get_elem0(result);
        Jac[i1*(alignN) + nSpecies] = get_elem1(result);
        Jac[i2*(alignN) + nSpecies] = get_elem2(result);
        Jac[i3*(alignN) + nSpecies] = get_elem3(result);
    }
}



void 
OptReaction::ddTdtdYT_Vec_0
(
    const double* __restrict__ Cp,
    const double* __restrict__ dCpdT,
    const double* __restrict__ Ha,      
    double* __restrict__ dPhidt,  
    double* __restrict__ Jac
) const noexcept
{
    const double& CpM = Cp[this->nSpecies];
    const double& dCpMdT = dCpdT[this->nSpecies];
    const double invCpM = 1.0/CpM;
    __m256d dTdtv = _mm256_setzero_pd();
    for (unsigned int i=0; i<this->nSpecies; i=i+4)
    {
        dTdtv = _mm256_fmadd_pd(-_mm256_loadu_pd(&Ha[i+0]),_mm256_loadu_pd(&dPhidt[i+0]),dTdtv);
    }

    double dTdt = hsum4(dTdtv);    
    dTdt *= invCpM;
    dPhidt[this->nSpecies] = dTdt;
    double& ddTdtdT = Jac[this->nSpecies *(alignN)+ this->nSpecies];
    ddTdtdT = 0;
    __m256d ddTdtdTv = _mm256_setzero_pd();
    for (unsigned int i=0; i<this->nSpecies; i=i+4)
    {
        __m256d ddTdtdYiv = _mm256_set1_pd(0.0);
        for (unsigned int j=0; j<this->nSpecies; j=j+4)
        {
            __m256d ddYj0dtdYiv = _mm256_loadu_pd(&Jac[(j+0) *(alignN)+ (i+0)]);
            __m256d ddYj1dtdYiv = _mm256_loadu_pd(&Jac[(j+1) *(alignN)+ (i+0)]);
            __m256d ddYj2dtdYiv = _mm256_loadu_pd(&Jac[(j+2) *(alignN)+ (i+0)]);
            __m256d ddYj3dtdYiv = _mm256_loadu_pd(&Jac[(j+3) *(alignN)+ (i+0)]);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj0dtdYiv,_mm256_set1_pd(-Ha[j+0]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj1dtdYiv,_mm256_set1_pd(-Ha[j+1]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj2dtdYiv,_mm256_set1_pd(-Ha[j+2]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj3dtdYiv,_mm256_set1_pd(-Ha[j+3]),ddTdtdYiv);                     
        }
        ddTdtdYiv = _mm256_fmadd_pd(_mm256_loadu_pd(&Cp[i+0]),_mm256_set1_pd(-dTdt),ddTdtdYiv);
        ddTdtdYiv =_mm256_mul_pd(ddTdtdYiv,_mm256_set1_pd(invCpM));
        _mm256_storeu_pd(&Jac[this->nSpecies *(alignN)+ i+0],ddTdtdYiv);
        __m256d dYi0dtv = _mm256_loadu_pd(&dPhidt[i+0]);
        const double ddYi0dtdT = Jac[(i+0) *(alignN)+ this->nSpecies];
        const double ddYi1dtdT = Jac[(i+1) *(alignN)+ this->nSpecies];
        const double ddYi2dtdT = Jac[(i+2) *(alignN)+ this->nSpecies];
        const double ddYi3dtdT = Jac[(i+3) *(alignN)+ this->nSpecies];
        __m256d ddYi0dtdTv = _mm256_setr_pd(ddYi0dtdT,ddYi1dtdT,ddYi2dtdT,ddYi3dtdT);
        ddTdtdTv = _mm256_fmadd_pd(dYi0dtv,_mm256_loadu_pd(&Cp[i+0]),ddTdtdTv);
        ddTdtdTv = _mm256_fmadd_pd(ddYi0dtdTv,_mm256_loadu_pd(&Ha[i+0]),ddTdtdTv);
    }

    ddTdtdT = ddTdtdT - hsum4(ddTdtdTv);
    ddTdtdT -= dTdt*dCpMdT; 
    ddTdtdT *= invCpM;
}

void 
OptReaction::ddTdtdYT_Vec_1
(
    const double* __restrict__ Cp,
    const double* __restrict__ dCpdT,
    const double* __restrict__ Ha,      
    double* __restrict__ dPhidt,  
    double* __restrict__ Jac
) const noexcept
{
    const double& CpM = Cp[this->nSpecies];
    const double& dCpMdT = dCpdT[this->nSpecies];
    const double invCpM = 1.0/CpM;
    __m256d dTdtv = _mm256_setzero_pd();
    for (unsigned int i=0; i<this->nSpecies-1; i=i+4)
    {
        dTdtv = _mm256_fmadd_pd(-_mm256_loadu_pd(&Ha[i+0]),_mm256_loadu_pd(&dPhidt[i+0]),dTdtv);
    }
    double dTdt = hsum4(dTdtv);
    {
        unsigned int i = this->nSpecies-1;
        dTdt -= dPhidt[i+0]*Ha[i+0];       
    }
    dTdt *= invCpM;

    dPhidt[this->nSpecies] = dTdt;

    double& ddTdtdT = Jac[this->nSpecies *(alignN)+ this->nSpecies];
    ddTdtdT = 0;

    __m256d ddTdtdTv = _mm256_setzero_pd();
    for (unsigned int i=0; i<this->nSpecies-1; i=i+4)
    {
        __m256d ddTdtdYiv = _mm256_set1_pd(0.0);
        for (unsigned int j=0; j<this->nSpecies-1; j=j+4)
        {
            __m256d ddYj0dtdYiv = _mm256_loadu_pd(&Jac[(j+0) *(alignN)+ (i+0)]);
            __m256d ddYj1dtdYiv = _mm256_loadu_pd(&Jac[(j+1) *(alignN)+ (i+0)]);
            __m256d ddYj2dtdYiv = _mm256_loadu_pd(&Jac[(j+2) *(alignN)+ (i+0)]);
            __m256d ddYj3dtdYiv = _mm256_loadu_pd(&Jac[(j+3) *(alignN)+ (i+0)]);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj0dtdYiv,_mm256_set1_pd(-Ha[j+0]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj1dtdYiv,_mm256_set1_pd(-Ha[j+1]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj2dtdYiv,_mm256_set1_pd(-Ha[j+2]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj3dtdYiv,_mm256_set1_pd(-Ha[j+3]),ddTdtdYiv);                     
        }
        {
            unsigned int j = this->nSpecies-1;
            __m256d ddYj0dtdYi = _mm256_loadu_pd(&Jac[(j+0) *(alignN)+ (i+0)]);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj0dtdYi,_mm256_set1_pd(-Ha[j+0]),ddTdtdYiv);             
        }
        ddTdtdYiv = _mm256_fmadd_pd(_mm256_loadu_pd(&Cp[i+0]),_mm256_set1_pd(-dTdt),ddTdtdYiv);
        ddTdtdYiv =_mm256_mul_pd(ddTdtdYiv,_mm256_set1_pd(invCpM));
        _mm256_storeu_pd(&Jac[this->nSpecies *(alignN)+ i+0],ddTdtdYiv);
        __m256d dYi0dtv = _mm256_loadu_pd(&dPhidt[i+0]);

        const double ddYi0dtdT = Jac[(i+0) *(alignN)+ this->nSpecies];
        const double ddYi1dtdT = Jac[(i+1) *(alignN)+ this->nSpecies];
        const double ddYi2dtdT = Jac[(i+2) *(alignN)+ this->nSpecies];
        const double ddYi3dtdT = Jac[(i+3) *(alignN)+ this->nSpecies];
        __m256d ddYi0dtdTv = _mm256_setr_pd(ddYi0dtdT,ddYi1dtdT,ddYi2dtdT,ddYi3dtdT);

        ddTdtdTv = _mm256_fmadd_pd(dYi0dtv,_mm256_loadu_pd(&Cp[i+0]),ddTdtdTv);
        ddTdtdTv = _mm256_fmadd_pd(ddYi0dtdTv,_mm256_loadu_pd(&Ha[i+0]),ddTdtdTv);
    }  
    ddTdtdT = ddTdtdT - hsum4(ddTdtdTv);

    {
        unsigned int i = this->nSpecies-1;
        double& ddTdtdYi = Jac[this->nSpecies *(alignN)+ i];
        ddTdtdYi = 0;
        __m256d ddTdtdYiv = _mm256_setzero_pd();
        for (unsigned int j=0; j<this->nSpecies-1; j=j+4)
        {
            __m256d Av = _mm256_setr_pd(Jac[(j+0) *(alignN)+ i],Jac[(j+1) *(alignN)+ i],Jac[(j+2) *(alignN)+ i],Jac[(j+3) *(alignN)+ i]);
            __m256d Hav = _mm256_loadu_pd(&Ha[j+0]);
            ddTdtdYiv = _mm256_fmadd_pd(Av,Hav,ddTdtdYiv);
        }
        {
            ddTdtdYi = ddTdtdYi - hsum4(ddTdtdYiv);
            unsigned int j = nSpecies-1;
            ddTdtdYi -= Jac[(j+0) *(alignN)+ i]*Ha[j+0];
        }

        ddTdtdYi -= Cp[i]*dTdt;
        ddTdtdYi *= invCpM;
        const double dYidt = dPhidt[i];
        const double ddYidtdT = Jac[i *(alignN)+ this->nSpecies];
        ddTdtdT -= dYidt*Cp[i] + ddYidtdT*Ha[i];
    }
    ddTdtdT -= dTdt*dCpMdT;
    ddTdtdT *= invCpM;   
}

void 
OptReaction::ddTdtdYT_Vec_2
(
    const double* __restrict__ Cp,
    const double* __restrict__ dCpdT,
    const double* __restrict__ Ha,      
    double* __restrict__ dPhidt,  
    double* __restrict__ Jac
) const noexcept
{

    const double& CpM = Cp[this->nSpecies];
    const double& dCpMdT = dCpdT[this->nSpecies];
    const double invCpM = 1.0/CpM;
    __m256d dTdtv = _mm256_setzero_pd();
    for (unsigned int i=0; i<this->nSpecies-2; i=i+4)
    {
        dTdtv = _mm256_fmadd_pd(-_mm256_loadu_pd(&Ha[i+0]),_mm256_loadu_pd(&dPhidt[i+0]),dTdtv);
    }
    double dTdt = hsum4(dTdtv);
    {
        unsigned int i0 = this->nSpecies-2;
        unsigned int i1 = this->nSpecies-1;
        dTdt -= dPhidt[i0]*Ha[i0];
        dTdt -= dPhidt[i1]*Ha[i1];
    }
    dTdt *= invCpM;
    dPhidt[this->nSpecies] = dTdt;
    double& ddTdtdT = Jac[this->nSpecies *(alignN)+ this->nSpecies];
    ddTdtdT = 0;
    __m256d ddTdtdTv = _mm256_setzero_pd();
    for (unsigned int i=0; i<this->nSpecies-2; i=i+4)
    {
        __m256d ddTdtdYiv = _mm256_set1_pd(0.0);
        for (unsigned int j=0; j<this->nSpecies-2; j=j+4)
        {
            __m256d ddYj0dtdYiv = _mm256_loadu_pd(&Jac[(j+0) *(alignN)+ (i+0)]);
            __m256d ddYj1dtdYiv = _mm256_loadu_pd(&Jac[(j+1) *(alignN)+ (i+0)]);
            __m256d ddYj2dtdYiv = _mm256_loadu_pd(&Jac[(j+2) *(alignN)+ (i+0)]);
            __m256d ddYj3dtdYiv = _mm256_loadu_pd(&Jac[(j+3) *(alignN)+ (i+0)]);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj0dtdYiv,_mm256_set1_pd(-Ha[j+0]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj1dtdYiv,_mm256_set1_pd(-Ha[j+1]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj2dtdYiv,_mm256_set1_pd(-Ha[j+2]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj3dtdYiv,_mm256_set1_pd(-Ha[j+3]),ddTdtdYiv);                     
        }
        {
            unsigned int j = this->nSpecies-2;
            __m256d ddYj0dtdYi = _mm256_loadu_pd(&Jac[(j+0) *(alignN)+ (i+0)]);
            __m256d ddYj1dtdYi = _mm256_loadu_pd(&Jac[(j+1) *(alignN)+ (i+0)]);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj0dtdYi,_mm256_set1_pd(-Ha[j+0]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj1dtdYi,_mm256_set1_pd(-Ha[j+1]),ddTdtdYiv);
        }
        ddTdtdYiv = _mm256_fmadd_pd(_mm256_loadu_pd(&Cp[i+0]),_mm256_set1_pd(-dTdt),ddTdtdYiv);
        ddTdtdYiv =_mm256_mul_pd(ddTdtdYiv,_mm256_set1_pd(invCpM));
        _mm256_storeu_pd(&Jac[this->nSpecies *(alignN)+ i+0],ddTdtdYiv);
        __m256d dYi0dtv = _mm256_loadu_pd(&dPhidt[i+0]);

        const double ddYi0dtdT = Jac[(i+0) *(alignN)+ this->nSpecies];
        const double ddYi1dtdT = Jac[(i+1) *(alignN)+ this->nSpecies];
        const double ddYi2dtdT = Jac[(i+2) *(alignN)+ this->nSpecies];
        const double ddYi3dtdT = Jac[(i+3) *(alignN)+ this->nSpecies];
        __m256d ddYi0dtdTv = _mm256_setr_pd(ddYi0dtdT,ddYi1dtdT,ddYi2dtdT,ddYi3dtdT);
        ddTdtdTv = _mm256_fmadd_pd(dYi0dtv,_mm256_loadu_pd(&Cp[i+0]),ddTdtdTv);
        ddTdtdTv = _mm256_fmadd_pd(ddYi0dtdTv,_mm256_loadu_pd(&Ha[i+0]),ddTdtdTv);
    }
    ddTdtdT = ddTdtdT - hsum4(ddTdtdTv);
    {
        unsigned int i0 = this->nSpecies-2;
        unsigned int i1 = this->nSpecies-1;
        double& ddTdtdYi0 = Jac[this->nSpecies *(alignN)+ i0];
        double& ddTdtdYi1 = Jac[this->nSpecies *(alignN)+ i1];        
        ddTdtdYi0 = 0;
        ddTdtdYi1 = 0;   
        __m256d result0 = _mm256_setzero_pd();
        __m256d result1 = _mm256_setzero_pd();
        for (unsigned int j=0; j<this->nSpecies-2; j=j+4)
        {
            __m256d A0 = _mm256_setr_pd(Jac[(j+0) *(alignN)+ i0],Jac[(j+1) *(alignN)+ i0],Jac[(j+2) *(alignN)+ i0],Jac[(j+3) *(alignN)+ i0]);
            __m256d A1 = _mm256_setr_pd(Jac[(j+0) *(alignN)+ i1],Jac[(j+1) *(alignN)+ i1],Jac[(j+2) *(alignN)+ i1],Jac[(j+3) *(alignN)+ i1]);
            __m256d Hav = _mm256_loadu_pd(&Ha[j+0]);
            result0 = _mm256_fmadd_pd(A0,Hav,result0);
            result1 = _mm256_fmadd_pd(A1,Hav,result1);
        }
        {
            ddTdtdYi0 = ddTdtdYi0 - hsum4(result0);
            ddTdtdYi1 = ddTdtdYi1 - hsum4(result1);
            unsigned int j = nSpecies-2;
            ddTdtdYi0 -= Jac[(j+0) *(alignN)+ i0]*Ha[j+0];
            ddTdtdYi0 -= Jac[(j+1) *(alignN)+ i0]*Ha[j+1];
            ddTdtdYi1 -= Jac[(j+0) *(alignN)+ i1]*Ha[j+0];
            ddTdtdYi1 -= Jac[(j+1) *(alignN)+ i1]*Ha[j+1];
        }
        ddTdtdYi0 -= Cp[i0]*dTdt;
        ddTdtdYi1 -= Cp[i1]*dTdt;

        ddTdtdYi0 = ddTdtdYi0*invCpM;
        ddTdtdYi1 = ddTdtdYi1*invCpM;

        const double dYi0dt = dPhidt[i0];
        const double dYi1dt = dPhidt[i1];

        const double ddYi0dtdT = Jac[i0 *(alignN)+ this->nSpecies];
        const double ddYi1dtdT = Jac[i1 *(alignN)+ this->nSpecies];

        ddTdtdT -= dYi0dt*Cp[i0] + ddYi0dtdT*Ha[i0];
        ddTdtdT -= dYi1dt*Cp[i1] + ddYi1dtdT*Ha[i1];        
    }
    ddTdtdT -= dTdt*dCpMdT;
    ddTdtdT = ddTdtdT*invCpM;
}

void 
OptReaction::ddTdtdYT_Vec_3
(
    const double* __restrict__ Cp,
    const double* __restrict__ dCpdT,
    const double* __restrict__ Ha,      
    double* __restrict__ dPhidt,  
    double* __restrict__ Jac
) const noexcept
{

    const double& CpM = Cp[this->nSpecies];
    const double& dCpMdT = dCpdT[this->nSpecies];
    const double invCpM = 1.0/CpM;
    __m256d dTdtv = _mm256_setzero_pd();
    for (unsigned int i=0; i<this->nSpecies-3; i=i+4)
    {
        dTdtv = _mm256_fmadd_pd(-_mm256_loadu_pd(&Ha[i+0]),_mm256_loadu_pd(&dPhidt[i+0]),dTdtv);
    }
    double dTdt = hsum4(dTdtv);
    {
        unsigned int i0 = this->nSpecies-3;
        unsigned int i1 = this->nSpecies-2;
        unsigned int i2 = this->nSpecies-1;
        dTdt -= dPhidt[i0]*Ha[i0];
        dTdt -= dPhidt[i1]*Ha[i1];
        dTdt -= dPhidt[i2]*Ha[i2];
    }
    dTdt *= invCpM;
    dPhidt[this->nSpecies] = dTdt;
    double& ddTdtdT = Jac[this->nSpecies *(alignN)+ this->nSpecies];
    ddTdtdT = 0;
    __m256d ddTdtdTv = _mm256_setzero_pd();
    for (unsigned int i=0; i<this->nSpecies-3; i=i+4)
    {
        __m256d ddTdtdYiv = _mm256_set1_pd(0.0);
        for (unsigned int j=0; j<this->nSpecies-3; j=j+4)
        {
            __m256d ddYj0dtdYiv = _mm256_loadu_pd(&Jac[(j+0) *(alignN)+ (i+0)]);
            __m256d ddYj1dtdYiv = _mm256_loadu_pd(&Jac[(j+1) *(alignN)+ (i+0)]);
            __m256d ddYj2dtdYiv = _mm256_loadu_pd(&Jac[(j+2) *(alignN)+ (i+0)]);
            __m256d ddYj3dtdYiv = _mm256_loadu_pd(&Jac[(j+3) *(alignN)+ (i+0)]);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj0dtdYiv,_mm256_set1_pd(-Ha[j+0]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj1dtdYiv,_mm256_set1_pd(-Ha[j+1]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj2dtdYiv,_mm256_set1_pd(-Ha[j+2]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj3dtdYiv,_mm256_set1_pd(-Ha[j+3]),ddTdtdYiv);                     
        }
        {
            unsigned int j = this->nSpecies-3;
            __m256d ddYj0dtdYi = _mm256_loadu_pd(&Jac[(j+0) *(alignN)+ (i+0)]);
            __m256d ddYj1dtdYi = _mm256_loadu_pd(&Jac[(j+1) *(alignN)+ (i+0)]);
            __m256d ddYj2dtdYi = _mm256_loadu_pd(&Jac[(j+2) *(alignN)+ (i+0)]);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj0dtdYi,_mm256_set1_pd(-Ha[j+0]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj1dtdYi,_mm256_set1_pd(-Ha[j+1]),ddTdtdYiv);
            ddTdtdYiv = _mm256_fmadd_pd(ddYj2dtdYi,_mm256_set1_pd(-Ha[j+2]),ddTdtdYiv);
        }
        ddTdtdYiv = _mm256_fmadd_pd(_mm256_loadu_pd(&Cp[i+0]),_mm256_set1_pd(-dTdt),ddTdtdYiv);
        ddTdtdYiv =_mm256_mul_pd(ddTdtdYiv,_mm256_set1_pd(invCpM));
        _mm256_storeu_pd(&Jac[this->nSpecies *(alignN)+ i+0],ddTdtdYiv);
        __m256d dYi0dtv = _mm256_loadu_pd(&dPhidt[i+0]);

        const double ddYi0dtdT = Jac[(i+0) *(alignN)+ this->nSpecies];
        const double ddYi1dtdT = Jac[(i+1) *(alignN)+ this->nSpecies];
        const double ddYi2dtdT = Jac[(i+2) *(alignN)+ this->nSpecies];
        const double ddYi3dtdT = Jac[(i+3) *(alignN)+ this->nSpecies];
        __m256d ddYi0dtdTv = _mm256_setr_pd(ddYi0dtdT,ddYi1dtdT,ddYi2dtdT,ddYi3dtdT);

        ddTdtdTv = _mm256_fmadd_pd(dYi0dtv,_mm256_loadu_pd(&Cp[i+0]),ddTdtdTv);
        ddTdtdTv = _mm256_fmadd_pd(ddYi0dtdTv,_mm256_loadu_pd(&Ha[i+0]),ddTdtdTv);
    }
   
    ddTdtdT = ddTdtdT - hsum4(ddTdtdTv);
    {
        unsigned int i0 = this->nSpecies-3;
        unsigned int i1 = this->nSpecies-2;
        unsigned int i2 = this->nSpecies-1;

        double& ddTdtdYi0 = Jac[this->nSpecies *(alignN)+ i0];
        double& ddTdtdYi1 = Jac[this->nSpecies *(alignN)+ i1];  
        double& ddTdtdYi2 = Jac[this->nSpecies *(alignN)+ i2];                
        ddTdtdYi0 = 0;
        ddTdtdYi1 = 0;   
        ddTdtdYi2 = 0; 
        __m256d result0 = _mm256_setzero_pd();
        __m256d result1 = _mm256_setzero_pd();
        __m256d result2 = _mm256_setzero_pd();
        for (unsigned int j=0; j<this->nSpecies-3; j=j+4)
        {
            __m256d A0 = _mm256_setr_pd(Jac[(j+0) *(alignN)+ i0],Jac[(j+1) *(alignN)+ i0],Jac[(j+2) *(alignN)+ i0],Jac[(j+3) *(alignN)+ i0]);
            __m256d A1 = _mm256_setr_pd(Jac[(j+0) *(alignN)+ i1],Jac[(j+1) *(alignN)+ i1],Jac[(j+2) *(alignN)+ i1],Jac[(j+3) *(alignN)+ i1]);
            __m256d A2 = _mm256_setr_pd(Jac[(j+0) *(alignN)+ i2],Jac[(j+1) *(alignN)+ i2],Jac[(j+2) *(alignN)+ i2],Jac[(j+3) *(alignN)+ i2]);
            __m256d Hav = _mm256_loadu_pd(&Ha[j+0]);
            result0 = _mm256_fmadd_pd(A0,Hav,result0);
            result1 = _mm256_fmadd_pd(A1,Hav,result1);
            result2 = _mm256_fmadd_pd(A2,Hav,result2);
        }
        {

            ddTdtdYi0 = ddTdtdYi0 - hsum4(result0);
            ddTdtdYi1 = ddTdtdYi1 - hsum4(result1);
            ddTdtdYi2 = ddTdtdYi2 - hsum4(result2);

            unsigned int j = nSpecies-3;
            ddTdtdYi0 -= Jac[(j+0) *(alignN)+ i0]*Ha[j+0];
            ddTdtdYi0 -= Jac[(j+1) *(alignN)+ i0]*Ha[j+1];
            ddTdtdYi0 -= Jac[(j+2) *(alignN)+ i0]*Ha[j+2];

            ddTdtdYi1 -= Jac[(j+0) *(alignN)+ i1]*Ha[j+0];
            ddTdtdYi1 -= Jac[(j+1) *(alignN)+ i1]*Ha[j+1];
            ddTdtdYi1 -= Jac[(j+2) *(alignN)+ i1]*Ha[j+2];

            ddTdtdYi2 -= Jac[(j+0) *(alignN)+ i2]*Ha[j+0];
            ddTdtdYi2 -= Jac[(j+1) *(alignN)+ i2]*Ha[j+1];
            ddTdtdYi2 -= Jac[(j+2) *(alignN)+ i2]*Ha[j+2];
        }
        ddTdtdYi0 -= Cp[i0]*dTdt;
        ddTdtdYi1 -= Cp[i1]*dTdt;
        ddTdtdYi2 -= Cp[i2]*dTdt;

        ddTdtdYi0 *= invCpM;
        ddTdtdYi1 *= invCpM;
        ddTdtdYi2 *= invCpM;

        const double dYi0dt = dPhidt[i0];
        const double dYi1dt = dPhidt[i1];
        const double dYi2dt = dPhidt[i2];

        const double ddYi0dtdT = Jac[i0 *(alignN)+ this->nSpecies];
        const double ddYi1dtdT = Jac[i1 *(alignN)+ this->nSpecies];
        const double ddYi2dtdT = Jac[i2 *(alignN)+ this->nSpecies];

        ddTdtdT -= dYi0dt*Cp[i0] + ddYi0dtdT*Ha[i0];
        ddTdtdT -= dYi1dt*Cp[i1] + ddYi1dtdT*Ha[i1];        
        ddTdtdT -= dYi2dt*Cp[i2] + ddYi2dtdT*Ha[i2]; 
    }
    ddTdtdT -= dTdt*dCpMdT;
    ddTdtdT = ddTdtdT*invCpM;   
}

void 
OptReaction::FastddYdtdY_Vec0
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ rhoMByRhoi,
    const double* __restrict__ WiByrhoM,  
    const double* __restrict__ dYTpdt, 
    double* __restrict__ Jac
) const noexcept
{ 
    __m256d rhoMv = _mm256_set1_pd(this->rhoM);
    for(unsigned int i=0; i<this->nSpecies; i=i+4)
    {
        const double Wi0ByrhoM_ = WiByrhoM[i+0];
        const double Wi1ByrhoM_ = WiByrhoM[i+1];
        const double Wi2ByrhoM_ = WiByrhoM[i+2];
        const double Wi3ByrhoM_ = WiByrhoM[i+3];
        __m256d Wi0ByrhoMv = _mm256_set1_pd(Wi0ByrhoM_);
        __m256d Wi1ByrhoMv = _mm256_set1_pd(Wi1ByrhoM_);
        __m256d Wi2ByrhoMv = _mm256_set1_pd(Wi2ByrhoM_);
        __m256d Wi3ByrhoMv = _mm256_set1_pd(Wi3ByrhoM_);

        double dYi0dt = dYTpdt[i+0]*WiByrhoM[i+0];
        double dYi1dt = dYTpdt[i+1]*WiByrhoM[i+1];
        double dYi2dt = dYTpdt[i+2]*WiByrhoM[i+2];
        double dYi3dt = dYTpdt[i+3]*WiByrhoM[i+3];
        __m256d dYi0dtv = _mm256_set1_pd(dYi0dt);
        __m256d dYi1dtv = _mm256_set1_pd(dYi1dt);
        __m256d dYi2dtv = _mm256_set1_pd(dYi2dt);
        __m256d dYi3dtv = _mm256_set1_pd(dYi3dt);        
        for (unsigned int j=0; j<this->nSpecies; j=j+4)
        {
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,_mm256_loadu_pd(&this->invW[j+0]));
            __m256d ddNi0dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+0)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+1)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi2dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+2)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi3dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+3)*(alignN)+j+0]),dCjdYj);
            __m256d rhoMvj_ = _mm256_loadu_pd(&rhoMByRhoi[j+0]);
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi0dtv));
            _mm256_storeu_pd(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi1dtv));
            _mm256_storeu_pd(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);
            __m256d ddYi2dtdYj = _mm256_fmadd_pd(Wi2ByrhoMv,ddNi2dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi2dtv));
            _mm256_storeu_pd(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYj);
            __m256d ddYi3dtdYj = _mm256_fmadd_pd(Wi3ByrhoMv,ddNi3dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi3dtv));
            _mm256_storeu_pd(&Jac[(i+3)*(alignN) + j+0],ddYi3dtdYj);
        }
    }
}

void 
OptReaction::FastddYdtdY_Vec1
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ rhoMByRhoi,
    const double* __restrict__ WiByrhoM,  
    const double* __restrict__ dYTpdt, 
    double* __restrict__ Jac
) const noexcept
{ 
    __m256d rhoMv = _mm256_set1_pd(this->rhoM);
    for(unsigned int i=0; i<this->nSpecies-1; i=i+4)
    {
        const double Wi0ByrhoM_ = WiByrhoM[i+0];
        const double Wi1ByrhoM_ = WiByrhoM[i+1];
        const double Wi2ByrhoM_ = WiByrhoM[i+2];
        const double Wi3ByrhoM_ = WiByrhoM[i+3];
        __m256d Wi0ByrhoMv = _mm256_set1_pd(Wi0ByrhoM_);
        __m256d Wi1ByrhoMv = _mm256_set1_pd(Wi1ByrhoM_);
        __m256d Wi2ByrhoMv = _mm256_set1_pd(Wi2ByrhoM_);
        __m256d Wi3ByrhoMv = _mm256_set1_pd(Wi3ByrhoM_);

        double dYi0dt = dYTpdt[i+0]*WiByrhoM[i+0];
        double dYi1dt = dYTpdt[i+1]*WiByrhoM[i+1];
        double dYi2dt = dYTpdt[i+2]*WiByrhoM[i+2];
        double dYi3dt = dYTpdt[i+3]*WiByrhoM[i+3];
        __m256d dYi0dtv = _mm256_set1_pd(dYi0dt);
        __m256d dYi1dtv = _mm256_set1_pd(dYi1dt);
        __m256d dYi2dtv = _mm256_set1_pd(dYi2dt);
        __m256d dYi3dtv = _mm256_set1_pd(dYi3dt);        
        for (unsigned int j=0; j<this->nSpecies-1; j=j+4)
        {
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,_mm256_loadu_pd(&this->invW[j+0]));
            __m256d ddNi0dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+0)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+1)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi2dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+2)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi3dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+3)*(alignN)+j+0]),dCjdYj);
            __m256d rhoMvj_ = _mm256_loadu_pd(&rhoMByRhoi[j+0]);
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi0dtv));
            _mm256_storeu_pd(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi1dtv));
            _mm256_storeu_pd(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);
            __m256d ddYi2dtdYj = _mm256_fmadd_pd(Wi2ByrhoMv,ddNi2dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi2dtv));
            _mm256_storeu_pd(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYj);
            __m256d ddYi3dtdYj = _mm256_fmadd_pd(Wi3ByrhoMv,ddNi3dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi3dtv));
            _mm256_storeu_pd(&Jac[(i+3)*(alignN) + j+0],ddYi3dtdYj);
        }
        {
            unsigned int j = this->nSpecies-1;
            const double dCj0dYj0 = rhoM*this->invW[j+0];
            const double ddNi0dtByVdYj0 = ddNdtByVdcT[(i+0)*(alignN) + j+0]*dCj0dYj0;
            const double ddNi1dtByVdYj0 = ddNdtByVdcT[(i+1)*(alignN) + j+0]*dCj0dYj0;
            const double ddNi2dtByVdYj0 = ddNdtByVdcT[(i+2)*(alignN) + j+0]*dCj0dYj0;
            const double ddNi3dtByVdYj0 = ddNdtByVdcT[(i+3)*(alignN) + j+0]*dCj0dYj0;
            double& ddYi0dtdYj0 = Jac[(i+0)*(alignN) + j+0];
            double& ddYi1dtdYj0 = Jac[(i+1)*(alignN) + j+0];
            double& ddYi2dtdYj0 = Jac[(i+2)*(alignN) + j+0];
            double& ddYi3dtdYj0 = Jac[(i+3)*(alignN) + j+0];
            ddYi0dtdYj0 = Wi0ByrhoM_*ddNi0dtByVdYj0 + rhoMByRhoi[j+0]*dYi0dt;
            ddYi1dtdYj0 = Wi1ByrhoM_*ddNi1dtByVdYj0 + rhoMByRhoi[j+0]*dYi1dt;
            ddYi2dtdYj0 = Wi2ByrhoM_*ddNi2dtByVdYj0 + rhoMByRhoi[j+0]*dYi2dt;
            ddYi3dtdYj0 = Wi3ByrhoM_*ddNi3dtByVdYj0 + rhoMByRhoi[j+0]*dYi3dt;
        }
    }

    {
        unsigned int i = this->nSpecies-1;
        const double Wi0ByrhoM_ = WiByrhoM[i+0];
        __m256d Wi0ByrhoMv = _mm256_set1_pd(Wi0ByrhoM_);
        double dYi0dt = dYTpdt[i+0]*WiByrhoM[i+0];
        __m256d dYi0dtv = _mm256_set1_pd(dYi0dt);

        for (unsigned int j=0; j<this->nSpecies-1; j=j+4)
        {
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,_mm256_loadu_pd(&this->invW[j+0]));
            __m256d rhoMvj_ = _mm256_loadu_pd(&rhoMByRhoi[j+0]);

            __m256d ddNi0dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+0)*(alignN) + j+0]),dCjdYj);
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi0dtv));
            _mm256_storeu_pd(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);
        }
        {
            unsigned int j = this->nSpecies-1;
            const double dCj0dYj0 = rhoM*this->invW[j+0];
            const double ddNi0dtByVdYj0 = ddNdtByVdcT[(i+0)*(alignN) + j+0]*dCj0dYj0;
            double& ddYi0dtdYj0 = Jac[(i+0)*(alignN) + j+0];
            ddYi0dtdYj0 = Wi0ByrhoM_*ddNi0dtByVdYj0 + rhoMByRhoi[j+0]*dYi0dt;
        }
    }
}

void 
OptReaction::FastddYdtdY_Vec2
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ rhoMByRhoi,
    const double* __restrict__ WiByrhoM,  
    const double* __restrict__ dYTpdt, 
    double* __restrict__ Jac
) const noexcept
{ 
    __m256d rhoMv = _mm256_set1_pd(this->rhoM);

    for(unsigned int i=0; i<this->nSpecies-2; i=i+4)
    {
        const double Wi0ByrhoM_ = WiByrhoM[i+0];
        const double Wi1ByrhoM_ = WiByrhoM[i+1];
        const double Wi2ByrhoM_ = WiByrhoM[i+2];
        const double Wi3ByrhoM_ = WiByrhoM[i+3];

        __m256d Wi0ByrhoMv = _mm256_set1_pd(Wi0ByrhoM_);
        __m256d Wi1ByrhoMv = _mm256_set1_pd(Wi1ByrhoM_);
        __m256d Wi2ByrhoMv = _mm256_set1_pd(Wi2ByrhoM_);
        __m256d Wi3ByrhoMv = _mm256_set1_pd(Wi3ByrhoM_);

        double dYi0dt = dYTpdt[i+0]*WiByrhoM[i+0];
        double dYi1dt = dYTpdt[i+1]*WiByrhoM[i+1];
        double dYi2dt = dYTpdt[i+2]*WiByrhoM[i+2];
        double dYi3dt = dYTpdt[i+3]*WiByrhoM[i+3];

        __m256d dYi0dtv = _mm256_set1_pd(dYi0dt);
        __m256d dYi1dtv = _mm256_set1_pd(dYi1dt);
        __m256d dYi2dtv = _mm256_set1_pd(dYi2dt);
        __m256d dYi3dtv = _mm256_set1_pd(dYi3dt);

        for (unsigned int j=0; j<this->nSpecies-2; j=j+4)
        {
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,_mm256_loadu_pd(&this->invW[j+0]));
            __m256d ddNi0dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+0)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+1)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi2dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+2)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi3dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+3)*(alignN)+j+0]),dCjdYj);
            __m256d rhoMvj_ = _mm256_loadu_pd(&rhoMByRhoi[j+0]);
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi0dtv));
            _mm256_storeu_pd(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi1dtv));
            _mm256_storeu_pd(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);
            __m256d ddYi2dtdYj = _mm256_fmadd_pd(Wi2ByrhoMv,ddNi2dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi2dtv));
            _mm256_storeu_pd(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYj);
            __m256d ddYi3dtdYj = _mm256_fmadd_pd(Wi3ByrhoMv,ddNi3dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi3dtv));
            _mm256_storeu_pd(&Jac[(i+3)*(alignN) + j+0],ddYi3dtdYj);
        }

        {
            unsigned int j = this->nSpecies-2;

            __m128d invWv = _mm_loadu_pd(&this->invW[j+0]);
            __m128d dCjdYjv =  _mm_mul_pd(invWv,_mm256_castpd256_pd128(rhoMv));
            __m128d ddNi0dtByVdYjv = _mm_mul_pd(dCjdYjv,_mm_loadu_pd(&ddNdtByVdcT[(i+0)*(alignN) + j+0]));
            __m128d ddNi1dtByVdYjv = _mm_mul_pd(dCjdYjv,_mm_loadu_pd(&ddNdtByVdcT[(i+1)*(alignN) + j+0]));
            __m128d ddNi2dtByVdYjv = _mm_mul_pd(dCjdYjv,_mm_loadu_pd(&ddNdtByVdcT[(i+2)*(alignN) + j+0]));
            __m128d ddNi3dtByVdYjv = _mm_mul_pd(dCjdYjv,_mm_loadu_pd(&ddNdtByVdcT[(i+3)*(alignN) + j+0]));
            __m128d rhoMByRhov = _mm_loadu_pd(&rhoMByRhoi[j+0]);
            __m128d ddYi0dtdYjv = _mm_mul_pd(_mm256_castpd256_pd128(Wi0ByrhoMv),ddNi0dtByVdYjv);
            __m128d ddYi1dtdYjv = _mm_mul_pd(_mm256_castpd256_pd128(Wi1ByrhoMv),ddNi1dtByVdYjv);
            __m128d ddYi2dtdYjv = _mm_mul_pd(_mm256_castpd256_pd128(Wi2ByrhoMv),ddNi2dtByVdYjv);
            __m128d ddYi3dtdYjv = _mm_mul_pd(_mm256_castpd256_pd128(Wi3ByrhoMv),ddNi3dtByVdYjv);
            ddYi0dtdYjv = _mm_fmadd_pd(rhoMByRhov,_mm256_castpd256_pd128(dYi0dtv),ddYi0dtdYjv);
            ddYi1dtdYjv = _mm_fmadd_pd(rhoMByRhov,_mm256_castpd256_pd128(dYi1dtv),ddYi1dtdYjv);
            ddYi2dtdYjv = _mm_fmadd_pd(rhoMByRhov,_mm256_castpd256_pd128(dYi2dtv),ddYi2dtdYjv);
            ddYi3dtdYjv = _mm_fmadd_pd(rhoMByRhov,_mm256_castpd256_pd128(dYi3dtv),ddYi3dtdYjv);
            _mm_storeu_pd(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYjv);
            _mm_storeu_pd(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYjv);
            _mm_storeu_pd(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYjv);
            _mm_storeu_pd(&Jac[(i+3)*(alignN) + j+0],ddYi3dtdYjv);
        }
    }
    {
        unsigned int i = this->nSpecies-2;
        const double Wi0ByrhoM_ = WiByrhoM[i+0];
        const double Wi1ByrhoM_ = WiByrhoM[i+1];

        __m256d Wi0ByrhoMv = _mm256_set1_pd(Wi0ByrhoM_);
        __m256d Wi1ByrhoMv = _mm256_set1_pd(Wi1ByrhoM_);

        double dYi0dt = dYTpdt[i+0]*WiByrhoM[i+0];
        double dYi1dt = dYTpdt[i+1]*WiByrhoM[i+1];

        __m256d dYi0dtv = _mm256_set1_pd(dYi0dt);
        __m256d dYi1dtv = _mm256_set1_pd(dYi1dt);

        for (unsigned int j=0; j<this->nSpecies-2; j=j+4)
        {
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,_mm256_loadu_pd(&this->invW[j+0]));
            __m256d rhoMByRhov = _mm256_loadu_pd(&rhoMByRhoi[j+0]);

            __m256d ddNi0dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+0)*(alignN) + j+0]),dCjdYj);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+1)*(alignN) + j+0]),dCjdYj);
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,_mm256_mul_pd(rhoMByRhov,dYi0dtv));
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,_mm256_mul_pd(rhoMByRhov,dYi1dtv));            
            _mm256_storeu_pd(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);
            _mm256_storeu_pd(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);
        }
        {
            unsigned int j = this->nSpecies-2;
            __m128d dCjdYj = _mm_mul_pd(_mm256_castpd256_pd128(rhoMv),_mm_loadu_pd(&this->invW[j+0]));
            __m128d rhoMByRhov = _mm_loadu_pd(&rhoMByRhoi[j+0]);
            __m128d ddNi0dtByVdYj = _mm_mul_pd(_mm_loadu_pd(&ddNdtByVdcT[(i+0)*(alignN) + j+0]),dCjdYj);
            __m128d ddNi1dtByVdYj = _mm_mul_pd(_mm_loadu_pd(&ddNdtByVdcT[(i+1)*(alignN) + j+0]),dCjdYj);
            __m128d ddYi0dtdYj = _mm_fmadd_pd(_mm256_castpd256_pd128(Wi0ByrhoMv),ddNi0dtByVdYj,_mm_mul_pd(rhoMByRhov,_mm256_castpd256_pd128(dYi0dtv)));
            __m128d ddYi1dtdYj = _mm_fmadd_pd(_mm256_castpd256_pd128(Wi1ByrhoMv),ddNi1dtByVdYj,_mm_mul_pd(rhoMByRhov,_mm256_castpd256_pd128(dYi1dtv)));            
            _mm_storeu_pd(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);
            _mm_storeu_pd(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);
        }
    }
}
void 
OptReaction::FastddYdtdY_Vec3
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ rhoMByRhoi,
    const double* __restrict__ WiByrhoM,  
    const double* __restrict__ dYTpdt, 
    double* __restrict__ Jac
) const noexcept
{ 
    __m256d rhoMv = _mm256_set1_pd(this->rhoM);    
    for(unsigned int i=0; i<this->nSpecies-3; i=i+4)
    {
        const double Wi0ByrhoM_ = WiByrhoM[i+0];
        const double Wi1ByrhoM_ = WiByrhoM[i+1];
        const double Wi2ByrhoM_ = WiByrhoM[i+2];
        const double Wi3ByrhoM_ = WiByrhoM[i+3];

        __m256d Wi0ByrhoMv = _mm256_set1_pd(Wi0ByrhoM_);
        __m256d Wi1ByrhoMv = _mm256_set1_pd(Wi1ByrhoM_);
        __m256d Wi2ByrhoMv = _mm256_set1_pd(Wi2ByrhoM_);
        __m256d Wi3ByrhoMv = _mm256_set1_pd(Wi3ByrhoM_);

        double dYi0dt = dYTpdt[i+0]*WiByrhoM[i+0];
        double dYi1dt = dYTpdt[i+1]*WiByrhoM[i+1];
        double dYi2dt = dYTpdt[i+2]*WiByrhoM[i+2];
        double dYi3dt = dYTpdt[i+3]*WiByrhoM[i+3];

        __m256d dYi0dtv = _mm256_set1_pd(dYi0dt);
        __m256d dYi1dtv = _mm256_set1_pd(dYi1dt);
        __m256d dYi2dtv = _mm256_set1_pd(dYi2dt);
        __m256d dYi3dtv = _mm256_set1_pd(dYi3dt);

        for (unsigned int j=0; j<this->nSpecies-3; j=j+4)
        {
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,_mm256_loadu_pd(&this->invW[j+0]));
            __m256d ddNi0dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+0)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+1)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi2dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+2)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi3dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+3)*(alignN)+j+0]),dCjdYj);
            __m256d rhoMvj_ = _mm256_loadu_pd(&rhoMByRhoi[j+0]);
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi0dtv));
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi1dtv));
            __m256d ddYi2dtdYj = _mm256_fmadd_pd(Wi2ByrhoMv,ddNi2dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi2dtv));
            __m256d ddYi3dtdYj = _mm256_fmadd_pd(Wi3ByrhoMv,ddNi3dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi3dtv));
            _mm256_storeu_pd(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);
            _mm256_storeu_pd(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);
            _mm256_storeu_pd(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYj);
            _mm256_storeu_pd(&Jac[(i+3)*(alignN) + j+0],ddYi3dtdYj);
        }
        {
            unsigned int j = this->nSpecies-3;
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,_mm256_setr_pd(this->invW[j+0],this->invW[j+1],this->invW[j+2],0));
            __m256d ddNi0dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+0)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+1)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi2dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+2)*(alignN)+j+0]),dCjdYj);
            __m256d ddNi3dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+3)*(alignN)+j+0]),dCjdYj);  
            __m256d rhoMvj_ = _mm256_setr_pd(rhoMByRhoi[j+0],rhoMByRhoi[j+1],rhoMByRhoi[j+2],0);
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi0dtv));
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi1dtv));
            __m256d ddYi2dtdYj = _mm256_fmadd_pd(Wi2ByrhoMv,ddNi2dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi2dtv));
            __m256d ddYi3dtdYj = _mm256_fmadd_pd(Wi3ByrhoMv,ddNi3dtByVdYj,_mm256_mul_pd(rhoMvj_,dYi3dtv)); 
            _mm256_storeu_pd(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);
            _mm256_storeu_pd(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);
            _mm256_storeu_pd(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYj);
            _mm256_storeu_pd(&Jac[(i+3)*(alignN) + j+0],ddYi3dtdYj);                                 
        }
    }
    {
        unsigned int i = this->nSpecies - 3;
        const double Wi0ByrhoM_ = WiByrhoM[i+0];
        const double Wi1ByrhoM_ = WiByrhoM[i+1];
        const double Wi2ByrhoM_ = WiByrhoM[i+2];
        __m256d Wi0ByrhoMv = _mm256_set1_pd(Wi0ByrhoM_);
        __m256d Wi1ByrhoMv = _mm256_set1_pd(Wi1ByrhoM_);
        __m256d Wi2ByrhoMv = _mm256_set1_pd(Wi2ByrhoM_);
        double dYi0dt = dYTpdt[i+0]*WiByrhoM[i+0];
        double dYi1dt = dYTpdt[i+1]*WiByrhoM[i+1];
        double dYi2dt = dYTpdt[i+2]*WiByrhoM[i+2];
        __m256d dYi0dtv = _mm256_set1_pd(dYi0dt);
        __m256d dYi1dtv = _mm256_set1_pd(dYi1dt);
        __m256d dYi2dtv = _mm256_set1_pd(dYi2dt);

        for (unsigned int j=0; j<this->nSpecies-3; j=j+4)
        {
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,_mm256_loadu_pd(&this->invW[j+0]));
            __m256d rhoMByRhov = _mm256_loadu_pd(&rhoMByRhoi[j+0]);
            __m256d ddNi0dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+0)*(alignN) + j+0]),dCjdYj);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+1)*(alignN) + j+0]),dCjdYj);
            __m256d ddNi2dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+2)*(alignN) + j+0]),dCjdYj);

            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,_mm256_mul_pd(rhoMByRhov,dYi0dtv));
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,_mm256_mul_pd(rhoMByRhov,dYi1dtv));
            __m256d ddYi2dtdYj = _mm256_fmadd_pd(Wi2ByrhoMv,ddNi2dtByVdYj,_mm256_mul_pd(rhoMByRhov,dYi2dtv));
            _mm256_storeu_pd(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);
            _mm256_storeu_pd(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);
            _mm256_storeu_pd(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYj);
        }
        {
            unsigned int j = this->nSpecies-3;
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,_mm256_setr_pd(this->invW[j+0],this->invW[j+1],this->invW[j+2],0)); 
            __m256d rhoMByRhov = _mm256_setr_pd(rhoMByRhoi[j+0],rhoMByRhoi[j+1],rhoMByRhoi[j+2],0);
            __m256d ddNi0dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+0)*(alignN) + j+0]),dCjdYj);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+1)*(alignN) + j+0]),dCjdYj);
            __m256d ddNi2dtByVdYj = _mm256_mul_pd(_mm256_loadu_pd(&ddNdtByVdcT[(i+2)*(alignN) + j+0]),dCjdYj); 
            
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,_mm256_mul_pd(rhoMByRhov,dYi0dtv));
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,_mm256_mul_pd(rhoMByRhov,dYi1dtv));
            __m256d ddYi2dtdYj = _mm256_fmadd_pd(Wi2ByrhoMv,ddNi2dtByVdYj,_mm256_mul_pd(rhoMByRhov,dYi2dtv));
            _mm256_storeu_pd(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);
            _mm256_storeu_pd(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);
            _mm256_storeu_pd(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYj);
        }
    }
}

void 
OptReaction::ddTdtdYT
(
    const double* __restrict__ Cp,
    const double* __restrict__ dCpdT,
    const double* __restrict__ Ha,      
    double* __restrict__ dPhidt,  
    double* __restrict__ Jac
) const noexcept
{

    //unsigned int remain = this->nSpecies%4;
    const double& CpM = Cp[this->nSpecies];
    const double& dCpMdT = dCpdT[this->nSpecies];

    // dT/dt
    double& dTdt = dPhidt[this->nSpecies];
    dTdt = 0;

    for (unsigned int i=0; i<this->nSpecies; i++)
    {
        dTdt -= dPhidt[i]*Ha[i];
    }

    dTdt /= CpM;


    // d(dTdt)/dY
    double& ddTdtdT = Jac[this->nSpecies *(n_)+ this->nSpecies];
    ddTdtdT = 0;
    for (unsigned int i=0; i<this->nSpecies; i++)
    {
        double& ddTdtdYi = Jac[this->nSpecies *(n_)+ i];
        ddTdtdYi = 0;

        for (unsigned int j=0; j<this->nSpecies; j++)
        {
            const double ddYjdtdYi = Jac[j *(n_)+ i];
            ddTdtdYi -= ddYjdtdYi*Ha[j];
        }
        ddTdtdYi -= Cp[i]*dTdt;
        ddTdtdYi /= CpM;

        const double dYidt = dPhidt[i];
        const double ddYidtdT = Jac[i *(n_)+ this->nSpecies];

        ddTdtdT -= dYidt*Cp[i] + ddYidtdT*Ha[i];
    }

    ddTdtdT -= dTdt*dCpMdT;
    ddTdtdT /= CpM;   
}

void 
OptReaction::ddYdtdTP
(
    const double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ WiByrhoM,
    const double* __restrict__ C,  
    double* __restrict__ dPhidt, 
    double* __restrict__ Jac
) const noexcept
{
    //double* const* __restrict__ ArrRes = this->ArrPtr.data();
    double alphavM = this->invT;

    for (unsigned int i=0; i<this->nSpecies; i=i+1)
    {
        unsigned int i0 = i;

        //const scalar WiByrhoM = reaction.W[i]/rhoM;
        const double Wi0ByrhoM_ = WiByrhoM[i0];
        double& dYi0dt = dPhidt[i0];
        dYi0dt *= Wi0ByrhoM_;
        double ddNi0dtByVdT = ddNdtByVdcTp[(i0)*(n_) + nSpecies];
        for (unsigned int j=0; j<this->nSpecies; j++)
        {
            const double ddNi0dtByVdcj = ddNdtByVdcTp[i0*(n_) + j];
            ddNi0dtByVdT -= ddNi0dtByVdcj*C[j]*alphavM;                                 
        }

        double& ddYi0dtdT = Jac[i0*(n_) + nSpecies];
        ddYi0dtdT = Wi0ByrhoM_*ddNi0dtByVdT + alphavM*dYi0dt;
     
    }
}

void 
OptReaction::ddYdtdY
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ rhoMByRhoi,
    const double* __restrict__ WiByrhoM,  
    const double* __restrict__ dPhidt,  
    const double* __restrict__ Phi, 
    double* __restrict__ Jac
) const noexcept
{
    std::vector<double> Array0(this->nSpecies,0);

    // Convert d(dCdt)dC into d(dYdt)dY
    for (unsigned int j=0; j<this->nSpecies; j++)
    {
        // Compute dcdY
        for(unsigned int i=0; i<this->nSpecies; i++)
        {
            const double rhoMByWi = rhoM*this->invW[i];
            const double rhoMvj_ = rhoMByRhoi[j];
            Array0[i] = rhoMByWi*((i == j) - rhoMvj_*Phi[i]);
        }

        // Convert d(dCdt)dC into d(dYdt)dY
        for(unsigned int i=0; i<this->nSpecies; i++)
        {

            //const scalar WiByrhoM = reaction.W[i]/rhoM;
            const double WiByrhoM_ = WiByrhoM[i];

            //const scalar rhoMvj = MW*reaction.invW[j];
            const double rhoMvj_ = rhoMByRhoi[j];
            double dYidt = dPhidt[i];
            dYidt *= WiByrhoM_;

            double ddNidtByVdYj = 0;

            for (unsigned int k=0; k<this->nSpecies; k++)
            {
                const double ddNidtByVdck = ddNdtByVdcT[i*(alignN) + k];
                ddNidtByVdYj += ddNidtByVdck*Array0[k];            
            }

            double& ddYidtdYj = Jac[i*(alignN) + j];
            ddYidtdYj = WiByrhoM_*ddNidtByVdYj + rhoMvj_*dYidt;
        }
    }


}

