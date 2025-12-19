
#include "OptReaction.H"

void 
OptReaction::Tc
(
    int celli,
    double p,
    double Temperature,
    double* C,
    double& sumW,
    double& sumWRateByCTot
) const noexcept
{

    this->update_Pow_pByRT_SumVki(Temperature);
    this->update_Pow_pByRT_SumVki2(Temperature);

    {
        for(unsigned i = 0; i <this->n_Troe;i++)
        {
            const unsigned j0 = i + this->nSpecies;
            const unsigned j1 = i + this->nSpecies + this->n_Troe;
            const unsigned j2 = i + this->nSpecies + this->n_Troe*2;         
            this->tmp_Exp[j0] = -Temperature*this->invTsss_[i];
            this->tmp_Exp[j1] = -this->Tss_[i]*invT;    
            this->tmp_Exp[j2] = -Temperature*this->invTs_[i];
        }
    }

    {
        for(unsigned i = 0; i <this->n_SRI;i++)
        {
            const unsigned j0 = i + this->nSpecies + this->n_Troe*3;
            const unsigned j1 = i + this->nSpecies + this->n_Troe*3 + this->n_SRI;
            this->tmp_Exp[j0] = -this->b_[i]*invT;
            this->tmp_Exp[j1] = -Temperature*this->invc_[i];
        }   

        unsigned int remain = this->tmp_ExpSize%4;
        for(unsigned int i = 0; i < this->tmp_ExpSize-remain;i=i+4)
        {
            __m256d tmp = _mm256_loadu_pd(&this->tmp_Exp[i]);
            tmp = vec256_expd(tmp);
            _mm256_storeu_pd(&this->tmp_Exp[i],tmp);
        }
        if(remain==1)
        {
            const unsigned i = this->tmp_ExpSize-1;
            this->tmp_Exp[i] = std::exp(this->tmp_Exp[i]);
        }
        else if(remain==2)
        {
            const unsigned i0 = this->tmp_ExpSize-2;
            const unsigned i1 = this->tmp_ExpSize-1;
            __m256d tmp = _mm256_setr_pd(tmp_Exp[i0],tmp_Exp[i1],0,0);
            tmp = vec256_expd(tmp);
            this->tmp_Exp[i0] = get_elem0(tmp);
            this->tmp_Exp[i1] = get_elem1(tmp);
        }
        else if(remain==3)
        {
            const unsigned i0 = this->tmp_ExpSize-3;
            const unsigned i1 = this->tmp_ExpSize-2;
            const unsigned i2 = this->tmp_ExpSize-1;

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

                this->A[i+this->Ikf[6]] = A0;
                this->A[i+this->Ikf[11]] = A0;
                this->beta[i+this->Ikf[6]] = beta0;
                this->beta[i+this->Ikf[11]] = beta0;
                this->Ta[i+this->Ikf[6]] = Ta0;
                this->Ta[i+this->Ikf[11]] = Ta0;
                this->Pindex[i] = 0;
            }
            else if(p>=this->Prange[i][length-1])
            {
                double A1 = this->APlog[i][length-1];
                double beta1 = this->betaPlog[i][length-1];
                double Ta1 = this->TaPlog[i][length-1];
                
                this->A[i+this->Ikf[6]] = A1;
                this->A[i+this->Ikf[11]] = A1;
                this->beta[i+this->Ikf[6]] = beta1;
                this->beta[i+this->Ikf[11]] = beta1;
                this->Ta[i+this->Ikf[6]] = Ta1;
                this->Ta[i+this->Ikf[11]] = Ta1;
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
                this->A[i+this->Ikf[6]] = A0;
                this->A[i+this->Ikf[11]] = A1;
                this->beta[i+this->Ikf[6]] = beta0;
                this->beta[i+this->Ikf[11]] = beta1;
                this->Ta[i+this->Ikf[6]] = Ta0;
                this->Ta[i+this->Ikf[11]] = Ta1;
            }
        }
    }


    {
    
        __m256d LogT = _mm256_set1_pd(logT);
        __m256d InvT = _mm256_set1_pd(-invT);
        unsigned int remain = (this->Ikf[12]-this->n_Temperature_Independent_Reaction)%4;
        unsigned int times = (this->Ikf[12]-this->n_Temperature_Independent_Reaction)/4;
        for(unsigned int z = 0; z <times;z=z+1)
        {
            unsigned int i = z*4 + this->n_Temperature_Independent_Reaction;
            __m256d A_ = _mm256_loadu_pd(&this->A[i]);
            __m256d beta_ = _mm256_loadu_pd(&this->beta[i]);
            __m256d Ta_ = _mm256_loadu_pd(&this->Ta[i]);
            __m256d Kf = _mm256_mul_pd(Ta_,InvT);
            Kf = _mm256_fmadd_pd(beta_,LogT,Kf);
            Kf = vec256_expd(Kf);
            Kf = _mm256_mul_pd(A_,Kf);
            _mm256_storeu_pd(&this->Kf_[i],Kf);
        }
        if(remain==1)
        {
            unsigned int i = this->Ikf[12]-1;
            this->Kf_[i] = this->A[i]*std::exp(this->beta[i+0]*logT-this->Ta[i+0]*invT);   
        }
        else if(remain==2)
        {
            unsigned int i0 = this->Ikf[12]-2;
            unsigned int i1 = this->Ikf[12]-1;
            __m256d A_ = _mm256_setr_pd(this->A[i0],this->A[i1],0,0);
            __m256d beta_ = _mm256_setr_pd(this->beta[i0],this->beta[i1],0,0);
            __m256d Ta_ = _mm256_setr_pd(this->Ta[i0],this->Ta[i1],0,0);
            __m256d Kf = _mm256_mul_pd(Ta_,InvT);
            Kf = _mm256_fmadd_pd(beta_,LogT,Kf);
            Kf = vec256_expd(Kf);
            Kf = _mm256_mul_pd(A_,Kf);
            this->Kf_[i0] = get_elem0(Kf);
            this->Kf_[i1] = get_elem1(Kf);
        }
        else if(remain==3)
        {
            unsigned int i0 = this->Ikf[12]-3;
            unsigned int i1 = this->Ikf[12]-2;
            unsigned int i2 = this->Ikf[12]-1;
            __m256d A_ = _mm256_setr_pd(this->A[i0],this->A[i1],this->A[i2],0);
            __m256d beta_ = _mm256_setr_pd(this->beta[i0],this->beta[i1],this->beta[i2],0);
            __m256d Ta_ = _mm256_setr_pd(this->Ta[i0],this->Ta[i1],this->Ta[i2],0);
            __m256d Kf = _mm256_mul_pd(Ta_,InvT);
            Kf = _mm256_fmadd_pd(beta_,LogT,Kf);
            Kf = vec256_expd(Kf);
            Kf = _mm256_mul_pd(A_,Kf);
            this->Kf_[i0] = get_elem0(Kf);
            this->Kf_[i1] = get_elem1(Kf); 
            this->Kf_[i2] = get_elem2(Kf); 
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
                __m256d C_ = _mm256_loadu_pd(&C[j+0]);
                
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
                __m256d C_ = _mm256_loadu_pd(&C[j+0]);
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
                __m256d C_ = _mm256_loadu_pd(&C[j+0]);
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
                __m256d C_ = _mm256_loadu_pd(&C[j+0]);
                arrM_0 = _mm256_fmadd_pd(Factor0,C_,arrM_0);
            }

            M0 = M0 + hsum4(arrM_0);

            this->tmp_M[i+0] = M0;
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
                unsigned index = this->Pindex[i];
                double weight = (this->logP - this->logPi[i][index])*this->rDeltaP_[i][index];
                double Kf0 = this->Kf_[i+this->Ikf[6]];
                double Kf1 = this->Kf_[i+this->Ikf[11]];
                this->Kf_[i+this->Ikf[6]] = Kf0*std::pow(Kf1/Kf0,weight);
            }
        }
    }


    {
        for(unsigned int i = 0; i < this->n_ThirdBodyReaction; i++)
        {
            const unsigned int j = i + this->Ikf[3];
            this->Kf_[j] = this->Kf_[j]*this->tmp_M[i+this->Itbr[1]];
        }
    }

    {
        for(unsigned int i = 0; i < this->n_NonEquilibriumThirdBodyReaction; i++)
        {
            double Mfwd = this->tmp_M[i];
            double Mrev = this->tmp_M[this->Itbr[4]+i];
            this->Kf_[Ikf[2]+i] = this->Kf_[Ikf[2]+i]*Mfwd;
            this->Kf_[Ikf[10]+i] = this->Kf_[Ikf[10]+i]*Mrev;
        } 
    }

    if(this->n_Lindemann>0)
    {
        this->evalLindemannRateConstant();
        /*unsigned int remain = this->n_LindemannFO%4;
        for (unsigned int i = 0;i<this->n_LindemannFO-remain;i=i+4)
        {
            const unsigned int j0 = this->LindemannFO[i+0]+0;
            const unsigned int j1 = this->LindemannFO[i+0]+1;
            const unsigned int j2 = this->LindemannFO[i+0]+2;
            const unsigned int j3 = this->LindemannFO[i+0]+3;

            const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
            const unsigned int k0 = j0 - this->Ikf[4];
            const unsigned int k1 = j1 - this->Ikf[4];
            const unsigned int k2 = j2 - this->Ikf[4];
            const unsigned int k3 = j3 - this->Ikf[4];
            __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
            __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
            __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);         
            __m256d Pr = _mm256_div_pd(_mm256_mul_pd(K0,M),Kinf);
            __m256d N = _mm256_div_pd(K0,_mm256_add_pd(Pr,_mm256_set1_pd(1.0)));
            __m256d k = _mm256_setr_pd(k0,k1,k2,k3);
            __m256d cmp = _mm256_cmp_pd(k,_mm256_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
            __m256d Kf = _mm256_blendv_pd(N,_mm256_mul_pd(M,N),cmp);
            _mm256_storeu_pd(&this->Kf_[j0],Kf);
        }
        if(remain==1)
        {
            const unsigned i = this->n_LindemannFO-1;
            const unsigned int j = this->LindemannFO[i];
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
            const unsigned i = this->n_LindemannFO-2;
            const unsigned int j0 = this->LindemannFO[i+0]+0;
            const unsigned int j1 = this->LindemannFO[i+0]+1;
            const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
            const unsigned int k0 = j0 - this->Ikf[4];
            const unsigned int k1 = j1 - this->Ikf[4];
            __m128d Kinf = _mm_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
            __m128d M = _mm_loadu_pd(&this->tmp_M[m0]);
            __m128d K0 = _mm_loadu_pd(&this->Kf_[j0]);
            __m128d Pr = _mm_div_pd(_mm_mul_pd(K0,M),Kinf);
            __m128d N = _mm_div_pd(K0,_mm_add_pd(Pr,_mm_set1_pd(1.0)));
            __m128d k = _mm_setr_pd(k0,k1);
            __m128d cmp = _mm_cmp_pd(k,_mm_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
            __m128d Kf = _mm_blendv_pd(N,_mm_mul_pd(M,N),cmp);
            _mm_storeu_pd(&this->Kf_[j0],Kf);
        }
        else if(remain==3)
        {
            const unsigned i = this->n_LindemannFO-3;
            const unsigned int j0 = this->LindemannFO[i+0];
            const unsigned int j1 = this->LindemannFO[i+1];
            const unsigned int j2 = this->LindemannFO[i+2];
            const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
            const unsigned int m1 = j1 - this->Ikf[4] + this->Itbr[2];
            const unsigned int m2 = j2 - this->Ikf[4] + this->Itbr[2];
            const unsigned int k0 = j0 - this->Ikf[4];
            const unsigned int k1 = j1 - this->Ikf[4];
            const unsigned int k2 = j2 - this->Ikf[4];
            const double Kinf0 = this->Kf_[j0+this->offset_kinf];
            const double Kinf1 = this->Kf_[j1+this->offset_kinf];
            const double Kinf2 = this->Kf_[j2+this->offset_kinf];
            double M0 = this->tmp_M[m0];
            double M1 = this->tmp_M[m1];
            double M2 = this->tmp_M[m2];
            const double K00 = this->Kf_[j0];
            const double K01 = this->Kf_[j1];
            const double K02 = this->Kf_[j2];
            __m256d Kinf = _mm256_setr_pd(Kinf0,Kinf1,Kinf2,1);
            __m256d M = _mm256_setr_pd(M0,M1,M2,0);
            __m256d K0 = _mm256_setr_pd(K00,K01,K02,0);         
            __m256d Pr = _mm256_div_pd(_mm256_mul_pd(K0,M),Kinf);
            __m256d N = _mm256_div_pd(K0,_mm256_add_pd(Pr,_mm256_set1_pd(1.0)));
            __m256d k = _mm256_setr_pd(k0,k1,k2,0);
            __m256d cmp = _mm256_cmp_pd(k,_mm256_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
            __m256d Kf = _mm256_blendv_pd(N,_mm256_mul_pd(M,N),cmp);
            this->Kf_[j0] = get_elem0(Kf);
            this->Kf_[j1] = get_elem1(Kf);
            this->Kf_[j2] = get_elem2(Kf);
        }*/

    }

    if(this->n_Troe>0)
    {
        this->evalTroeRateConstant();
        /*unsigned int remain = this->n_Troe%4;
        for(unsigned int i = 0; i < this->n_Troe-remain;i=i+4)
        {
            const unsigned int j0 = this->TroeFO[i+0];
            const unsigned int j1 = j0 + 1;
            const unsigned int j2 = j0 + 2;
            const unsigned int j3 = j0 + 3;
            const unsigned int k0 = j0 - this->Ikf[4];
            const unsigned int k1 = j1 - this->Ikf[4];
            const unsigned int k2 = j2 - this->Ikf[4];
            const unsigned int k3 = j3 - this->Ikf[4];
            const unsigned int m0 = j0 - this->Ikf[4] + Itbr[2];
            __m256d Kinf = _mm256_loadu_pd(&this->Kf_[j0+this->offset_kinf]);
            __m256d M = _mm256_loadu_pd(&this->tmp_M[m0]);
            __m256d K0 = _mm256_loadu_pd(&this->Kf_[j0]);    
            __m256d Pr_ = _mm256_div_pd(_mm256_mul_pd(K0,M),Kinf);
            __m256d small = _mm256_set1_pd(2.2e-16);
            Pr_ = _mm256_max_pd(small,Pr_);
            const double invLog10 = 0.43429448190325182765112891891661;
            __m256d logPr_ = _mm256_mul_pd(vec256_logd(Pr_),_mm256_set1_pd(invLog10));
            __m256d alpha = _mm256_loadu_pd(&this->alpha_[i]);
            __m256d one = _mm256_set1_pd(1.0);
            __m256d expTTsss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies]);
            __m256d expTTss = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe]);
            __m256d expTTs = _mm256_loadu_pd(&this->tmp_Exp[i+this->nSpecies+this->n_Troe*2]);
            __m256d Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha), expTTsss);
            Fcent = _mm256_fmadd_pd(alpha, expTTs,Fcent);
            Fcent = _mm256_add_pd(expTTss,Fcent);
            __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),_mm256_set1_pd(invLog10));
            __m256d cc = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(0.67),_mm256_set1_pd(0.4));
            __m256d n = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(-1.27),_mm256_set1_pd(0.75));
            __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(cc,logPr_),_mm256_set1_pd(0.14),n);
            __m256d x2 = _mm256_div_pd(_mm256_sub_pd(logPr_,cc),x1);
            __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
            __m256d x4 = _mm256_div_pd(logFcent,x3);
            __m256d F_ = vec256_powd(_mm256_set1_pd(10),x4);
            __m256d N = _mm256_div_pd(_mm256_mul_pd(K0,F_),_mm256_add_pd(_mm256_set1_pd(1.0),Pr_));
            __m256d k = _mm256_setr_pd(k0,k1,k2,k3);
            __m256d cmp = _mm256_cmp_pd(k,_mm256_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
            __m256d Kf = _mm256_blendv_pd(N,_mm256_mul_pd(M,N),cmp);
            _mm256_storeu_pd(&this->Kf_[j0],Kf);         
        }        
        if(remain==1)       {this->Troe_F_1();}
        else if(remain==2)  {this->Troe_F_2();}
        else if(remain==3)  {this->Troe_F_3();}*/
    }
    if(this->n_SRI>0)
    {
        this->evalSRIRateConstant(T);
        /*for (unsigned int i = 0;i<this->n_SRI;i++)
        {
            const unsigned int j = this->SRIFO[i];

            const unsigned int m = j - this->Ikf[4] + this->Itbr[2];
            const unsigned int k = j - this->Ikf[4];
            const double Kinf = this->Kf_[j+this->offset_kinf];
            double M = this->tmp_M[m]; 
            const double K0 = this->Kf_[j];
            const double Pr = K0*M/Kinf;   
            const double F  = this->SRI_F(Temperature,Pr,i);
            const double N  = 1/(1+Pr)*F*K0;
            this->Kf_[j] = k<this->n_Fall_Off_Reaction ? M*N : N;   
        }*/
    }



    for (unsigned int i = 0; i < this->Ikf[7]; i++) 
    {
        if(this->isGlobal[i]==1)
        {
            double CR = 1.0;
            double CF = 1.0;
            double Kc_ = 1.0;
            double Kr = 0;
            if(this->isIrreversible[i]==0) // Reversible reaction
            {
                double Kp = 0;
                double sumVki = 0;
                for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
                {
                    const unsigned int si = this->lhsSpeciesIndex[i][j];
                    const double sl = this->lhsStoichCoeff[i][j];
                    const double el = this->lhsReactionOrder[i][j];
                    Kp += sl*this->negGstdByRT[si];
                    sumVki = sumVki - sl;
                    CF = CF * (C[si] >= small || el >= 1 ? std::pow(std::max(C[si], 0.0), el) : 0.0);
                }
                
                for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
                {
                    const unsigned int si = this->rhsSpeciesIndex[i][j];
                    const double sr = this->rhsStoichCoeff[i][j];
                    const double er = this->rhsReactionOrder[i][j];
                    Kp -= sr*this->negGstdByRT[si];
                    sumVki = sumVki + sr;
                    CR = CR * (C[si] >= small || er >= 1 ? std::pow(std::max(C[si], 0.0), er) : 0.0);            
                }
                Kp = std::exp(Kp);
                Kc_ = Kp*std::pow(this->Pstd/(this->Ru*this->T),sumVki);
                Kc_ = std::max(Kc_,1.49011611938476E-08);
                Kr = this->Kf_[i]/Kc_;
            }
            else if(this->isIrreversible[i]==2)
            {
                unsigned int J = i - this->Ikf[1] + this->Ikf[9];
                Kr = this->Kf_[J];

                for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
                {
                    const unsigned int si = this->lhsSpeciesIndex[i][j];
                    const double el = this->lhsReactionOrder[i][j];
                    CF = CF * (C[si] >= small || el >= 1 ? std::pow(std::max(C[si], 0.0), el) : 0.0);
                }
                
                for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
                {
                    const unsigned int si = this->rhsSpeciesIndex[i][j];
                    const double er = this->rhsReactionOrder[i][j];
                    CR = CR * (C[si] >= small || er >= 1 ? std::pow(std::max(C[si], 0.0), er) : 0.0);            
                }
            }
            else
            {
                for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
                {
                    const unsigned int si = this->lhsSpeciesIndex[i][j];
                    const double el = this->lhsReactionOrder[i][j];
                    CF = CF * (C[si] >= small || el >= 1 ? std::pow(std::max(C[si], 0.0), el) : 0.0);
                }
                
                for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
                {
                    const unsigned int si = this->rhsSpeciesIndex[i][j];
                    const double er = this->rhsReactionOrder[i][j];
                    CR = CR * (C[si] >= small || er >= 1 ? std::pow(std::max(C[si], 0.0), er) : 0.0);            
                }
            }


            const double omegaf = this->Kf_[i]*CF;
            const double omegar = Kr*CR;

            double wf = 0;
            for(unsigned int j = 0; j < lhsSpeciesIndex[i].size();j++)
            {
                const double sl = lhsStoichCoeff[i][j];
                wf += sl*omegaf;
            }
            sumW += wf;
            sumWRateByCTot += (wf*wf);

            double wr = 0;
            for(unsigned int j = 0; j < rhsSpeciesIndex[i].size();j++)
            {
                const double sr = rhsStoichCoeff[i][j];
                wr += sr*omegar   ;
            }
            sumW += wr;
            sumWRateByCTot += (wr*wr);            
            continue;
        }
   
        auto J = lhsOffset[i+1]-lhsOffset[i];
        auto K = rhsOffset[i+1]-rhsOffset[i];
        const double* __restrict__ const ExpNegGbyRT = &tmp_Exp[0];
        if(J==2)
        {
            if(K==2)
            {
                const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
                const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
                const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
                const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];
                double Kr = 0;
                if(this->isIrreversible[i]==0)
                {
                    const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
                    double Kc = Kp;
                    Kc = std::max(Kc,1.4901171103413047e-8);  
                    Kr = this->Kf_[i]/Kc;         
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int l = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[l];
                }

                const double CF = C[sl0]*C[sl1];
                const double CR = C[sr0]*C[sr1];

                const double omegaf = (this->Kf_[i]*CF);
                const double omegar = (Kr*CR);

                double wf = omegaf + omegaf;
                double wr = omegar + omegar;
                sumW += wf;
                sumWRateByCTot += (wf*wf);
                sumW += wr;
                sumWRateByCTot += (wr*wr);
            }
            else if(K==1)
            {
                const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
                const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
                const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
                double Kr = 0;
                if(this->isIrreversible[i]==0)
                {
                    const double Kp = ExpNegGbyRT[sr0]/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
                    double Kc = Kp*this->Pow_pByRT_SumVki[1];
                    Kc = std::max(Kc,1.4901171103413047e-8);
                    Kr = this->Kf_[i]/Kc;        
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }
                const double CF = C[sl0]*C[sl1];
                const double CR = C[sr0];

                const double omegaf = (this->Kf_[i]*CF);
                const double omegar = (Kr*CR);

                double wf = omegaf + omegaf;
                double wr = omegar;
                sumW += wf;
                sumWRateByCTot += (wf*wf);
                sumW += wr;
                sumWRateByCTot += (wr*wr);
 
            }
            else if(K==3)
            {
                const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
                const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
                const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
                const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];
                const unsigned int sr2 = rhsSpeciesIndex1D[rhsOffset[i]+2];

                double Kr = 0;

                if(this->isIrreversible[i]==0)
                {
                    const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
                    double Kc = Kp*this->Pow_pByRT_SumVki[3];
                    Kc = std::max(Kc,1.4901171103413047e-8);
                    Kr = this->Kf_[i]/Kc;    
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }

                const double CF = C[sl0]*C[sl1];
                const double CR = C[sr0]*C[sr1]*C[sr2];

                const double omegaf = (this->Kf_[i]*CF);
                const double omegar = (Kr*CR);

                double wf = omegaf + omegaf;
                double wr = omegar + omegar + omegar;
                sumW += wf;
                sumWRateByCTot += (wf*wf);
                sumW += wr;
                sumWRateByCTot += (wr*wr);
            }
        }
        else if(J==1)
        {
            if(K==2)
            {
                const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
                const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
                const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];
                double Kr = 0;

                if(this->isIrreversible[i]==0)
                {
                    const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1])/(ExpNegGbyRT[sl0]);
                    double Kc = Kp*this->Pow_pByRT_SumVki[3];
                    Kc = std::max(Kc,1.4901171103413047e-8);
                    Kr = this->Kf_[i]/Kc;       
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }

                const double CF = C[sl0];
                const double CR = C[sr0]*C[sr1];

                const double omegaf = (this->Kf_[i]*CF);
                const double omegar = (Kr*CR);

                double wf = omegaf;
                double wr = omegar + omegar;
                sumW += wf;
                sumWRateByCTot += (wf*wf);
                sumW += wr;
                sumWRateByCTot += (wr*wr);
            }
            else if(K==1)
            {
                const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
                const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
                double Kr = 0;
                if(this->isIrreversible[i]==0)
                {
                    const double Kp = ExpNegGbyRT[sr0]/ExpNegGbyRT[sl0];
                    double Kc = Kp;
                    Kc = std::max(Kc,1.4901171103413047e-8);
                    Kr = this->Kf_[i]/Kc;       
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }
                const double CF = C[sl0];
                const double CR = C[sr0];


                const double omegaf = (this->Kf_[i]*CF);
                const double omegar = (Kr*CR);

                double wf = omegaf;
                double wr = omegar;
                sumW += wf;
                sumWRateByCTot += (wf*wf);
                sumW += wr;
                sumWRateByCTot += (wr*wr);
            }
            else if(K==3)
            {
                const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
                const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
                const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];
                const unsigned int sr2 = rhsSpeciesIndex1D[rhsOffset[i]+2];

                double Kr = 0;
                if(this->isIrreversible[i]==0)
                {
                    const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]);
                    double Kc = Kp*this->Pow_pByRT_SumVki[4];
                    Kc = std::max(Kc,1.4901171103413047e-8);
                    Kr = this->Kf_[i]/Kc;         
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }

                const double CF = C[sl0];
                const double CR = C[sr0]*C[sr1]*C[sr2];

                const double omegaf = (this->Kf_[i]*CF);
                const double omegar = (Kr*CR);

                double wf = omegaf;
                double wr = omegar + omegar + omegar;
                sumW += wf;
                sumWRateByCTot += (wf*wf);
                sumW += wr;
                sumWRateByCTot += (wr*wr);
            }
        }
        else if(J==3)
        {
            if(K==2)
            {
                const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
                const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
                const unsigned int sl2 = lhsSpeciesIndex1D[lhsOffset[i]+2];
                const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
                const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];

                double Kr = 0;

                if(this->isIrreversible[i]==0)
                {
                    const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
                    double Kc = Kp*this->Pow_pByRT_SumVki[1];
                    Kc = std::max(Kc,1.4901171103413047e-8);
                    Kr = this->Kf_[i]/Kc;
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }

                const double CF = C[sl0]*C[sl1]*C[sl2];
                const double CR = C[sr0]*C[sr1];

                const double omegaf = (this->Kf_[i]*CF);
                const double omegar = (Kr*CR);

                double wf = omegaf + omegaf + omegaf;
                double wr = omegar + omegar;
                sumW += wf;
                sumWRateByCTot += (wf*wf);
                sumW += wr;
                sumWRateByCTot += (wr*wr);
            }
            else if(K==1)
            {
                const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
                const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
                const unsigned int sl2 = lhsSpeciesIndex1D[lhsOffset[i]+2];
                const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];

                double Kr = 0;

                if(this->isIrreversible[i]==0)
                {
                    const double Kp = (ExpNegGbyRT[sr0])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
                    double Kc = Kp*this->Pow_pByRT_SumVki[0];
                    Kc = std::max(Kc,1.4901171103413047e-8);
                    Kr = this->Kf_[i]/Kc;         
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }

                const double CF = C[sl0]*C[sl1]*C[sl2];
                const double CR = C[sr0];

                const double omegaf = (this->Kf_[i]*CF);
                const double omegar = (Kr*CR);

                double wf = omegaf + omegaf + omegaf;
                double wr = omegar;
                sumW += wf;
                sumWRateByCTot += (wf*wf);
                sumW += wr;
                sumWRateByCTot += (wr*wr);
            }
            else if(K==3)
            {
                const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
                const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
                const unsigned int sl2 = lhsSpeciesIndex1D[lhsOffset[i]+2];
                const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
                const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];
                const unsigned int sr2 = rhsSpeciesIndex1D[rhsOffset[i]+2];
                double Kr = 0;

                if(this->isIrreversible[i]==0)
                {
                    const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
                    double Kc = Kp;
                    Kc = std::max(Kc,1.4901171103413047e-8);
                    Kr = this->Kf_[i]/Kc;   
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }

                const double CF = C[sl0]*C[sl1]*C[sl2];
                const double CR = C[sr0]*C[sr1];

                const double omegaf = (this->Kf_[i]*CF);
                const double omegar = (Kr*CR);

                double wf = omegaf + omegaf + omegaf;
                double wr = omegar + omegar + omegar;
                sumW += wf;
                sumWRateByCTot += (wf*wf);
                sumW += wr;
                sumWRateByCTot += (wr*wr);
            }
        }
        if(J>3 || K>3)
        {
            //this->RFGI(i,this->Kf_[i],c,dNdtByV,&tmp_Exp[0]);
            double Kf = this->Kf_[i];
            int sumVki = 0;
            double Kp = 1.0;
            double CR = 1.0;
            double CF = 1.0;
            double Kc_ = 1.0;
            double Kr = 0;
            if(this->isIrreversible[i]==0)
            {
                for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
                {
                    const unsigned int si = this->lhsSpeciesIndex[i][j];
                    Kp = Kp / this->tmp_Exp[si];
                    sumVki = sumVki - 1;
                    CF = CF * C[si];
                }
                
                for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
                {
                    const unsigned int si = this->rhsSpeciesIndex[i][j];
                    Kp = Kp * this->tmp_Exp[si];
                    sumVki = sumVki + 1;
                    CR = CR * C[si];
                }

                Kc_ = Kp*this->Pow_pByRT_SumVki_I[sumVki];
                Kc_ = std::max(Kc_,1.49011611938476E-08);
                Kr = Kf/Kc_;
            }
            else if(this->isIrreversible[i]==2)
            {
                unsigned int JJ = i - this->Ikf[1] + this->Ikf[9];

                for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
                {
                    const unsigned int si = this->lhsSpeciesIndex[i][j];
                    CF = CF * C[si];
                }
                
                for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
                {
                    const unsigned int si = this->rhsSpeciesIndex[i][j];
                    CR = CR * C[si];
                }
                Kr = this->Kf_[JJ];
            }
            else
            {
                for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
                {
                    const unsigned int si = this->lhsSpeciesIndex[i][j];
                    CF = CF * C[si];
                }
                
                for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
                {
                    const unsigned int si = this->rhsSpeciesIndex[i][j];
                    CR = CR * C[si];
                }
            }

            const double omegaf = (Kf*CF);
            const double omegar = (Kr*CR);

            double wf = 0;
            double wr = 0;

            for(unsigned int j = 0; j < lhsSpeciesIndex[i].size();j++)
            {
                wf = wf + omegaf;
            }
            
            for(unsigned int j = 0; j < rhsSpeciesIndex[i].size();j++)
            {
                wr = wr + omegar;
            }       
            sumW += wf;
            sumWRateByCTot += (wf*wf);
            sumW += wr;
            sumWRateByCTot += (wr*wr);
        }
    }
}
