
#include "OptReaction.H"

void 
OptReaction::dNdtByV
(
    double p,
    double Temperature,
    double* __restrict__ c,
    double* __restrict__ dNdtByV
) const noexcept
{

    this->update_Pow_pByRT_SumVki(Temperature);
    this->update_Pow_pByRT_SumVki2(Temperature);
    
    {
        for(size_t i = 0; i <this->n_Troe;i++)
        {
            size_t j0 = i + this->nSpecies;
            size_t j1 = i + this->nSpecies + this->n_Troe;
            size_t j2 = i + this->nSpecies + this->n_Troe*2;         
            this->tmp_Exp[j0] = -Temperature*this->invTsss_[i];
            this->tmp_Exp[j1] = -this->Tss_[i]*invT;    
            this->tmp_Exp[j2] = -Temperature*this->invTs_[i];
        }
    }

    {
        for(size_t i = 0; i <this->n_SRI;i++)
        {
            size_t j0 = i + this->nSpecies + this->n_Troe*3;
            size_t j1 = i + this->nSpecies + this->n_Troe*3 + this->n_SRI;
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

    if(this->n_Lindemann)
    {
        this->evalLindemannRateConstant();
    }

    if(this->n_Troe)
    {
        this->evalTroeRateConstant();
    }

    if(this->n_SRI)
    {
        for (unsigned int i = 0;i<this->n_SRI;i++)
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
        }
    }



    this->update11Reaction(c,dNdtByV,tmp_Exp);
    this->update12Reaction(c,dNdtByV,tmp_Exp);
    this->update13Reaction(c,dNdtByV,tmp_Exp);

    //this->update21Reaction(c,dNdtByV,tmp_Exp);

    this->update21ReversibleReaction(c,dNdtByV,tmp_Exp);
    this->update21IrreversibleReaction(c,dNdtByV,tmp_Exp);    
    this->update21NonEquilibriumReaction(c,dNdtByV,tmp_Exp);    

    this->update22ReversibleReaction(c,dNdtByV,tmp_Exp);
    this->update22IrreversibleReaction(c,dNdtByV,tmp_Exp);    
    this->update22NonEquilibriumReaction(c,dNdtByV,tmp_Exp);

    this->update23ReversibleReaction(c,dNdtByV,tmp_Exp);
    this->update23IrreversibleReaction(c,dNdtByV,tmp_Exp);    
    this->update23NonEquilibriumReaction(c,dNdtByV,tmp_Exp);

    this->update31Reaction(c,dNdtByV,tmp_Exp);
    this->update32Reaction(c,dNdtByV,tmp_Exp);
    this->update33Reaction(c,dNdtByV,tmp_Exp);
    this->updateGlobalIntegerReaction(c,dNdtByV,tmp_Exp);
    this->updateGlobalNonIntegerReaction(c,dNdtByV,tmp_Exp);



    
    return;

    for (unsigned int i = 0; i < this->Ikf[7]; i++) 
    {
        if(this->isGlobal[i]==1)
        {
            this->RFGNI(i,this->Kf_[i],c,dNdtByV,&tmp_Exp[0]);            
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
                    //Kc = std::max(Kc,1.4901171103413047e-8);
                    Kc = std::max(Kc,KcLimiter);
                    Kr = this->Kf_[i]/Kc;         
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int l = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[l];
                }

                const double CF = c[sl0]*c[sl1];
                const double CR = c[sr0]*c[sr1];
                const double q = (this->Kf_[i]*CF) - (Kr*CR);
//std::cout<<this->Kf_[i]<<" "<<CF<<" "<<Kr<<" "<<CR<<std::endl;
                dNdtByV[sl0] = dNdtByV[sl0] - q; 
                dNdtByV[sl1] = dNdtByV[sl1] - q;
                dNdtByV[sr0] = dNdtByV[sr0] + q;
                dNdtByV[sr1] = dNdtByV[sr1] + q; 
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
                    //Kc = std::max(Kc,1.4901171103413047e-8);
                    Kc = std::max(Kc,KcLimiter);
                    Kr = this->Kf_[i]/Kc;        
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }
                const double CF = c[sl0]*c[sl1];
                const double CR = c[sr0];
                const double q = (this->Kf_[i]*CF) - (Kr*CR);
                dNdtByV[sl0] = dNdtByV[sl0] - q; 
                dNdtByV[sl1] = dNdtByV[sl1] - q;
                dNdtByV[sr0] = dNdtByV[sr0] + q;   
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
                    //Kc = std::max(Kc,1.4901171103413047e-8);
                    Kc = std::max(Kc,KcLimiter);
                    Kr = this->Kf_[i]/Kc;    
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }

                const double CF = c[sl0]*c[sl1];
                const double CR = c[sr0]*c[sr1]*c[sr2];
                const double q = (this->Kf_[i]*CF) - (Kr*CR);


                dNdtByV[sl0] = dNdtByV[sl0] - q;
                dNdtByV[sl1] = dNdtByV[sl1] - q;
                dNdtByV[sr0] = dNdtByV[sr0] + q;
                dNdtByV[sr1] = dNdtByV[sr1] + q;
                dNdtByV[sr2] = dNdtByV[sr2] + q;
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
                    //Kc = std::max(Kc,1.4901171103413047e-8);
                    Kc = std::max(Kc,KcLimiter);
                    Kr = this->Kf_[i]/Kc;       
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }

                const double CF = c[sl0];
                const double CR = c[sr0]*c[sr1];

                const double q = (this->Kf_[i]*CF) - (Kr*CR);

            
                dNdtByV[sl0] = dNdtByV[sl0] - q; 
                dNdtByV[sr0] = dNdtByV[sr0] + q;
                dNdtByV[sr1] = dNdtByV[sr1] + q;
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
                    //Kc = std::max(Kc,1.4901171103413047e-8);
                    Kc = std::max(Kc,KcLimiter);
                    Kr = this->Kf_[i]/Kc;       
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }
                const double CF = c[sl0];
                const double CR = c[sr0];
                const double q = (this->Kf_[i]*CF) - (Kr*CR);
                dNdtByV[sl0] = dNdtByV[sl0] - q; 
                dNdtByV[sr0] = dNdtByV[sr0] + q;
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
                    //Kc = std::max(Kc,1.4901171103413047e-8);
                    
                    Kc = std::max(Kc,KcLimiter);
                    Kr = this->Kf_[i]/Kc;         
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }

                const double CF = c[sl0];
                const double CR = c[sr0]*c[sr1]*c[sr2];

                const double q = (this->Kf_[i]*CF) - (Kr*CR);

                dNdtByV[sl0] = dNdtByV[sl0] - q;
                dNdtByV[sr0] = dNdtByV[sr0] + q;
                dNdtByV[sr1] = dNdtByV[sr1] + q;
                dNdtByV[sr2] = dNdtByV[sr2] + q; 
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
                    //Kc = std::max(Kc,1.4901171103413047e-8);
                    Kc = std::max(Kc,KcLimiter);
                    Kr = this->Kf_[i]/Kc;
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }

                const double CF = c[sl0]*c[sl1]*c[sl2];
                const double CR = c[sr0]*c[sr1];
                const double q = (this->Kf_[i]*CF) - (Kr*CR);

                dNdtByV[sl0] = dNdtByV[sl0] - q; 
                dNdtByV[sl1] = dNdtByV[sl1] - q; 
                dNdtByV[sl2] = dNdtByV[sl2] - q; 
                dNdtByV[sr0] = dNdtByV[sr0] + q;
                dNdtByV[sr1] = dNdtByV[sr1] + q;
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
                    //Kc = std::max(Kc,1.4901171103413047e-8);
                    Kc = std::max(Kc,KcLimiter);
                    Kr = this->Kf_[i]/Kc;         
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }

                const double CF = c[sl0]*c[sl1]*c[sl2];
                const double CR = c[sr0];

                const double q = (this->Kf_[i]*CF) - (Kr*CR);

                dNdtByV[sl0] = dNdtByV[sl0] - q; 
                dNdtByV[sl1] = dNdtByV[sl1] - q; 
                dNdtByV[sl2] = dNdtByV[sl2] - q; 
                dNdtByV[sr0] = dNdtByV[sr0] + q;
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
                    //Kc = std::max(Kc,1.4901171103413047e-8);
                    Kc = std::max(Kc,KcLimiter);
                    Kr = this->Kf_[i]/Kc;   
                }
                else if(this->isIrreversible[i]==2)
                {
                    unsigned int j = i - this->Ikf[1] + this->Ikf[9];
                    Kr = this->Kf_[j];
                }

                const double CF = c[sl0]*c[sl1]*c[sl2];
                const double CR = c[sr0]*c[sr1];

                const double q = (this->Kf_[i]*CF) - (Kr*CR);

                dNdtByV[sl0] = dNdtByV[sl0] - q; 
                dNdtByV[sl1] = dNdtByV[sl1] - q; 
                dNdtByV[sl2] = dNdtByV[sl2] - q; 
                dNdtByV[sr0] = dNdtByV[sr0] + q;
                dNdtByV[sr1] = dNdtByV[sr1] + q;
                dNdtByV[sr2] = dNdtByV[sr2] + q;
            }
        }
        if(J>3 || K>3)
        {
            this->RFGI(i,this->Kf_[i],c,dNdtByV,&tmp_Exp[0]);
        }
    }
//std::exit(0);
}

