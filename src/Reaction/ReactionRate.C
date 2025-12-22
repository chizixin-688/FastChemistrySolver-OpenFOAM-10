
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
    {
        __m256d onev = _mm256_set1_pd(1);
        unsigned int remain = this->nSpecies%4;
        for(unsigned int i=0; i<this->nSpecies-remain; i=i+4)
        {
            __m256d r = load256d(&this->tmp_Exp[i]);
            __m256d invr = _mm256_div_pd(onev,r);
            store256d(&this->invNegGstdByRT[i],invr);
        }
        for(unsigned int i=this->nSpecies-remain; i<this->nSpecies; i=i+1)
        {
            this->invNegGstdByRT[i] = 1.0/this->tmp_Exp[i];
        }
    }
    //for(int i =0;i<this->nSpecies;i++)
    //{
    //    std::cout<<this->tmp_Exp[i]<<" "<<this->invNegGstdByRT[i]<<std::endl;
    //}
    //std::exit(0);
    if(this->n_PlogReaction>0)
    {
        this->findPlogPressureRange(p);
    }


    {
    
        __m256d LogT = _mm256_set1_pd(logT);
        __m256d InvT = _mm256_set1_pd(-invT);
        const unsigned int end = this->Ikf[11];
        unsigned int remain = (end-this->n_Temperature_Independent_Reaction)%4;
        unsigned int times = (end-this->n_Temperature_Independent_Reaction)/4;
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
            unsigned int i = end-1;
            this->Kf_[i] = this->A[i]*std::exp(this->beta[i+0]*logT-this->Ta[i+0]*invT);   
        }
        else if(remain==2)
        {
            unsigned int i0 = end-2;
            unsigned int i1 = end-1;
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
            unsigned int i0 = end-3;
            unsigned int i1 = end-2;
            unsigned int i2 = end-1;
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
        this->evalPlogRateConstant();
        /*for(unsigned int i = 0; i< this->n_PlogReaction; i ++)
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
        }*/

            /*double logT = this->logT;
            double invT = this->invT;

            for(unsigned int i=0; i<this->ActivePlogReactionNumber; i=i+1)
            {
                unsigned j = ActivePlogReactionIndex[i];

                unsigned index = this->Pindex[j];
                double weight = (this->logP - this->logPi[j][index])*this->rDeltaP_[j][index];
                double Kf0 = this->Kf_[j+this->Ikf[6]];
                double Kf1 = this->Kf_[j+this->Ikf[11]];

                double logA1A0 = this->logAPlog[j][index+1] - this->logAPlog[j][index];
                double beta1beta0 = this->betaPlog[j][index+1] - this->betaPlog[j][index];
                double Ta0Ta1 = this->TaPlog[j][index] - this->TaPlog[j][index+1];
                double logk1k0 = logA1A0 + beta1beta0*logT + Ta0Ta1*invT;
                double k = Kf0*std::exp(weight*logk1k0);
                //this->Kf_[j+this->Ikf[6]] = k;
                this->Kf_[0] = k;
            }*/

        //std::setprecision(20);
        //for(unsigned int i = 0; i< this->n_PlogReaction; i ++)
        //{
        //    std::cout<<std::setprecision(20)<<this->Kf_[i+this->Ikf[6]]<<std::endl;
        //}
        //std::cout<<"debug"<<std::endl;
    }
    //std::exit(0);

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
        this->evalSRIRateConstant();
    }

//auto GlobalTimeStart = std::chrono::high_resolution_clock::now();

//for(unsigned int i =0;i<10000000;i++)
//{
    for(auto funcPtr : RFptr)
    {
        (this->*funcPtr)(c,dNdtByV,tmp_Exp);
    }
    
    //this->update11Reaction(c,dNdtByV,tmp_Exp);
    /*this->RF11RR(c,dNdtByV,tmp_Exp);
    this->RF11IR(c,dNdtByV,tmp_Exp);
    this->RF11NER(c,dNdtByV,tmp_Exp);

    //this->update12Reaction(c,dNdtByV,tmp_Exp);
    this->RF12RR(c,dNdtByV,tmp_Exp);
    this->RF12IR(c,dNdtByV,tmp_Exp);
    this->RF12NER(c,dNdtByV,tmp_Exp);

    //this->update13Reaction(c,dNdtByV,tmp_Exp);

    this->RF13RR(c,dNdtByV,tmp_Exp);
    this->RF13IR(c,dNdtByV,tmp_Exp);
    this->RF13NER(c,dNdtByV,tmp_Exp);
    //std::exit(0);
    //this->update21Reaction(c,dNdtByV,tmp_Exp);

    this->RF21RR(c,dNdtByV,tmp_Exp);
    this->RF21IR(c,dNdtByV,tmp_Exp);
    this->RF21NER(c,dNdtByV,tmp_Exp);
    //this->update22Reaction(c,dNdtByV,tmp_Exp);

    
    this->RF22RR(c,dNdtByV,tmp_Exp);
    this->RF22IR(c,dNdtByV,tmp_Exp);    
    this->RF22NER(c,dNdtByV,tmp_Exp);

    this->RF23RR(c,dNdtByV,tmp_Exp);
    this->RF23IR(c,dNdtByV,tmp_Exp);    
    this->RF23NER(c,dNdtByV,tmp_Exp);

    //this->update23Reaction(c,dNdtByV,tmp_Exp);

    //this->update31Reaction(c,dNdtByV,tmp_Exp);
    this->RF31RR(c,dNdtByV,tmp_Exp);
    this->RF31IR(c,dNdtByV,tmp_Exp);    
    this->RF31NER(c,dNdtByV,tmp_Exp);

    //this->update32Reaction(c,dNdtByV,tmp_Exp);
    this->RF32RR(c,dNdtByV,tmp_Exp);
    this->RF32IR(c,dNdtByV,tmp_Exp);    
    this->RF32NER(c,dNdtByV,tmp_Exp);


    //this->update33Reaction(c,dNdtByV,tmp_Exp);
    this->RF33RR(c,dNdtByV,tmp_Exp);
    this->RF33IR(c,dNdtByV,tmp_Exp);    
    this->RF33NER(c,dNdtByV,tmp_Exp);
    this->updateGlobalIntegerReaction(c,dNdtByV,tmp_Exp);
    this->updateGlobalNonIntegerReaction(c,dNdtByV,tmp_Exp);*/
    
//}
//auto duration = (std::chrono::duration_cast<std::chrono::microseconds>
//(std::chrono::high_resolution_clock::now()-GlobalTimeStart));
//    std::cout<<duration.count()<<std::endl;
//    std::exit(0);




    
    //return;

    /*for (unsigned int i = 0; i < this->Ikf[7]; i++) 
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
    }*/
//std::exit(0);
}

