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
    double* __restrict__ ddNdtByVdcT
) const noexcept
{
    this->update_Pow_pByRT_SumVki(Temperature);
    this->update_Pow_pByRT_SumVki2(Temperature);


    for(unsigned int i = 0; i <this->n_Troe;i++)
    {
        unsigned int j0 = i + this->nSpecies;
        unsigned int j1 = i + this->nSpecies + this->n_Troe;
        unsigned int j2 = i + this->nSpecies + this->n_Troe*2;         
        this->tmp_Exp[j0] = -Temperature*this->invTsss_[i];            
        this->tmp_Exp[j1] = -this->Tss_[i]*invT; 
        this->tmp_Exp[j2] = -Temperature*this->invTs_[i];            
    }
    
    for(unsigned int i = 0; i <this->n_SRI;i++)
    {
        unsigned int j0 = i + this->nSpecies + this->n_Troe*3;
        unsigned int j1 = i + this->nSpecies + this->n_Troe*3 + this->n_SRI;
        this->tmp_Exp[j0] = -this->b_[i]*invT;
        this->tmp_Exp[j1] = -Temperature*this->invc_[i];            
    }   



    {
        unsigned int remain = this->tmp_ExpSize%4;
        for(unsigned int i = 0; i < this->tmp_ExpSize-remain;i=i+4)
        {
            __m256d tmp = _mm256_loadu_pd(&this->tmp_Exp[i]);
            tmp = vec256_expd(tmp);
            _mm256_storeu_pd(&this->tmp_Exp[i],tmp);
        }
        if(remain==1)
        {
            unsigned int i = this->tmp_ExpSize-1;
            this->tmp_Exp[i] = std::exp(this->tmp_Exp[i]);
        }
        else if(remain==2)
        {
            unsigned int i0 = this->tmp_ExpSize-2;
            unsigned int i1 = this->tmp_ExpSize-1;
            __m256d tmp = _mm256_setr_pd(tmp_Exp[i0],tmp_Exp[i1],0,0);
            tmp = vec256_expd(tmp);
            this->tmp_Exp[i0] = get_elem0(tmp);
            this->tmp_Exp[i1] = get_elem1(tmp);
        }
        else if(remain==3)
        {
            unsigned int i0 = this->tmp_ExpSize-3;
            unsigned int i1 = this->tmp_ExpSize-2;
            unsigned int i2 = this->tmp_ExpSize-1;
            __m256d tmp = _mm256_setr_pd(tmp_Exp[i0],tmp_Exp[i1],tmp_Exp[i2],0);
            tmp = vec256_expd(tmp);
            this->tmp_Exp[i0] = get_elem0(tmp);
            this->tmp_Exp[i1] = get_elem1(tmp);
            this->tmp_Exp[i2] = get_elem2(tmp);
        }
    }

    if(this->n_PlogReaction>0)
    {
        this->findPlogPressureRange(p);
        /*for(unsigned int i = 0; i< this->n_PlogReaction; i ++)
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
        }*/
    }

    {
        __m256d LogTv = _mm256_set1_pd(logT);
        __m256d InvTv = _mm256_set1_pd(invT);  
        const unsigned int end =       this->Ikf[11];
        unsigned int remain = (end-this->n_Temperature_Independent_Reaction)%4;
        unsigned int times = (end-this->n_Temperature_Independent_Reaction)/4;
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
            unsigned int i = end-1;
            this->Kf_[i] = this->A[i]*std::exp(this->beta[i+0]*logT-this->Ta[i+0]*invT);   
            this->dKfdT_[i+0] = this->Kf_[i+0]*(this->beta[i+0]+this->Ta[i+0]*invT)*invT;  
        }
        else if(remain==2)
        {
            unsigned int i0 = end-2;
            unsigned int i1 = end-1;    
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
            unsigned int i0 = end-3;
            unsigned int i1 = end-2;
            unsigned int i2 = end-1;    
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
        this->evalPlogPartialDerivative();
        /*for(unsigned int i = 0; i< this->n_PlogReaction; i ++)
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
        }*/

        //for(unsigned int i = 0; i< this->n_PlogReaction; i ++)
        //{
        //    std::cout<<this->Kf_[i+this->Ikf[6]]<<" "<<this->dKfdT_[i+this->Ikf[6]]<<std::endl;
        //}
        //std::exit(0);
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

    if(this->n_Lindemann>0)
    {
        this->evalLindemannPartialDerivative();
    }

    if(this->n_Troe>0)
    {
        this->evalTroePartialDerivative();
    }        


    if(this->n_SRI>0)
    {
        this->evalSRIPartialDerivative();
        /*for (unsigned int i = 0;i<this->n_SRI;i++)
        {
            const unsigned int j = this->SRIFO[i];
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
        }*/
    }

    this->updateJacobian11(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);
    this->updateJacobian12(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);
    this->updateJacobian13(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);
    this->updateJacobian21(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);
    this->updateJacobian22ReversibleReaction(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);
    this->updateJacobian22IrreversibleReaction(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);
    this->updateJacobian22NonEquilibriumReaction(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);
    this->updateJacobian23(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);
    this->updateJacobian31(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);
    this->updateJacobian32(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);
    this->updateJacobian33(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);
    this->updateJacobianGlobalNonIntegerReaction(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);
    this->updateJacobianGlobalIntegerReaction(c,dNdtByV,ddNdtByVdcT,tmp_Exp,dBdT);




    /*for(unsigned int z = 0; z < this->Ikf[7];z++)
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
            else 
            if(k==1)   {this->JF21(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else
            if(k==3)   {this->JF23(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
        }
        else 
        if(j==1)
        {
            if(k==2)        {this->JF12(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else 
            if(k==1)   {this->JF11(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else 
            if(k==3)   {this->JF13(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
        }
        else 
        if(j==3)
        {
            if(k==2)        {this->JF32(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else 
            if(k==1)   {this->JF31(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
            else 
            if(k==3)   {this->JF33(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
        }
        if(j>3 || k>3){this->JFGI(i,this->Kf_[z],this->dKfdT_[z],c,dNdtByV,ddNdtByVdcT,&this->tmp_Exp[0],dBdT);}
        //std::cout<<j<<" "<<k<<std::endl;
    }*/

/*for(int i = 0; i < this->nSpecies;i++)
{
    for(int j = 0; j < this->nSpecies-1;j++)
    {
        std::cout<<ddNdtByVdcT[i*(this->alignN)+j]<<" ";
    }
    int j = this->nSpecies-1;
    std::cout<<ddNdtByVdcT[i*(this->alignN)+j]<<std::endl;
}

for(int i = 0; i < this->nSpecies;i++)
{
    std::cout<<dNdtByV[i]<<std::endl;
}
std::exit(0);*/

}

