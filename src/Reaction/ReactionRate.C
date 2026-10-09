/*---------------------------------------------------------------------------*\
  Description
      Computing molar reaction rate (kmol/m^3/s)

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/


//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"
#include <iostream>
#include <set>
//=============================================================================//
void 
FastChemistry::OptReaction::dNdtByV
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
            this->tmp_Exp[i0] = get0(tmp);
            this->tmp_Exp[i1] = get1(tmp);
        }
        else if(remain==3)
        {
            size_t i0 = this->tmp_ExpSize-3;
            size_t i1 = this->tmp_ExpSize-2;
            size_t i2 = this->tmp_ExpSize-1;

            __m256d tmp = _mm256_setr_pd(tmp_Exp[i0],tmp_Exp[i1],tmp_Exp[i2],0);
            tmp = vec256_expd(tmp);
            this->tmp_Exp[i0] = get0(tmp);
            this->tmp_Exp[i1] = get1(tmp);
            this->tmp_Exp[i2] = get2(tmp);
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

    if(this->n_PlogReaction>0)
    {
        this->findPlogPressureRange(p);
    }


    {
        double* __restrict__ Kf_ptr = &this->Kf_[0];
        const double* __restrict__ A_ptr = &this->A[0];
        const double* __restrict__ beta_ptr = &this->beta[0];
        const double* __restrict__ Ta_ptr = &this->Ta[0];

        __m256d LogT = _mm256_set1_pd(logT);
        __m256d InvT = _mm256_set1_pd(-invT);
        const unsigned int end = this->Ikf[11];
        unsigned int remain = (end-this->n_Temperature_Independent_Reaction)%4;
        unsigned int times = (end-this->n_Temperature_Independent_Reaction)/4;
        for(unsigned int z = 0; z <times;z=z+1)
        {
            unsigned int i = z*4 + this->n_Temperature_Independent_Reaction;
            __m256d A_ = _mm256_loadu_pd(&A_ptr[i]);
            __m256d beta_ = _mm256_loadu_pd(&beta_ptr[i]);
            __m256d Ta_ = _mm256_loadu_pd(&Ta_ptr[i]);
            __m256d Kf = _mm256_mul_pd(Ta_,InvT);
            Kf = fmadd256d(beta_,LogT,Kf);
            Kf = vec256_expd(Kf);
            Kf = _mm256_mul_pd(A_,Kf);
            _mm256_storeu_pd(&Kf_ptr[i],Kf);
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
            Kf = fmadd256d(beta_,LogT,Kf);
            Kf = vec256_expd(Kf);
            Kf = _mm256_mul_pd(A_,Kf);
            this->Kf_[i0] = get0(Kf);
            this->Kf_[i1] = get1(Kf);
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
            Kf = fmadd256d(beta_,LogT,Kf);
            Kf = vec256_expd(Kf);
            Kf = _mm256_mul_pd(A_,Kf);
            this->Kf_[i0] = get0(Kf);
            this->Kf_[i1] = get1(Kf); 
            this->Kf_[i2] = get2(Kf); 
        }
    }

    {
        unsigned int Tremain = (this->Itbr[4])%8;
        unsigned int Spremain = this->AlignSpecies%8;
        for(unsigned int i = 0; i < this->Itbr[4]-Tremain; i=i+8)
        {
            __m256d M0 = _mm256_setzero_pd();
            __m256d M1 = _mm256_setzero_pd();
            __m256d M2 = _mm256_setzero_pd();
            __m256d M3 = _mm256_setzero_pd();
            __m256d M4 = _mm256_setzero_pd();
            __m256d M5 = _mm256_setzero_pd();
            __m256d M6 = _mm256_setzero_pd();
            __m256d M7 = _mm256_setzero_pd();
            double* __restrict__ TBF1DRowi0 = &ThirdBodyFactor1D[(i+0)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi1 = &ThirdBodyFactor1D[(i+1)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi2 = &ThirdBodyFactor1D[(i+2)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi3 = &ThirdBodyFactor1D[(i+3)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi4 = &ThirdBodyFactor1D[(i+4)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi5 = &ThirdBodyFactor1D[(i+5)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi6 = &ThirdBodyFactor1D[(i+6)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi7 = &ThirdBodyFactor1D[(i+7)*this->AlignSpecies];
            for(unsigned int j  = 0;j<this->AlignSpecies-Spremain;j=j+8)
            {
                __m256d C03 = load256d(&c[j+0]);
                M0 = fmadd256d(load256d(&TBF1DRowi0[j+0]),C03,M0);
                M1 = fmadd256d(load256d(&TBF1DRowi1[j+0]),C03,M1);
                M2 = fmadd256d(load256d(&TBF1DRowi2[j+0]),C03,M2);
                M3 = fmadd256d(load256d(&TBF1DRowi3[j+0]),C03,M3);
                M4 = fmadd256d(load256d(&TBF1DRowi4[j+0]),C03,M4);
                M5 = fmadd256d(load256d(&TBF1DRowi5[j+0]),C03,M5);
                M6 = fmadd256d(load256d(&TBF1DRowi6[j+0]),C03,M6);
                M7 = fmadd256d(load256d(&TBF1DRowi7[j+0]),C03,M7);

                __m256d C47 = load256d(&c[j+4]);
                M0 = fmadd256d(load256d(&TBF1DRowi0[j+4]),C47,M0);
                M1 = fmadd256d(load256d(&TBF1DRowi1[j+4]),C47,M1);
                M2 = fmadd256d(load256d(&TBF1DRowi2[j+4]),C47,M2);
                M3 = fmadd256d(load256d(&TBF1DRowi3[j+4]),C47,M3);
                M4 = fmadd256d(load256d(&TBF1DRowi4[j+4]),C47,M4);
                M5 = fmadd256d(load256d(&TBF1DRowi5[j+4]),C47,M5);
                M6 = fmadd256d(load256d(&TBF1DRowi6[j+4]),C47,M6);
                M7 = fmadd256d(load256d(&TBF1DRowi7[j+4]),C47,M7);
            }
            if(Spremain==4)
            {
                unsigned int j = this->AlignSpecies-4;
                __m256d C03 = load256d(&c[j+0]);
                M0 = fmadd256d(load256d(&TBF1DRowi0[j+0]),C03,M0);
                M1 = fmadd256d(load256d(&TBF1DRowi1[j+0]),C03,M1);
                M2 = fmadd256d(load256d(&TBF1DRowi2[j+0]),C03,M2);
                M3 = fmadd256d(load256d(&TBF1DRowi3[j+0]),C03,M3);
                M4 = fmadd256d(load256d(&TBF1DRowi4[j+0]),C03,M4);
                M5 = fmadd256d(load256d(&TBF1DRowi5[j+0]),C03,M5);
                M6 = fmadd256d(load256d(&TBF1DRowi6[j+0]),C03,M6);
                M7 = fmadd256d(load256d(&TBF1DRowi7[j+0]),C03,M7);
            }
            __m256d M03 = hsum4x4(M0,M1,M2,M3);
            __m256d M47 = hsum4x4(M4,M5,M6,M7);
           _mm256_storeu_pd(&this->tmp_M[i+0],M03);
           _mm256_storeu_pd(&this->tmp_M[i+4],M47);
        }
        if(Tremain>=4)
        {
            unsigned int i =(this->Itbr[4]) -Tremain;
            Tremain = Tremain - 4;
            __m256d M0 = _mm256_setzero_pd();
            __m256d M1 = _mm256_setzero_pd();
            __m256d M2 = _mm256_setzero_pd();
            __m256d M3 = _mm256_setzero_pd();
            double* __restrict__ TBF1DRowi0 = &ThirdBodyFactor1D[(i+0)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi1 = &ThirdBodyFactor1D[(i+1)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi2 = &ThirdBodyFactor1D[(i+2)*this->AlignSpecies];
            double* __restrict__ TBF1DRowi3 = &ThirdBodyFactor1D[(i+3)*this->AlignSpecies];
            for(unsigned int j  = 0;j<this->AlignSpecies-Spremain;j=j+8)
            {
                __m256d C03 = load256d(&c[j+0]);
                M0 = fmadd256d(load256d(&TBF1DRowi0[j+0]),C03,M0);
                M1 = fmadd256d(load256d(&TBF1DRowi1[j+0]),C03,M1);
                M2 = fmadd256d(load256d(&TBF1DRowi2[j+0]),C03,M2);
                M3 = fmadd256d(load256d(&TBF1DRowi3[j+0]),C03,M3);

                __m256d C47 = load256d(&c[j+4]);
                M0 = fmadd256d(load256d(&TBF1DRowi0[j+4]),C47,M0);
                M1 = fmadd256d(load256d(&TBF1DRowi1[j+4]),C47,M1);
                M2 = fmadd256d(load256d(&TBF1DRowi2[j+4]),C47,M2);
                M3 = fmadd256d(load256d(&TBF1DRowi3[j+4]),C47,M3);
            }
            if(Spremain==4)
            {
                unsigned int j = this->AlignSpecies-4;
                __m256d C03 = load256d(&c[j+0]);
                M0 = fmadd256d(load256d(&TBF1DRowi0[j+0]),C03,M0);
                M1 = fmadd256d(load256d(&TBF1DRowi1[j+0]),C03,M1);
                M2 = fmadd256d(load256d(&TBF1DRowi2[j+0]),C03,M2);
                M3 = fmadd256d(load256d(&TBF1DRowi3[j+0]),C03,M3);
            }
            __m256d M03 = hsum4x4(M0,M1,M2,M3);
            _mm256_storeu_pd(&this->tmp_M[i+0],M03);
            
        }
        if(Tremain==3)
        {
            unsigned int i =(this->Itbr[4]) -3;
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
                arrM_0 = fmadd256d(Factor0,C_,arrM_0);
                arrM_1 = fmadd256d(Factor1,C_,arrM_1);
                arrM_2 = fmadd256d(Factor2,C_,arrM_2);
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
            unsigned int i =(this->Itbr[4]) -2;
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
                arrM_0 = fmadd256d(Factor0,C_,arrM_0);
                arrM_1 = fmadd256d(Factor1,C_,arrM_1);
            }

            M0 = M0 + hsum4(arrM_0);
            M1 = M1 + hsum4(arrM_1);

            this->tmp_M[i+0] = M0;
            this->tmp_M[i+1] = M1;
        }
        else if(Tremain==1)
        {
            unsigned int i =(this->Itbr[4]) -1;
            double* __restrict__ TBF1DRowi0 = &ThirdBodyFactor1D[(i+0)*this->AlignSpecies];            
            double M0 = 0;
            __m256d arrM_0 = _mm256_setzero_pd();
            for(unsigned int j  = 0;j<this->AlignSpecies;j=j+4)
            {
                __m256d Factor0 = _mm256_loadu_pd(&TBF1DRowi0[j+0]);
                __m256d C_ = _mm256_loadu_pd(&c[j+0]);
                arrM_0 = fmadd256d(Factor0,C_,arrM_0);
            }

            M0 = M0 + hsum4(arrM_0);

            this->tmp_M[i+0] = M0;
        }
    }

    if(this->n_PlogReaction>0)
    {
        this->evalPlogRateConstant();
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
            this->Kf_[Ikf[2]+i] = this->Kf_[Ikf[2]+i]*Mfwd;
            this->Kf_[Ikf[10]+i] = this->Kf_[Ikf[10]+i]*Mfwd;
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

    for(auto funcPtr: RFptr)
    {
        (this->*funcPtr)(c,dNdtByV,tmp_Exp);
    }
}

