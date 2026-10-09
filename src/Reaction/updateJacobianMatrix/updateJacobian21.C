/*---------------------------------------------------------------------------*\
  Description
      Update the net rates of production (dNdtByV) and the molar-concentration
      based Jacobian matrix (ddNdtByVdcTp) for two-reactant / one-product
      (2-1) reactions, e.g. A + A -> B and A + B -> C.

      The Jacobian matrix is stored row-wise with one column per species
      plus one extra column for the temperature derivative (column index
      nSpecies).

      RR : reversible reaction; the reverse rate constant is obtained from
           the equilibrium constant
      IR : irreversible reaction; the reverse rate constant is zero
      NER: non-equilibrium reaction; the reverse rate constant is computed
           with the Arrhenius form instead of the equilibrium constant

  Functions
      JF21RR : reversible 2-1 reactions, incl. duplicate reactants (A+A=B)
      JF21IR : irreversible 2-1 reactions, incl. duplicate reactants (A+A=B)
      JF21NER: non-equilibrium 2-1 reactions

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"
#include <iostream>
#include <algorithm>
//=============================================================================//
void FastChemistry::OptReaction::JF21RR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
) const noexcept
{
    const double PowPByRTSumVki = this->Pow_pByRT_SumVki[3];
    const double invKcLimiter = FastChemistry::invKcLimiter;
    std::size_t endRR21 = this->RR21Size;
    std::size_t remainRR21 = endRR21 % 4;
    for(std::size_t k=0; k<endRR21-remainRR21; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->RR21AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->RR21AllIndex[(k+0)*4+1];
        const unsigned int sl1a = this->RR21AllIndex[(k+0)*4+2];
        const unsigned int sr0a = this->RR21AllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double C_sr0a = C[sr0a];

        const double Kf0 = this->Kf_[i0];
        const double Kf0dCfdC0 = Kf0*C_sl1a;
        const double Kf0dCfdC1 = Kf0*C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0;
        const double J_sl1a_sr0a_increment = J_sl1a_ptr[sr0a] + Kr0;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sl1a_ptr[sr0a] = J_sl1a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a]) - (dBdT[sl0a] + dBdT[sl1a]);
        const double dKcdTByKc0 = sumVdBdT0 + invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C_sl0a*C_sl1a;
        const double CR0 = C_sr0a;
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF0-invKc0*CR0);
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d WdMdC4 = _mm256_mul_pd(Wv,dMdC4);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);

                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1a_ptr[j+0]),WdMdC0);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1a_ptr[j+4]),WdMdC4);
                store256d(&J_sl1a_ptr[j+0],sl1v0);
                store256d(&J_sl1a_ptr[j+4],sl1v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1a_ptr[j+0]),WdMdC0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl1a_ptr[j+0],sl1v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR21AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->RR21AllIndex[(k+1)*4+1];
        const unsigned int sl1b = this->RR21AllIndex[(k+1)*4+2];
        const unsigned int sr0b = this->RR21AllIndex[(k+1)*4+3];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sl1b_ptr = &ddNdtByVdcTp[sl1b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        const double C_sl1b = C[sl1b];
        const double C_sl0b = C[sl0b];
        const double C_sr0b = C[sr0b];

        const double Kf1 = this->Kf_[i1];
        const double Kf1dCfdC0 = Kf1*C_sl1b;
        const double Kf1dCfdC1 = Kf1*C_sl0b;
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1dCfdC0;
        const double J_sl1b_sl0b_increment = J_sl1b_ptr[sl0b] - Kf1dCfdC0;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sl1b_ptr[sl0b] = J_sl1b_sl0b_increment;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1dCfdC0;
        const double J_sl0b_sl1b_increment = J_sl0b_ptr[sl1b] - Kf1dCfdC1;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sl0b_ptr[sl1b] = J_sl0b_sl1b_increment;
        const double J_sl1b_sl1b_increment = J_sl1b_ptr[sl1b] - Kf1dCfdC1;
        const double J_sr0b_sl1b_increment = J_sr0b_ptr[sl1b] + Kf1dCfdC1;
        J_sl1b_ptr[sl1b] = J_sl1b_sl1b_increment;
        J_sr0b_ptr[sl1b] = J_sr0b_sl1b_increment;

        const double invKp1 = (invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1;
        const double J_sl1b_sr0b_increment = J_sl1b_ptr[sr0b] + Kr1;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sl1b_ptr[sr0b] = J_sl1b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double sumVdBdT1 = (dBdT[sr0b]) - (dBdT[sl0b] + dBdT[sl1b]);
        const double dKcdTByKc1 = sumVdBdT1 + invT;
        const double dKrdT1 = (dKfdT1*invKc1 - (invKc1 < invKcLimiter ? Kr1*dKcdTByKc1 : 0));
        const double CF1 = C_sl0b*C_sl1b;
        const double CR1 = C_sr0b;
        const double dqdT1 = (dKfdT1*CF1) - (dKrdT1*CR1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sl1b_T_increment = J_sl1b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sl1b_ptr[this->nSpecies] = J_sl1b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF1-invKc1*CR1);
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d WdMdC4 = _mm256_mul_pd(Wv,dMdC4);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1b_ptr[j+0]),WdMdC0);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl1b_ptr[j+0],sl1v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+4]),WdMdC4);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1b_ptr[j+4]),WdMdC4);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+4]),WdMdC4);
                store256d(&J_sl1b_ptr[j+4],sl1v4);
                store256d(&J_sr0b_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1b_ptr[j+0]),WdMdC0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl1b_ptr[j+0],sl1v0);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
            }
        }
        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR21AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->RR21AllIndex[(k+2)*4+1];
        const unsigned int sl1c = this->RR21AllIndex[(k+2)*4+2];
        const unsigned int sr0c = this->RR21AllIndex[(k+2)*4+3];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sl1c_ptr = &ddNdtByVdcTp[sl1c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        const double C_sl1c = C[sl1c];
        const double C_sl0c = C[sl0c];
        const double C_sr0c = C[sr0c];

        const double Kf2 = this->Kf_[i2];
        const double Kf2dCfdC0 = Kf2*C_sl1c;
        const double Kf2dCfdC1 = Kf2*C_sl0c;
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2dCfdC0;
        const double J_sl1c_sl0c_increment = J_sl1c_ptr[sl0c] - Kf2dCfdC0;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sl1c_ptr[sl0c] = J_sl1c_sl0c_increment;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2dCfdC0;
        const double J_sl0c_sl1c_increment = J_sl0c_ptr[sl1c] - Kf2dCfdC1;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sl0c_ptr[sl1c] = J_sl0c_sl1c_increment;
        const double J_sl1c_sl1c_increment = J_sl1c_ptr[sl1c] - Kf2dCfdC1;
        const double J_sr0c_sl1c_increment = J_sr0c_ptr[sl1c] + Kf2dCfdC1;
        J_sl1c_ptr[sl1c] = J_sl1c_sl1c_increment;
        J_sr0c_ptr[sl1c] = J_sr0c_sl1c_increment;

        const double invKp2 = (invNegGstdByRT[sr0c])*(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]);
        double invKc2 = invKp2*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2;
        const double J_sl1c_sr0c_increment = J_sl1c_ptr[sr0c] + Kr2;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sl1c_ptr[sr0c] = J_sl1c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double sumVdBdT2 = (dBdT[sr0c]) - (dBdT[sl0c] + dBdT[sl1c]);
        const double dKcdTByKc2 = sumVdBdT2 + invT;
        const double dKrdT2 = (dKfdT2*invKc2 - (invKc2 < invKcLimiter ? Kr2*dKcdTByKc2 : 0));
        const double CF2 = C_sl0c*C_sl1c;
        const double CR2 = C_sr0c;
        const double dqdT2 = (dKfdT2*CF2) - (dKrdT2*CR2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sl1c_T_increment = J_sl1c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sl1c_ptr[this->nSpecies] = J_sl1c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF2-invKc2*CR2);
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d WdMdC4 = _mm256_mul_pd(Wv,dMdC4);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1c_ptr[j+0]),WdMdC0);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl1c_ptr[j+0],sl1v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+4]),WdMdC4);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1c_ptr[j+4]),WdMdC4);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+4]),WdMdC4);
                store256d(&J_sl1c_ptr[j+4],sl1v4);
                store256d(&J_sr0c_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1c_ptr[j+0]),WdMdC0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl1c_ptr[j+0],sl1v0);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
            }
        }
        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR21AllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->RR21AllIndex[(k+3)*4+1];
        const unsigned int sl1d = this->RR21AllIndex[(k+3)*4+2];
        const unsigned int sr0d = this->RR21AllIndex[(k+3)*4+3];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sl1d_ptr = &ddNdtByVdcTp[sl1d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        const double C_sl1d = C[sl1d];
        const double C_sl0d = C[sl0d];
        const double C_sr0d = C[sr0d];

        const double Kf3 = this->Kf_[i3];
        const double Kf3dCfdC0 = Kf3*C_sl1d;
        const double Kf3dCfdC1 = Kf3*C_sl0d;
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3dCfdC0;
        const double J_sl1d_sl0d_increment = J_sl1d_ptr[sl0d] - Kf3dCfdC0;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sl1d_ptr[sl0d] = J_sl1d_sl0d_increment;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3dCfdC0;
        const double J_sl0d_sl1d_increment = J_sl0d_ptr[sl1d] - Kf3dCfdC1;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sl0d_ptr[sl1d] = J_sl0d_sl1d_increment;
        const double J_sl1d_sl1d_increment = J_sl1d_ptr[sl1d] - Kf3dCfdC1;
        const double J_sr0d_sl1d_increment = J_sr0d_ptr[sl1d] + Kf3dCfdC1;
        J_sl1d_ptr[sl1d] = J_sl1d_sl1d_increment;
        J_sr0d_ptr[sl1d] = J_sr0d_sl1d_increment;

        const double invKp3 = (invNegGstdByRT[sr0d])*(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]);
        double invKc3 = invKp3*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3;
        const double J_sl1d_sr0d_increment = J_sl1d_ptr[sr0d] + Kr3;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sl1d_ptr[sr0d] = J_sl1d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double sumVdBdT3 = (dBdT[sr0d]) - (dBdT[sl0d] + dBdT[sl1d]);
        const double dKcdTByKc3 = sumVdBdT3 + invT;
        const double dKrdT3 = (dKfdT3*invKc3 - (invKc3 < invKcLimiter ? Kr3*dKcdTByKc3 : 0));
        const double CF3 = C_sl0d*C_sl1d;
        const double CR3 = C_sr0d;
        const double dqdT3 = (dKfdT3*CF3) - (dKrdT3*CR3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sl1d_T_increment = J_sl1d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sl1d_ptr[this->nSpecies] = J_sl1d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF3-invKc3*CR3);
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d WdMdC4 = _mm256_mul_pd(Wv,dMdC4);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1d_ptr[j+0]),WdMdC0);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl1d_ptr[j+0],sl1v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+4]),WdMdC4);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1d_ptr[j+4]),WdMdC4);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+4]),WdMdC4);
                store256d(&J_sl1d_ptr[j+4],sl1v4);
                store256d(&J_sr0d_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1d_ptr[j+0]),WdMdC0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl1d_ptr[j+0],sl1v0);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
            }
        }
    }

    for(std::size_t k=endRR21-remainRR21; k<endRR21; k++)
    {
        const unsigned int i0 = this->RR21AllIndex[k*4+0];
        const unsigned int sl0a = this->RR21AllIndex[k*4+1];
        const unsigned int sl1a = this->RR21AllIndex[k*4+2];
        const unsigned int sr0a = this->RR21AllIndex[k*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double C_sr0a = C[sr0a];

        const double Kf0 = this->Kf_[i0];
        const double Kf0dCfdC0 = Kf0*C_sl1a;
        const double Kf0dCfdC1 = Kf0*C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0;
        const double J_sl1a_sr0a_increment = J_sl1a_ptr[sr0a] + Kr0;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sl1a_ptr[sr0a] = J_sl1a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a]) - (dBdT[sl0a] + dBdT[sl1a]);
        const double dKcdTByKc0 = sumVdBdT0 + invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C_sl0a*C_sl1a;
        const double CR0 = C_sr0a;
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF0-invKc0*CR0);
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d WdMdC4 = _mm256_mul_pd(Wv,dMdC4);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1a_ptr[j+0]),WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl1a_ptr[j+0],sl1v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+4]),WdMdC4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1a_ptr[j+4]),WdMdC4);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl1a_ptr[j+4],sl1v4);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1a_ptr[j+0]),WdMdC0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl1a_ptr[j+0],sl1v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
            }
        }
    }



    std::size_t endRR21Dup = this->RR21DupSize;
    std::size_t remainRR21Dup = endRR21Dup % 4;
    for(std::size_t k=0; k<endRR21Dup-remainRR21Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->RR21DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR21DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR21DupAllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];

        const double C_sl0a = C[sl0a];
        const double C_sr0a = C[sr0a];


        const double Kf0 = this->Kf_[i0];
        const double Kf0dCfdC0a = 2*Kf0*C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - 2*Kf0dCfdC0a;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0 + Kr0;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a]) - (dBdT[sl0a] + dBdT[sl0a]);
        const double dKcdTByKc0 = sumVdBdT0 + invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C_sl0a*C_sl0a;
        const double CR0 = C_sr0a;
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - 2*dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - 2*q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF0-invKc0*CR0);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);

                __m256d sl0v = load256d(&J_sl0a_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);

                __m256d sr0v = load256d(&J_sr0a_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR21DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->RR21DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->RR21DupAllIndex[(k+1)*3+2];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];

        const double C_sl0b = C[sl0b];
        const double C_sr0b = C[sr0b];


        const double Kf1 = this->Kf_[i1];
        const double Kf1dCfdC0b = 2*Kf1*C_sl0b;
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - 2*Kf1dCfdC0b;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1dCfdC0b;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;

        const double invKp1 = (invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1 + Kr1;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double sumVdBdT1 = (dBdT[sr0b]) - (dBdT[sl0b] + dBdT[sl0b]);
        const double dKcdTByKc1 = sumVdBdT1 + invT;
        const double dKrdT1 = (dKfdT1*invKc1 - (invKc1 < invKcLimiter ? Kr1*dKcdTByKc1 : 0));
        const double CF1 = C_sl0b*C_sl0b;
        const double CR1 = C_sr0b;
        const double dqdT1 = (dKfdT1*CF1) - (dKrdT1*CR1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - 2*dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - 2*q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF1-invKc1*CR1);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);

                __m256d sl0v = load256d(&J_sl0b_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);

                __m256d sr0v = load256d(&J_sr0b_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0b_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR21DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->RR21DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->RR21DupAllIndex[(k+2)*3+2];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];

        const double C_sl0c = C[sl0c];
        const double C_sr0c = C[sr0c];


        const double Kf2 = this->Kf_[i2];
        const double Kf2dCfdC0c = 2*Kf2*C_sl0c;
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - 2*Kf2dCfdC0c;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2dCfdC0c;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;

        const double invKp2 = (invNegGstdByRT[sr0c])*(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2 + Kr2;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double sumVdBdT2 = (dBdT[sr0c]) - (dBdT[sl0c] + dBdT[sl0c]);
        const double dKcdTByKc2 = sumVdBdT2 + invT;
        const double dKrdT2 = (dKfdT2*invKc2 - (invKc2 < invKcLimiter ? Kr2*dKcdTByKc2 : 0));
        const double CF2 = C_sl0c*C_sl0c;
        const double CR2 = C_sr0c;
        const double dqdT2 = (dKfdT2*CF2) - (dKrdT2*CR2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - 2*dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - 2*q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF2-invKc2*CR2);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);

                __m256d sl0v = load256d(&J_sl0c_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);

                __m256d sr0v = load256d(&J_sr0c_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0c_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR21DupAllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->RR21DupAllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->RR21DupAllIndex[(k+3)*3+2];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];

        const double C_sl0d = C[sl0d];
        const double C_sr0d = C[sr0d];


        const double Kf3 = this->Kf_[i3];
        const double Kf3dCfdC0d = 2*Kf3*C_sl0d;
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - 2*Kf3dCfdC0d;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3dCfdC0d;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;

        const double invKp3 = (invNegGstdByRT[sr0d])*(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl0d]);
        double invKc3 = invKp3*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3 + Kr3;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double sumVdBdT3 = (dBdT[sr0d]) - (dBdT[sl0d] + dBdT[sl0d]);
        const double dKcdTByKc3 = sumVdBdT3 + invT;
        const double dKrdT3 = (dKfdT3*invKc3 - (invKc3 < invKcLimiter ? Kr3*dKcdTByKc3 : 0));
        const double CF3 = C_sl0d*C_sl0d;
        const double CR3 = C_sr0d;
        const double dqdT3 = (dKfdT3*CF3) - (dKrdT3*CR3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - 2*dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - 2*q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF3-invKc3*CR3);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);

                __m256d sl0v = load256d(&J_sl0d_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);

                __m256d sr0v = load256d(&J_sr0d_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0d_ptr[j],sr0v);
            }
        }
    }

    for(std::size_t k=endRR21Dup-remainRR21Dup; k<endRR21Dup; k++)
    {
        const unsigned int i0 = this->RR21DupAllIndex[k*3+0];
        const unsigned int sl0a = this->RR21DupAllIndex[k*3+1];
        const unsigned int sr0a = this->RR21DupAllIndex[k*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];

        const double C_sl0a = C[sl0a];
        const double C_sr0a = C[sr0a];


        const double Kf0 = this->Kf_[i0];
        const double Kf0dCfdC0a = 2*Kf0*C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - 2*Kf0dCfdC0a;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0 + Kr0;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a]) - (dBdT[sl0a] + dBdT[sl0a]);
        const double dKcdTByKc0 = sumVdBdT0 + invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C_sl0a*C_sl0a;
        const double CR0 = C_sr0a;
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - 2*dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - 2*q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF0-invKc0*CR0);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);

                __m256d sl0v = load256d(&J_sl0a_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);

                __m256d sr0v = load256d(&J_sr0a_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }
    }

}


void FastChemistry::OptReaction::JF21IR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
) const noexcept
{

    std::size_t end21IR = this->IR21Size;
    std::size_t remain21IR = end21IR % 4;

    for(std::size_t k=0; k<end21IR-remain21IR; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->IR21AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->IR21AllIndex[(k+0)*4+1];
        const unsigned int sl1a = this->IR21AllIndex[(k+0)*4+2];
        const unsigned int sr0a = this->IR21AllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];

        const double Kf0 = this->Kf_[i0];
        const double Kf0dCfdC0 = Kf0*C_sl1a;
        const double Kf0dCfdC1 = Kf0*C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C_sl0a*C_sl1a;
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF0;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&J_sl1a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sl1a_ptr[j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR21AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->IR21AllIndex[(k+1)*4+1];
        const unsigned int sl1b = this->IR21AllIndex[(k+1)*4+2];
        const unsigned int sr0b = this->IR21AllIndex[(k+1)*4+3];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sl1b_ptr = &ddNdtByVdcTp[sl1b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        const double C_sl1b = C[sl1b];
        const double C_sl0b = C[sl0b];

        const double Kf1 = this->Kf_[i1];
        const double Kf1dCfdC0 = Kf1*C_sl1b;
        const double Kf1dCfdC1 = Kf1*C_sl0b;
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1dCfdC0;
        const double J_sl1b_sl0b_increment = J_sl1b_ptr[sl0b] - Kf1dCfdC0;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sl1b_ptr[sl0b] = J_sl1b_sl0b_increment;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1dCfdC0;
        const double J_sl0b_sl1b_increment = J_sl0b_ptr[sl1b] - Kf1dCfdC1;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sl0b_ptr[sl1b] = J_sl0b_sl1b_increment;
        const double J_sl1b_sl1b_increment = J_sl1b_ptr[sl1b] - Kf1dCfdC1;
        const double J_sr0b_sl1b_increment = J_sr0b_ptr[sl1b] + Kf1dCfdC1;
        J_sl1b_ptr[sl1b] = J_sl1b_sl1b_increment;
        J_sr0b_ptr[sl1b] = J_sr0b_sl1b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double CF1 = C_sl0b*C_sl1b;
        const double dqdT1 = (dKfdT1*CF1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sl1b_T_increment = J_sl1b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sl1b_ptr[this->nSpecies] = J_sl1b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;

        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF1;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0b_ptr[j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&J_sl1b_ptr[j]),WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sl1b_ptr[j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0b_ptr[j]),WdMdC);
                store256d(&J_sr0b_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR21AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->IR21AllIndex[(k+2)*4+1];
        const unsigned int sl1c = this->IR21AllIndex[(k+2)*4+2];
        const unsigned int sr0c = this->IR21AllIndex[(k+2)*4+3];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sl1c_ptr = &ddNdtByVdcTp[sl1c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        const double C_sl1c = C[sl1c];
        const double C_sl0c = C[sl0c];

        const double Kf2 = this->Kf_[i2];
        const double Kf2dCfdC0 = Kf2*C_sl1c;
        const double Kf2dCfdC1 = Kf2*C_sl0c;
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2dCfdC0;
        const double J_sl1c_sl0c_increment = J_sl1c_ptr[sl0c] - Kf2dCfdC0;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sl1c_ptr[sl0c] = J_sl1c_sl0c_increment;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2dCfdC0;
        const double J_sl0c_sl1c_increment = J_sl0c_ptr[sl1c] - Kf2dCfdC1;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sl0c_ptr[sl1c] = J_sl0c_sl1c_increment;
        const double J_sl1c_sl1c_increment = J_sl1c_ptr[sl1c] - Kf2dCfdC1;
        const double J_sr0c_sl1c_increment = J_sr0c_ptr[sl1c] + Kf2dCfdC1;
        J_sl1c_ptr[sl1c] = J_sl1c_sl1c_increment;
        J_sr0c_ptr[sl1c] = J_sr0c_sl1c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double CF2 = C_sl0c*C_sl1c;
        const double dqdT2 = (dKfdT2*CF2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sl1c_T_increment = J_sl1c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sl1c_ptr[this->nSpecies] = J_sl1c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;

        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF2;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0c_ptr[j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&J_sl1c_ptr[j]),WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sl1c_ptr[j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0c_ptr[j]),WdMdC);
                store256d(&J_sr0c_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR21AllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->IR21AllIndex[(k+3)*4+1];
        const unsigned int sl1d = this->IR21AllIndex[(k+3)*4+2];
        const unsigned int sr0d = this->IR21AllIndex[(k+3)*4+3];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sl1d_ptr = &ddNdtByVdcTp[sl1d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        const double C_sl1d = C[sl1d];
        const double C_sl0d = C[sl0d];

        const double Kf3 = this->Kf_[i3];
        const double Kf3dCfdC0 = Kf3*C_sl1d;
        const double Kf3dCfdC1 = Kf3*C_sl0d;
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3dCfdC0;
        const double J_sl1d_sl0d_increment = J_sl1d_ptr[sl0d] - Kf3dCfdC0;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sl1d_ptr[sl0d] = J_sl1d_sl0d_increment;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3dCfdC0;
        const double J_sl0d_sl1d_increment = J_sl0d_ptr[sl1d] - Kf3dCfdC1;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sl0d_ptr[sl1d] = J_sl0d_sl1d_increment;
        const double J_sl1d_sl1d_increment = J_sl1d_ptr[sl1d] - Kf3dCfdC1;
        const double J_sr0d_sl1d_increment = J_sr0d_ptr[sl1d] + Kf3dCfdC1;
        J_sl1d_ptr[sl1d] = J_sl1d_sl1d_increment;
        J_sr0d_ptr[sl1d] = J_sr0d_sl1d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double CF3 = C_sl0d*C_sl1d;
        const double dqdT3 = (dKfdT3*CF3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sl1d_T_increment = J_sl1d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sl1d_ptr[this->nSpecies] = J_sl1d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;

        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF3;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0d_ptr[j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&J_sl1d_ptr[j]),WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sl1d_ptr[j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0d_ptr[j]),WdMdC);
                store256d(&J_sr0d_ptr[j],sr0v);
            }
        }
    }

    for(std::size_t k=end21IR-remain21IR; k<end21IR; k++)
    {
        const unsigned int i0 = this->IR21AllIndex[k*4+0];
        const unsigned int sl0a = this->IR21AllIndex[k*4+1];
        const unsigned int sl1a = this->IR21AllIndex[k*4+2];
        const unsigned int sr0a = this->IR21AllIndex[k*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];

        const double Kf0 = this->Kf_[i0];
        const double Kf0dCfdC0 = Kf0*C_sl1a;
        const double Kf0dCfdC1 = Kf0*C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C_sl0a*C_sl1a;
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF0;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&J_sl1a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sl1a_ptr[j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }
    }

    {
    std::size_t end21IRDup = this->IR21DupSize;
    std::size_t remain21IRDup = end21IRDup%4;

    for(std::size_t k=0; k<end21IRDup-remain21IRDup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->IR21DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR21DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR21DupAllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        const double C_sl0a = C[sl0a];


        const double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C_sl0a;

        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - 4*Kf0*dCfdC0a;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0*(dCfdC0a+dCfdC0a);
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C_sl0a*C_sl0a;
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - 2*dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - 2*q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF0;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&J_sl0a_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                __m256d sr0v = load256d(&J_sr0a_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR21DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->IR21DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->IR21DupAllIndex[(k+1)*3+2];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        const double C_sl0b = C[sl0b];


        const double Kf1 = this->Kf_[i1];
        const double dCfdC0b = C_sl0b;
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - 4*Kf1*dCfdC0b;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1*(dCfdC0b+dCfdC0b);
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double CF1 = C_sl0b*C_sl0b;
        const double dqdT1 = (dKfdT1*CF1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - 2*dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;

        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - 2*q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF1;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&J_sl0b_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                __m256d sr0v = load256d(&J_sr0b_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0b_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR21DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->IR21DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->IR21DupAllIndex[(k+2)*3+2];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        const double C_sl0c = C[sl0c];


        const double Kf2 = this->Kf_[i2];
        const double dCfdC0c = C_sl0c;
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - 4*Kf2*dCfdC0c;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2*(dCfdC0c+dCfdC0c);
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double CF2 = C_sl0c*C_sl0c;
        const double dqdT2 = (dKfdT2*CF2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - 2*dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;

        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - 2*q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF2;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&J_sl0c_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                __m256d sr0v = load256d(&J_sr0c_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0c_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR21DupAllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->IR21DupAllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->IR21DupAllIndex[(k+3)*3+2];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        const double C_sl0d = C[sl0d];


        const double Kf3 = this->Kf_[i3];
        const double dCfdC0d = C_sl0d;
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - 4*Kf3*dCfdC0d;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3*(dCfdC0d+dCfdC0d);
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double CF3 = C_sl0d*C_sl0d;
        const double dqdT3 = (dKfdT3*CF3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - 2*dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;

        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - 2*q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF3;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&J_sl0d_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                __m256d sr0v = load256d(&J_sr0d_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0d_ptr[j],sr0v);
            }
        }
    }

    for(std::size_t k=end21IRDup-remain21IRDup; k<end21IRDup; k++)
    {
        const unsigned int i0 = this->IR21DupAllIndex[k*3+0];
        const unsigned int sl0a = this->IR21DupAllIndex[k*3+1];
        const unsigned int sr0a = this->IR21DupAllIndex[k*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        const double C_sl0a = C[sl0a];


        const double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - 4*Kf0*dCfdC0a;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0*(dCfdC0a+dCfdC0a);
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C_sl0a*C_sl0a;
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - 2*dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - 2*q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF0;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&J_sl0a_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                __m256d sr0v = load256d(&J_sr0a_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }
    }
    }
}


void  FastChemistry::OptReaction::JF21NER
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t endNER21 = this->NER21AllIndex.size()/4;
    std::size_t remainNER21 = endNER21 % 4;

    for(std::size_t k=0; k<endNER21-remainNER21; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->NER21AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->NER21AllIndex[(k+0)*4+1];
        const unsigned int sl1a = this->NER21AllIndex[(k+0)*4+2];
        const unsigned int sr0a = this->NER21AllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double C_sr0a = C[sr0a];

        const double Kf0 = this->Kf_[i0];
        const double Kf0dCfdC0 = Kf0*C_sl1a;
        const double Kf0dCfdC1 = Kf0*C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;

        const double Kr0 = this->Kf_[i0 - this->Ikf[1] + this->Ikf[9]];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0;
        const double J_sl1a_sr0a_increment = J_sl1a_ptr[sr0a] + Kr0;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sl1a_ptr[sr0a] = J_sl1a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - this->Ikf[1] + this->Ikf[9]];
        const double CF0 = C_sl0a*C_sl1a;
        const double CR0 = C_sr0a;
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF0 - this->dKfdC_[idx+this->Itbr[4]]*CR0;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&J_sl1a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sl1a_ptr[j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->NER21AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->NER21AllIndex[(k+1)*4+1];
        const unsigned int sl1b = this->NER21AllIndex[(k+1)*4+2];
        const unsigned int sr0b = this->NER21AllIndex[(k+1)*4+3];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sl1b_ptr = &ddNdtByVdcTp[sl1b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        const double C_sl1b = C[sl1b];
        const double C_sl0b = C[sl0b];
        const double C_sr0b = C[sr0b];

        const double Kf1 = this->Kf_[i1];
        const double Kf1dCfdC0 = Kf1*C_sl1b;
        const double Kf1dCfdC1 = Kf1*C_sl0b;
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1dCfdC0;
        const double J_sl1b_sl0b_increment = J_sl1b_ptr[sl0b] - Kf1dCfdC0;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sl1b_ptr[sl0b] = J_sl1b_sl0b_increment;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1dCfdC0;
        const double J_sl0b_sl1b_increment = J_sl0b_ptr[sl1b] - Kf1dCfdC1;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sl0b_ptr[sl1b] = J_sl0b_sl1b_increment;
        const double J_sl1b_sl1b_increment = J_sl1b_ptr[sl1b] - Kf1dCfdC1;
        const double J_sr0b_sl1b_increment = J_sr0b_ptr[sl1b] + Kf1dCfdC1;
        J_sl1b_ptr[sl1b] = J_sl1b_sl1b_increment;
        J_sr0b_ptr[sl1b] = J_sr0b_sl1b_increment;

        const double Kr1 = this->Kf_[i1 - this->Ikf[1] + this->Ikf[9]];
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1;
        const double J_sl1b_sr0b_increment = J_sl1b_ptr[sr0b] + Kr1;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sl1b_ptr[sr0b] = J_sl1b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double dKrdT1 = this->dKfdT_[i1 - this->Ikf[1] + this->Ikf[9]];
        const double CF1 = C_sl0b*C_sl1b;
        const double CR1 = C_sr0b;
        const double dqdT1 = (dKfdT1*CF1) - (dKrdT1*CR1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sl1b_T_increment = J_sl1b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sl1b_ptr[this->nSpecies] = J_sl1b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF1 - this->dKfdC_[idx+this->Itbr[4]]*CR1;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0b_ptr[j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&J_sl1b_ptr[j]),WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sl1b_ptr[j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0b_ptr[j]),WdMdC);
                store256d(&J_sr0b_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->NER21AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->NER21AllIndex[(k+2)*4+1];
        const unsigned int sl1c = this->NER21AllIndex[(k+2)*4+2];
        const unsigned int sr0c = this->NER21AllIndex[(k+2)*4+3];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sl1c_ptr = &ddNdtByVdcTp[sl1c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        const double C_sl1c = C[sl1c];
        const double C_sl0c = C[sl0c];
        const double C_sr0c = C[sr0c];

        const double Kf2 = this->Kf_[i2];
        const double Kf2dCfdC0 = Kf2*C_sl1c;
        const double Kf2dCfdC1 = Kf2*C_sl0c;
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2dCfdC0;
        const double J_sl1c_sl0c_increment = J_sl1c_ptr[sl0c] - Kf2dCfdC0;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sl1c_ptr[sl0c] = J_sl1c_sl0c_increment;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2dCfdC0;
        const double J_sl0c_sl1c_increment = J_sl0c_ptr[sl1c] - Kf2dCfdC1;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sl0c_ptr[sl1c] = J_sl0c_sl1c_increment;
        const double J_sl1c_sl1c_increment = J_sl1c_ptr[sl1c] - Kf2dCfdC1;
        const double J_sr0c_sl1c_increment = J_sr0c_ptr[sl1c] + Kf2dCfdC1;
        J_sl1c_ptr[sl1c] = J_sl1c_sl1c_increment;
        J_sr0c_ptr[sl1c] = J_sr0c_sl1c_increment;

        const double Kr2 = this->Kf_[i2 - this->Ikf[1] + this->Ikf[9]];
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2;
        const double J_sl1c_sr0c_increment = J_sl1c_ptr[sr0c] + Kr2;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sl1c_ptr[sr0c] = J_sl1c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double dKrdT2 = this->dKfdT_[i2 - this->Ikf[1] + this->Ikf[9]];
        const double CF2 = C_sl0c*C_sl1c;
        const double CR2 = C_sr0c;
        const double dqdT2 = (dKfdT2*CF2) - (dKrdT2*CR2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sl1c_T_increment = J_sl1c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sl1c_ptr[this->nSpecies] = J_sl1c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF2 - this->dKfdC_[idx+this->Itbr[4]]*CR2;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0c_ptr[j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&J_sl1c_ptr[j]),WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sl1c_ptr[j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0c_ptr[j]),WdMdC);
                store256d(&J_sr0c_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->NER21AllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->NER21AllIndex[(k+3)*4+1];
        const unsigned int sl1d = this->NER21AllIndex[(k+3)*4+2];
        const unsigned int sr0d = this->NER21AllIndex[(k+3)*4+3];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sl1d_ptr = &ddNdtByVdcTp[sl1d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        const double C_sl1d = C[sl1d];
        const double C_sl0d = C[sl0d];
        const double C_sr0d = C[sr0d];

        const double Kf3 = this->Kf_[i3];
        const double Kf3dCfdC0 = Kf3*C_sl1d;
        const double Kf3dCfdC1 = Kf3*C_sl0d;
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3dCfdC0;
        const double J_sl1d_sl0d_increment = J_sl1d_ptr[sl0d] - Kf3dCfdC0;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sl1d_ptr[sl0d] = J_sl1d_sl0d_increment;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3dCfdC0;
        const double J_sl0d_sl1d_increment = J_sl0d_ptr[sl1d] - Kf3dCfdC1;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sl0d_ptr[sl1d] = J_sl0d_sl1d_increment;
        const double J_sl1d_sl1d_increment = J_sl1d_ptr[sl1d] - Kf3dCfdC1;
        const double J_sr0d_sl1d_increment = J_sr0d_ptr[sl1d] + Kf3dCfdC1;
        J_sl1d_ptr[sl1d] = J_sl1d_sl1d_increment;
        J_sr0d_ptr[sl1d] = J_sr0d_sl1d_increment;

        const double Kr3 = this->Kf_[i3 - this->Ikf[1] + this->Ikf[9]];
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3;
        const double J_sl1d_sr0d_increment = J_sl1d_ptr[sr0d] + Kr3;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sl1d_ptr[sr0d] = J_sl1d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double dKrdT3 = this->dKfdT_[i3 - this->Ikf[1] + this->Ikf[9]];
        const double CF3 = C_sl0d*C_sl1d;
        const double CR3 = C_sr0d;
        const double dqdT3 = (dKfdT3*CF3) - (dKrdT3*CR3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sl1d_T_increment = J_sl1d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sl1d_ptr[this->nSpecies] = J_sl1d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF3 - this->dKfdC_[idx+this->Itbr[4]]*CR3;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0d_ptr[j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&J_sl1d_ptr[j]),WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sl1d_ptr[j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0d_ptr[j]),WdMdC);
                store256d(&J_sr0d_ptr[j],sr0v);
            }
        }

    }

    for(std::size_t k=endNER21-remainNER21; k<endNER21; k++)
    {
        const unsigned int i0 = this->NER21AllIndex[k*4+0];
        const unsigned int sl0a = this->NER21AllIndex[k*4+1];
        const unsigned int sl1a = this->NER21AllIndex[k*4+2];
        const unsigned int sr0a = this->NER21AllIndex[k*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double C_sr0a = C[sr0a];

        const double Kf0 = this->Kf_[i0];
        const double Kf0dCfdC0 = Kf0*C_sl1a;
        const double Kf0dCfdC1 = Kf0*C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;

        const double Kr0 = this->Kf_[i0 - this->Ikf[1] + this->Ikf[9]];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0;
        const double J_sl1a_sr0a_increment = J_sl1a_ptr[sr0a] + Kr0;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sl1a_ptr[sr0a] = J_sl1a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - this->Ikf[1] + this->Ikf[9]];
        const double CF0 = C_sl0a*C_sl1a;
        const double CR0 = C_sr0a;
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF0 - this->dKfdC_[idx+this->Itbr[4]]*CR0;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&J_sl1a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sl1a_ptr[j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }
    }

    std::size_t endNER21Dup = this->NER21DupAllIndex.size()/3;
    std::size_t remainNER21Dup = endNER21Dup % 4;

    for(std::size_t k=0; k<endNER21Dup-remainNER21Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->NER21DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER21DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER21DupAllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        const double C_sl0a = C[sl0a];
        const double C_sr0a = C[sr0a];

        const double Kf0 = this->Kf_[i0];
        const double Kf0dCfdC0a = 2*Kf0*C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - 2*Kf0dCfdC0a;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double Kr0 = this->Kf_[i0 - this->Ikf[1] + this->Ikf[9]];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0 + Kr0;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - this->Ikf[1] + this->Ikf[9]];
        const double CF0 = C_sl0a*C_sl0a;
        const double CR0 = C_sr0a;
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - 2*dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - 2*q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF0 - this->dKfdC_[idx+this->Itbr[4]]*CR0;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&J_sl0a_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                __m256d sr0v = load256d(&J_sr0a_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->NER21DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->NER21DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->NER21DupAllIndex[(k+1)*3+2];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        const double C_sl0b = C[sl0b];
        const double C_sr0b = C[sr0b];

        const double Kf1 = this->Kf_[i1];
        const double Kf1dCfdC0b = 2*Kf1*C_sl0b;
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - 2*Kf1dCfdC0b;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1dCfdC0b;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;

        const double Kr1 = this->Kf_[i1 - this->Ikf[1] + this->Ikf[9]];
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1 + Kr1;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double dKrdT1 = this->dKfdT_[i1 - this->Ikf[1] + this->Ikf[9]];
        const double CF1 = C_sl0b*C_sl0b;
        const double CR1 = C_sr0b;
        const double dqdT1 = (dKfdT1*CF1) - (dKrdT1*CR1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - 2*dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - 2*q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF1 - this->dKfdC_[idx+this->Itbr[4]]*CR1;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&J_sl0b_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                __m256d sr0v = load256d(&J_sr0b_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0b_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->NER21DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->NER21DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->NER21DupAllIndex[(k+2)*3+2];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        const double C_sl0c = C[sl0c];
        const double C_sr0c = C[sr0c];

        const double Kf2 = this->Kf_[i2];
        const double Kf2dCfdC0c = 2*Kf2*C_sl0c;
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - 2*Kf2dCfdC0c;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2dCfdC0c;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;

        const double Kr2 = this->Kf_[i2 - this->Ikf[1] + this->Ikf[9]];
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2 + Kr2;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double dKrdT2 = this->dKfdT_[i2 - this->Ikf[1] + this->Ikf[9]];
        const double CF2 = C_sl0c*C_sl0c;
        const double CR2 = C_sr0c;
        const double dqdT2 = (dKfdT2*CF2) - (dKrdT2*CR2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - 2*dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - 2*q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF2 - this->dKfdC_[idx+this->Itbr[4]]*CR2;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&J_sl0c_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                __m256d sr0v = load256d(&J_sr0c_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0c_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->NER21DupAllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->NER21DupAllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->NER21DupAllIndex[(k+3)*3+2];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        const double C_sl0d = C[sl0d];
        const double C_sr0d = C[sr0d];

        const double Kf3 = this->Kf_[i3];
        const double Kf3dCfdC0d = 2*Kf3*C_sl0d;
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - 2*Kf3dCfdC0d;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3dCfdC0d;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;

        const double Kr3 = this->Kf_[i3 - this->Ikf[1] + this->Ikf[9]];
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3 + Kr3;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double dKrdT3 = this->dKfdT_[i3 - this->Ikf[1] + this->Ikf[9]];
        const double CF3 = C_sl0d*C_sl0d;
        const double CR3 = C_sr0d;
        const double dqdT3 = (dKfdT3*CF3) - (dKrdT3*CR3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - 2*dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - 2*q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF3 - this->dKfdC_[idx+this->Itbr[4]]*CR3;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&J_sl0d_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                __m256d sr0v = load256d(&J_sr0d_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0d_ptr[j],sr0v);
            }
        }

    }

    for(std::size_t k=endNER21Dup-remainNER21Dup; k<endNER21Dup; k++)
    {
        const unsigned int i0 = this->NER21DupAllIndex[k*3+0];
        const unsigned int sl0a = this->NER21DupAllIndex[k*3+1];
        const unsigned int sr0a = this->NER21DupAllIndex[k*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        const double C_sl0a = C[sl0a];
        const double C_sr0a = C[sr0a];

        const double Kf0 = this->Kf_[i0];
        const double Kf0dCfdC0a = 2*Kf0*C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - 2*Kf0dCfdC0a;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double Kr0 = this->Kf_[i0 - this->Ikf[1] + this->Ikf[9]];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0 + Kr0;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - this->Ikf[1] + this->Ikf[9]];
        const double CF0 = C_sl0a*C_sl0a;
        const double CR0 = C_sr0a;
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - 2*dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - 2*q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF0 - this->dKfdC_[idx+this->Itbr[4]]*CR0;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&J_sl0a_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                __m256d sr0v = load256d(&J_sr0a_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }
    }

}