/*---------------------------------------------------------------------------*\
  Description
      Update the net rates of production (dNdtByV) and the molar-concentration
      based Jacobian matrix (ddNdtByVdcTp) for two-reactant / two-product
      (2-2) reactions, e.g. A + A -> B + B and A + B -> C + D.

      The Jacobian matrix is stored row-wise with one column per species
      plus one extra column for the temperature derivative (column index
      nSpecies).

      RR : reversible reaction; the reverse rate constant is obtained from
           the equilibrium constant
      IR : irreversible reaction; the reverse rate constant is zero
      NER: non-equilibrium reaction; the reverse rate constant is computed
           with the Arrhenius form instead of the equilibrium constant

  Functions
      JF22RR : reversible 2-2 reactions, incl. duplicate reactants/products
      JF22IR : irreversible 2-2 reactions, incl. duplicate reactants/products
      JF22NER: non-equilibrium 2-2 reactions

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"
#include <algorithm>

//=============================================================================//

void FastChemistry::OptReaction::JF22RR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
) const noexcept
{

    std::size_t endRR22 = this->RR22Size;
    std::size_t remainRR22 = endRR22 % 4;

    const double invKcLimiter = FastChemistry::invKcLimiter;
    for(std::size_t k=0; k<endRR22-remainRR22; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->RR22AllIndex[k*5+0];
        const unsigned int sl0a = this->RR22AllIndex[k*5+1];
        const unsigned int sl1a = this->RR22AllIndex[k*5+2];
        const unsigned int sr0a = this->RR22AllIndex[k*5+3];
        const unsigned int sr1a = this->RR22AllIndex[k*5+4];

        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double Kf0dCfdC0a = Kf0*C_sl1a;
        const double Kf0dCfdC1a = Kf0*C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0a;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1a;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0dCfdC0a;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1a;
        const double J_sr1a_sl1a_increment = J_sr1a_ptr[sl1a] + Kf0dCfdC1a;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;
        J_sr1a_ptr[sl1a] = J_sr1a_sl1a_increment;

        double invKc0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double C_sr1a = C[sr1a];
        const double C_sr0a = C[sr0a];
        const double Kr0dCrdC0a = -Kr0*C_sr1a;
        const double Kr0dCrdC1a = -Kr0*C_sr0a;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sl1a_sr0a_increment = J_sl1a_ptr[sr0a] - Kr0dCrdC0a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sl1a_ptr[sr0a] = J_sl1a_sr0a_increment;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] - Kr0dCrdC1a;
        const double J_sl1a_sr1a_increment = J_sl1a_ptr[sr1a] - Kr0dCrdC1a;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sl1a_ptr[sr1a] = J_sl1a_sr1a_increment;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] + Kr0dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] + Kr0dCrdC0a;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] + Kr0dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] + Kr0dCrdC1a;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C_sl0a*C_sl1a;
        const double CR0 = C_sr0a*C_sr1a;
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a]) - (dBdT[sl0a] + dBdT[sl1a]);
        const double dKcdTByKc0 = sumVdBdT0;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF0-dKcdTByKc0*CR0);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sl1a_ptr[j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]),WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR22AllIndex[k*5+5];
        const unsigned int sl0b = this->RR22AllIndex[k*5+6];
        const unsigned int sl1b = this->RR22AllIndex[k*5+7];
        const unsigned int sr0b = this->RR22AllIndex[k*5+8];
        const unsigned int sr1b = this->RR22AllIndex[k*5+9];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sl1b_ptr = &ddNdtByVdcTp[sl1b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double C_sl1b = C[sl1b];
        const double C_sl0b = C[sl0b];
        const double Kf1dCfdC0b = Kf1*C_sl1b;
        const double Kf1dCfdC1b = Kf1*C_sl0b;
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1dCfdC0b;
        const double J_sl1b_sl0b_increment = J_sl1b_ptr[sl0b] - Kf1dCfdC0b;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sl1b_ptr[sl0b] = J_sl1b_sl0b_increment;
        const double J_sl0b_sl1b_increment = J_sl0b_ptr[sl1b] - Kf1dCfdC1b;
        const double J_sl1b_sl1b_increment = J_sl1b_ptr[sl1b] - Kf1dCfdC1b;
        J_sl0b_ptr[sl1b] = J_sl0b_sl1b_increment;
        J_sl1b_ptr[sl1b] = J_sl1b_sl1b_increment;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1dCfdC0b;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1dCfdC0b;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;
        const double J_sr0b_sl1b_increment = J_sr0b_ptr[sl1b] + Kf1dCfdC1b;
        const double J_sr1b_sl1b_increment = J_sr1b_ptr[sl1b] + Kf1dCfdC1b;
        J_sr0b_ptr[sl1b] = J_sr0b_sl1b_increment;
        J_sr1b_ptr[sl1b] = J_sr1b_sl1b_increment;

        double invKc1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]);
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double C_sr1b = C[sr1b];
        const double C_sr0b = C[sr0b];
        const double Kr1dCrdC0b = -Kr1*C_sr1b;
        const double Kr1dCrdC1b = -Kr1*C_sr0b;
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] - Kr1dCrdC0b;
        const double J_sl1b_sr0b_increment = J_sl1b_ptr[sr0b] - Kr1dCrdC0b;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sl1b_ptr[sr0b] = J_sl1b_sr0b_increment;
        const double J_sl0b_sr1b_increment = J_sl0b_ptr[sr1b] - Kr1dCrdC1b;
        const double J_sl1b_sr1b_increment = J_sl1b_ptr[sr1b] - Kr1dCrdC1b;
        J_sl0b_ptr[sr1b] = J_sl0b_sr1b_increment;
        J_sl1b_ptr[sr1b] = J_sl1b_sr1b_increment;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] + Kr1dCrdC0b;
        const double J_sr1b_sr0b_increment = J_sr1b_ptr[sr0b] + Kr1dCrdC0b;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;
        J_sr1b_ptr[sr0b] = J_sr1b_sr0b_increment;
        const double J_sr0b_sr1b_increment = J_sr0b_ptr[sr1b] + Kr1dCrdC1b;
        const double J_sr1b_sr1b_increment = J_sr1b_ptr[sr1b] + Kr1dCrdC1b;
        J_sr0b_ptr[sr1b] = J_sr0b_sr1b_increment;
        J_sr1b_ptr[sr1b] = J_sr1b_sr1b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double sumVdBdT1 = (dBdT[sr0b] + dBdT[sr1b]) - (dBdT[sl0b] + dBdT[sl1b]);
        const double dKcdTByKc1 = sumVdBdT1;
        const double dKrdT1 = (dKfdT1*invKc1 - (invKc1 < invKcLimiter ? Kr1*dKcdTByKc1 : 0));
        const double CF1 = C_sl0b*C_sl1b;
        const double CR1 = C_sr0b*C_sr1b;
        const double dqdT1 = (dKfdT1*CF1) - (dKrdT1*CR1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sl1b_T_increment = J_sl1b_ptr[this->nSpecies] - dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sl1b_ptr[this->nSpecies] = J_sl1b_T_increment;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF1-invKc1*CR1);
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
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1b_ptr[j]),WdMdC);
                store256d(&J_sr0b_ptr[j],sr0v);
                store256d(&J_sr1b_ptr[j],sr1v);
            }
        }
        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR22AllIndex[k*5+10];
        const unsigned int sl0c = this->RR22AllIndex[k*5+11];
        const unsigned int sl1c = this->RR22AllIndex[k*5+12];
        const unsigned int sr0c = this->RR22AllIndex[k*5+13];
        const unsigned int sr1c = this->RR22AllIndex[k*5+14];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sl1c_ptr = &ddNdtByVdcTp[sl1c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        double C_sl1c = C[sl1c];
        double C_sl0c = C[sl0c];
        const double Kf2dCfdC0c = Kf2*C_sl1c;
        const double Kf2dCfdC1c = Kf2*C_sl0c;
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2dCfdC0c;
        const double J_sl1c_sl0c_increment = J_sl1c_ptr[sl0c] - Kf2dCfdC0c;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sl1c_ptr[sl0c] = J_sl1c_sl0c_increment;
        const double J_sl0c_sl1c_increment = J_sl0c_ptr[sl1c] - Kf2dCfdC1c;
        const double J_sl1c_sl1c_increment = J_sl1c_ptr[sl1c] - Kf2dCfdC1c;
        J_sl0c_ptr[sl1c] = J_sl0c_sl1c_increment;
        J_sl1c_ptr[sl1c] = J_sl1c_sl1c_increment;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2dCfdC0c;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2dCfdC0c;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;
        const double J_sr0c_sl1c_increment = J_sr0c_ptr[sl1c] + Kf2dCfdC1c;
        const double J_sr1c_sl1c_increment = J_sr1c_ptr[sl1c] + Kf2dCfdC1c;
        J_sr0c_ptr[sl1c] = J_sr0c_sl1c_increment;
        J_sr1c_ptr[sl1c] = J_sr1c_sl1c_increment;

        const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]);
        const double invKc2 = std::min(invKp2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double C_sr0c = C[sr0c];
        const double C_sr1c = C[sr1c];
        const double Kr2dCrdC0c = -Kr2*C_sr1c;
        const double Kr2dCrdC1c = -Kr2*C_sr0c;
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] - Kr2dCrdC0c;
        const double J_sl1c_sr0c_increment = J_sl1c_ptr[sr0c] - Kr2dCrdC0c;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sl1c_ptr[sr0c] = J_sl1c_sr0c_increment;
        const double J_sl0c_sr1c_increment = J_sl0c_ptr[sr1c] - Kr2dCrdC1c;
        const double J_sl1c_sr1c_increment = J_sl1c_ptr[sr1c] - Kr2dCrdC1c;
        J_sl0c_ptr[sr1c] = J_sl0c_sr1c_increment;
        J_sl1c_ptr[sr1c] = J_sl1c_sr1c_increment;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] + Kr2dCrdC0c;
        const double J_sr1c_sr0c_increment = J_sr1c_ptr[sr0c] + Kr2dCrdC0c;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;
        J_sr1c_ptr[sr0c] = J_sr1c_sr0c_increment;
        const double J_sr0c_sr1c_increment = J_sr0c_ptr[sr1c] + Kr2dCrdC1c;
        const double J_sr1c_sr1c_increment = J_sr1c_ptr[sr1c] + Kr2dCrdC1c;
        J_sr0c_ptr[sr1c] = J_sr0c_sr1c_increment;
        J_sr1c_ptr[sr1c] = J_sr1c_sr1c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double sumVdBdT2 = (dBdT[sr0c] + dBdT[sr1c]) - (dBdT[sl0c] + dBdT[sl1c]);
        const double dKcdTByKc2 = sumVdBdT2;
        const double dKrdT2 = (dKfdT2*invKc2 - (invKc2 < invKcLimiter ? Kr2*dKcdTByKc2 : 0));
        const double CF2 = C_sl0c*C_sl1c;
        const double CR2 = C_sr0c*C_sr1c;
        const double dqdT2 = (dKfdT2*CF2) - (dKrdT2*CR2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sl1c_T_increment = J_sl1c_ptr[this->nSpecies] - dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sl1c_ptr[this->nSpecies] = J_sl1c_T_increment;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF2-invKc2*CR2);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[j]),WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sl1c_ptr[j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[j]),WdMdC);
                store256d(&J_sr0c_ptr[j],sr0v);
                store256d(&J_sr1c_ptr[j],sr1v);
            }
        }
        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR22AllIndex[k*5+15];
        const unsigned int sl0d = this->RR22AllIndex[k*5+16];
        const unsigned int sl1d = this->RR22AllIndex[k*5+17];
        const unsigned int sr0d = this->RR22AllIndex[k*5+18];
        const unsigned int sr1d = this->RR22AllIndex[k*5+19];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sl1d_ptr = &ddNdtByVdcTp[sl1d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];
        const double Kf3 = this->Kf_[i3];
        const double C_sl1d = C[sl1d];
        const double C_sl0d = C[sl0d];
        const double Kf3dCfdC0d = Kf3*C_sl1d;
        const double Kf3dCfdC1d = Kf3*C_sl0d;
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3dCfdC0d;
        const double J_sl1d_sl0d_increment = J_sl1d_ptr[sl0d] - Kf3dCfdC0d;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sl1d_ptr[sl0d] = J_sl1d_sl0d_increment;
        const double J_sl0d_sl1d_increment = J_sl0d_ptr[sl1d] - Kf3dCfdC1d;
        const double J_sl1d_sl1d_increment = J_sl1d_ptr[sl1d] - Kf3dCfdC1d;
        J_sl0d_ptr[sl1d] = J_sl0d_sl1d_increment;
        J_sl1d_ptr[sl1d] = J_sl1d_sl1d_increment;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3dCfdC0d;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3dCfdC0d;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;
        const double J_sr0d_sl1d_increment = J_sr0d_ptr[sl1d] + Kf3dCfdC1d;
        const double J_sr1d_sl1d_increment = J_sr1d_ptr[sl1d] + Kf3dCfdC1d;
        J_sr0d_ptr[sl1d] = J_sr0d_sl1d_increment;
        J_sr1d_ptr[sl1d] = J_sr1d_sl1d_increment;

        const double invKp3 = (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d])*(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]);
        const double invKc3 = std::min(invKp3, invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double C_sr1d = C[sr1d];
        const double C_sr0d = C[sr0d];
        const double Kr3dCrdC0d = -Kr3*C_sr1d;
        const double Kr3dCrdC1d = -Kr3*C_sr0d;
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] - Kr3dCrdC0d;
        const double J_sl1d_sr0d_increment = J_sl1d_ptr[sr0d] - Kr3dCrdC0d;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sl1d_ptr[sr0d] = J_sl1d_sr0d_increment;
        const double J_sl0d_sr1d_increment = J_sl0d_ptr[sr1d] - Kr3dCrdC1d;
        const double J_sl1d_sr1d_increment = J_sl1d_ptr[sr1d] - Kr3dCrdC1d;
        J_sl0d_ptr[sr1d] = J_sl0d_sr1d_increment;
        J_sl1d_ptr[sr1d] = J_sl1d_sr1d_increment;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] + Kr3dCrdC0d;
        const double J_sr1d_sr0d_increment = J_sr1d_ptr[sr0d] + Kr3dCrdC0d;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;
        J_sr1d_ptr[sr0d] = J_sr1d_sr0d_increment;
        const double J_sr0d_sr1d_increment = J_sr0d_ptr[sr1d] + Kr3dCrdC1d;
        const double J_sr1d_sr1d_increment = J_sr1d_ptr[sr1d] + Kr3dCrdC1d;
        J_sr0d_ptr[sr1d] = J_sr0d_sr1d_increment;
        J_sr1d_ptr[sr1d] = J_sr1d_sr1d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double sumVdBdT3 = (dBdT[sr0d] + dBdT[sr1d]) - (dBdT[sl0d] + dBdT[sl1d]);
        const double dKcdTByKc3 = sumVdBdT3;
        const double dKrdT3 = (dKfdT3*invKc3 - (invKc3 < invKcLimiter ? Kr3*dKcdTByKc3 : 0));
        const double CF3 = C[sl0d]*C[sl1d];
        const double CR3 = C[sr0d]*C[sr1d];
        const double dqdT3 = (dKfdT3*CF3) - (dKrdT3*CR3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sl1d_T_increment = J_sl1d_ptr[this->nSpecies] - dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sl1d_ptr[this->nSpecies] = J_sl1d_T_increment;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF3-invKc3*CR3);
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
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1d_ptr[j]),WdMdC);
                store256d(&J_sr0d_ptr[j],sr0v);
                store256d(&J_sr1d_ptr[j],sr1v);
            }
        }
    }
    for(std::size_t k=endRR22-remainRR22; k<endRR22; k=k+1)
    {
        const unsigned int i0 = this->RR22AllIndex[k*5+0];
        const unsigned int sl0a = this->RR22AllIndex[k*5+1];
        const unsigned int sl1a = this->RR22AllIndex[k*5+2];
        const unsigned int sr0a = this->RR22AllIndex[k*5+3];
        const unsigned int sr1a = this->RR22AllIndex[k*5+4];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double Kf0dCfdC0a = Kf0*C_sl1a;
        const double Kf0dCfdC1a = Kf0*C_sl0a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0a;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1a;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0dCfdC0a;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1a;
        const double J_sr1a_sl1a_increment = J_sr1a_ptr[sl1a] + Kf0dCfdC1a;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;
        J_sr1a_ptr[sl1a] = J_sr1a_sl1a_increment;

        double invKc0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double C_sr1a = C[sr1a];
        const double C_sr0a = C[sr0a];
        const double Kr0dCrdC0a = -Kr0*C_sr1a;
        const double Kr0dCrdC1a = -Kr0*C_sr0a;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sl1a_sr0a_increment = J_sl1a_ptr[sr0a] - Kr0dCrdC0a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sl1a_ptr[sr0a] = J_sl1a_sr0a_increment;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] - Kr0dCrdC1a;
        const double J_sl1a_sr1a_increment = J_sl1a_ptr[sr1a] - Kr0dCrdC1a;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sl1a_ptr[sr1a] = J_sl1a_sr1a_increment;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] + Kr0dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] + Kr0dCrdC0a;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] + Kr0dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] + Kr0dCrdC1a;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C_sl0a*C_sl1a;
        const double CR0 = C_sr0a*C_sr1a;
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a]) - (dBdT[sl0a] + dBdT[sl1a]);
        const double dKcdTByKc0 = sumVdBdT0;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF0-dKcdTByKc0*CR0);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sl1a_ptr[j],sl1v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }
    }
    {
    std::size_t endRRDup22 = this->RR22DupSize;
    std::size_t remainRRDup22 = endRRDup22 % 4;

    for(std::size_t k=0; k<endRRDup22-remainRRDup22; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->RR22DupAllIndex[k*5+0];
        const unsigned int sl0a = this->RR22DupAllIndex[k*5+1];
        const unsigned int sl1a = this->RR22DupAllIndex[k*5+2];
        const unsigned int sr0a = this->RR22DupAllIndex[k*5+3];
        const unsigned int sr1a = this->RR22DupAllIndex[k*5+4];

        const double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C[sl1a];
        const double dCfdC1a = C[sl0a];
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl1a] += Kf0*dCfdC1a;

        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
        double invKc0 = std::min(invKp0,FastChemistry::invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a];
        ddNdtByVdcTp[sl0a*(this->alignN)+ sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl1a*(this->alignN)+ sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl0a*(this->alignN)+ sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sl1a*(this->alignN)+ sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sr0a*(this->alignN)+ sr0a] += (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sr1a*(this->alignN)+ sr0a] += (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sr0a*(this->alignN)+ sr1a] += (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sr1a*(this->alignN)+ sr1a] += (-Kr0*dCrdC1a);

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a]) - (dBdT[sl0a] + dBdT[sl1a]);
        const double dKcdTByKc0 = sumVdBdT0;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < FastChemistry::invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C[sl0a]*C[sl1a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        ddNdtByVdcTp[sl0a*(this->alignN)+(this->nSpecies)] -= dqdT0;
        ddNdtByVdcTp[sl1a*(this->alignN)+(this->nSpecies)] -= dqdT0;
        ddNdtByVdcTp[sr0a*(this->alignN)+(this->nSpecies)] += dqdT0;
        ddNdtByVdcTp[sr1a*(this->alignN)+(this->nSpecies)] += dqdT0;


        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF0-invKc0*CR0);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j],sr1v);
            }
        }
        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1   = this->RR22DupAllIndex[k*5+5];
        const unsigned int sl0b = this->RR22DupAllIndex[k*5+6];
        const unsigned int sl1b = this->RR22DupAllIndex[k*5+7];
        const unsigned int sr0b = this->RR22DupAllIndex[k*5+8];
        const unsigned int sr1b = this->RR22DupAllIndex[k*5+9];
        const double Kf1 = this->Kf_[i1];
        const double dCfdC0b = C[sl1b];
        const double dCfdC1b = C[sl0b];
        ddNdtByVdcTp[sl0b*(this->alignN)+ sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl1b*(this->alignN)+ sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl0b*(this->alignN)+ sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sl1b*(this->alignN)+ sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sr0b*(this->alignN)+ sl0b] += Kf1*dCfdC0b;
        ddNdtByVdcTp[sr1b*(this->alignN)+ sl0b] += Kf1*dCfdC0b;
        ddNdtByVdcTp[sr0b*(this->alignN)+ sl1b] += Kf1*dCfdC1b;
        ddNdtByVdcTp[sr1b*(this->alignN)+ sl1b] += Kf1*dCfdC1b;

        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]);
        const double invKc1 = std::min(invKp1,FastChemistry::invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double dCrdC0b = C[sr1b];
        const double dCrdC1b = C[sr0b];
        ddNdtByVdcTp[sl0b*(this->alignN)+ sr0b] -= (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sl1b*(this->alignN)+ sr0b] -= (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sl0b*(this->alignN)+ sr1b] -= (-Kr1*dCrdC1b);
        ddNdtByVdcTp[sl1b*(this->alignN)+ sr1b] -= (-Kr1*dCrdC1b);
        ddNdtByVdcTp[sr0b*(this->alignN)+ sr0b] += (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sr1b*(this->alignN)+ sr0b] += (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sr0b*(this->alignN)+ sr1b] += (-Kr1*dCrdC1b);
        ddNdtByVdcTp[sr1b*(this->alignN)+ sr1b] += (-Kr1*dCrdC1b);

        const double dKfdT1 = this->dKfdT_[i1];
        const double sumVdBdT1 = (dBdT[sr0b] + dBdT[sr1b]) - (dBdT[sl0b] + dBdT[sl1b]);
        const double dKcdTByKc1 = sumVdBdT1;
        const double dKrdT1 = (dKfdT1*invKc1 - (invKc1 < FastChemistry::invKcLimiter ? Kr1*dKcdTByKc1 : 0));
        const double CF1 = C[sl0b]*C[sl1b];
        const double CR1 = C[sr0b]*C[sr1b];
        const double dqdT1 = (dKfdT1*CF1) - (dKrdT1*CR1);
        ddNdtByVdcTp[sl0b*(this->alignN)+(this->nSpecies)] -= dqdT1;
        ddNdtByVdcTp[sl1b*(this->alignN)+(this->nSpecies)] -= dqdT1;
        ddNdtByVdcTp[sr0b*(this->alignN)+(this->nSpecies)] += dqdT1;
        ddNdtByVdcTp[sr1b*(this->alignN)+(this->nSpecies)] += dqdT1;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF1-invKc1*CR1);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j],sl1v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j],sr1v);
            }
        }
        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2   = this->RR22DupAllIndex[k*5+10];
        const unsigned int sl0c = this->RR22DupAllIndex[k*5+11];
        const unsigned int sl1c = this->RR22DupAllIndex[k*5+12];
        const unsigned int sr0c = this->RR22DupAllIndex[k*5+13];
        const unsigned int sr1c = this->RR22DupAllIndex[k*5+14];
        const double Kf2 = this->Kf_[i2];
        const double dCfdC0c = C[sl1c];
        const double dCfdC1c = C[sl0c];
        ddNdtByVdcTp[sl0c*(this->alignN)+ sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl1c*(this->alignN)+ sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl0c*(this->alignN)+ sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sl1c*(this->alignN)+ sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sr0c*(this->alignN)+ sl0c] += Kf2*dCfdC0c;
        ddNdtByVdcTp[sr1c*(this->alignN)+ sl0c] += Kf2*dCfdC0c;
        ddNdtByVdcTp[sr0c*(this->alignN)+ sl1c] += Kf2*dCfdC1c;
        ddNdtByVdcTp[sr1c*(this->alignN)+ sl1c] += Kf2*dCfdC1c;

        const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]);
        const double invKc2 = std::min(invKp2,FastChemistry::invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double dCrdC0c = C[sr1c];
        const double dCrdC1c = C[sr0c];
        ddNdtByVdcTp[sl0c*(this->alignN)+ sr0c] -= (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sl1c*(this->alignN)+ sr0c] -= (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sl0c*(this->alignN)+ sr1c] -= (-Kr2*dCrdC1c);
        ddNdtByVdcTp[sl1c*(this->alignN)+ sr1c] -= (-Kr2*dCrdC1c);
        ddNdtByVdcTp[sr0c*(this->alignN)+ sr0c] += (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sr1c*(this->alignN)+ sr0c] += (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sr0c*(this->alignN)+ sr1c] += (-Kr2*dCrdC1c);
        ddNdtByVdcTp[sr1c*(this->alignN)+ sr1c] += (-Kr2*dCrdC1c);

        const double dKfdT2 = this->dKfdT_[i2];
        const double sumVdBdT2 = (dBdT[sr0c] + dBdT[sr1c]) - (dBdT[sl0c] + dBdT[sl1c]);
        const double dKcdTByKc2 = sumVdBdT2;
        const double dKrdT2 = (dKfdT2*invKc2 - (invKc2 < FastChemistry::invKcLimiter ? Kr2*dKcdTByKc2 : 0));
        const double CF2 = C[sl0c]*C[sl1c];
        const double CR2 = C[sr0c]*C[sr1c];
        const double dqdT2 = (dKfdT2*CF2) - (dKrdT2*CR2);
        ddNdtByVdcTp[sl0c*(this->alignN)+(this->nSpecies)] -= dqdT2;
        ddNdtByVdcTp[sl1c*(this->alignN)+(this->nSpecies)] -= dqdT2;
        ddNdtByVdcTp[sr0c*(this->alignN)+(this->nSpecies)] += dqdT2;
        ddNdtByVdcTp[sr1c*(this->alignN)+(this->nSpecies)] += dqdT2;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF2-invKc2*CR2);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j],sl1v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j],sr1v);
            }
        }
        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3   = this->RR22DupAllIndex[k*5+15];
        const unsigned int sl0d = this->RR22DupAllIndex[k*5+16];
        const unsigned int sl1d = this->RR22DupAllIndex[k*5+17];
        const unsigned int sr0d = this->RR22DupAllIndex[k*5+18];
        const unsigned int sr1d = this->RR22DupAllIndex[k*5+19];
        const double Kf3 = this->Kf_[i3];
        const double dCfdC0d = C[sl1d];
        const double dCfdC1d = C[sl0d];
        ddNdtByVdcTp[sl0d*(this->alignN)+ sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl1d*(this->alignN)+ sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl0d*(this->alignN)+ sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sl1d*(this->alignN)+ sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sr0d*(this->alignN)+ sl0d] += Kf3*dCfdC0d;
        ddNdtByVdcTp[sr1d*(this->alignN)+ sl0d] += Kf3*dCfdC0d;
        ddNdtByVdcTp[sr0d*(this->alignN)+ sl1d] += Kf3*dCfdC1d;
        ddNdtByVdcTp[sr1d*(this->alignN)+ sl1d] += Kf3*dCfdC1d;

        const double invKp3 = (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d])*(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]);
        const double invKc3 = std::min(invKp3,FastChemistry::invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double dCrdC0d = C[sr1d];
        const double dCrdC1d = C[sr0d];
        ddNdtByVdcTp[sl0d*(this->alignN)+ sr0d] -= (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sl1d*(this->alignN)+ sr0d] -= (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sl0d*(this->alignN)+ sr1d] -= (-Kr3*dCrdC1d);
        ddNdtByVdcTp[sl1d*(this->alignN)+ sr1d] -= (-Kr3*dCrdC1d);
        ddNdtByVdcTp[sr0d*(this->alignN)+ sr0d] += (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sr1d*(this->alignN)+ sr0d] += (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sr0d*(this->alignN)+ sr1d] += (-Kr3*dCrdC1d);
        ddNdtByVdcTp[sr1d*(this->alignN)+ sr1d] += (-Kr3*dCrdC1d);

        const double dKfdT3 = this->dKfdT_[i3];
        const double sumVdBdT3 = (dBdT[sr0d] + dBdT[sr1d]) - (dBdT[sl0d] + dBdT[sl1d]);
        const double dKcdTByKc3 = sumVdBdT3;
        const double dKrdT3 = (dKfdT3*invKc3 - (invKc3 < FastChemistry::invKcLimiter ? Kr3*dKcdTByKc3 : 0));
        const double CF3 = C[sl0d]*C[sl1d];
        const double CR3 = C[sr0d]*C[sr1d];
        const double dqdT3 = (dKfdT3*CF3) - (dKrdT3*CR3);
        ddNdtByVdcTp[sl0d*(this->alignN)+(this->nSpecies)] -= dqdT3;
        ddNdtByVdcTp[sl1d*(this->alignN)+(this->nSpecies)] -= dqdT3;
        ddNdtByVdcTp[sr0d*(this->alignN)+(this->nSpecies)] += dqdT3;
        ddNdtByVdcTp[sr1d*(this->alignN)+(this->nSpecies)] += dqdT3;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF3-invKc3*CR3);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j],sl1v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j],sr1v);
            }
        }

    }
    for(std::size_t k=endRRDup22-remainRRDup22; k<endRRDup22; k=k+1)
    {
        const unsigned int i0 = this->RR22DupAllIndex[k*5+0];
        const unsigned int sl0a = this->RR22DupAllIndex[k*5+1];
        const unsigned int sl1a = this->RR22DupAllIndex[k*5+2];
        const unsigned int sr0a = this->RR22DupAllIndex[k*5+3];
        const unsigned int sr1a = this->RR22DupAllIndex[k*5+4];
        const double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C[sl1a];
        const double dCfdC1a = C[sl0a];
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl1a] += Kf0*dCfdC1a;

        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
        const double invKc0 = std::min(invKp0,FastChemistry::invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a];
        ddNdtByVdcTp[sl0a*(this->alignN)+ sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl1a*(this->alignN)+ sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl0a*(this->alignN)+ sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sl1a*(this->alignN)+ sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sr0a*(this->alignN)+ sr0a] += (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sr1a*(this->alignN)+ sr0a] += (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sr0a*(this->alignN)+ sr1a] += (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sr1a*(this->alignN)+ sr1a] += (-Kr0*dCrdC1a);

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a]) - (dBdT[sl0a] + dBdT[sl1a]);
        const double dKcdTByKc0 = sumVdBdT0;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < FastChemistry::invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C[sl0a]*C[sl1a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        ddNdtByVdcTp[sl0a*(this->alignN)+(this->nSpecies)] -= dqdT0;
        ddNdtByVdcTp[sl1a*(this->alignN)+(this->nSpecies)] -= dqdT0;
        ddNdtByVdcTp[sr0a*(this->alignN)+(this->nSpecies)] += dqdT0;
        ddNdtByVdcTp[sr1a*(this->alignN)+(this->nSpecies)] += dqdT0;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF0-invKc0*CR0);
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j],sr1v);
            }
        }
    }
    }
}

void FastChemistry::OptReaction::JF22IR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
) const noexcept
{

    std::size_t endIR22 = this->IR22Size;
    std::size_t remainIR22 = endIR22 % 4;

    for(std::size_t k=0; k<endIR22-remainIR22; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->IR22AllIndex[k*5+0];
        const unsigned int sl0a = this->IR22AllIndex[k*5+1];
        const unsigned int sl1a = this->IR22AllIndex[k*5+2];
        const unsigned int sr0a = this->IR22AllIndex[k*5+3];
        const unsigned int sr1a = this->IR22AllIndex[k*5+4];

        const double Kf0 = this->Kf_[i0];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double Kf0dCfdC0a = Kf0*C_sl1a;
        const double Kf0dCfdC1a = Kf0*C_sl0a;
        const double CF0 = C_sl0a*C_sl1a;
        const double dqdT0 = this->dKfdT_[i0]*CF0;
        ddNdtByVdcTp[sl0a*(this->alignN)+(this->nSpecies)] -= dqdT0;
        ddNdtByVdcTp[sl1a*(this->alignN)+(this->nSpecies)] -= dqdT0;
        ddNdtByVdcTp[sr0a*(this->alignN)+(this->nSpecies)] += dqdT0;
        ddNdtByVdcTp[sr1a*(this->alignN)+(this->nSpecies)] += dqdT0;

        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        ddNdtByVdcTp[sl0a*(this->alignN)+ sl0a] -= Kf0dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl0a] -= Kf0dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl1a] -= Kf0dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl1a] -= Kf0dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl0a] += Kf0dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl0a] += Kf0dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl1a] += Kf0dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl1a] += Kf0dCfdC1a;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF0;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j],sr0v);
                store256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR22AllIndex[k*5+5];
        const unsigned int sl0b = this->IR22AllIndex[k*5+6];
        const unsigned int sl1b = this->IR22AllIndex[k*5+7];
        const unsigned int sr0b = this->IR22AllIndex[k*5+8];
        const unsigned int sr1b = this->IR22AllIndex[k*5+9];

        const double Kf1 = this->Kf_[i1];
        const double C_sl1b = C[sl1b];
        const double C_sl0b = C[sl0b];
        const double Kf1dCfdC0b = Kf1*C_sl1b;
        const double Kf1dCfdC1b = Kf1*C_sl0b;
        const double CF1 = C_sl0b*C_sl1b;
        const double dqdT1 = this->dKfdT_[i1]*CF1;
        ddNdtByVdcTp[sl0b*(this->alignN)+(this->nSpecies)] -= dqdT1;
        ddNdtByVdcTp[sl1b*(this->alignN)+(this->nSpecies)] -= dqdT1;
        ddNdtByVdcTp[sr0b*(this->alignN)+(this->nSpecies)] += dqdT1;
        ddNdtByVdcTp[sr1b*(this->alignN)+(this->nSpecies)] += dqdT1;

        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        ddNdtByVdcTp[sl0b*(this->alignN)+ sl0b] -= Kf1dCfdC0b;
        ddNdtByVdcTp[sl1b*(this->alignN)+ sl0b] -= Kf1dCfdC0b;
        ddNdtByVdcTp[sl0b*(this->alignN)+ sl1b] -= Kf1dCfdC1b;
        ddNdtByVdcTp[sl1b*(this->alignN)+ sl1b] -= Kf1dCfdC1b;
        ddNdtByVdcTp[sr0b*(this->alignN)+ sl0b] += Kf1dCfdC0b;
        ddNdtByVdcTp[sr1b*(this->alignN)+ sl0b] += Kf1dCfdC0b;
        ddNdtByVdcTp[sr0b*(this->alignN)+ sl1b] += Kf1dCfdC1b;
        ddNdtByVdcTp[sr1b*(this->alignN)+ sl1b] += Kf1dCfdC1b;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF1;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j],sl0v);
                store256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j],sr0v);
                store256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR22AllIndex[k*5+10];
        const unsigned int sl0c = this->IR22AllIndex[k*5+11];
        const unsigned int sl1c = this->IR22AllIndex[k*5+12];
        const unsigned int sr0c = this->IR22AllIndex[k*5+13];
        const unsigned int sr1c = this->IR22AllIndex[k*5+14];

        const double Kf2 = this->Kf_[i2];
        const double C_sl1c = C[sl1c];
        const double C_sl0c = C[sl0c];
        const double Kf2dCfdC0c = Kf2*C_sl1c;
        const double Kf2dCfdC1c = Kf2*C_sl0c;
        const double CF2 = C_sl0c*C_sl1c;
        const double dqdT2 = this->dKfdT_[i2]*CF2;
        ddNdtByVdcTp[sl0c*(this->alignN)+(this->nSpecies)] -= dqdT2;
        ddNdtByVdcTp[sl1c*(this->alignN)+(this->nSpecies)] -= dqdT2;
        ddNdtByVdcTp[sr0c*(this->alignN)+(this->nSpecies)] += dqdT2;
        ddNdtByVdcTp[sr1c*(this->alignN)+(this->nSpecies)] += dqdT2;

        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        ddNdtByVdcTp[sl0c*(this->alignN)+ sl0c] -= Kf2dCfdC0c;
        ddNdtByVdcTp[sl1c*(this->alignN)+ sl0c] -= Kf2dCfdC0c;
        ddNdtByVdcTp[sl0c*(this->alignN)+ sl1c] -= Kf2dCfdC1c;
        ddNdtByVdcTp[sl1c*(this->alignN)+ sl1c] -= Kf2dCfdC1c;
        ddNdtByVdcTp[sr0c*(this->alignN)+ sl0c] += Kf2dCfdC0c;
        ddNdtByVdcTp[sr1c*(this->alignN)+ sl0c] += Kf2dCfdC0c;
        ddNdtByVdcTp[sr0c*(this->alignN)+ sl1c] += Kf2dCfdC1c;
        ddNdtByVdcTp[sr1c*(this->alignN)+ sl1c] += Kf2dCfdC1c;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF2;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j],sl0v);
                store256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j],sr0v);
                store256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR22AllIndex[k*5+15];
        const unsigned int sl0d = this->IR22AllIndex[k*5+16];
        const unsigned int sl1d = this->IR22AllIndex[k*5+17];
        const unsigned int sr0d = this->IR22AllIndex[k*5+18];
        const unsigned int sr1d = this->IR22AllIndex[k*5+19];

        const double Kf3 = this->Kf_[i3];
        const double C_sl1d = C[sl1d];
        const double C_sl0d = C[sl0d];
        const double Kf3dCfdC0d = Kf3*C_sl1d;
        const double Kf3dCfdC1d = Kf3*C_sl0d;
        const double CF3 = C_sl0d*C_sl1d;
        const double dqdT3 = this->dKfdT_[i3]*CF3;
        ddNdtByVdcTp[sl0d*(this->alignN)+(this->nSpecies)] -= dqdT3;
        ddNdtByVdcTp[sl1d*(this->alignN)+(this->nSpecies)] -= dqdT3;
        ddNdtByVdcTp[sr0d*(this->alignN)+(this->nSpecies)] += dqdT3;
        ddNdtByVdcTp[sr1d*(this->alignN)+(this->nSpecies)] += dqdT3;

        const double q3 = (Kf3*CF3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

        ddNdtByVdcTp[sl0d*(this->alignN)+ sl0d] -= Kf3dCfdC0d;
        ddNdtByVdcTp[sl1d*(this->alignN)+ sl0d] -= Kf3dCfdC0d;
        ddNdtByVdcTp[sl0d*(this->alignN)+ sl1d] -= Kf3dCfdC1d;
        ddNdtByVdcTp[sl1d*(this->alignN)+ sl1d] -= Kf3dCfdC1d;
        ddNdtByVdcTp[sr0d*(this->alignN)+ sl0d] += Kf3dCfdC0d;
        ddNdtByVdcTp[sr1d*(this->alignN)+ sl0d] += Kf3dCfdC0d;
        ddNdtByVdcTp[sr0d*(this->alignN)+ sl1d] += Kf3dCfdC1d;
        ddNdtByVdcTp[sr1d*(this->alignN)+ sl1d] += Kf3dCfdC1d;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF3;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j],sl0v);
                store256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j],sr0v);
                store256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j],sr1v);
            }
        }
    }

    for(std::size_t k=endIR22-remainIR22; k<endIR22; k=k+1)
    {
        const unsigned int i0 = this->IR22AllIndex[k*5+0];
        const unsigned int sl0a = this->IR22AllIndex[k*5+1];
        const unsigned int sl1a = this->IR22AllIndex[k*5+2];
        const unsigned int sr0a = this->IR22AllIndex[k*5+3];
        const unsigned int sr1a = this->IR22AllIndex[k*5+4];

        const double Kf0 = this->Kf_[i0];
        const double CF0 = C[sl0a]*C[sl1a];
        const double dqdT0 = this->dKfdT_[i0]*CF0;
        ddNdtByVdcTp[sl0a*(this->alignN)+(this->nSpecies)] -= dqdT0;
        ddNdtByVdcTp[sl1a*(this->alignN)+(this->nSpecies)] -= dqdT0;
        ddNdtByVdcTp[sr0a*(this->alignN)+(this->nSpecies)] += dqdT0;
        ddNdtByVdcTp[sr1a*(this->alignN)+(this->nSpecies)] += dqdT0;

        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const double dCfdC0a = C[sl1a];
        const double dCfdC1a = C[sl0a];
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl1a] += Kf0*dCfdC1a;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF0;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j],sr0v);
                store256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j],sr1v);
            }
        }
    }

    std::size_t endIR22Dup = this->IR22DupSize;
    std::size_t remainIR22Dup = endIR22Dup % 4;
    for(std::size_t k=0; k<endIR22Dup-remainIR22Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->IR22DupAllIndex[k*5+0];
        const unsigned int sl0a = this->IR22DupAllIndex[k*5+1];
        const unsigned int sl1a = this->IR22DupAllIndex[k*5+2];
        const unsigned int sr0a = this->IR22DupAllIndex[k*5+3];
        const unsigned int sr1a = this->IR22DupAllIndex[k*5+4];

        const double Kf0 = this->Kf_[i0];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double Kf0dCfdC0a = Kf0*C_sl1a;
        const double Kf0dCfdC1a = Kf0*C_sl0a;
        const double CF0 = C_sl0a*C_sl1a;
        const double dqdT0 = this->dKfdT_[i0]*CF0;
        ddNdtByVdcTp[sl0a*(this->alignN)+(this->nSpecies)] -= dqdT0;
        ddNdtByVdcTp[sl1a*(this->alignN)+(this->nSpecies)] -= dqdT0;
        ddNdtByVdcTp[sr0a*(this->alignN)+(this->nSpecies)] += dqdT0;
        ddNdtByVdcTp[sr1a*(this->alignN)+(this->nSpecies)] += dqdT0;

        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        ddNdtByVdcTp[sl0a*(this->alignN)+ sl0a] -= Kf0dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl0a] -= Kf0dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl1a] -= Kf0dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl1a] -= Kf0dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl0a] += Kf0dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl0a] += Kf0dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl1a] += Kf0dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl1a] += Kf0dCfdC1a;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF0;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j],sr0v);
                store256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR22DupAllIndex[k*5+5];
        const unsigned int sl0b = this->IR22DupAllIndex[k*5+6];
        const unsigned int sl1b = this->IR22DupAllIndex[k*5+7];
        const unsigned int sr0b = this->IR22DupAllIndex[k*5+8];
        const unsigned int sr1b = this->IR22DupAllIndex[k*5+9];

        const double Kf1 = this->Kf_[i1];
        const double C_sl1b = C[sl1b];
        const double C_sl0b = C[sl0b];
        const double Kf1dCfdC0b = Kf1*C_sl1b;
        const double Kf1dCfdC1b = Kf1*C_sl0b;
        const double CF1 = C_sl0b*C_sl1b;
        const double dqdT1 = this->dKfdT_[i1]*CF1;
        ddNdtByVdcTp[sl0b*(this->alignN)+(this->nSpecies)] -= dqdT1;
        ddNdtByVdcTp[sl1b*(this->alignN)+(this->nSpecies)] -= dqdT1;
        ddNdtByVdcTp[sr0b*(this->alignN)+(this->nSpecies)] += dqdT1;
        ddNdtByVdcTp[sr1b*(this->alignN)+(this->nSpecies)] += dqdT1;

        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        ddNdtByVdcTp[sl0b*(this->alignN)+ sl0b] -= Kf1dCfdC0b;
        ddNdtByVdcTp[sl1b*(this->alignN)+ sl0b] -= Kf1dCfdC0b;
        ddNdtByVdcTp[sl0b*(this->alignN)+ sl1b] -= Kf1dCfdC1b;
        ddNdtByVdcTp[sl1b*(this->alignN)+ sl1b] -= Kf1dCfdC1b;
        ddNdtByVdcTp[sr0b*(this->alignN)+ sl0b] += Kf1dCfdC0b;
        ddNdtByVdcTp[sr1b*(this->alignN)+ sl0b] += Kf1dCfdC0b;
        ddNdtByVdcTp[sr0b*(this->alignN)+ sl1b] += Kf1dCfdC1b;
        ddNdtByVdcTp[sr1b*(this->alignN)+ sl1b] += Kf1dCfdC1b;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF1;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j],sl0v);
                store256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j],sr0v);
                store256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR22DupAllIndex[k*5+10];
        const unsigned int sl0c = this->IR22DupAllIndex[k*5+11];
        const unsigned int sl1c = this->IR22DupAllIndex[k*5+12];
        const unsigned int sr0c = this->IR22DupAllIndex[k*5+13];
        const unsigned int sr1c = this->IR22DupAllIndex[k*5+14];

        const double Kf2 = this->Kf_[i2];
        const double C_sl1c = C[sl1c];
        const double C_sl0c = C[sl0c];
        const double Kf2dCfdC0c = Kf2*C_sl1c;
        const double Kf2dCfdC1c = Kf2*C_sl0c;
        const double CF2 = C_sl0c*C_sl1c;
        const double dqdT2 = this->dKfdT_[i2]*CF2;
        ddNdtByVdcTp[sl0c*(this->alignN)+(this->nSpecies)] -= dqdT2;
        ddNdtByVdcTp[sl1c*(this->alignN)+(this->nSpecies)] -= dqdT2;
        ddNdtByVdcTp[sr0c*(this->alignN)+(this->nSpecies)] += dqdT2;
        ddNdtByVdcTp[sr1c*(this->alignN)+(this->nSpecies)] += dqdT2;

        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        ddNdtByVdcTp[sl0c*(this->alignN)+ sl0c] -= Kf2dCfdC0c;
        ddNdtByVdcTp[sl1c*(this->alignN)+ sl0c] -= Kf2dCfdC0c;
        ddNdtByVdcTp[sl0c*(this->alignN)+ sl1c] -= Kf2dCfdC1c;
        ddNdtByVdcTp[sl1c*(this->alignN)+ sl1c] -= Kf2dCfdC1c;
        ddNdtByVdcTp[sr0c*(this->alignN)+ sl0c] += Kf2dCfdC0c;
        ddNdtByVdcTp[sr1c*(this->alignN)+ sl0c] += Kf2dCfdC0c;
        ddNdtByVdcTp[sr0c*(this->alignN)+ sl1c] += Kf2dCfdC1c;
        ddNdtByVdcTp[sr1c*(this->alignN)+ sl1c] += Kf2dCfdC1c;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF2;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j],sl0v);
                store256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j],sr0v);
                store256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR22DupAllIndex[k*5+15];
        const unsigned int sl0d = this->IR22DupAllIndex[k*5+16];
        const unsigned int sl1d = this->IR22DupAllIndex[k*5+17];
        const unsigned int sr0d = this->IR22DupAllIndex[k*5+18];
        const unsigned int sr1d = this->IR22DupAllIndex[k*5+19];

        const double Kf3 = this->Kf_[i3];
        const double C_sl1d = C[sl1d];
        const double C_sl0d = C[sl0d];
        const double Kf3dCfdC0d = Kf3*C_sl1d;
        const double Kf3dCfdC1d = Kf3*C_sl0d;
        const double CF3 = C_sl0d*C_sl1d;
        const double dqdT3 = this->dKfdT_[i3]*CF3;
        ddNdtByVdcTp[sl0d*(this->alignN)+(this->nSpecies)] -= dqdT3;
        ddNdtByVdcTp[sl1d*(this->alignN)+(this->nSpecies)] -= dqdT3;
        ddNdtByVdcTp[sr0d*(this->alignN)+(this->nSpecies)] += dqdT3;
        ddNdtByVdcTp[sr1d*(this->alignN)+(this->nSpecies)] += dqdT3;

        const double q3 = (Kf3*CF3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

        ddNdtByVdcTp[sl0d*(this->alignN)+ sl0d] -= Kf3dCfdC0d;
        ddNdtByVdcTp[sl1d*(this->alignN)+ sl0d] -= Kf3dCfdC0d;
        ddNdtByVdcTp[sl0d*(this->alignN)+ sl1d] -= Kf3dCfdC1d;
        ddNdtByVdcTp[sl1d*(this->alignN)+ sl1d] -= Kf3dCfdC1d;
        ddNdtByVdcTp[sr0d*(this->alignN)+ sl0d] += Kf3dCfdC0d;
        ddNdtByVdcTp[sr1d*(this->alignN)+ sl0d] += Kf3dCfdC0d;
        ddNdtByVdcTp[sr0d*(this->alignN)+ sl1d] += Kf3dCfdC1d;
        ddNdtByVdcTp[sr1d*(this->alignN)+ sl1d] += Kf3dCfdC1d;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF3;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j],sl0v);
                store256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j],sr0v);
                store256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j],sr1v);
            }
        }
    }

    for(std::size_t k=endIR22Dup-remainIR22Dup; k<endIR22Dup; k=k+1)
    {
        const unsigned int i0 = this->IR22DupAllIndex[k*5+0];
        const unsigned int sl0a = this->IR22DupAllIndex[k*5+1];
        const unsigned int sl1a = this->IR22DupAllIndex[k*5+2];
        const unsigned int sr0a = this->IR22DupAllIndex[k*5+3];
        const unsigned int sr1a = this->IR22DupAllIndex[k*5+4];

        const double Kf0 = this->Kf_[i0];
        const double CF0 = C[sl0a]*C[sl1a];
        const double dqdT0 = this->dKfdT_[i0]*CF0;
        ddNdtByVdcTp[sl0a*(this->alignN)+(this->nSpecies)] -= dqdT0;
        ddNdtByVdcTp[sl1a*(this->alignN)+(this->nSpecies)] -= dqdT0;
        ddNdtByVdcTp[sr0a*(this->alignN)+(this->nSpecies)] += dqdT0;
        ddNdtByVdcTp[sr1a*(this->alignN)+(this->nSpecies)] += dqdT0;

        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const double dCfdC0a = C[sl1a];
        const double dCfdC1a = C[sl0a];
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl1a] += Kf0*dCfdC1a;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CF0;
            __m256d Wv = _mm256_set1_pd(W);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv,dMdC);
                __m256d sl0v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]),WdMdC);
                __m256d sl1v = _mm256_sub_pd(load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);
                __m256d sr0v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]),WdMdC);
                store256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j],sr0v);
                store256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j],sr1v);
            }
        }
    }
}

void  FastChemistry::OptReaction::JF22NER
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t endNER22 = this->NER22Size;
    std::size_t remainNER22 = endNER22%4;
    for(std::size_t k=0; k<endNER22-remainNER22; k=k+4)
    {

        const unsigned int i0 = this->NER22AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER22AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER22AllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->NER22AllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->NER22AllIndex[(k+0)*5+4];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C[sl1a];
        const double dCfdC1a = C[sl0a];
        const double Kf0dCfdC0a = Kf0*dCfdC0a;
        const double Kf0dCfdC1a = Kf0*dCfdC1a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0dCfdC0a;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1a;
        const double J_sr1a_sl1a_increment = J_sr1a_ptr[sl1a] + Kf0dCfdC1a;

        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;
        J_sr1a_ptr[sl1a] = J_sr1a_sl1a_increment;

        double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a];
        const double Kr0dCrdC0a = Kr0*dCrdC0a;
        const double Kr0dCrdC1a = Kr0*dCrdC1a;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0dCrdC0a;
        const double J_sl1a_sr0a_increment = J_sl1a_ptr[sr0a] + Kr0dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0dCrdC1a;
        const double J_sl1a_sr1a_increment = J_sl1a_ptr[sr1a] + Kr0dCrdC1a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0dCrdC1a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sl1a_ptr[sr0a] = J_sl1a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sl1a_ptr[sr1a] = J_sl1a_sr1a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a]*C[sl1a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;

        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;


        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;

        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[2] && i0 <this->Ikf[3])
        { 
            const unsigned int idx = i0 - this->Ikf[2];
            double W = this->dKfdC_[idx]*CF0 - this->dKfdC_[idx+this->Itbr[4]]*CR0;
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
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1a_ptr[j+0]),WdMdC0);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1a_ptr[j+4]),WdMdC4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl1a_ptr[j+0],sl1v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sl1a_ptr[j+4],sl1v4);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
            } 
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1a_ptr[j+0]),WdMdC0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl1a_ptr[j+0],sl1v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
            }
        }

        const unsigned int i1 = this->NER22AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER22AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER22AllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->NER22AllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->NER22AllIndex[(k+1)*5+4];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sl1b_ptr = &ddNdtByVdcTp[sl1b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];

        double Kf1 = this->Kf_[i1];
        const double dCfdC0b = C[sl1b];
        const double dCfdC1b = C[sl0b];
        const double Kf1dCfdC0b = Kf1*dCfdC0b;
        const double Kf1dCfdC1b = Kf1*dCfdC1b;
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1dCfdC0b;
        const double J_sl1b_sl0b_increment = J_sl1b_ptr[sl0b] - Kf1dCfdC0b;
        const double J_sl0b_sl1b_increment = J_sl0b_ptr[sl1b] - Kf1dCfdC1b;
        const double J_sl1b_sl1b_increment = J_sl1b_ptr[sl1b] - Kf1dCfdC1b;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1dCfdC0b;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1dCfdC0b;
        const double J_sr0b_sl1b_increment = J_sr0b_ptr[sl1b] + Kf1dCfdC1b;
        const double J_sr1b_sl1b_increment = J_sr1b_ptr[sl1b] + Kf1dCfdC1b;

        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sl1b_ptr[sl0b] = J_sl1b_sl0b_increment;
        J_sl0b_ptr[sl1b] = J_sl0b_sl1b_increment;
        J_sl1b_ptr[sl1b] = J_sl1b_sl1b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;
        J_sr0b_ptr[sl1b] = J_sr0b_sl1b_increment;
        J_sr1b_ptr[sl1b] = J_sr1b_sl1b_increment;

        double Kr1 = this->Kf_[i1 - Ikf[1] + Ikf[9]];
        const double dCrdC0b = C[sr1b];
        const double dCrdC1b = C[sr0b];
        const double Kr1dCrdC0b = Kr1*dCrdC0b;
        const double Kr1dCrdC1b = Kr1*dCrdC1b;
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1dCrdC0b;
        const double J_sl1b_sr0b_increment = J_sl1b_ptr[sr0b] + Kr1dCrdC0b;
        const double J_sl0b_sr1b_increment = J_sl0b_ptr[sr1b] + Kr1dCrdC1b;
        const double J_sl1b_sr1b_increment = J_sl1b_ptr[sr1b] + Kr1dCrdC1b;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1dCrdC0b;
        const double J_sr1b_sr0b_increment = J_sr1b_ptr[sr0b] - Kr1dCrdC0b;
        const double J_sr0b_sr1b_increment = J_sr0b_ptr[sr1b] - Kr1dCrdC1b;
        const double J_sr1b_sr1b_increment = J_sr1b_ptr[sr1b] - Kr1dCrdC1b;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sl1b_ptr[sr0b] = J_sl1b_sr0b_increment;
        J_sl0b_ptr[sr1b] = J_sl0b_sr1b_increment;
        J_sl1b_ptr[sr1b] = J_sl1b_sr1b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;
        J_sr1b_ptr[sr0b] = J_sr1b_sr0b_increment;
        J_sr0b_ptr[sr1b] = J_sr0b_sr1b_increment;
        J_sr1b_ptr[sr1b] = J_sr1b_sr1b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double dKrdT1 = this->dKfdT_[i1 - Ikf[1] + Ikf[9]];
        const double CF1 = C[sl0b]*C[sl1b];
        const double CR1 = C[sr0b]*C[sr1b];
        const double dqdT1 = (dKfdT1*CF1)-(dKrdT1*CR1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sl1b_T_increment = J_sl1b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1;

        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sl1b_ptr[this->nSpecies] = J_sl1b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;


        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;

        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        if(i1>=this->Ikf[2] && i1 <this->Ikf[3])
        { 
            const unsigned int idx = i1 - this->Ikf[2];
            double W = this->dKfdC_[idx]*CF1 - this->dKfdC_[idx+this->Itbr[4]]*CR1;
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);

                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d WdMdC4 = _mm256_mul_pd(Wv,dMdC4);

                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+4]),WdMdC4);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1b_ptr[j+0]),WdMdC0);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1b_ptr[j+4]),WdMdC4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+4]),WdMdC4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+4]),WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl1b_ptr[j+0],sl1v0);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sr1b_ptr[j+0],sr1v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
                store256d(&J_sl1b_ptr[j+4],sl1v4);
                store256d(&J_sr0b_ptr[j+4],sr0v4);
                store256d(&J_sr1b_ptr[j+4],sr1v4);
            } 
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1b_ptr[j+0]),WdMdC0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl1b_ptr[j+0],sl1v0);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sr1b_ptr[j+0],sr1v0);
            }
        }

        const unsigned int i2 = this->NER22AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER22AllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->NER22AllIndex[(k+2)*5+2];
        const unsigned int sr0c = this->NER22AllIndex[(k+2)*5+3];
        const unsigned int sr1c = this->NER22AllIndex[(k+2)*5+4];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sl1c_ptr = &ddNdtByVdcTp[sl1c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];

        double Kf2 = this->Kf_[i2];
        const double dCfdC0c = C[sl1c];
        const double dCfdC1c = C[sl0c];
        const double Kf2dCfdC0c = Kf2*dCfdC0c;
        const double Kf2dCfdC1c = Kf2*dCfdC1c;
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2dCfdC0c;
        const double J_sl1c_sl0c_increment = J_sl1c_ptr[sl0c] - Kf2dCfdC0c;
        const double J_sl0c_sl1c_increment = J_sl0c_ptr[sl1c] - Kf2dCfdC1c;
        const double J_sl1c_sl1c_increment = J_sl1c_ptr[sl1c] - Kf2dCfdC1c;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2dCfdC0c;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2dCfdC0c;
        const double J_sr0c_sl1c_increment = J_sr0c_ptr[sl1c] + Kf2dCfdC1c;
        const double J_sr1c_sl1c_increment = J_sr1c_ptr[sl1c] + Kf2dCfdC1c;

        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sl1c_ptr[sl0c] = J_sl1c_sl0c_increment;
        J_sl0c_ptr[sl1c] = J_sl0c_sl1c_increment;
        J_sl1c_ptr[sl1c] = J_sl1c_sl1c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;
        J_sr0c_ptr[sl1c] = J_sr0c_sl1c_increment;
        J_sr1c_ptr[sl1c] = J_sr1c_sl1c_increment;

        double Kr2 = this->Kf_[i2 - Ikf[1] + Ikf[9]];
        const double dCrdC0c = C[sr1c];
        const double dCrdC1c = C[sr0c];
        const double Kr2dCrdC0c = Kr2*dCrdC0c;
        const double Kr2dCrdC1c = Kr2*dCrdC1c;
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2dCrdC0c;
        const double J_sl1c_sr0c_increment = J_sl1c_ptr[sr0c] + Kr2dCrdC0c;
        const double J_sl0c_sr1c_increment = J_sl0c_ptr[sr1c] + Kr2dCrdC1c;
        const double J_sl1c_sr1c_increment = J_sl1c_ptr[sr1c] + Kr2dCrdC1c;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2dCrdC0c;
        const double J_sr1c_sr0c_increment = J_sr1c_ptr[sr0c] - Kr2dCrdC0c;
        const double J_sr0c_sr1c_increment = J_sr0c_ptr[sr1c] - Kr2dCrdC1c;
        const double J_sr1c_sr1c_increment = J_sr1c_ptr[sr1c] - Kr2dCrdC1c;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sl1c_ptr[sr0c] = J_sl1c_sr0c_increment;
        J_sl0c_ptr[sr1c] = J_sl0c_sr1c_increment;
        J_sl1c_ptr[sr1c] = J_sl1c_sr1c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;
        J_sr1c_ptr[sr0c] = J_sr1c_sr0c_increment;
        J_sr0c_ptr[sr1c] = J_sr0c_sr1c_increment;
        J_sr1c_ptr[sr1c] = J_sr1c_sr1c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double dKrdT2 = this->dKfdT_[i2 - Ikf[1] + Ikf[9]];
        const double CF2 = C[sl0c]*C[sl1c];
        const double CR2 = C[sr0c]*C[sr1c];
        const double dqdT2 = (dKfdT2*CF2)-(dKrdT2*CR2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sl1c_T_increment = J_sl1c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2;

        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sl1c_ptr[this->nSpecies] = J_sl1c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;


        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;

        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        if(i2>=this->Ikf[2] && i2 <this->Ikf[3])
        { 
            const unsigned int idx = i2 - this->Ikf[2];
            double W = this->dKfdC_[idx]*CF2 - this->dKfdC_[idx+this->Itbr[4]]*CR2;
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);

                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d WdMdC4 = _mm256_mul_pd(Wv,dMdC4);

                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+4]),WdMdC4);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1c_ptr[j+0]),WdMdC0);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1c_ptr[j+4]),WdMdC4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+4]),WdMdC4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+4]),WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl1c_ptr[j+0],sl1v0);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sr1c_ptr[j+0],sr1v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
                store256d(&J_sl1c_ptr[j+4],sl1v4);
                store256d(&J_sr0c_ptr[j+4],sr0v4);
                store256d(&J_sr1c_ptr[j+4],sr1v4);
            } 
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1c_ptr[j+0]),WdMdC0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl1c_ptr[j+0],sl1v0);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sr1c_ptr[j+0],sr1v0);
            }
        }

        const unsigned int i3 = this->NER22AllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->NER22AllIndex[(k+3)*5+1];
        const unsigned int sl1d = this->NER22AllIndex[(k+3)*5+2];
        const unsigned int sr0d = this->NER22AllIndex[(k+3)*5+3];
        const unsigned int sr1d = this->NER22AllIndex[(k+3)*5+4];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sl1d_ptr = &ddNdtByVdcTp[sl1d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];

        double Kf3 = this->Kf_[i3];
        const double dCfdC0d = C[sl1d];
        const double dCfdC1d = C[sl0d];
        const double Kf3dCfdC0d = Kf3*dCfdC0d;
        const double Kf3dCfdC1d = Kf3*dCfdC1d;
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3dCfdC0d;
        const double J_sl1d_sl0d_increment = J_sl1d_ptr[sl0d] - Kf3dCfdC0d;
        const double J_sl0d_sl1d_increment = J_sl0d_ptr[sl1d] - Kf3dCfdC1d;
        const double J_sl1d_sl1d_increment = J_sl1d_ptr[sl1d] - Kf3dCfdC1d;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3dCfdC0d;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3dCfdC0d;
        const double J_sr0d_sl1d_increment = J_sr0d_ptr[sl1d] + Kf3dCfdC1d;
        const double J_sr1d_sl1d_increment = J_sr1d_ptr[sl1d] + Kf3dCfdC1d;

        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sl1d_ptr[sl0d] = J_sl1d_sl0d_increment;
        J_sl0d_ptr[sl1d] = J_sl0d_sl1d_increment;
        J_sl1d_ptr[sl1d] = J_sl1d_sl1d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;
        J_sr0d_ptr[sl1d] = J_sr0d_sl1d_increment;
        J_sr1d_ptr[sl1d] = J_sr1d_sl1d_increment;

        double Kr3 = this->Kf_[i3 - Ikf[1] + Ikf[9]];
        const double dCrdC0d = C[sr1d];
        const double dCrdC1d = C[sr0d];
        const double Kr3dCrdC0d = Kr3*dCrdC0d;
        const double Kr3dCrdC1d = Kr3*dCrdC1d;
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3dCrdC0d;
        const double J_sl1d_sr0d_increment = J_sl1d_ptr[sr0d] + Kr3dCrdC0d;
        const double J_sl0d_sr1d_increment = J_sl0d_ptr[sr1d] + Kr3dCrdC1d;
        const double J_sl1d_sr1d_increment = J_sl1d_ptr[sr1d] + Kr3dCrdC1d;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3dCrdC0d;
        const double J_sr1d_sr0d_increment = J_sr1d_ptr[sr0d] - Kr3dCrdC0d;
        const double J_sr0d_sr1d_increment = J_sr0d_ptr[sr1d] - Kr3dCrdC1d;
        const double J_sr1d_sr1d_increment = J_sr1d_ptr[sr1d] - Kr3dCrdC1d;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sl1d_ptr[sr0d] = J_sl1d_sr0d_increment;
        J_sl0d_ptr[sr1d] = J_sl0d_sr1d_increment;
        J_sl1d_ptr[sr1d] = J_sl1d_sr1d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;
        J_sr1d_ptr[sr0d] = J_sr1d_sr0d_increment;
        J_sr0d_ptr[sr1d] = J_sr0d_sr1d_increment;
        J_sr1d_ptr[sr1d] = J_sr1d_sr1d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double dKrdT3 = this->dKfdT_[i3 - Ikf[1] + Ikf[9]];
        const double CF3 = C[sl0d]*C[sl1d];
        const double CR3 = C[sr0d]*C[sr1d];
        const double dqdT3 = (dKfdT3*CF3)-(dKrdT3*CR3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sl1d_T_increment = J_sl1d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3;

        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sl1d_ptr[this->nSpecies] = J_sl1d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;


        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;

        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;

        if(i3>=this->Ikf[2] && i3 <this->Ikf[3])
        { 
            const unsigned int idx = i3 - this->Ikf[2];
            double W = this->dKfdC_[idx]*CF3 - this->dKfdC_[idx+this->Itbr[4]]*CR3;
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);

                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d WdMdC4 = _mm256_mul_pd(Wv,dMdC4);

                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+4]),WdMdC4);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1d_ptr[j+0]),WdMdC0);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1d_ptr[j+4]),WdMdC4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+4]),WdMdC4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+4]),WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl1d_ptr[j+0],sl1v0);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sr1d_ptr[j+0],sr1v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
                store256d(&J_sl1d_ptr[j+4],sl1v4);
                store256d(&J_sr0d_ptr[j+4],sr0v4);
                store256d(&J_sr1d_ptr[j+4],sr1v4);
            } 
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1d_ptr[j+0]),WdMdC0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl1d_ptr[j+0],sl1v0);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sr1d_ptr[j+0],sr1v0);
            }
        }
    }
    for(std::size_t k=endNER22-remainNER22; k<endNER22; k=k+1)
    {

        const unsigned int i0 = this->NER22AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER22AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER22AllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->NER22AllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->NER22AllIndex[(k+0)*5+4];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C[sl1a];
        const double dCfdC1a = C[sl0a];
        const double Kf0dCfdC0a = Kf0*dCfdC0a;
        const double Kf0dCfdC1a = Kf0*dCfdC1a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0dCfdC0a;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1a;
        const double J_sr1a_sl1a_increment = J_sr1a_ptr[sl1a] + Kf0dCfdC1a;

        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;
        J_sr1a_ptr[sl1a] = J_sr1a_sl1a_increment;

        double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a];
        const double Kr0dCrdC0a = Kr0*dCrdC0a;
        const double Kr0dCrdC1a = Kr0*dCrdC1a;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0dCrdC0a;
        const double J_sl1a_sr0a_increment = J_sl1a_ptr[sr0a] + Kr0dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0dCrdC1a;
        const double J_sl1a_sr1a_increment = J_sl1a_ptr[sr1a] + Kr0dCrdC1a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0dCrdC1a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sl1a_ptr[sr0a] = J_sl1a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sl1a_ptr[sr1a] = J_sl1a_sr1a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a]*C[sl1a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;

        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;


        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;

        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[2] && i0 <this->Ikf[3])
        { 
            const unsigned int idx = i0 - this->Ikf[2];
            double W = this->dKfdC_[idx]*CF0 - this->dKfdC_[idx+this->Itbr[4]]*CR0;
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
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1a_ptr[j+0]),WdMdC0);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1a_ptr[j+4]),WdMdC4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl1a_ptr[j+0],sl1v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sl1a_ptr[j+4],sl1v4);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
            } 
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1a_ptr[j+0]),WdMdC0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl1a_ptr[j+0],sl1v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
            }
        }
    }
    
    std::size_t endNER22Dup = this->NER22DupSize;
    std::size_t remainNER22Dup = endNER22Dup%4;
    for(std::size_t k=0; k<endNER22Dup-remainNER22Dup; k=k+4)
    {

        const unsigned int i0 = this->NER22DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER22DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER22DupAllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->NER22DupAllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->NER22DupAllIndex[(k+0)*5+4];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C[sl1a];
        const double dCfdC1a = C[sl0a];
        const double Kf0dCfdC0a = Kf0*dCfdC0a;
        const double Kf0dCfdC1a = Kf0*dCfdC1a;
        J_sl0a_ptr[sl0a] -= Kf0dCfdC0a;
        J_sl1a_ptr[sl0a] -= Kf0dCfdC0a;
        J_sl0a_ptr[sl1a] -= Kf0dCfdC1a;
        J_sl1a_ptr[sl1a] -= Kf0dCfdC1a;
        J_sr0a_ptr[sl0a] += Kf0dCfdC0a;
        J_sr1a_ptr[sl0a] += Kf0dCfdC0a;
        J_sr0a_ptr[sl1a] += Kf0dCfdC1a;
        J_sr1a_ptr[sl1a] += Kf0dCfdC1a;

        double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a]; 
        J_sl0a_ptr[sr0a] -= (-Kr0*dCrdC0a);
        J_sl1a_ptr[sr0a] -= (-Kr0*dCrdC0a);
        J_sl0a_ptr[sr1a] -= (-Kr0*dCrdC1a);
        J_sl1a_ptr[sr1a] -= (-Kr0*dCrdC1a);
        J_sr0a_ptr[sr0a] += (-Kr0*dCrdC0a);
        J_sr1a_ptr[sr0a] += (-Kr0*dCrdC0a);
        J_sr0a_ptr[sr1a] += (-Kr0*dCrdC1a);
        J_sr1a_ptr[sr1a] += (-Kr0*dCrdC1a);

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a]*C[sl1a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        J_sl0a_ptr[this->nSpecies] -= dqdT0;
        J_sl1a_ptr[this->nSpecies] -= dqdT0;
        J_sr0a_ptr[this->nSpecies] += dqdT0;
        J_sr1a_ptr[this->nSpecies] += dqdT0;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        if(i0>=this->Ikf[2] && i0 <this->Ikf[3])
        { 
            const unsigned int idx = i0 - this->Ikf[2];
            double W = this->dKfdC_[idx]*CF0 - this->dKfdC_[idx+this->Itbr[4]]*CR0;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
            } 
        }

        const unsigned int i1 = this->NER22DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER22DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER22DupAllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->NER22DupAllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->NER22DupAllIndex[(k+1)*5+4];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sl1b_ptr = &ddNdtByVdcTp[sl1b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];

        double Kf1 = this->Kf_[i1];
        const double dCfdC0b = C[sl1b];
        const double dCfdC1b = C[sl0b];
        const double Kf1dCfdC0b = Kf1*dCfdC0b;
        const double Kf1dCfdC1b = Kf1*dCfdC1b;
        J_sl0b_ptr[sl0b] -= Kf1dCfdC0b;
        J_sl1b_ptr[sl0b] -= Kf1dCfdC0b;
        J_sl0b_ptr[sl1b] -= Kf1dCfdC1b;
        J_sl1b_ptr[sl1b] -= Kf1dCfdC1b;
        J_sr0b_ptr[sl0b] += Kf1dCfdC0b;
        J_sr1b_ptr[sl0b] += Kf1dCfdC0b;
        J_sr0b_ptr[sl1b] += Kf1dCfdC1b;
        J_sr1b_ptr[sl1b] += Kf1dCfdC1b;

        double Kr1 = this->Kf_[i1 - Ikf[1] + Ikf[9]];
        const double dCrdC0b = C[sr1b];
        const double dCrdC1b = C[sr0b]; 
        J_sl0b_ptr[sr0b] -= (-Kr1*dCrdC0b);
        J_sl1b_ptr[sr0b] -= (-Kr1*dCrdC0b);
        J_sl0b_ptr[sr1b] -= (-Kr1*dCrdC1b);
        J_sl1b_ptr[sr1b] -= (-Kr1*dCrdC1b);
        J_sr0b_ptr[sr0b] += (-Kr1*dCrdC0b);
        J_sr1b_ptr[sr0b] += (-Kr1*dCrdC0b);
        J_sr0b_ptr[sr1b] += (-Kr1*dCrdC1b);
        J_sr1b_ptr[sr1b] += (-Kr1*dCrdC1b);

        const double dKfdT1 = this->dKfdT_[i1];
        const double dKrdT1 = this->dKfdT_[i1 - Ikf[1] + Ikf[9]];
        const double CF1 = C[sl0b]*C[sl1b];
        const double CR1 = C[sr0b]*C[sr1b];
        const double dqdT1 = (dKfdT1*CF1)-(dKrdT1*CR1);
        J_sl0b_ptr[this->nSpecies] -= dqdT1;
        J_sl1b_ptr[this->nSpecies] -= dqdT1;
        J_sr0b_ptr[this->nSpecies] += dqdT1;
        J_sr1b_ptr[this->nSpecies] += dqdT1;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        if(i1>=this->Ikf[2] && i1 <this->Ikf[3])
        { 
            const unsigned int idx = i1 - this->Ikf[2];
            double W = this->dKfdC_[idx]*CF1 - this->dKfdC_[idx+this->Itbr[4]]*CR1;
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d WdMdC4 = _mm256_mul_pd(Wv,dMdC4);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+4]),WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1b_ptr[j+0]),WdMdC0);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1b_ptr[j+4]),WdMdC4);
                store256d(&J_sl1b_ptr[j+0],sl1v0);
                store256d(&J_sl1b_ptr[j+4],sl1v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+4]),WdMdC4);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sr0b_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+4]),WdMdC4);
                store256d(&J_sr1b_ptr[j+0],sr1v0);
                store256d(&J_sr1b_ptr[j+4],sr1v4);
            } 
        }

        const unsigned int i2 = this->NER22DupAllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER22DupAllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->NER22DupAllIndex[(k+2)*5+2];
        const unsigned int sr0c = this->NER22DupAllIndex[(k+2)*5+3];
        const unsigned int sr1c = this->NER22DupAllIndex[(k+2)*5+4];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sl1c_ptr = &ddNdtByVdcTp[sl1c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];

        double Kf2 = this->Kf_[i2];
        const double dCfdC0c = C[sl1c];
        const double dCfdC1c = C[sl0c];
        const double Kf2dCfdC0c = Kf2*dCfdC0c;
        const double Kf2dCfdC1c = Kf2*dCfdC1c;
        J_sl0c_ptr[sl0c] -= Kf2dCfdC0c;
        J_sl1c_ptr[sl0c] -= Kf2dCfdC0c;
        J_sl0c_ptr[sl1c] -= Kf2dCfdC1c;
        J_sl1c_ptr[sl1c] -= Kf2dCfdC1c;
        J_sr0c_ptr[sl0c] += Kf2dCfdC0c;
        J_sr1c_ptr[sl0c] += Kf2dCfdC0c;
        J_sr0c_ptr[sl1c] += Kf2dCfdC1c;
        J_sr1c_ptr[sl1c] += Kf2dCfdC1c;

        double Kr2 = this->Kf_[i2 - Ikf[1] + Ikf[9]];
        const double dCrdC0c = C[sr1c];
        const double dCrdC1c = C[sr0c]; 
        J_sl0c_ptr[sr0c] -= (-Kr2*dCrdC0c);
        J_sl1c_ptr[sr0c] -= (-Kr2*dCrdC0c);
        J_sl0c_ptr[sr1c] -= (-Kr2*dCrdC1c);
        J_sl1c_ptr[sr1c] -= (-Kr2*dCrdC1c);
        J_sr0c_ptr[sr0c] += (-Kr2*dCrdC0c);
        J_sr1c_ptr[sr0c] += (-Kr2*dCrdC0c);
        J_sr0c_ptr[sr1c] += (-Kr2*dCrdC1c);
        J_sr1c_ptr[sr1c] += (-Kr2*dCrdC1c);

        const double dKfdT2 = this->dKfdT_[i2];
        const double dKrdT2 = this->dKfdT_[i2 - Ikf[1] + Ikf[9]];
        const double CF2 = C[sl0c]*C[sl1c];
        const double CR2 = C[sr0c]*C[sr1c];
        const double dqdT2 = (dKfdT2*CF2)-(dKrdT2*CR2);
        J_sl0c_ptr[this->nSpecies] -= dqdT2;
        J_sl1c_ptr[this->nSpecies] -= dqdT2;
        J_sr0c_ptr[this->nSpecies] += dqdT2;
        J_sr1c_ptr[this->nSpecies] += dqdT2;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        if(i2>=this->Ikf[2] && i2 <this->Ikf[3])
        { 
            const unsigned int idx = i2 - this->Ikf[2];
            double W = this->dKfdC_[idx]*CF2 - this->dKfdC_[idx+this->Itbr[4]]*CR2;
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d WdMdC4 = _mm256_mul_pd(Wv,dMdC4);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+4]),WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1c_ptr[j+0]),WdMdC0);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1c_ptr[j+4]),WdMdC4);
                store256d(&J_sl1c_ptr[j+0],sl1v0);
                store256d(&J_sl1c_ptr[j+4],sl1v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+4]),WdMdC4);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sr0c_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+4]),WdMdC4);
                store256d(&J_sr1c_ptr[j+0],sr1v0);
                store256d(&J_sr1c_ptr[j+4],sr1v4);
            } 
        }

        const unsigned int i3 = this->NER22DupAllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->NER22DupAllIndex[(k+3)*5+1];
        const unsigned int sl1d = this->NER22DupAllIndex[(k+3)*5+2];
        const unsigned int sr0d = this->NER22DupAllIndex[(k+3)*5+3];
        const unsigned int sr1d = this->NER22DupAllIndex[(k+3)*5+4];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sl1d_ptr = &ddNdtByVdcTp[sl1d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];

        double Kf3 = this->Kf_[i3];
        const double dCfdC0d = C[sl1d];
        const double dCfdC1d = C[sl0d];
        const double Kf3dCfdC0d = Kf3*dCfdC0d;
        const double Kf3dCfdC1d = Kf3*dCfdC1d;
        J_sl0d_ptr[sl0d] -= Kf3dCfdC0d;
        J_sl1d_ptr[sl0d] -= Kf3dCfdC0d;
        J_sl0d_ptr[sl1d] -= Kf3dCfdC1d;
        J_sl1d_ptr[sl1d] -= Kf3dCfdC1d;
        J_sr0d_ptr[sl0d] += Kf3dCfdC0d;
        J_sr1d_ptr[sl0d] += Kf3dCfdC0d;
        J_sr0d_ptr[sl1d] += Kf3dCfdC1d;
        J_sr1d_ptr[sl1d] += Kf3dCfdC1d;

        double Kr3 = this->Kf_[i3 - Ikf[1] + Ikf[9]];
        const double dCrdC0d = C[sr1d];
        const double dCrdC1d = C[sr0d]; 
        J_sl0d_ptr[sr0d] -= (-Kr3*dCrdC0d);
        J_sl1d_ptr[sr0d] -= (-Kr3*dCrdC0d);
        J_sl0d_ptr[sr1d] -= (-Kr3*dCrdC1d);
        J_sl1d_ptr[sr1d] -= (-Kr3*dCrdC1d);
        J_sr0d_ptr[sr0d] += (-Kr3*dCrdC0d);
        J_sr1d_ptr[sr0d] += (-Kr3*dCrdC0d);
        J_sr0d_ptr[sr1d] += (-Kr3*dCrdC1d);
        J_sr1d_ptr[sr1d] += (-Kr3*dCrdC1d);

        const double dKfdT3 = this->dKfdT_[i3];
        const double dKrdT3 = this->dKfdT_[i3 - Ikf[1] + Ikf[9]];
        const double CF3 = C[sl0d]*C[sl1d];
        const double CR3 = C[sr0d]*C[sr1d];
        const double dqdT3 = (dKfdT3*CF3)-(dKrdT3*CR3);
        J_sl0d_ptr[this->nSpecies] -= dqdT3;
        J_sl1d_ptr[this->nSpecies] -= dqdT3;
        J_sr0d_ptr[this->nSpecies] += dqdT3;
        J_sr1d_ptr[this->nSpecies] += dqdT3;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

        if(i3>=this->Ikf[2] && i3 <this->Ikf[3])
        { 
            const unsigned int idx = i3 - this->Ikf[2];
            double W = this->dKfdC_[idx]*CF3 - this->dKfdC_[idx+this->Itbr[4]]*CR3;
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(Wv,dMdC0);
                __m256d WdMdC4 = _mm256_mul_pd(Wv,dMdC4);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+4]),WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
                __m256d sl1v0 = _mm256_sub_pd(load256d(&J_sl1d_ptr[j+0]),WdMdC0);
                __m256d sl1v4 = _mm256_sub_pd(load256d(&J_sl1d_ptr[j+4]),WdMdC4);
                store256d(&J_sl1d_ptr[j+0],sl1v0);
                store256d(&J_sl1d_ptr[j+4],sl1v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+4]),WdMdC4);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sr0d_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+4]),WdMdC4);
                store256d(&J_sr1d_ptr[j+0],sr1v0);
                store256d(&J_sr1d_ptr[j+4],sr1v4);
            } 
        }
    }
    for(std::size_t k=endNER22Dup-remainNER22Dup; k<endNER22Dup; k=k+1)
    {

        const unsigned int i0 = this->NER22DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER22DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER22DupAllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->NER22DupAllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->NER22DupAllIndex[(k+0)*5+4];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C[sl1a];
        const double dCfdC1a = C[sl0a];
        const double Kf0dCfdC0a = Kf0*dCfdC0a;
        const double Kf0dCfdC1a = Kf0*dCfdC1a;
        J_sl0a_ptr[sl0a] -= Kf0dCfdC0a;
        J_sl1a_ptr[sl0a] -= Kf0dCfdC0a;
        J_sl0a_ptr[sl1a] -= Kf0dCfdC1a;
        J_sl1a_ptr[sl1a] -= Kf0dCfdC1a;
        J_sr0a_ptr[sl0a] += Kf0dCfdC0a;
        J_sr1a_ptr[sl0a] += Kf0dCfdC0a;
        J_sr0a_ptr[sl1a] += Kf0dCfdC1a;
        J_sr1a_ptr[sl1a] += Kf0dCfdC1a;

        double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a]; 
        J_sl0a_ptr[sr0a] -= (-Kr0*dCrdC0a);
        J_sl1a_ptr[sr0a] -= (-Kr0*dCrdC0a);
        J_sl0a_ptr[sr1a] -= (-Kr0*dCrdC1a);
        J_sl1a_ptr[sr1a] -= (-Kr0*dCrdC1a);
        J_sr0a_ptr[sr0a] += (-Kr0*dCrdC0a);
        J_sr1a_ptr[sr0a] += (-Kr0*dCrdC0a);
        J_sr0a_ptr[sr1a] += (-Kr0*dCrdC1a);
        J_sr1a_ptr[sr1a] += (-Kr0*dCrdC1a);

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a]*C[sl1a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        J_sl0a_ptr[this->nSpecies] -= dqdT0;
        J_sl1a_ptr[this->nSpecies] -= dqdT0;
        J_sr0a_ptr[this->nSpecies] += dqdT0;
        J_sr1a_ptr[this->nSpecies] += dqdT0;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        if(i0>=this->Ikf[2] && i0 <this->Ikf[3])
        { 
            const unsigned int idx = i0 - this->Ikf[2];
            double W = this->dKfdC_[idx]*CF0 - this->dKfdC_[idx+this->Itbr[4]]*CR0;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
            } 
        }
    }
}   
