/*---------------------------------------------------------------------------*\
  Description
      Computing the molar concentration based jacobian matrix. The function 
      is used for three-two reaction, e.g. A+A+A=B+B. A+B+C=D+E 
      
      RR:  reversible reaction
      IR:  irreversible reaction reverse rate constant is zero
      NER: non-equilibrium reaction, reverse rate constant is computed using
           Arrhenius form instead of equilibrium rate constant
  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"
#include <iostream>
//=============================================================================//

void  FastChemistry::OptReaction::JF32RR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{
    const double PowPByRTSumVki = this->Pow_pByRT_SumVki[3];
    const double invKcLimiter = FastChemistry::invKcLimiter;
    std::size_t endRR32 = this->RR32Size;
    std::size_t remainRR32 = endRR32%4;
    for(std::size_t k=0; k<endRR32-remainRR32; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->RR32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->RR32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->RR32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->RR32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->RR32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->RR32AllIndex[(k+0)*6+5];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sl2a_ptr = &ddNdtByVdcTp[sl2a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double C_sl2a = C[sl2a];
        const double Kf0dCfdC0a = Kf0*C_sl1a*C_sl2a;
        const double Kf0dCfdC1a = Kf0*C_sl0a*C_sl2a;
        const double Kf0dCfdC2a = Kf0*C_sl0a*C_sl1a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl2a_sl0a_increment = J_sl2a_ptr[sl0a] - Kf0dCfdC0a;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        J_sl2a_ptr[sl0a] = J_sl2a_sl0a_increment;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl2a_sl1a_increment = J_sl2a_ptr[sl1a] - Kf0dCfdC1a;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sl2a_ptr[sl1a] = J_sl2a_sl1a_increment;
        const double J_sl0a_sl2a_increment = J_sl0a_ptr[sl2a] - Kf0dCfdC2a;
        const double J_sl1a_sl2a_increment = J_sl1a_ptr[sl2a] - Kf0dCfdC2a;
        const double J_sl2a_sl2a_increment = J_sl2a_ptr[sl2a] - Kf0dCfdC2a;
        J_sl0a_ptr[sl2a] = J_sl0a_sl2a_increment;
        J_sl1a_ptr[sl2a] = J_sl1a_sl2a_increment;
        J_sl2a_ptr[sl2a] = J_sl2a_sl2a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0dCfdC0a;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1a;
        const double J_sr1a_sl1a_increment = J_sr1a_ptr[sl1a] + Kf0dCfdC1a;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;
        J_sr1a_ptr[sl1a] = J_sr1a_sl1a_increment;
        const double J_sr0a_sl2a_increment = J_sr0a_ptr[sl2a] + Kf0dCfdC2a;
        const double J_sr1a_sl2a_increment = J_sr1a_ptr[sl2a] + Kf0dCfdC2a;
        J_sr0a_ptr[sl2a] = J_sr0a_sl2a_increment;
        J_sr1a_ptr[sl2a] = J_sr1a_sl2a_increment;

        double invKc0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
            (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);

        const double Kr0 = Kf0*invKc0;
        const double C_sr1a = C[sr1a];
        const double C_sr0a = C[sr0a];
        const double Kr0dCrdC0a = -Kr0*C_sr1a;
        const double Kr0dCrdC1a = -Kr0*C_sr0a;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sl1a_sr0a_increment = J_sl1a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sl2a_sr0a_increment = J_sl2a_ptr[sr0a] - Kr0dCrdC0a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sl1a_ptr[sr0a] = J_sl1a_sr0a_increment;
        J_sl2a_ptr[sr0a] = J_sl2a_sr0a_increment;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] - Kr0dCrdC1a;
        const double J_sl1a_sr1a_increment = J_sl1a_ptr[sr1a] - Kr0dCrdC1a;
        const double J_sl2a_sr1a_increment = J_sl2a_ptr[sr1a] - Kr0dCrdC1a;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sl1a_ptr[sr1a] = J_sl1a_sr1a_increment;
        J_sl2a_ptr[sr1a] = J_sl2a_sr1a_increment;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] + Kr0dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] + Kr0dCrdC0a;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] + Kr0dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] + Kr0dCrdC1a;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a]) - (dBdT[sl0a] + dBdT[sl1a] + dBdT[sl2a]);
        const double dKcdTByKc0 = sumVdBdT0 + invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C[sl0a]*C[sl1a]*C[sl2a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sl2a_T_increment = J_sl2a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sl2a_ptr[this->nSpecies] = J_sl2a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*(CF0-invKc0*CR0);
            __m256d Wv0 = _mm256_set1_pd(W0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv0,dMdC);
                __m256d sl0v = load256d(&J_sl0a_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1a_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1a_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2a_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2a_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0a_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1a_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->RR32AllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->RR32AllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->RR32AllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->RR32AllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->RR32AllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->RR32AllIndex[(k+1)*6+5];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sl1b_ptr = &ddNdtByVdcTp[sl1b*(this->alignN)];
        double* __restrict__ J_sl2b_ptr = &ddNdtByVdcTp[sl2b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double C_sl1b = C[sl1b];
        const double C_sl0b = C[sl0b];
        const double C_sl2b = C[sl2b];
        const double Kf1dCfdC0b = Kf1*C_sl1b*C_sl2b;
        const double Kf1dCfdC1b = Kf1*C_sl0b*C_sl2b;
        const double Kf1dCfdC2b = Kf1*C_sl0b*C_sl1b;
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1dCfdC0b;
        const double J_sl1b_sl0b_increment = J_sl1b_ptr[sl0b] - Kf1dCfdC0b;
        const double J_sl2b_sl0b_increment = J_sl2b_ptr[sl0b] - Kf1dCfdC0b;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sl1b_ptr[sl0b] = J_sl1b_sl0b_increment;
        J_sl2b_ptr[sl0b] = J_sl2b_sl0b_increment;
        const double J_sl0b_sl1b_increment = J_sl0b_ptr[sl1b] - Kf1dCfdC1b;
        const double J_sl1b_sl1b_increment = J_sl1b_ptr[sl1b] - Kf1dCfdC1b;
        const double J_sl2b_sl1b_increment = J_sl2b_ptr[sl1b] - Kf1dCfdC1b;
        J_sl0b_ptr[sl1b] = J_sl0b_sl1b_increment;
        J_sl1b_ptr[sl1b] = J_sl1b_sl1b_increment;
        J_sl2b_ptr[sl1b] = J_sl2b_sl1b_increment;
        const double J_sl0b_sl2b_increment = J_sl0b_ptr[sl2b] - Kf1dCfdC2b;
        const double J_sl1b_sl2b_increment = J_sl1b_ptr[sl2b] - Kf1dCfdC2b;
        const double J_sl2b_sl2b_increment = J_sl2b_ptr[sl2b] - Kf1dCfdC2b;
        J_sl0b_ptr[sl2b] = J_sl0b_sl2b_increment;
        J_sl1b_ptr[sl2b] = J_sl1b_sl2b_increment;
        J_sl2b_ptr[sl2b] = J_sl2b_sl2b_increment;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1dCfdC0b;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1dCfdC0b;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;
        const double J_sr0b_sl1b_increment = J_sr0b_ptr[sl1b] + Kf1dCfdC1b;
        const double J_sr1b_sl1b_increment = J_sr1b_ptr[sl1b] + Kf1dCfdC1b;
        J_sr0b_ptr[sl1b] = J_sr0b_sl1b_increment;
        J_sr1b_ptr[sl1b] = J_sr1b_sl1b_increment;
        const double J_sr0b_sl2b_increment = J_sr0b_ptr[sl2b] + Kf1dCfdC2b;
        const double J_sr1b_sl2b_increment = J_sr1b_ptr[sl2b] + Kf1dCfdC2b;
        J_sr0b_ptr[sl2b] = J_sr0b_sl2b_increment;
        J_sr1b_ptr[sl2b] = J_sr1b_sl2b_increment;

        double invKc1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*
            (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);

        const double Kr1 = Kf1*invKc1;
        const double C_sr1b = C[sr1b];
        const double C_sr0b = C[sr0b];
        const double Kr1dCrdC0b = -Kr1*C_sr1b;
        const double Kr1dCrdC1b = -Kr1*C_sr0b;
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] - Kr1dCrdC0b;
        const double J_sl1b_sr0b_increment = J_sl1b_ptr[sr0b] - Kr1dCrdC0b;
        const double J_sl2b_sr0b_increment = J_sl2b_ptr[sr0b] - Kr1dCrdC0b;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sl1b_ptr[sr0b] = J_sl1b_sr0b_increment;
        J_sl2b_ptr[sr0b] = J_sl2b_sr0b_increment;
        const double J_sl0b_sr1b_increment = J_sl0b_ptr[sr1b] - Kr1dCrdC1b;
        const double J_sl1b_sr1b_increment = J_sl1b_ptr[sr1b] - Kr1dCrdC1b;
        const double J_sl2b_sr1b_increment = J_sl2b_ptr[sr1b] - Kr1dCrdC1b;
        J_sl0b_ptr[sr1b] = J_sl0b_sr1b_increment;
        J_sl1b_ptr[sr1b] = J_sl1b_sr1b_increment;
        J_sl2b_ptr[sr1b] = J_sl2b_sr1b_increment;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] + Kr1dCrdC0b;
        const double J_sr1b_sr0b_increment = J_sr1b_ptr[sr0b] + Kr1dCrdC0b;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;
        J_sr1b_ptr[sr0b] = J_sr1b_sr0b_increment;
        const double J_sr0b_sr1b_increment = J_sr0b_ptr[sr1b] + Kr1dCrdC1b;
        const double J_sr1b_sr1b_increment = J_sr1b_ptr[sr1b] + Kr1dCrdC1b;
        J_sr0b_ptr[sr1b] = J_sr0b_sr1b_increment;
        J_sr1b_ptr[sr1b] = J_sr1b_sr1b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double sumVdBdT1 = (dBdT[sr0b] + dBdT[sr1b]) - (dBdT[sl0b] + dBdT[sl1b] + dBdT[sl2b]);
        const double dKcdTByKc1 = sumVdBdT1 + invT;
        const double dKrdT1 = (dKfdT1*invKc1 - (invKc1 < invKcLimiter ? Kr1*dKcdTByKc1 : 0));
        const double CF1 = C[sl0b]*C[sl1b]*C[sl2b];
        const double CR1 = C[sr0b]*C[sr1b];
        const double dqdT1 = (dKfdT1*CF1) - (dKrdT1*CR1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sl1b_T_increment = J_sl1b_ptr[this->nSpecies] - dqdT1;
        const double J_sl2b_T_increment = J_sl2b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sl1b_ptr[this->nSpecies] = J_sl1b_T_increment;
        J_sl2b_ptr[this->nSpecies] = J_sl2b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*(CF1-invKc1*CR1);
            __m256d Wv1 = _mm256_set1_pd(W1);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv1,dMdC);
                __m256d sl0v = load256d(&J_sl0b_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1b_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1b_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2b_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2b_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0b_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0b_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1b_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1b_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->RR32AllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->RR32AllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->RR32AllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->RR32AllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->RR32AllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->RR32AllIndex[(k+2)*6+5];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sl1c_ptr = &ddNdtByVdcTp[sl1c*(this->alignN)];
        double* __restrict__ J_sl2c_ptr = &ddNdtByVdcTp[sl2c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double C_sl1c = C[sl1c];
        const double C_sl0c = C[sl0c];
        const double C_sl2c = C[sl2c];
        const double Kf2dCfdC0c = Kf2*C_sl1c*C_sl2c;
        const double Kf2dCfdC1c = Kf2*C_sl0c*C_sl2c;
        const double Kf2dCfdC2c = Kf2*C_sl0c*C_sl1c;
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2dCfdC0c;
        const double J_sl1c_sl0c_increment = J_sl1c_ptr[sl0c] - Kf2dCfdC0c;
        const double J_sl2c_sl0c_increment = J_sl2c_ptr[sl0c] - Kf2dCfdC0c;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sl1c_ptr[sl0c] = J_sl1c_sl0c_increment;
        J_sl2c_ptr[sl0c] = J_sl2c_sl0c_increment;
        const double J_sl0c_sl1c_increment = J_sl0c_ptr[sl1c] - Kf2dCfdC1c;
        const double J_sl1c_sl1c_increment = J_sl1c_ptr[sl1c] - Kf2dCfdC1c;
        const double J_sl2c_sl1c_increment = J_sl2c_ptr[sl1c] - Kf2dCfdC1c;
        J_sl0c_ptr[sl1c] = J_sl0c_sl1c_increment;
        J_sl1c_ptr[sl1c] = J_sl1c_sl1c_increment;
        J_sl2c_ptr[sl1c] = J_sl2c_sl1c_increment;
        const double J_sl0c_sl2c_increment = J_sl0c_ptr[sl2c] - Kf2dCfdC2c;
        const double J_sl1c_sl2c_increment = J_sl1c_ptr[sl2c] - Kf2dCfdC2c;
        const double J_sl2c_sl2c_increment = J_sl2c_ptr[sl2c] - Kf2dCfdC2c;
        J_sl0c_ptr[sl2c] = J_sl0c_sl2c_increment;
        J_sl1c_ptr[sl2c] = J_sl1c_sl2c_increment;
        J_sl2c_ptr[sl2c] = J_sl2c_sl2c_increment;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2dCfdC0c;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2dCfdC0c;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;
        const double J_sr0c_sl1c_increment = J_sr0c_ptr[sl1c] + Kf2dCfdC1c;
        const double J_sr1c_sl1c_increment = J_sr1c_ptr[sl1c] + Kf2dCfdC1c;
        J_sr0c_ptr[sl1c] = J_sr0c_sl1c_increment;
        J_sr1c_ptr[sl1c] = J_sr1c_sl1c_increment;
        const double J_sr0c_sl2c_increment = J_sr0c_ptr[sl2c] + Kf2dCfdC2c;
        const double J_sr1c_sl2c_increment = J_sr1c_ptr[sl2c] + Kf2dCfdC2c;
        J_sr0c_ptr[sl2c] = J_sr0c_sl2c_increment;
        J_sr1c_ptr[sl2c] = J_sr1c_sl2c_increment;

        double invKc2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*
            (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);

        const double Kr2 = Kf2*invKc2;
        const double C_sr1c = C[sr1c];
        const double C_sr0c = C[sr0c];
        const double Kr2dCrdC0c = -Kr2*C_sr1c;
        const double Kr2dCrdC1c = -Kr2*C_sr0c;
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] - Kr2dCrdC0c;
        const double J_sl1c_sr0c_increment = J_sl1c_ptr[sr0c] - Kr2dCrdC0c;
        const double J_sl2c_sr0c_increment = J_sl2c_ptr[sr0c] - Kr2dCrdC0c;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sl1c_ptr[sr0c] = J_sl1c_sr0c_increment;
        J_sl2c_ptr[sr0c] = J_sl2c_sr0c_increment;
        const double J_sl0c_sr1c_increment = J_sl0c_ptr[sr1c] - Kr2dCrdC1c;
        const double J_sl1c_sr1c_increment = J_sl1c_ptr[sr1c] - Kr2dCrdC1c;
        const double J_sl2c_sr1c_increment = J_sl2c_ptr[sr1c] - Kr2dCrdC1c;
        J_sl0c_ptr[sr1c] = J_sl0c_sr1c_increment;
        J_sl1c_ptr[sr1c] = J_sl1c_sr1c_increment;
        J_sl2c_ptr[sr1c] = J_sl2c_sr1c_increment;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] + Kr2dCrdC0c;
        const double J_sr1c_sr0c_increment = J_sr1c_ptr[sr0c] + Kr2dCrdC0c;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;
        J_sr1c_ptr[sr0c] = J_sr1c_sr0c_increment;
        const double J_sr0c_sr1c_increment = J_sr0c_ptr[sr1c] + Kr2dCrdC1c;
        const double J_sr1c_sr1c_increment = J_sr1c_ptr[sr1c] + Kr2dCrdC1c;
        J_sr0c_ptr[sr1c] = J_sr0c_sr1c_increment;
        J_sr1c_ptr[sr1c] = J_sr1c_sr1c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double sumVdBdT2 = (dBdT[sr0c] + dBdT[sr1c]) - (dBdT[sl0c] + dBdT[sl1c] + dBdT[sl2c]);
        const double dKcdTByKc2 = sumVdBdT2 + invT;
        const double dKrdT2 = (dKfdT2*invKc2 - (invKc2 < invKcLimiter ? Kr2*dKcdTByKc2 : 0));
        const double CF2 = C[sl0c]*C[sl1c]*C[sl2c];
        const double CR2 = C[sr0c]*C[sr1c];
        const double dqdT2 = (dKfdT2*CF2) - (dKrdT2*CR2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sl1c_T_increment = J_sl1c_ptr[this->nSpecies] - dqdT2;
        const double J_sl2c_T_increment = J_sl2c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sl1c_ptr[this->nSpecies] = J_sl1c_T_increment;
        J_sl2c_ptr[this->nSpecies] = J_sl2c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*(CF2-invKc2*CR2);
            __m256d Wv2 = _mm256_set1_pd(W2);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv2,dMdC);
                __m256d sl0v = load256d(&J_sl0c_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1c_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1c_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2c_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2c_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0c_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0c_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1c_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1c_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->RR32AllIndex[(k+3)*6+0];
        const unsigned int sl0d = this->RR32AllIndex[(k+3)*6+1];
        const unsigned int sl1d = this->RR32AllIndex[(k+3)*6+2];
        const unsigned int sl2d = this->RR32AllIndex[(k+3)*6+3];
        const unsigned int sr0d = this->RR32AllIndex[(k+3)*6+4];
        const unsigned int sr1d = this->RR32AllIndex[(k+3)*6+5];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sl1d_ptr = &ddNdtByVdcTp[sl1d*(this->alignN)];
        double* __restrict__ J_sl2d_ptr = &ddNdtByVdcTp[sl2d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double C_sl1d = C[sl1d];
        const double C_sl0d = C[sl0d];
        const double C_sl2d = C[sl2d];
        const double Kf3dCfdC0d = Kf3*C_sl1d*C_sl2d;
        const double Kf3dCfdC1d = Kf3*C_sl0d*C_sl2d;
        const double Kf3dCfdC2d = Kf3*C_sl0d*C_sl1d;
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3dCfdC0d;
        const double J_sl1d_sl0d_increment = J_sl1d_ptr[sl0d] - Kf3dCfdC0d;
        const double J_sl2d_sl0d_increment = J_sl2d_ptr[sl0d] - Kf3dCfdC0d;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sl1d_ptr[sl0d] = J_sl1d_sl0d_increment;
        J_sl2d_ptr[sl0d] = J_sl2d_sl0d_increment;
        const double J_sl0d_sl1d_increment = J_sl0d_ptr[sl1d] - Kf3dCfdC1d;
        const double J_sl1d_sl1d_increment = J_sl1d_ptr[sl1d] - Kf3dCfdC1d;
        const double J_sl2d_sl1d_increment = J_sl2d_ptr[sl1d] - Kf3dCfdC1d;
        J_sl0d_ptr[sl1d] = J_sl0d_sl1d_increment;
        J_sl1d_ptr[sl1d] = J_sl1d_sl1d_increment;
        J_sl2d_ptr[sl1d] = J_sl2d_sl1d_increment;
        const double J_sl0d_sl2d_increment = J_sl0d_ptr[sl2d] - Kf3dCfdC2d;
        const double J_sl1d_sl2d_increment = J_sl1d_ptr[sl2d] - Kf3dCfdC2d;
        const double J_sl2d_sl2d_increment = J_sl2d_ptr[sl2d] - Kf3dCfdC2d;
        J_sl0d_ptr[sl2d] = J_sl0d_sl2d_increment;
        J_sl1d_ptr[sl2d] = J_sl1d_sl2d_increment;
        J_sl2d_ptr[sl2d] = J_sl2d_sl2d_increment;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3dCfdC0d;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3dCfdC0d;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;
        const double J_sr0d_sl1d_increment = J_sr0d_ptr[sl1d] + Kf3dCfdC1d;
        const double J_sr1d_sl1d_increment = J_sr1d_ptr[sl1d] + Kf3dCfdC1d;
        J_sr0d_ptr[sl1d] = J_sr0d_sl1d_increment;
        J_sr1d_ptr[sl1d] = J_sr1d_sl1d_increment;
        const double J_sr0d_sl2d_increment = J_sr0d_ptr[sl2d] + Kf3dCfdC2d;
        const double J_sr1d_sl2d_increment = J_sr1d_ptr[sl2d] + Kf3dCfdC2d;
        J_sr0d_ptr[sl2d] = J_sr0d_sl2d_increment;
        J_sr1d_ptr[sl2d] = J_sr1d_sl2d_increment;

        double invKc3 = (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d])*
            (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]*ExpNegGbyRT[sl2d])*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);

        const double Kr3 = Kf3*invKc3;
        const double C_sr1d = C[sr1d];
        const double C_sr0d = C[sr0d];
        const double Kr3dCrdC0d = -Kr3*C_sr1d;
        const double Kr3dCrdC1d = -Kr3*C_sr0d;
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] - Kr3dCrdC0d;
        const double J_sl1d_sr0d_increment = J_sl1d_ptr[sr0d] - Kr3dCrdC0d;
        const double J_sl2d_sr0d_increment = J_sl2d_ptr[sr0d] - Kr3dCrdC0d;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sl1d_ptr[sr0d] = J_sl1d_sr0d_increment;
        J_sl2d_ptr[sr0d] = J_sl2d_sr0d_increment;
        const double J_sl0d_sr1d_increment = J_sl0d_ptr[sr1d] - Kr3dCrdC1d;
        const double J_sl1d_sr1d_increment = J_sl1d_ptr[sr1d] - Kr3dCrdC1d;
        const double J_sl2d_sr1d_increment = J_sl2d_ptr[sr1d] - Kr3dCrdC1d;
        J_sl0d_ptr[sr1d] = J_sl0d_sr1d_increment;
        J_sl1d_ptr[sr1d] = J_sl1d_sr1d_increment;
        J_sl2d_ptr[sr1d] = J_sl2d_sr1d_increment;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] + Kr3dCrdC0d;
        const double J_sr1d_sr0d_increment = J_sr1d_ptr[sr0d] + Kr3dCrdC0d;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;
        J_sr1d_ptr[sr0d] = J_sr1d_sr0d_increment;
        const double J_sr0d_sr1d_increment = J_sr0d_ptr[sr1d] + Kr3dCrdC1d;
        const double J_sr1d_sr1d_increment = J_sr1d_ptr[sr1d] + Kr3dCrdC1d;
        J_sr0d_ptr[sr1d] = J_sr0d_sr1d_increment;
        J_sr1d_ptr[sr1d] = J_sr1d_sr1d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double sumVdBdT3 = (dBdT[sr0d] + dBdT[sr1d]) - (dBdT[sl0d] + dBdT[sl1d] + dBdT[sl2d]);
        const double dKcdTByKc3 = sumVdBdT3 + invT;
        const double dKrdT3 = (dKfdT3*invKc3 - (invKc3 < invKcLimiter ? Kr3*dKcdTByKc3 : 0));
        const double CF3 = C[sl0d]*C[sl1d]*C[sl2d];
        const double CR3 = C[sr0d]*C[sr1d];
        const double dqdT3 = (dKfdT3*CF3) - (dKrdT3*CR3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sl1d_T_increment = J_sl1d_ptr[this->nSpecies] - dqdT3;
        const double J_sl2d_T_increment = J_sl2d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sl1d_ptr[this->nSpecies] = J_sl1d_T_increment;
        J_sl2d_ptr[this->nSpecies] = J_sl2d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sl2d_increment = dNdtByV[sl2d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sl2d] = omega_sl2d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*(CF3-invKc3*CR3);
            __m256d Wv3 = _mm256_set1_pd(W3);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv3,dMdC);
                __m256d sl0v = load256d(&J_sl0d_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1d_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1d_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2d_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2d_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0d_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0d_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1d_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1d_ptr[j],sr1v);
            }
        }
    }
    for(std::size_t k=endRR32-remainRR32; k<endRR32; k=k+1)
    {
        const unsigned int i0 = this->RR32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->RR32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->RR32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->RR32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->RR32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->RR32AllIndex[(k+0)*6+5];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sl2a_ptr = &ddNdtByVdcTp[sl2a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double C_sl2a = C[sl2a];
        const double Kf0dCfdC0a = Kf0*C_sl1a*C_sl2a;
        const double Kf0dCfdC1a = Kf0*C_sl0a*C_sl2a;
        const double Kf0dCfdC2a = Kf0*C_sl0a*C_sl1a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl2a_sl0a_increment = J_sl2a_ptr[sl0a] - Kf0dCfdC0a;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        J_sl2a_ptr[sl0a] = J_sl2a_sl0a_increment;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl2a_sl1a_increment = J_sl2a_ptr[sl1a] - Kf0dCfdC1a;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sl2a_ptr[sl1a] = J_sl2a_sl1a_increment;
        const double J_sl0a_sl2a_increment = J_sl0a_ptr[sl2a] - Kf0dCfdC2a;
        const double J_sl1a_sl2a_increment = J_sl1a_ptr[sl2a] - Kf0dCfdC2a;
        const double J_sl2a_sl2a_increment = J_sl2a_ptr[sl2a] - Kf0dCfdC2a;
        J_sl0a_ptr[sl2a] = J_sl0a_sl2a_increment;
        J_sl1a_ptr[sl2a] = J_sl1a_sl2a_increment;
        J_sl2a_ptr[sl2a] = J_sl2a_sl2a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0dCfdC0a;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1a;
        const double J_sr1a_sl1a_increment = J_sr1a_ptr[sl1a] + Kf0dCfdC1a;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;
        J_sr1a_ptr[sl1a] = J_sr1a_sl1a_increment;
        const double J_sr0a_sl2a_increment = J_sr0a_ptr[sl2a] + Kf0dCfdC2a;
        const double J_sr1a_sl2a_increment = J_sr1a_ptr[sl2a] + Kf0dCfdC2a;
        J_sr0a_ptr[sl2a] = J_sr0a_sl2a_increment;
        J_sr1a_ptr[sl2a] = J_sr1a_sl2a_increment;

        double invKc0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
            (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);

        const double Kr0 = Kf0*invKc0;
        const double C_sr1a = C[sr1a];
        const double C_sr0a = C[sr0a];
        const double Kr0dCrdC0a = -Kr0*C_sr1a;
        const double Kr0dCrdC1a = -Kr0*C_sr0a;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sl1a_sr0a_increment = J_sl1a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sl2a_sr0a_increment = J_sl2a_ptr[sr0a] - Kr0dCrdC0a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sl1a_ptr[sr0a] = J_sl1a_sr0a_increment;
        J_sl2a_ptr[sr0a] = J_sl2a_sr0a_increment;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] - Kr0dCrdC1a;
        const double J_sl1a_sr1a_increment = J_sl1a_ptr[sr1a] - Kr0dCrdC1a;
        const double J_sl2a_sr1a_increment = J_sl2a_ptr[sr1a] - Kr0dCrdC1a;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sl1a_ptr[sr1a] = J_sl1a_sr1a_increment;
        J_sl2a_ptr[sr1a] = J_sl2a_sr1a_increment;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] + Kr0dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] + Kr0dCrdC0a;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] + Kr0dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] + Kr0dCrdC1a;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a]) - (dBdT[sl0a] + dBdT[sl1a] + dBdT[sl2a]);
        const double dKcdTByKc0 = sumVdBdT0 + invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C[sl0a]*C[sl1a]*C[sl2a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sl2a_T_increment = J_sl2a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sl2a_ptr[this->nSpecies] = J_sl2a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*(CF0-invKc0*CR0);
            __m256d Wv0 = _mm256_set1_pd(W0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv0,dMdC);
                __m256d sl0v = load256d(&J_sl0a_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1a_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1a_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2a_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2a_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0a_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1a_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }
    }

    std::size_t endRR32Dup = this->RR32DupSize;
    std::size_t remainRR32Dup = endRR32Dup%4;
    for(std::size_t k=0; k<endRR32Dup-remainRR32Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->RR32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->RR32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->RR32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->RR32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->RR32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->RR32DupAllIndex[(k+0)*6+5];

        const double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C[sl1a]*C[sl2a];
        const double dCfdC1a = C[sl0a]*C[sl2a];
        const double dCfdC2a = C[sl0a]*C[sl1a];
        ddNdtByVdcTp[sl0a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl0a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl2a] += Kf0*dCfdC2a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl2a] += Kf0*dCfdC2a;

        double invKc0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
            (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a];
        ddNdtByVdcTp[sl0a*(this->alignN)+sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl1a*(this->alignN)+sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl2a*(this->alignN)+sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl0a*(this->alignN)+sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sl1a*(this->alignN)+sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sl2a*(this->alignN)+sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sr0a*(this->alignN)+sr0a] += (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sr1a*(this->alignN)+sr0a] += (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sr0a*(this->alignN)+sr1a] += (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sr1a*(this->alignN)+sr1a] += (-Kr0*dCrdC1a);

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a]) - (dBdT[sl0a] + dBdT[sl1a] + dBdT[sl2a]);
        const double dKcdTByKc0 = sumVdBdT0 + invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C[sl0a]*C[sl1a]*C[sl2a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        ddNdtByVdcTp[sl0a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sl1a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sl2a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sr0a*(this->alignN)+this->nSpecies] += dqdT0;
        ddNdtByVdcTp[sr1a*(this->alignN)+this->nSpecies] += dqdT0;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*(CF0-invKc0*CR0);
            __m256d Wv0 = _mm256_set1_pd(W0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv0,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2a*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2a*(this->alignN)+j],sl2v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->RR32DupAllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->RR32DupAllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->RR32DupAllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->RR32DupAllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->RR32DupAllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->RR32DupAllIndex[(k+1)*6+5];

        const double Kf1 = this->Kf_[i1];
        const double dCfdC0b = C[sl1b]*C[sl2b];
        const double dCfdC1b = C[sl0b]*C[sl2b];
        const double dCfdC2b = C[sl0b]*C[sl1b];
        ddNdtByVdcTp[sl0b*(this->alignN)+sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl1b*(this->alignN)+sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl2b*(this->alignN)+sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl0b*(this->alignN)+sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sl1b*(this->alignN)+sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sl2b*(this->alignN)+sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sl0b*(this->alignN)+sl2b] -= Kf1*dCfdC2b;
        ddNdtByVdcTp[sl1b*(this->alignN)+sl2b] -= Kf1*dCfdC2b;
        ddNdtByVdcTp[sl2b*(this->alignN)+sl2b] -= Kf1*dCfdC2b;
        ddNdtByVdcTp[sr0b*(this->alignN)+sl0b] += Kf1*dCfdC0b;
        ddNdtByVdcTp[sr1b*(this->alignN)+sl0b] += Kf1*dCfdC0b;
        ddNdtByVdcTp[sr0b*(this->alignN)+sl1b] += Kf1*dCfdC1b;
        ddNdtByVdcTp[sr1b*(this->alignN)+sl1b] += Kf1*dCfdC1b;
        ddNdtByVdcTp[sr0b*(this->alignN)+sl2b] += Kf1*dCfdC2b;
        ddNdtByVdcTp[sr1b*(this->alignN)+sl2b] += Kf1*dCfdC2b;

        double invKc1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*
            (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double dCrdC0b = C[sr1b];
        const double dCrdC1b = C[sr0b];
        ddNdtByVdcTp[sl0b*(this->alignN)+sr0b] -= (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sl1b*(this->alignN)+sr0b] -= (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sl2b*(this->alignN)+sr0b] -= (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sl0b*(this->alignN)+sr1b] -= (-Kr1*dCrdC1b);
        ddNdtByVdcTp[sl1b*(this->alignN)+sr1b] -= (-Kr1*dCrdC1b);
        ddNdtByVdcTp[sl2b*(this->alignN)+sr1b] -= (-Kr1*dCrdC1b);
        ddNdtByVdcTp[sr0b*(this->alignN)+sr0b] += (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sr1b*(this->alignN)+sr0b] += (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sr0b*(this->alignN)+sr1b] += (-Kr1*dCrdC1b);
        ddNdtByVdcTp[sr1b*(this->alignN)+sr1b] += (-Kr1*dCrdC1b);

        const double dKfdT1 = this->dKfdT_[i1];
        const double sumVdBdT1 = (dBdT[sr0b] + dBdT[sr1b]) - (dBdT[sl0b] + dBdT[sl1b] + dBdT[sl2b]);
        const double dKcdTByKc1 = sumVdBdT1 + invT;
        const double dKrdT1 = (dKfdT1*invKc1 - (invKc1 < invKcLimiter ? Kr1*dKcdTByKc1 : 0));
        const double CF1 = C[sl0b]*C[sl1b]*C[sl2b];
        const double CR1 = C[sr0b]*C[sr1b];
        const double dqdT1 = (dKfdT1*CF1) - (dKrdT1*CR1);
        ddNdtByVdcTp[sl0b*(this->alignN)+this->nSpecies] -= dqdT1;
        ddNdtByVdcTp[sl1b*(this->alignN)+this->nSpecies] -= dqdT1;
        ddNdtByVdcTp[sl2b*(this->alignN)+this->nSpecies] -= dqdT1;
        ddNdtByVdcTp[sr0b*(this->alignN)+this->nSpecies] += dqdT1;
        ddNdtByVdcTp[sr1b*(this->alignN)+this->nSpecies] += dqdT1;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*(CF1-invKc1*CR1);
            __m256d Wv1 = _mm256_set1_pd(W1);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv1,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2b*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2b*(this->alignN)+j],sl2v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->RR32DupAllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->RR32DupAllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->RR32DupAllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->RR32DupAllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->RR32DupAllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->RR32DupAllIndex[(k+2)*6+5];

        const double Kf2 = this->Kf_[i2];
        const double dCfdC0c = C[sl1c]*C[sl2c];
        const double dCfdC1c = C[sl0c]*C[sl2c];
        const double dCfdC2c = C[sl0c]*C[sl1c];
        ddNdtByVdcTp[sl0c*(this->alignN)+sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl1c*(this->alignN)+sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl2c*(this->alignN)+sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl0c*(this->alignN)+sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sl1c*(this->alignN)+sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sl2c*(this->alignN)+sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sl0c*(this->alignN)+sl2c] -= Kf2*dCfdC2c;
        ddNdtByVdcTp[sl1c*(this->alignN)+sl2c] -= Kf2*dCfdC2c;
        ddNdtByVdcTp[sl2c*(this->alignN)+sl2c] -= Kf2*dCfdC2c;
        ddNdtByVdcTp[sr0c*(this->alignN)+sl0c] += Kf2*dCfdC0c;
        ddNdtByVdcTp[sr1c*(this->alignN)+sl0c] += Kf2*dCfdC0c;
        ddNdtByVdcTp[sr0c*(this->alignN)+sl1c] += Kf2*dCfdC1c;
        ddNdtByVdcTp[sr1c*(this->alignN)+sl1c] += Kf2*dCfdC1c;
        ddNdtByVdcTp[sr0c*(this->alignN)+sl2c] += Kf2*dCfdC2c;
        ddNdtByVdcTp[sr1c*(this->alignN)+sl2c] += Kf2*dCfdC2c;

        double invKc2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*
            (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double dCrdC0c = C[sr1c];
        const double dCrdC1c = C[sr0c];
        ddNdtByVdcTp[sl0c*(this->alignN)+sr0c] -= (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sl1c*(this->alignN)+sr0c] -= (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sl2c*(this->alignN)+sr0c] -= (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sl0c*(this->alignN)+sr1c] -= (-Kr2*dCrdC1c);
        ddNdtByVdcTp[sl1c*(this->alignN)+sr1c] -= (-Kr2*dCrdC1c);
        ddNdtByVdcTp[sl2c*(this->alignN)+sr1c] -= (-Kr2*dCrdC1c);
        ddNdtByVdcTp[sr0c*(this->alignN)+sr0c] += (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sr1c*(this->alignN)+sr0c] += (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sr0c*(this->alignN)+sr1c] += (-Kr2*dCrdC1c);
        ddNdtByVdcTp[sr1c*(this->alignN)+sr1c] += (-Kr2*dCrdC1c);

        const double dKfdT2 = this->dKfdT_[i2];
        const double sumVdBdT2 = (dBdT[sr0c] + dBdT[sr1c]) - (dBdT[sl0c] + dBdT[sl1c] + dBdT[sl2c]);
        const double dKcdTByKc2 = sumVdBdT2 + invT;
        const double dKrdT2 = (dKfdT2*invKc2 - (invKc2 < invKcLimiter ? Kr2*dKcdTByKc2 : 0));
        const double CF2 = C[sl0c]*C[sl1c]*C[sl2c];
        const double CR2 = C[sr0c]*C[sr1c];
        const double dqdT2 = (dKfdT2*CF2) - (dKrdT2*CR2);
        ddNdtByVdcTp[sl0c*(this->alignN)+this->nSpecies] -= dqdT2;
        ddNdtByVdcTp[sl1c*(this->alignN)+this->nSpecies] -= dqdT2;
        ddNdtByVdcTp[sl2c*(this->alignN)+this->nSpecies] -= dqdT2;
        ddNdtByVdcTp[sr0c*(this->alignN)+this->nSpecies] += dqdT2;
        ddNdtByVdcTp[sr1c*(this->alignN)+this->nSpecies] += dqdT2;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*(CF2-invKc2*CR2);
            __m256d Wv2 = _mm256_set1_pd(W2);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv2,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2c*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2c*(this->alignN)+j],sl2v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->RR32DupAllIndex[(k+3)*6+0];
        const unsigned int sl0d = this->RR32DupAllIndex[(k+3)*6+1];
        const unsigned int sl1d = this->RR32DupAllIndex[(k+3)*6+2];
        const unsigned int sl2d = this->RR32DupAllIndex[(k+3)*6+3];
        const unsigned int sr0d = this->RR32DupAllIndex[(k+3)*6+4];
        const unsigned int sr1d = this->RR32DupAllIndex[(k+3)*6+5];

        const double Kf3 = this->Kf_[i3];
        const double dCfdC0d = C[sl1d]*C[sl2d];
        const double dCfdC1d = C[sl0d]*C[sl2d];
        const double dCfdC2d = C[sl0d]*C[sl1d];
        ddNdtByVdcTp[sl0d*(this->alignN)+sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl1d*(this->alignN)+sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl2d*(this->alignN)+sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl0d*(this->alignN)+sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sl1d*(this->alignN)+sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sl2d*(this->alignN)+sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sl0d*(this->alignN)+sl2d] -= Kf3*dCfdC2d;
        ddNdtByVdcTp[sl1d*(this->alignN)+sl2d] -= Kf3*dCfdC2d;
        ddNdtByVdcTp[sl2d*(this->alignN)+sl2d] -= Kf3*dCfdC2d;
        ddNdtByVdcTp[sr0d*(this->alignN)+sl0d] += Kf3*dCfdC0d;
        ddNdtByVdcTp[sr1d*(this->alignN)+sl0d] += Kf3*dCfdC0d;
        ddNdtByVdcTp[sr0d*(this->alignN)+sl1d] += Kf3*dCfdC1d;
        ddNdtByVdcTp[sr1d*(this->alignN)+sl1d] += Kf3*dCfdC1d;
        ddNdtByVdcTp[sr0d*(this->alignN)+sl2d] += Kf3*dCfdC2d;
        ddNdtByVdcTp[sr1d*(this->alignN)+sl2d] += Kf3*dCfdC2d;

        double invKc3 = (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d])*
            (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]*ExpNegGbyRT[sl2d])*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double dCrdC0d = C[sr1d];
        const double dCrdC1d = C[sr0d];
        ddNdtByVdcTp[sl0d*(this->alignN)+sr0d] -= (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sl1d*(this->alignN)+sr0d] -= (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sl2d*(this->alignN)+sr0d] -= (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sl0d*(this->alignN)+sr1d] -= (-Kr3*dCrdC1d);
        ddNdtByVdcTp[sl1d*(this->alignN)+sr1d] -= (-Kr3*dCrdC1d);
        ddNdtByVdcTp[sl2d*(this->alignN)+sr1d] -= (-Kr3*dCrdC1d);
        ddNdtByVdcTp[sr0d*(this->alignN)+sr0d] += (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sr1d*(this->alignN)+sr0d] += (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sr0d*(this->alignN)+sr1d] += (-Kr3*dCrdC1d);
        ddNdtByVdcTp[sr1d*(this->alignN)+sr1d] += (-Kr3*dCrdC1d);

        const double dKfdT3 = this->dKfdT_[i3];
        const double sumVdBdT3 = (dBdT[sr0d] + dBdT[sr1d]) - (dBdT[sl0d] + dBdT[sl1d] + dBdT[sl2d]);
        const double dKcdTByKc3 = sumVdBdT3 + invT;
        const double dKrdT3 = (dKfdT3*invKc3 - (invKc3 < invKcLimiter ? Kr3*dKcdTByKc3 : 0));
        const double CF3 = C[sl0d]*C[sl1d]*C[sl2d];
        const double CR3 = C[sr0d]*C[sr1d];
        const double dqdT3 = (dKfdT3*CF3) - (dKrdT3*CR3);
        ddNdtByVdcTp[sl0d*(this->alignN)+this->nSpecies] -= dqdT3;
        ddNdtByVdcTp[sl1d*(this->alignN)+this->nSpecies] -= dqdT3;
        ddNdtByVdcTp[sl2d*(this->alignN)+this->nSpecies] -= dqdT3;
        ddNdtByVdcTp[sr0d*(this->alignN)+this->nSpecies] += dqdT3;
        ddNdtByVdcTp[sr1d*(this->alignN)+this->nSpecies] += dqdT3;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*(CF3-invKc3*CR3);
            __m256d Wv3 = _mm256_set1_pd(W3);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv3,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2d*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2d*(this->alignN)+j],sl2v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j],sr1v);
            }
        }
    }
    for(std::size_t k=endRR32Dup-remainRR32Dup; k<endRR32Dup; k=k+1)
    {
        const unsigned int i0 = this->RR32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->RR32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->RR32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->RR32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->RR32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->RR32DupAllIndex[(k+0)*6+5];

        const double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C[sl1a]*C[sl2a];
        const double dCfdC1a = C[sl0a]*C[sl2a];
        const double dCfdC2a = C[sl0a]*C[sl1a];
        ddNdtByVdcTp[sl0a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl0a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl2a] += Kf0*dCfdC2a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl2a] += Kf0*dCfdC2a;

        double invKc0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
            (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a];
        ddNdtByVdcTp[sl0a*(this->alignN)+sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl1a*(this->alignN)+sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl2a*(this->alignN)+sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl0a*(this->alignN)+sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sl1a*(this->alignN)+sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sl2a*(this->alignN)+sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sr0a*(this->alignN)+sr0a] += (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sr1a*(this->alignN)+sr0a] += (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sr0a*(this->alignN)+sr1a] += (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sr1a*(this->alignN)+sr1a] += (-Kr0*dCrdC1a);

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a]) - (dBdT[sl0a] + dBdT[sl1a] + dBdT[sl2a]);
        const double dKcdTByKc0 = sumVdBdT0 + invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C[sl0a]*C[sl1a]*C[sl2a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        ddNdtByVdcTp[sl0a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sl1a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sl2a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sr0a*(this->alignN)+this->nSpecies] += dqdT0;
        ddNdtByVdcTp[sr1a*(this->alignN)+this->nSpecies] += dqdT0;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*(CF0-invKc0*CR0);
            __m256d Wv0 = _mm256_set1_pd(W0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv0,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2a*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2a*(this->alignN)+j],sl2v);
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


void  FastChemistry::OptReaction::JF32IR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{
    std::size_t endIR32 = this->IR32Size;
    std::size_t remainIR32 = endIR32%4;
    for(std::size_t k=0; k<endIR32-remainIR32; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->IR32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->IR32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->IR32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->IR32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->IR32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->IR32AllIndex[(k+0)*6+5];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sl2a_ptr = &ddNdtByVdcTp[sl2a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double C_sl2a = C[sl2a];
        const double Kf0dCfdC0a = Kf0*C_sl1a*C_sl2a;
        const double Kf0dCfdC1a = Kf0*C_sl0a*C_sl2a;
        const double Kf0dCfdC2a = Kf0*C_sl0a*C_sl1a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl2a_sl0a_increment = J_sl2a_ptr[sl0a] - Kf0dCfdC0a;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        J_sl2a_ptr[sl0a] = J_sl2a_sl0a_increment;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl2a_sl1a_increment = J_sl2a_ptr[sl1a] - Kf0dCfdC1a;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sl2a_ptr[sl1a] = J_sl2a_sl1a_increment;
        const double J_sl0a_sl2a_increment = J_sl0a_ptr[sl2a] - Kf0dCfdC2a;
        const double J_sl1a_sl2a_increment = J_sl1a_ptr[sl2a] - Kf0dCfdC2a;
        const double J_sl2a_sl2a_increment = J_sl2a_ptr[sl2a] - Kf0dCfdC2a;
        J_sl0a_ptr[sl2a] = J_sl0a_sl2a_increment;
        J_sl1a_ptr[sl2a] = J_sl1a_sl2a_increment;
        J_sl2a_ptr[sl2a] = J_sl2a_sl2a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0dCfdC0a;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1a;
        const double J_sr1a_sl1a_increment = J_sr1a_ptr[sl1a] + Kf0dCfdC1a;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;
        J_sr1a_ptr[sl1a] = J_sr1a_sl1a_increment;
        const double J_sr0a_sl2a_increment = J_sr0a_ptr[sl2a] + Kf0dCfdC2a;
        const double J_sr1a_sl2a_increment = J_sr1a_ptr[sl2a] + Kf0dCfdC2a;
        J_sr0a_ptr[sl2a] = J_sr0a_sl2a_increment;
        J_sr1a_ptr[sl2a] = J_sr1a_sl2a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C[sl0a]*C[sl1a]*C[sl2a];
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sl2a_T_increment = J_sl2a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sl2a_ptr[this->nSpecies] = J_sl2a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0;
            __m256d Wv0 = _mm256_set1_pd(W0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv0,dMdC);
                __m256d sl0v = load256d(&J_sl0a_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1a_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1a_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2a_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2a_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0a_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1a_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->IR32AllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->IR32AllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->IR32AllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->IR32AllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->IR32AllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->IR32AllIndex[(k+1)*6+5];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sl1b_ptr = &ddNdtByVdcTp[sl1b*(this->alignN)];
        double* __restrict__ J_sl2b_ptr = &ddNdtByVdcTp[sl2b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double C_sl1b = C[sl1b];
        const double C_sl0b = C[sl0b];
        const double C_sl2b = C[sl2b];
        const double Kf1dCfdC0b = Kf1*C_sl1b*C_sl2b;
        const double Kf1dCfdC1b = Kf1*C_sl0b*C_sl2b;
        const double Kf1dCfdC2b = Kf1*C_sl0b*C_sl1b;
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1dCfdC0b;
        const double J_sl1b_sl0b_increment = J_sl1b_ptr[sl0b] - Kf1dCfdC0b;
        const double J_sl2b_sl0b_increment = J_sl2b_ptr[sl0b] - Kf1dCfdC0b;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sl1b_ptr[sl0b] = J_sl1b_sl0b_increment;
        J_sl2b_ptr[sl0b] = J_sl2b_sl0b_increment;
        const double J_sl0b_sl1b_increment = J_sl0b_ptr[sl1b] - Kf1dCfdC1b;
        const double J_sl1b_sl1b_increment = J_sl1b_ptr[sl1b] - Kf1dCfdC1b;
        const double J_sl2b_sl1b_increment = J_sl2b_ptr[sl1b] - Kf1dCfdC1b;
        J_sl0b_ptr[sl1b] = J_sl0b_sl1b_increment;
        J_sl1b_ptr[sl1b] = J_sl1b_sl1b_increment;
        J_sl2b_ptr[sl1b] = J_sl2b_sl1b_increment;
        const double J_sl0b_sl2b_increment = J_sl0b_ptr[sl2b] - Kf1dCfdC2b;
        const double J_sl1b_sl2b_increment = J_sl1b_ptr[sl2b] - Kf1dCfdC2b;
        const double J_sl2b_sl2b_increment = J_sl2b_ptr[sl2b] - Kf1dCfdC2b;
        J_sl0b_ptr[sl2b] = J_sl0b_sl2b_increment;
        J_sl1b_ptr[sl2b] = J_sl1b_sl2b_increment;
        J_sl2b_ptr[sl2b] = J_sl2b_sl2b_increment;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1dCfdC0b;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1dCfdC0b;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;
        const double J_sr0b_sl1b_increment = J_sr0b_ptr[sl1b] + Kf1dCfdC1b;
        const double J_sr1b_sl1b_increment = J_sr1b_ptr[sl1b] + Kf1dCfdC1b;
        J_sr0b_ptr[sl1b] = J_sr0b_sl1b_increment;
        J_sr1b_ptr[sl1b] = J_sr1b_sl1b_increment;
        const double J_sr0b_sl2b_increment = J_sr0b_ptr[sl2b] + Kf1dCfdC2b;
        const double J_sr1b_sl2b_increment = J_sr1b_ptr[sl2b] + Kf1dCfdC2b;
        J_sr0b_ptr[sl2b] = J_sr0b_sl2b_increment;
        J_sr1b_ptr[sl2b] = J_sr1b_sl2b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double CF1 = C[sl0b]*C[sl1b]*C[sl2b];
        const double dqdT1 = (dKfdT1*CF1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sl1b_T_increment = J_sl1b_ptr[this->nSpecies] - dqdT1;
        const double J_sl2b_T_increment = J_sl2b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sl1b_ptr[this->nSpecies] = J_sl1b_T_increment;
        J_sl2b_ptr[this->nSpecies] = J_sl2b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;

        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*CF1;
            __m256d Wv1 = _mm256_set1_pd(W1);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv1,dMdC);
                __m256d sl0v = load256d(&J_sl0b_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1b_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1b_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2b_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2b_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0b_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0b_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1b_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1b_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->IR32AllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->IR32AllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->IR32AllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->IR32AllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->IR32AllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->IR32AllIndex[(k+2)*6+5];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sl1c_ptr = &ddNdtByVdcTp[sl1c*(this->alignN)];
        double* __restrict__ J_sl2c_ptr = &ddNdtByVdcTp[sl2c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double C_sl1c = C[sl1c];
        const double C_sl0c = C[sl0c];
        const double C_sl2c = C[sl2c];
        const double Kf2dCfdC0c = Kf2*C_sl1c*C_sl2c;
        const double Kf2dCfdC1c = Kf2*C_sl0c*C_sl2c;
        const double Kf2dCfdC2c = Kf2*C_sl0c*C_sl1c;
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2dCfdC0c;
        const double J_sl1c_sl0c_increment = J_sl1c_ptr[sl0c] - Kf2dCfdC0c;
        const double J_sl2c_sl0c_increment = J_sl2c_ptr[sl0c] - Kf2dCfdC0c;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sl1c_ptr[sl0c] = J_sl1c_sl0c_increment;
        J_sl2c_ptr[sl0c] = J_sl2c_sl0c_increment;
        const double J_sl0c_sl1c_increment = J_sl0c_ptr[sl1c] - Kf2dCfdC1c;
        const double J_sl1c_sl1c_increment = J_sl1c_ptr[sl1c] - Kf2dCfdC1c;
        const double J_sl2c_sl1c_increment = J_sl2c_ptr[sl1c] - Kf2dCfdC1c;
        J_sl0c_ptr[sl1c] = J_sl0c_sl1c_increment;
        J_sl1c_ptr[sl1c] = J_sl1c_sl1c_increment;
        J_sl2c_ptr[sl1c] = J_sl2c_sl1c_increment;
        const double J_sl0c_sl2c_increment = J_sl0c_ptr[sl2c] - Kf2dCfdC2c;
        const double J_sl1c_sl2c_increment = J_sl1c_ptr[sl2c] - Kf2dCfdC2c;
        const double J_sl2c_sl2c_increment = J_sl2c_ptr[sl2c] - Kf2dCfdC2c;
        J_sl0c_ptr[sl2c] = J_sl0c_sl2c_increment;
        J_sl1c_ptr[sl2c] = J_sl1c_sl2c_increment;
        J_sl2c_ptr[sl2c] = J_sl2c_sl2c_increment;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2dCfdC0c;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2dCfdC0c;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;
        const double J_sr0c_sl1c_increment = J_sr0c_ptr[sl1c] + Kf2dCfdC1c;
        const double J_sr1c_sl1c_increment = J_sr1c_ptr[sl1c] + Kf2dCfdC1c;
        J_sr0c_ptr[sl1c] = J_sr0c_sl1c_increment;
        J_sr1c_ptr[sl1c] = J_sr1c_sl1c_increment;
        const double J_sr0c_sl2c_increment = J_sr0c_ptr[sl2c] + Kf2dCfdC2c;
        const double J_sr1c_sl2c_increment = J_sr1c_ptr[sl2c] + Kf2dCfdC2c;
        J_sr0c_ptr[sl2c] = J_sr0c_sl2c_increment;
        J_sr1c_ptr[sl2c] = J_sr1c_sl2c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double CF2 = C[sl0c]*C[sl1c]*C[sl2c];
        const double dqdT2 = (dKfdT2*CF2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sl1c_T_increment = J_sl1c_ptr[this->nSpecies] - dqdT2;
        const double J_sl2c_T_increment = J_sl2c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sl1c_ptr[this->nSpecies] = J_sl1c_T_increment;
        J_sl2c_ptr[this->nSpecies] = J_sl2c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;

        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*CF2;
            __m256d Wv2 = _mm256_set1_pd(W2);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv2,dMdC);
                __m256d sl0v = load256d(&J_sl0c_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1c_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1c_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2c_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2c_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0c_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0c_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1c_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1c_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->IR32AllIndex[(k+3)*6+0];
        const unsigned int sl0d = this->IR32AllIndex[(k+3)*6+1];
        const unsigned int sl1d = this->IR32AllIndex[(k+3)*6+2];
        const unsigned int sl2d = this->IR32AllIndex[(k+3)*6+3];
        const unsigned int sr0d = this->IR32AllIndex[(k+3)*6+4];
        const unsigned int sr1d = this->IR32AllIndex[(k+3)*6+5];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sl1d_ptr = &ddNdtByVdcTp[sl1d*(this->alignN)];
        double* __restrict__ J_sl2d_ptr = &ddNdtByVdcTp[sl2d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double C_sl1d = C[sl1d];
        const double C_sl0d = C[sl0d];
        const double C_sl2d = C[sl2d];
        const double Kf3dCfdC0d = Kf3*C_sl1d*C_sl2d;
        const double Kf3dCfdC1d = Kf3*C_sl0d*C_sl2d;
        const double Kf3dCfdC2d = Kf3*C_sl0d*C_sl1d;
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3dCfdC0d;
        const double J_sl1d_sl0d_increment = J_sl1d_ptr[sl0d] - Kf3dCfdC0d;
        const double J_sl2d_sl0d_increment = J_sl2d_ptr[sl0d] - Kf3dCfdC0d;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sl1d_ptr[sl0d] = J_sl1d_sl0d_increment;
        J_sl2d_ptr[sl0d] = J_sl2d_sl0d_increment;
        const double J_sl0d_sl1d_increment = J_sl0d_ptr[sl1d] - Kf3dCfdC1d;
        const double J_sl1d_sl1d_increment = J_sl1d_ptr[sl1d] - Kf3dCfdC1d;
        const double J_sl2d_sl1d_increment = J_sl2d_ptr[sl1d] - Kf3dCfdC1d;
        J_sl0d_ptr[sl1d] = J_sl0d_sl1d_increment;
        J_sl1d_ptr[sl1d] = J_sl1d_sl1d_increment;
        J_sl2d_ptr[sl1d] = J_sl2d_sl1d_increment;
        const double J_sl0d_sl2d_increment = J_sl0d_ptr[sl2d] - Kf3dCfdC2d;
        const double J_sl1d_sl2d_increment = J_sl1d_ptr[sl2d] - Kf3dCfdC2d;
        const double J_sl2d_sl2d_increment = J_sl2d_ptr[sl2d] - Kf3dCfdC2d;
        J_sl0d_ptr[sl2d] = J_sl0d_sl2d_increment;
        J_sl1d_ptr[sl2d] = J_sl1d_sl2d_increment;
        J_sl2d_ptr[sl2d] = J_sl2d_sl2d_increment;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3dCfdC0d;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3dCfdC0d;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;
        const double J_sr0d_sl1d_increment = J_sr0d_ptr[sl1d] + Kf3dCfdC1d;
        const double J_sr1d_sl1d_increment = J_sr1d_ptr[sl1d] + Kf3dCfdC1d;
        J_sr0d_ptr[sl1d] = J_sr0d_sl1d_increment;
        J_sr1d_ptr[sl1d] = J_sr1d_sl1d_increment;
        const double J_sr0d_sl2d_increment = J_sr0d_ptr[sl2d] + Kf3dCfdC2d;
        const double J_sr1d_sl2d_increment = J_sr1d_ptr[sl2d] + Kf3dCfdC2d;
        J_sr0d_ptr[sl2d] = J_sr0d_sl2d_increment;
        J_sr1d_ptr[sl2d] = J_sr1d_sl2d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double CF3 = C[sl0d]*C[sl1d]*C[sl2d];
        const double dqdT3 = (dKfdT3*CF3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sl1d_T_increment = J_sl1d_ptr[this->nSpecies] - dqdT3;
        const double J_sl2d_T_increment = J_sl2d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sl1d_ptr[this->nSpecies] = J_sl1d_T_increment;
        J_sl2d_ptr[this->nSpecies] = J_sl2d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;

        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sl2d_increment = dNdtByV[sl2d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sl2d] = omega_sl2d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*CF3;
            __m256d Wv3 = _mm256_set1_pd(W3);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv3,dMdC);
                __m256d sl0v = load256d(&J_sl0d_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1d_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1d_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2d_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2d_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0d_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0d_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1d_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1d_ptr[j],sr1v);
            }
        }
    }
    for(std::size_t k=endIR32-remainIR32; k<endIR32; k=k+1)
    {
        const unsigned int i0 = this->IR32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->IR32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->IR32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->IR32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->IR32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->IR32AllIndex[(k+0)*6+5];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sl2a_ptr = &ddNdtByVdcTp[sl2a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double C_sl2a = C[sl2a];
        const double Kf0dCfdC0a = Kf0*C_sl1a*C_sl2a;
        const double Kf0dCfdC1a = Kf0*C_sl0a*C_sl2a;
        const double Kf0dCfdC2a = Kf0*C_sl0a*C_sl1a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl2a_sl0a_increment = J_sl2a_ptr[sl0a] - Kf0dCfdC0a;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        J_sl2a_ptr[sl0a] = J_sl2a_sl0a_increment;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl2a_sl1a_increment = J_sl2a_ptr[sl1a] - Kf0dCfdC1a;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sl2a_ptr[sl1a] = J_sl2a_sl1a_increment;
        const double J_sl0a_sl2a_increment = J_sl0a_ptr[sl2a] - Kf0dCfdC2a;
        const double J_sl1a_sl2a_increment = J_sl1a_ptr[sl2a] - Kf0dCfdC2a;
        const double J_sl2a_sl2a_increment = J_sl2a_ptr[sl2a] - Kf0dCfdC2a;
        J_sl0a_ptr[sl2a] = J_sl0a_sl2a_increment;
        J_sl1a_ptr[sl2a] = J_sl1a_sl2a_increment;
        J_sl2a_ptr[sl2a] = J_sl2a_sl2a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0dCfdC0a;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1a;
        const double J_sr1a_sl1a_increment = J_sr1a_ptr[sl1a] + Kf0dCfdC1a;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;
        J_sr1a_ptr[sl1a] = J_sr1a_sl1a_increment;
        const double J_sr0a_sl2a_increment = J_sr0a_ptr[sl2a] + Kf0dCfdC2a;
        const double J_sr1a_sl2a_increment = J_sr1a_ptr[sl2a] + Kf0dCfdC2a;
        J_sr0a_ptr[sl2a] = J_sr0a_sl2a_increment;
        J_sr1a_ptr[sl2a] = J_sr1a_sl2a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C[sl0a]*C[sl1a]*C[sl2a];
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sl2a_T_increment = J_sl2a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sl2a_ptr[this->nSpecies] = J_sl2a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0;
            __m256d Wv0 = _mm256_set1_pd(W0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv0,dMdC);
                __m256d sl0v = load256d(&J_sl0a_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1a_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1a_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2a_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2a_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0a_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1a_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }
    }

    std::size_t endIR32Dup = this->IR32DupSize;
    std::size_t remainIR32Dup = endIR32Dup%4;
    for(std::size_t k=0; k<endIR32Dup-remainIR32Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->IR32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->IR32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->IR32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->IR32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->IR32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->IR32DupAllIndex[(k+0)*6+5];

        const double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C[sl1a]*C[sl2a];
        const double dCfdC1a = C[sl0a]*C[sl2a];
        const double dCfdC2a = C[sl0a]*C[sl1a];
        ddNdtByVdcTp[sl0a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl0a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl2a] += Kf0*dCfdC2a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl2a] += Kf0*dCfdC2a;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C[sl0a]*C[sl1a]*C[sl2a];
        const double dqdT0 = (dKfdT0*CF0);
        ddNdtByVdcTp[sl0a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sl1a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sl2a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sr0a*(this->alignN)+this->nSpecies] += dqdT0;
        ddNdtByVdcTp[sr1a*(this->alignN)+this->nSpecies] += dqdT0;

        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0;
            __m256d Wv0 = _mm256_set1_pd(W0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv0,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2a*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2a*(this->alignN)+j],sl2v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->IR32DupAllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->IR32DupAllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->IR32DupAllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->IR32DupAllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->IR32DupAllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->IR32DupAllIndex[(k+1)*6+5];

        const double Kf1 = this->Kf_[i1];
        const double dCfdC0b = C[sl1b]*C[sl2b];
        const double dCfdC1b = C[sl0b]*C[sl2b];
        const double dCfdC2b = C[sl0b]*C[sl1b];
        ddNdtByVdcTp[sl0b*(this->alignN)+sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl1b*(this->alignN)+sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl2b*(this->alignN)+sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl0b*(this->alignN)+sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sl1b*(this->alignN)+sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sl2b*(this->alignN)+sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sl0b*(this->alignN)+sl2b] -= Kf1*dCfdC2b;
        ddNdtByVdcTp[sl1b*(this->alignN)+sl2b] -= Kf1*dCfdC2b;
        ddNdtByVdcTp[sl2b*(this->alignN)+sl2b] -= Kf1*dCfdC2b;
        ddNdtByVdcTp[sr0b*(this->alignN)+sl0b] += Kf1*dCfdC0b;
        ddNdtByVdcTp[sr1b*(this->alignN)+sl0b] += Kf1*dCfdC0b;
        ddNdtByVdcTp[sr0b*(this->alignN)+sl1b] += Kf1*dCfdC1b;
        ddNdtByVdcTp[sr1b*(this->alignN)+sl1b] += Kf1*dCfdC1b;
        ddNdtByVdcTp[sr0b*(this->alignN)+sl2b] += Kf1*dCfdC2b;
        ddNdtByVdcTp[sr1b*(this->alignN)+sl2b] += Kf1*dCfdC2b;

        const double dKfdT1 = this->dKfdT_[i1];
        const double CF1 = C[sl0b]*C[sl1b]*C[sl2b];
        const double dqdT1 = (dKfdT1*CF1);
        ddNdtByVdcTp[sl0b*(this->alignN)+this->nSpecies] -= dqdT1;
        ddNdtByVdcTp[sl1b*(this->alignN)+this->nSpecies] -= dqdT1;
        ddNdtByVdcTp[sl2b*(this->alignN)+this->nSpecies] -= dqdT1;
        ddNdtByVdcTp[sr0b*(this->alignN)+this->nSpecies] += dqdT1;
        ddNdtByVdcTp[sr1b*(this->alignN)+this->nSpecies] += dqdT1;

        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*CF1;
            __m256d Wv1 = _mm256_set1_pd(W1);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv1,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2b*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2b*(this->alignN)+j],sl2v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->IR32DupAllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->IR32DupAllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->IR32DupAllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->IR32DupAllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->IR32DupAllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->IR32DupAllIndex[(k+2)*6+5];

        const double Kf2 = this->Kf_[i2];
        const double dCfdC0c = C[sl1c]*C[sl2c];
        const double dCfdC1c = C[sl0c]*C[sl2c];
        const double dCfdC2c = C[sl0c]*C[sl1c];
        ddNdtByVdcTp[sl0c*(this->alignN)+sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl1c*(this->alignN)+sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl2c*(this->alignN)+sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl0c*(this->alignN)+sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sl1c*(this->alignN)+sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sl2c*(this->alignN)+sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sl0c*(this->alignN)+sl2c] -= Kf2*dCfdC2c;
        ddNdtByVdcTp[sl1c*(this->alignN)+sl2c] -= Kf2*dCfdC2c;
        ddNdtByVdcTp[sl2c*(this->alignN)+sl2c] -= Kf2*dCfdC2c;
        ddNdtByVdcTp[sr0c*(this->alignN)+sl0c] += Kf2*dCfdC0c;
        ddNdtByVdcTp[sr1c*(this->alignN)+sl0c] += Kf2*dCfdC0c;
        ddNdtByVdcTp[sr0c*(this->alignN)+sl1c] += Kf2*dCfdC1c;
        ddNdtByVdcTp[sr1c*(this->alignN)+sl1c] += Kf2*dCfdC1c;
        ddNdtByVdcTp[sr0c*(this->alignN)+sl2c] += Kf2*dCfdC2c;
        ddNdtByVdcTp[sr1c*(this->alignN)+sl2c] += Kf2*dCfdC2c;

        const double dKfdT2 = this->dKfdT_[i2];
        const double CF2 = C[sl0c]*C[sl1c]*C[sl2c];
        const double dqdT2 = (dKfdT2*CF2);
        ddNdtByVdcTp[sl0c*(this->alignN)+this->nSpecies] -= dqdT2;
        ddNdtByVdcTp[sl1c*(this->alignN)+this->nSpecies] -= dqdT2;
        ddNdtByVdcTp[sl2c*(this->alignN)+this->nSpecies] -= dqdT2;
        ddNdtByVdcTp[sr0c*(this->alignN)+this->nSpecies] += dqdT2;
        ddNdtByVdcTp[sr1c*(this->alignN)+this->nSpecies] += dqdT2;

        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*CF2;
            __m256d Wv2 = _mm256_set1_pd(W2);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv2,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2c*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2c*(this->alignN)+j],sl2v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->IR32DupAllIndex[(k+3)*6+0];
        const unsigned int sl0d = this->IR32DupAllIndex[(k+3)*6+1];
        const unsigned int sl1d = this->IR32DupAllIndex[(k+3)*6+2];
        const unsigned int sl2d = this->IR32DupAllIndex[(k+3)*6+3];
        const unsigned int sr0d = this->IR32DupAllIndex[(k+3)*6+4];
        const unsigned int sr1d = this->IR32DupAllIndex[(k+3)*6+5];

        const double Kf3 = this->Kf_[i3];
        const double dCfdC0d = C[sl1d]*C[sl2d];
        const double dCfdC1d = C[sl0d]*C[sl2d];
        const double dCfdC2d = C[sl0d]*C[sl1d];
        ddNdtByVdcTp[sl0d*(this->alignN)+sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl1d*(this->alignN)+sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl2d*(this->alignN)+sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl0d*(this->alignN)+sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sl1d*(this->alignN)+sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sl2d*(this->alignN)+sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sl0d*(this->alignN)+sl2d] -= Kf3*dCfdC2d;
        ddNdtByVdcTp[sl1d*(this->alignN)+sl2d] -= Kf3*dCfdC2d;
        ddNdtByVdcTp[sl2d*(this->alignN)+sl2d] -= Kf3*dCfdC2d;
        ddNdtByVdcTp[sr0d*(this->alignN)+sl0d] += Kf3*dCfdC0d;
        ddNdtByVdcTp[sr1d*(this->alignN)+sl0d] += Kf3*dCfdC0d;
        ddNdtByVdcTp[sr0d*(this->alignN)+sl1d] += Kf3*dCfdC1d;
        ddNdtByVdcTp[sr1d*(this->alignN)+sl1d] += Kf3*dCfdC1d;
        ddNdtByVdcTp[sr0d*(this->alignN)+sl2d] += Kf3*dCfdC2d;
        ddNdtByVdcTp[sr1d*(this->alignN)+sl2d] += Kf3*dCfdC2d;

        const double dKfdT3 = this->dKfdT_[i3];
        const double CF3 = C[sl0d]*C[sl1d]*C[sl2d];
        const double dqdT3 = (dKfdT3*CF3);
        ddNdtByVdcTp[sl0d*(this->alignN)+this->nSpecies] -= dqdT3;
        ddNdtByVdcTp[sl1d*(this->alignN)+this->nSpecies] -= dqdT3;
        ddNdtByVdcTp[sl2d*(this->alignN)+this->nSpecies] -= dqdT3;
        ddNdtByVdcTp[sr0d*(this->alignN)+this->nSpecies] += dqdT3;
        ddNdtByVdcTp[sr1d*(this->alignN)+this->nSpecies] += dqdT3;

        const double q3 = (Kf3*CF3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*CF3;
            __m256d Wv3 = _mm256_set1_pd(W3);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv3,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2d*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2d*(this->alignN)+j],sl2v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j],sr1v);
            }
        }
    }
    for(std::size_t k=endIR32Dup-remainIR32Dup; k<endIR32Dup; k=k+1)
    {
        const unsigned int i0 = this->IR32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->IR32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->IR32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->IR32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->IR32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->IR32DupAllIndex[(k+0)*6+5];

        const double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C[sl1a]*C[sl2a];
        const double dCfdC1a = C[sl0a]*C[sl2a];
        const double dCfdC2a = C[sl0a]*C[sl1a];
        ddNdtByVdcTp[sl0a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl0a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl2a] += Kf0*dCfdC2a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl2a] += Kf0*dCfdC2a;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C[sl0a]*C[sl1a]*C[sl2a];
        const double dqdT0 = (dKfdT0*CF0);
        ddNdtByVdcTp[sl0a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sl1a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sl2a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sr0a*(this->alignN)+this->nSpecies] += dqdT0;
        ddNdtByVdcTp[sr1a*(this->alignN)+this->nSpecies] += dqdT0;

        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0;
            __m256d Wv0 = _mm256_set1_pd(W0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv0,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2a*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2a*(this->alignN)+j],sl2v);
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


void  FastChemistry::OptReaction::JF32NER
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{
    std::size_t endNER32 = this->NER32Size;
    std::size_t remainNER32 = endNER32%4;
    for(std::size_t k=0; k<endNER32-remainNER32; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->NER32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->NER32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->NER32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->NER32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->NER32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->NER32AllIndex[(k+0)*6+5];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sl2a_ptr = &ddNdtByVdcTp[sl2a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double C_sl2a = C[sl2a];
        const double Kf0dCfdC0a = Kf0*C_sl1a*C_sl2a;
        const double Kf0dCfdC1a = Kf0*C_sl0a*C_sl2a;
        const double Kf0dCfdC2a = Kf0*C_sl0a*C_sl1a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl2a_sl0a_increment = J_sl2a_ptr[sl0a] - Kf0dCfdC0a;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        J_sl2a_ptr[sl0a] = J_sl2a_sl0a_increment;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl2a_sl1a_increment = J_sl2a_ptr[sl1a] - Kf0dCfdC1a;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sl2a_ptr[sl1a] = J_sl2a_sl1a_increment;
        const double J_sl0a_sl2a_increment = J_sl0a_ptr[sl2a] - Kf0dCfdC2a;
        const double J_sl1a_sl2a_increment = J_sl1a_ptr[sl2a] - Kf0dCfdC2a;
        const double J_sl2a_sl2a_increment = J_sl2a_ptr[sl2a] - Kf0dCfdC2a;
        J_sl0a_ptr[sl2a] = J_sl0a_sl2a_increment;
        J_sl1a_ptr[sl2a] = J_sl1a_sl2a_increment;
        J_sl2a_ptr[sl2a] = J_sl2a_sl2a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0dCfdC0a;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1a;
        const double J_sr1a_sl1a_increment = J_sr1a_ptr[sl1a] + Kf0dCfdC1a;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;
        J_sr1a_ptr[sl1a] = J_sr1a_sl1a_increment;
        const double J_sr0a_sl2a_increment = J_sr0a_ptr[sl2a] + Kf0dCfdC2a;
        const double J_sr1a_sl2a_increment = J_sr1a_ptr[sl2a] + Kf0dCfdC2a;
        J_sr0a_ptr[sl2a] = J_sr0a_sl2a_increment;
        J_sr1a_ptr[sl2a] = J_sr1a_sl2a_increment;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double C_sr1a = C[sr1a];
        const double C_sr0a = C[sr0a];
        const double Kr0dCrdC0a = -Kr0*C_sr1a;
        const double Kr0dCrdC1a = -Kr0*C_sr0a;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sl1a_sr0a_increment = J_sl1a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sl2a_sr0a_increment = J_sl2a_ptr[sr0a] - Kr0dCrdC0a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sl1a_ptr[sr0a] = J_sl1a_sr0a_increment;
        J_sl2a_ptr[sr0a] = J_sl2a_sr0a_increment;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] - Kr0dCrdC1a;
        const double J_sl1a_sr1a_increment = J_sl1a_ptr[sr1a] - Kr0dCrdC1a;
        const double J_sl2a_sr1a_increment = J_sl2a_ptr[sr1a] - Kr0dCrdC1a;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sl1a_ptr[sr1a] = J_sl1a_sr1a_increment;
        J_sl2a_ptr[sr1a] = J_sl2a_sr1a_increment;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] + Kr0dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] + Kr0dCrdC0a;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] + Kr0dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] + Kr0dCrdC1a;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a]*C[sl1a]*C[sl2a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sl2a_T_increment = J_sl2a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sl2a_ptr[this->nSpecies] = J_sl2a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[2] && i0<this->Ikf[3])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0 - this->dKfdC_[idx0+this->Itbr[4]]*CR0;
            __m256d Wv0 = _mm256_set1_pd(W0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv0,dMdC);
                __m256d sl0v = load256d(&J_sl0a_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1a_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1a_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2a_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2a_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0a_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1a_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->NER32AllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->NER32AllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->NER32AllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->NER32AllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->NER32AllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->NER32AllIndex[(k+1)*6+5];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sl1b_ptr = &ddNdtByVdcTp[sl1b*(this->alignN)];
        double* __restrict__ J_sl2b_ptr = &ddNdtByVdcTp[sl2b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double C_sl1b = C[sl1b];
        const double C_sl0b = C[sl0b];
        const double C_sl2b = C[sl2b];
        const double Kf1dCfdC0b = Kf1*C_sl1b*C_sl2b;
        const double Kf1dCfdC1b = Kf1*C_sl0b*C_sl2b;
        const double Kf1dCfdC2b = Kf1*C_sl0b*C_sl1b;
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1dCfdC0b;
        const double J_sl1b_sl0b_increment = J_sl1b_ptr[sl0b] - Kf1dCfdC0b;
        const double J_sl2b_sl0b_increment = J_sl2b_ptr[sl0b] - Kf1dCfdC0b;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sl1b_ptr[sl0b] = J_sl1b_sl0b_increment;
        J_sl2b_ptr[sl0b] = J_sl2b_sl0b_increment;
        const double J_sl0b_sl1b_increment = J_sl0b_ptr[sl1b] - Kf1dCfdC1b;
        const double J_sl1b_sl1b_increment = J_sl1b_ptr[sl1b] - Kf1dCfdC1b;
        const double J_sl2b_sl1b_increment = J_sl2b_ptr[sl1b] - Kf1dCfdC1b;
        J_sl0b_ptr[sl1b] = J_sl0b_sl1b_increment;
        J_sl1b_ptr[sl1b] = J_sl1b_sl1b_increment;
        J_sl2b_ptr[sl1b] = J_sl2b_sl1b_increment;
        const double J_sl0b_sl2b_increment = J_sl0b_ptr[sl2b] - Kf1dCfdC2b;
        const double J_sl1b_sl2b_increment = J_sl1b_ptr[sl2b] - Kf1dCfdC2b;
        const double J_sl2b_sl2b_increment = J_sl2b_ptr[sl2b] - Kf1dCfdC2b;
        J_sl0b_ptr[sl2b] = J_sl0b_sl2b_increment;
        J_sl1b_ptr[sl2b] = J_sl1b_sl2b_increment;
        J_sl2b_ptr[sl2b] = J_sl2b_sl2b_increment;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1dCfdC0b;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1dCfdC0b;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;
        const double J_sr0b_sl1b_increment = J_sr0b_ptr[sl1b] + Kf1dCfdC1b;
        const double J_sr1b_sl1b_increment = J_sr1b_ptr[sl1b] + Kf1dCfdC1b;
        J_sr0b_ptr[sl1b] = J_sr0b_sl1b_increment;
        J_sr1b_ptr[sl1b] = J_sr1b_sl1b_increment;
        const double J_sr0b_sl2b_increment = J_sr0b_ptr[sl2b] + Kf1dCfdC2b;
        const double J_sr1b_sl2b_increment = J_sr1b_ptr[sl2b] + Kf1dCfdC2b;
        J_sr0b_ptr[sl2b] = J_sr0b_sl2b_increment;
        J_sr1b_ptr[sl2b] = J_sr1b_sl2b_increment;

        const double Kr1 = this->Kf_[i1 - Ikf[1] + Ikf[9]];
        const double C_sr1b = C[sr1b];
        const double C_sr0b = C[sr0b];
        const double Kr1dCrdC0b = -Kr1*C_sr1b;
        const double Kr1dCrdC1b = -Kr1*C_sr0b;
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] - Kr1dCrdC0b;
        const double J_sl1b_sr0b_increment = J_sl1b_ptr[sr0b] - Kr1dCrdC0b;
        const double J_sl2b_sr0b_increment = J_sl2b_ptr[sr0b] - Kr1dCrdC0b;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sl1b_ptr[sr0b] = J_sl1b_sr0b_increment;
        J_sl2b_ptr[sr0b] = J_sl2b_sr0b_increment;
        const double J_sl0b_sr1b_increment = J_sl0b_ptr[sr1b] - Kr1dCrdC1b;
        const double J_sl1b_sr1b_increment = J_sl1b_ptr[sr1b] - Kr1dCrdC1b;
        const double J_sl2b_sr1b_increment = J_sl2b_ptr[sr1b] - Kr1dCrdC1b;
        J_sl0b_ptr[sr1b] = J_sl0b_sr1b_increment;
        J_sl1b_ptr[sr1b] = J_sl1b_sr1b_increment;
        J_sl2b_ptr[sr1b] = J_sl2b_sr1b_increment;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] + Kr1dCrdC0b;
        const double J_sr1b_sr0b_increment = J_sr1b_ptr[sr0b] + Kr1dCrdC0b;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;
        J_sr1b_ptr[sr0b] = J_sr1b_sr0b_increment;
        const double J_sr0b_sr1b_increment = J_sr0b_ptr[sr1b] + Kr1dCrdC1b;
        const double J_sr1b_sr1b_increment = J_sr1b_ptr[sr1b] + Kr1dCrdC1b;
        J_sr0b_ptr[sr1b] = J_sr0b_sr1b_increment;
        J_sr1b_ptr[sr1b] = J_sr1b_sr1b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double dKrdT1 = this->dKfdT_[i1 - Ikf[1] + Ikf[9]];
        const double CF1 = C[sl0b]*C[sl1b]*C[sl2b];
        const double CR1 = C[sr0b]*C[sr1b];
        const double dqdT1 = (dKfdT1*CF1) - (dKrdT1*CR1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sl1b_T_increment = J_sl1b_ptr[this->nSpecies] - dqdT1;
        const double J_sl2b_T_increment = J_sl2b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sl1b_ptr[this->nSpecies] = J_sl1b_T_increment;
        J_sl2b_ptr[this->nSpecies] = J_sl2b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        if(i1>=this->Ikf[2] && i1<this->Ikf[3])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*CF1 - this->dKfdC_[idx1+this->Itbr[4]]*CR1;
            __m256d Wv1 = _mm256_set1_pd(W1);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv1,dMdC);
                __m256d sl0v = load256d(&J_sl0b_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1b_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1b_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2b_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2b_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0b_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0b_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1b_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1b_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->NER32AllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->NER32AllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->NER32AllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->NER32AllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->NER32AllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->NER32AllIndex[(k+2)*6+5];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sl1c_ptr = &ddNdtByVdcTp[sl1c*(this->alignN)];
        double* __restrict__ J_sl2c_ptr = &ddNdtByVdcTp[sl2c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double C_sl1c = C[sl1c];
        const double C_sl0c = C[sl0c];
        const double C_sl2c = C[sl2c];
        const double Kf2dCfdC0c = Kf2*C_sl1c*C_sl2c;
        const double Kf2dCfdC1c = Kf2*C_sl0c*C_sl2c;
        const double Kf2dCfdC2c = Kf2*C_sl0c*C_sl1c;
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2dCfdC0c;
        const double J_sl1c_sl0c_increment = J_sl1c_ptr[sl0c] - Kf2dCfdC0c;
        const double J_sl2c_sl0c_increment = J_sl2c_ptr[sl0c] - Kf2dCfdC0c;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sl1c_ptr[sl0c] = J_sl1c_sl0c_increment;
        J_sl2c_ptr[sl0c] = J_sl2c_sl0c_increment;
        const double J_sl0c_sl1c_increment = J_sl0c_ptr[sl1c] - Kf2dCfdC1c;
        const double J_sl1c_sl1c_increment = J_sl1c_ptr[sl1c] - Kf2dCfdC1c;
        const double J_sl2c_sl1c_increment = J_sl2c_ptr[sl1c] - Kf2dCfdC1c;
        J_sl0c_ptr[sl1c] = J_sl0c_sl1c_increment;
        J_sl1c_ptr[sl1c] = J_sl1c_sl1c_increment;
        J_sl2c_ptr[sl1c] = J_sl2c_sl1c_increment;
        const double J_sl0c_sl2c_increment = J_sl0c_ptr[sl2c] - Kf2dCfdC2c;
        const double J_sl1c_sl2c_increment = J_sl1c_ptr[sl2c] - Kf2dCfdC2c;
        const double J_sl2c_sl2c_increment = J_sl2c_ptr[sl2c] - Kf2dCfdC2c;
        J_sl0c_ptr[sl2c] = J_sl0c_sl2c_increment;
        J_sl1c_ptr[sl2c] = J_sl1c_sl2c_increment;
        J_sl2c_ptr[sl2c] = J_sl2c_sl2c_increment;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2dCfdC0c;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2dCfdC0c;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;
        const double J_sr0c_sl1c_increment = J_sr0c_ptr[sl1c] + Kf2dCfdC1c;
        const double J_sr1c_sl1c_increment = J_sr1c_ptr[sl1c] + Kf2dCfdC1c;
        J_sr0c_ptr[sl1c] = J_sr0c_sl1c_increment;
        J_sr1c_ptr[sl1c] = J_sr1c_sl1c_increment;
        const double J_sr0c_sl2c_increment = J_sr0c_ptr[sl2c] + Kf2dCfdC2c;
        const double J_sr1c_sl2c_increment = J_sr1c_ptr[sl2c] + Kf2dCfdC2c;
        J_sr0c_ptr[sl2c] = J_sr0c_sl2c_increment;
        J_sr1c_ptr[sl2c] = J_sr1c_sl2c_increment;

        const double Kr2 = this->Kf_[i2 - Ikf[1] + Ikf[9]];
        const double C_sr1c = C[sr1c];
        const double C_sr0c = C[sr0c];
        const double Kr2dCrdC0c = -Kr2*C_sr1c;
        const double Kr2dCrdC1c = -Kr2*C_sr0c;
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] - Kr2dCrdC0c;
        const double J_sl1c_sr0c_increment = J_sl1c_ptr[sr0c] - Kr2dCrdC0c;
        const double J_sl2c_sr0c_increment = J_sl2c_ptr[sr0c] - Kr2dCrdC0c;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sl1c_ptr[sr0c] = J_sl1c_sr0c_increment;
        J_sl2c_ptr[sr0c] = J_sl2c_sr0c_increment;
        const double J_sl0c_sr1c_increment = J_sl0c_ptr[sr1c] - Kr2dCrdC1c;
        const double J_sl1c_sr1c_increment = J_sl1c_ptr[sr1c] - Kr2dCrdC1c;
        const double J_sl2c_sr1c_increment = J_sl2c_ptr[sr1c] - Kr2dCrdC1c;
        J_sl0c_ptr[sr1c] = J_sl0c_sr1c_increment;
        J_sl1c_ptr[sr1c] = J_sl1c_sr1c_increment;
        J_sl2c_ptr[sr1c] = J_sl2c_sr1c_increment;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] + Kr2dCrdC0c;
        const double J_sr1c_sr0c_increment = J_sr1c_ptr[sr0c] + Kr2dCrdC0c;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;
        J_sr1c_ptr[sr0c] = J_sr1c_sr0c_increment;
        const double J_sr0c_sr1c_increment = J_sr0c_ptr[sr1c] + Kr2dCrdC1c;
        const double J_sr1c_sr1c_increment = J_sr1c_ptr[sr1c] + Kr2dCrdC1c;
        J_sr0c_ptr[sr1c] = J_sr0c_sr1c_increment;
        J_sr1c_ptr[sr1c] = J_sr1c_sr1c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double dKrdT2 = this->dKfdT_[i2 - Ikf[1] + Ikf[9]];
        const double CF2 = C[sl0c]*C[sl1c]*C[sl2c];
        const double CR2 = C[sr0c]*C[sr1c];
        const double dqdT2 = (dKfdT2*CF2) - (dKrdT2*CR2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sl1c_T_increment = J_sl1c_ptr[this->nSpecies] - dqdT2;
        const double J_sl2c_T_increment = J_sl2c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sl1c_ptr[this->nSpecies] = J_sl1c_T_increment;
        J_sl2c_ptr[this->nSpecies] = J_sl2c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        if(i2>=this->Ikf[2] && i2<this->Ikf[3])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*CF2 - this->dKfdC_[idx2+this->Itbr[4]]*CR2;
            __m256d Wv2 = _mm256_set1_pd(W2);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv2,dMdC);
                __m256d sl0v = load256d(&J_sl0c_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1c_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1c_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2c_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2c_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0c_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0c_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1c_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1c_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->NER32AllIndex[(k+3)*6+0];
        const unsigned int sl0d = this->NER32AllIndex[(k+3)*6+1];
        const unsigned int sl1d = this->NER32AllIndex[(k+3)*6+2];
        const unsigned int sl2d = this->NER32AllIndex[(k+3)*6+3];
        const unsigned int sr0d = this->NER32AllIndex[(k+3)*6+4];
        const unsigned int sr1d = this->NER32AllIndex[(k+3)*6+5];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sl1d_ptr = &ddNdtByVdcTp[sl1d*(this->alignN)];
        double* __restrict__ J_sl2d_ptr = &ddNdtByVdcTp[sl2d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double C_sl1d = C[sl1d];
        const double C_sl0d = C[sl0d];
        const double C_sl2d = C[sl2d];
        const double Kf3dCfdC0d = Kf3*C_sl1d*C_sl2d;
        const double Kf3dCfdC1d = Kf3*C_sl0d*C_sl2d;
        const double Kf3dCfdC2d = Kf3*C_sl0d*C_sl1d;
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3dCfdC0d;
        const double J_sl1d_sl0d_increment = J_sl1d_ptr[sl0d] - Kf3dCfdC0d;
        const double J_sl2d_sl0d_increment = J_sl2d_ptr[sl0d] - Kf3dCfdC0d;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sl1d_ptr[sl0d] = J_sl1d_sl0d_increment;
        J_sl2d_ptr[sl0d] = J_sl2d_sl0d_increment;
        const double J_sl0d_sl1d_increment = J_sl0d_ptr[sl1d] - Kf3dCfdC1d;
        const double J_sl1d_sl1d_increment = J_sl1d_ptr[sl1d] - Kf3dCfdC1d;
        const double J_sl2d_sl1d_increment = J_sl2d_ptr[sl1d] - Kf3dCfdC1d;
        J_sl0d_ptr[sl1d] = J_sl0d_sl1d_increment;
        J_sl1d_ptr[sl1d] = J_sl1d_sl1d_increment;
        J_sl2d_ptr[sl1d] = J_sl2d_sl1d_increment;
        const double J_sl0d_sl2d_increment = J_sl0d_ptr[sl2d] - Kf3dCfdC2d;
        const double J_sl1d_sl2d_increment = J_sl1d_ptr[sl2d] - Kf3dCfdC2d;
        const double J_sl2d_sl2d_increment = J_sl2d_ptr[sl2d] - Kf3dCfdC2d;
        J_sl0d_ptr[sl2d] = J_sl0d_sl2d_increment;
        J_sl1d_ptr[sl2d] = J_sl1d_sl2d_increment;
        J_sl2d_ptr[sl2d] = J_sl2d_sl2d_increment;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3dCfdC0d;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3dCfdC0d;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;
        const double J_sr0d_sl1d_increment = J_sr0d_ptr[sl1d] + Kf3dCfdC1d;
        const double J_sr1d_sl1d_increment = J_sr1d_ptr[sl1d] + Kf3dCfdC1d;
        J_sr0d_ptr[sl1d] = J_sr0d_sl1d_increment;
        J_sr1d_ptr[sl1d] = J_sr1d_sl1d_increment;
        const double J_sr0d_sl2d_increment = J_sr0d_ptr[sl2d] + Kf3dCfdC2d;
        const double J_sr1d_sl2d_increment = J_sr1d_ptr[sl2d] + Kf3dCfdC2d;
        J_sr0d_ptr[sl2d] = J_sr0d_sl2d_increment;
        J_sr1d_ptr[sl2d] = J_sr1d_sl2d_increment;

        const double Kr3 = this->Kf_[i3 - Ikf[1] + Ikf[9]];
        const double C_sr1d = C[sr1d];
        const double C_sr0d = C[sr0d];
        const double Kr3dCrdC0d = -Kr3*C_sr1d;
        const double Kr3dCrdC1d = -Kr3*C_sr0d;
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] - Kr3dCrdC0d;
        const double J_sl1d_sr0d_increment = J_sl1d_ptr[sr0d] - Kr3dCrdC0d;
        const double J_sl2d_sr0d_increment = J_sl2d_ptr[sr0d] - Kr3dCrdC0d;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sl1d_ptr[sr0d] = J_sl1d_sr0d_increment;
        J_sl2d_ptr[sr0d] = J_sl2d_sr0d_increment;
        const double J_sl0d_sr1d_increment = J_sl0d_ptr[sr1d] - Kr3dCrdC1d;
        const double J_sl1d_sr1d_increment = J_sl1d_ptr[sr1d] - Kr3dCrdC1d;
        const double J_sl2d_sr1d_increment = J_sl2d_ptr[sr1d] - Kr3dCrdC1d;
        J_sl0d_ptr[sr1d] = J_sl0d_sr1d_increment;
        J_sl1d_ptr[sr1d] = J_sl1d_sr1d_increment;
        J_sl2d_ptr[sr1d] = J_sl2d_sr1d_increment;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] + Kr3dCrdC0d;
        const double J_sr1d_sr0d_increment = J_sr1d_ptr[sr0d] + Kr3dCrdC0d;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;
        J_sr1d_ptr[sr0d] = J_sr1d_sr0d_increment;
        const double J_sr0d_sr1d_increment = J_sr0d_ptr[sr1d] + Kr3dCrdC1d;
        const double J_sr1d_sr1d_increment = J_sr1d_ptr[sr1d] + Kr3dCrdC1d;
        J_sr0d_ptr[sr1d] = J_sr0d_sr1d_increment;
        J_sr1d_ptr[sr1d] = J_sr1d_sr1d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double dKrdT3 = this->dKfdT_[i3 - Ikf[1] + Ikf[9]];
        const double CF3 = C[sl0d]*C[sl1d]*C[sl2d];
        const double CR3 = C[sr0d]*C[sr1d];
        const double dqdT3 = (dKfdT3*CF3) - (dKrdT3*CR3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sl1d_T_increment = J_sl1d_ptr[this->nSpecies] - dqdT3;
        const double J_sl2d_T_increment = J_sl2d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sl1d_ptr[this->nSpecies] = J_sl1d_T_increment;
        J_sl2d_ptr[this->nSpecies] = J_sl2d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sl2d_increment = dNdtByV[sl2d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sl2d] = omega_sl2d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;

        if(i3>=this->Ikf[2] && i3<this->Ikf[3])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*CF3 - this->dKfdC_[idx3+this->Itbr[4]]*CR3;
            __m256d Wv3 = _mm256_set1_pd(W3);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv3,dMdC);
                __m256d sl0v = load256d(&J_sl0d_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1d_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1d_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2d_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2d_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0d_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0d_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1d_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1d_ptr[j],sr1v);
            }
        }
    }
    for(std::size_t k=endNER32-remainNER32; k<endNER32; k=k+1)
    {
        const unsigned int i0 = this->NER32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->NER32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->NER32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->NER32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->NER32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->NER32AllIndex[(k+0)*6+5];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sl1a_ptr = &ddNdtByVdcTp[sl1a*(this->alignN)];
        double* __restrict__ J_sl2a_ptr = &ddNdtByVdcTp[sl2a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double C_sl1a = C[sl1a];
        const double C_sl0a = C[sl0a];
        const double C_sl2a = C[sl2a];
        const double Kf0dCfdC0a = Kf0*C_sl1a*C_sl2a;
        const double Kf0dCfdC1a = Kf0*C_sl0a*C_sl2a;
        const double Kf0dCfdC2a = Kf0*C_sl0a*C_sl1a;
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl1a_sl0a_increment = J_sl1a_ptr[sl0a] - Kf0dCfdC0a;
        const double J_sl2a_sl0a_increment = J_sl2a_ptr[sl0a] - Kf0dCfdC0a;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sl1a_ptr[sl0a] = J_sl1a_sl0a_increment;
        J_sl2a_ptr[sl0a] = J_sl2a_sl0a_increment;
        const double J_sl0a_sl1a_increment = J_sl0a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl1a_sl1a_increment = J_sl1a_ptr[sl1a] - Kf0dCfdC1a;
        const double J_sl2a_sl1a_increment = J_sl2a_ptr[sl1a] - Kf0dCfdC1a;
        J_sl0a_ptr[sl1a] = J_sl0a_sl1a_increment;
        J_sl1a_ptr[sl1a] = J_sl1a_sl1a_increment;
        J_sl2a_ptr[sl1a] = J_sl2a_sl1a_increment;
        const double J_sl0a_sl2a_increment = J_sl0a_ptr[sl2a] - Kf0dCfdC2a;
        const double J_sl1a_sl2a_increment = J_sl1a_ptr[sl2a] - Kf0dCfdC2a;
        const double J_sl2a_sl2a_increment = J_sl2a_ptr[sl2a] - Kf0dCfdC2a;
        J_sl0a_ptr[sl2a] = J_sl0a_sl2a_increment;
        J_sl1a_ptr[sl2a] = J_sl1a_sl2a_increment;
        J_sl2a_ptr[sl2a] = J_sl2a_sl2a_increment;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0dCfdC0a;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0dCfdC0a;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        const double J_sr0a_sl1a_increment = J_sr0a_ptr[sl1a] + Kf0dCfdC1a;
        const double J_sr1a_sl1a_increment = J_sr1a_ptr[sl1a] + Kf0dCfdC1a;
        J_sr0a_ptr[sl1a] = J_sr0a_sl1a_increment;
        J_sr1a_ptr[sl1a] = J_sr1a_sl1a_increment;
        const double J_sr0a_sl2a_increment = J_sr0a_ptr[sl2a] + Kf0dCfdC2a;
        const double J_sr1a_sl2a_increment = J_sr1a_ptr[sl2a] + Kf0dCfdC2a;
        J_sr0a_ptr[sl2a] = J_sr0a_sl2a_increment;
        J_sr1a_ptr[sl2a] = J_sr1a_sl2a_increment;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double C_sr1a = C[sr1a];
        const double C_sr0a = C[sr0a];
        const double Kr0dCrdC0a = -Kr0*C_sr1a;
        const double Kr0dCrdC1a = -Kr0*C_sr0a;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sl1a_sr0a_increment = J_sl1a_ptr[sr0a] - Kr0dCrdC0a;
        const double J_sl2a_sr0a_increment = J_sl2a_ptr[sr0a] - Kr0dCrdC0a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sl1a_ptr[sr0a] = J_sl1a_sr0a_increment;
        J_sl2a_ptr[sr0a] = J_sl2a_sr0a_increment;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] - Kr0dCrdC1a;
        const double J_sl1a_sr1a_increment = J_sl1a_ptr[sr1a] - Kr0dCrdC1a;
        const double J_sl2a_sr1a_increment = J_sl2a_ptr[sr1a] - Kr0dCrdC1a;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sl1a_ptr[sr1a] = J_sl1a_sr1a_increment;
        J_sl2a_ptr[sr1a] = J_sl2a_sr1a_increment;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] + Kr0dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] + Kr0dCrdC0a;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] + Kr0dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] + Kr0dCrdC1a;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a]*C[sl1a]*C[sl2a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sl1a_T_increment = J_sl1a_ptr[this->nSpecies] - dqdT0;
        const double J_sl2a_T_increment = J_sl2a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sl1a_ptr[this->nSpecies] = J_sl1a_T_increment;
        J_sl2a_ptr[this->nSpecies] = J_sl2a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[2] && i0<this->Ikf[3])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0 - this->dKfdC_[idx0+this->Itbr[4]]*CR0;
            __m256d Wv0 = _mm256_set1_pd(W0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv0,dMdC);
                __m256d sl0v = load256d(&J_sl0a_ptr[j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                __m256d sl1v = load256d(&J_sl1a_ptr[j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&J_sl1a_ptr[j],sl1v);
                __m256d sl2v = load256d(&J_sl2a_ptr[j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&J_sl2a_ptr[j],sl2v);
                __m256d sr0v = load256d(&J_sr0a_ptr[j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sr0a_ptr[j],sr0v);
                __m256d sr1v = load256d(&J_sr1a_ptr[j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }
    }

    std::size_t endNER32Dup = this->NER32DupSize;
    std::size_t remainNER32Dup = endNER32Dup%4;
    for(std::size_t k=0; k<endNER32Dup-remainNER32Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->NER32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->NER32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->NER32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->NER32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->NER32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->NER32DupAllIndex[(k+0)*6+5];

        const double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C[sl1a]*C[sl2a];
        const double dCfdC1a = C[sl0a]*C[sl2a];
        const double dCfdC2a = C[sl0a]*C[sl1a];
        ddNdtByVdcTp[sl0a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl0a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl2a] += Kf0*dCfdC2a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl2a] += Kf0*dCfdC2a;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a];
        ddNdtByVdcTp[sl0a*(this->alignN)+sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl1a*(this->alignN)+sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl2a*(this->alignN)+sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl0a*(this->alignN)+sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sl1a*(this->alignN)+sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sl2a*(this->alignN)+sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sr0a*(this->alignN)+sr0a] += (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sr1a*(this->alignN)+sr0a] += (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sr0a*(this->alignN)+sr1a] += (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sr1a*(this->alignN)+sr1a] += (-Kr0*dCrdC1a);

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a]*C[sl1a]*C[sl2a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        ddNdtByVdcTp[sl0a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sl1a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sl2a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sr0a*(this->alignN)+this->nSpecies] += dqdT0;
        ddNdtByVdcTp[sr1a*(this->alignN)+this->nSpecies] += dqdT0;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        if(i0>=this->Ikf[2] && i0<this->Ikf[3])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0 - this->dKfdC_[idx0+this->Itbr[4]]*CR0;
            __m256d Wv0 = _mm256_set1_pd(W0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv0,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2a*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2a*(this->alignN)+j],sl2v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->NER32DupAllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->NER32DupAllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->NER32DupAllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->NER32DupAllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->NER32DupAllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->NER32DupAllIndex[(k+1)*6+5];

        const double Kf1 = this->Kf_[i1];
        const double dCfdC0b = C[sl1b]*C[sl2b];
        const double dCfdC1b = C[sl0b]*C[sl2b];
        const double dCfdC2b = C[sl0b]*C[sl1b];
        ddNdtByVdcTp[sl0b*(this->alignN)+sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl1b*(this->alignN)+sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl2b*(this->alignN)+sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl0b*(this->alignN)+sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sl1b*(this->alignN)+sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sl2b*(this->alignN)+sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sl0b*(this->alignN)+sl2b] -= Kf1*dCfdC2b;
        ddNdtByVdcTp[sl1b*(this->alignN)+sl2b] -= Kf1*dCfdC2b;
        ddNdtByVdcTp[sl2b*(this->alignN)+sl2b] -= Kf1*dCfdC2b;
        ddNdtByVdcTp[sr0b*(this->alignN)+sl0b] += Kf1*dCfdC0b;
        ddNdtByVdcTp[sr1b*(this->alignN)+sl0b] += Kf1*dCfdC0b;
        ddNdtByVdcTp[sr0b*(this->alignN)+sl1b] += Kf1*dCfdC1b;
        ddNdtByVdcTp[sr1b*(this->alignN)+sl1b] += Kf1*dCfdC1b;
        ddNdtByVdcTp[sr0b*(this->alignN)+sl2b] += Kf1*dCfdC2b;
        ddNdtByVdcTp[sr1b*(this->alignN)+sl2b] += Kf1*dCfdC2b;

        const double Kr1 = this->Kf_[i1 - Ikf[1] + Ikf[9]];
        const double dCrdC0b = C[sr1b];
        const double dCrdC1b = C[sr0b];
        ddNdtByVdcTp[sl0b*(this->alignN)+sr0b] -= (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sl1b*(this->alignN)+sr0b] -= (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sl2b*(this->alignN)+sr0b] -= (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sl0b*(this->alignN)+sr1b] -= (-Kr1*dCrdC1b);
        ddNdtByVdcTp[sl1b*(this->alignN)+sr1b] -= (-Kr1*dCrdC1b);
        ddNdtByVdcTp[sl2b*(this->alignN)+sr1b] -= (-Kr1*dCrdC1b);
        ddNdtByVdcTp[sr0b*(this->alignN)+sr0b] += (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sr1b*(this->alignN)+sr0b] += (-Kr1*dCrdC0b);
        ddNdtByVdcTp[sr0b*(this->alignN)+sr1b] += (-Kr1*dCrdC1b);
        ddNdtByVdcTp[sr1b*(this->alignN)+sr1b] += (-Kr1*dCrdC1b);

        const double dKfdT1 = this->dKfdT_[i1];
        const double dKrdT1 = this->dKfdT_[i1 - Ikf[1] + Ikf[9]];
        const double CF1 = C[sl0b]*C[sl1b]*C[sl2b];
        const double CR1 = C[sr0b]*C[sr1b];
        const double dqdT1 = (dKfdT1*CF1)-(dKrdT1*CR1);
        ddNdtByVdcTp[sl0b*(this->alignN)+this->nSpecies] -= dqdT1;
        ddNdtByVdcTp[sl1b*(this->alignN)+this->nSpecies] -= dqdT1;
        ddNdtByVdcTp[sl2b*(this->alignN)+this->nSpecies] -= dqdT1;
        ddNdtByVdcTp[sr0b*(this->alignN)+this->nSpecies] += dqdT1;
        ddNdtByVdcTp[sr1b*(this->alignN)+this->nSpecies] += dqdT1;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        if(i1>=this->Ikf[2] && i1<this->Ikf[3])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*CF1 - this->dKfdC_[idx1+this->Itbr[4]]*CR1;
            __m256d Wv1 = _mm256_set1_pd(W1);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv1,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2b*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2b*(this->alignN)+j],sl2v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->NER32DupAllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->NER32DupAllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->NER32DupAllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->NER32DupAllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->NER32DupAllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->NER32DupAllIndex[(k+2)*6+5];

        const double Kf2 = this->Kf_[i2];
        const double dCfdC0c = C[sl1c]*C[sl2c];
        const double dCfdC1c = C[sl0c]*C[sl2c];
        const double dCfdC2c = C[sl0c]*C[sl1c];
        ddNdtByVdcTp[sl0c*(this->alignN)+sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl1c*(this->alignN)+sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl2c*(this->alignN)+sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl0c*(this->alignN)+sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sl1c*(this->alignN)+sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sl2c*(this->alignN)+sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sl0c*(this->alignN)+sl2c] -= Kf2*dCfdC2c;
        ddNdtByVdcTp[sl1c*(this->alignN)+sl2c] -= Kf2*dCfdC2c;
        ddNdtByVdcTp[sl2c*(this->alignN)+sl2c] -= Kf2*dCfdC2c;
        ddNdtByVdcTp[sr0c*(this->alignN)+sl0c] += Kf2*dCfdC0c;
        ddNdtByVdcTp[sr1c*(this->alignN)+sl0c] += Kf2*dCfdC0c;
        ddNdtByVdcTp[sr0c*(this->alignN)+sl1c] += Kf2*dCfdC1c;
        ddNdtByVdcTp[sr1c*(this->alignN)+sl1c] += Kf2*dCfdC1c;
        ddNdtByVdcTp[sr0c*(this->alignN)+sl2c] += Kf2*dCfdC2c;
        ddNdtByVdcTp[sr1c*(this->alignN)+sl2c] += Kf2*dCfdC2c;

        const double Kr2 = this->Kf_[i2 - Ikf[1] + Ikf[9]];
        const double dCrdC0c = C[sr1c];
        const double dCrdC1c = C[sr0c];
        ddNdtByVdcTp[sl0c*(this->alignN)+sr0c] -= (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sl1c*(this->alignN)+sr0c] -= (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sl2c*(this->alignN)+sr0c] -= (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sl0c*(this->alignN)+sr1c] -= (-Kr2*dCrdC1c);
        ddNdtByVdcTp[sl1c*(this->alignN)+sr1c] -= (-Kr2*dCrdC1c);
        ddNdtByVdcTp[sl2c*(this->alignN)+sr1c] -= (-Kr2*dCrdC1c);
        ddNdtByVdcTp[sr0c*(this->alignN)+sr0c] += (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sr1c*(this->alignN)+sr0c] += (-Kr2*dCrdC0c);
        ddNdtByVdcTp[sr0c*(this->alignN)+sr1c] += (-Kr2*dCrdC1c);
        ddNdtByVdcTp[sr1c*(this->alignN)+sr1c] += (-Kr2*dCrdC1c);

        const double dKfdT2 = this->dKfdT_[i2];
        const double dKrdT2 = this->dKfdT_[i2 - Ikf[1] + Ikf[9]];
        const double CF2 = C[sl0c]*C[sl1c]*C[sl2c];
        const double CR2 = C[sr0c]*C[sr1c];
        const double dqdT2 = (dKfdT2*CF2)-(dKrdT2*CR2);
        ddNdtByVdcTp[sl0c*(this->alignN)+this->nSpecies] -= dqdT2;
        ddNdtByVdcTp[sl1c*(this->alignN)+this->nSpecies] -= dqdT2;
        ddNdtByVdcTp[sl2c*(this->alignN)+this->nSpecies] -= dqdT2;
        ddNdtByVdcTp[sr0c*(this->alignN)+this->nSpecies] += dqdT2;
        ddNdtByVdcTp[sr1c*(this->alignN)+this->nSpecies] += dqdT2;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        if(i2>=this->Ikf[2] && i2<this->Ikf[3])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*CF2 - this->dKfdC_[idx2+this->Itbr[4]]*CR2;
            __m256d Wv2 = _mm256_set1_pd(W2);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv2,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2c*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2c*(this->alignN)+j],sl2v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j],sr1v);
            }
        }

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->NER32DupAllIndex[(k+3)*6+0];
        const unsigned int sl0d = this->NER32DupAllIndex[(k+3)*6+1];
        const unsigned int sl1d = this->NER32DupAllIndex[(k+3)*6+2];
        const unsigned int sl2d = this->NER32DupAllIndex[(k+3)*6+3];
        const unsigned int sr0d = this->NER32DupAllIndex[(k+3)*6+4];
        const unsigned int sr1d = this->NER32DupAllIndex[(k+3)*6+5];

        const double Kf3 = this->Kf_[i3];
        const double dCfdC0d = C[sl1d]*C[sl2d];
        const double dCfdC1d = C[sl0d]*C[sl2d];
        const double dCfdC2d = C[sl0d]*C[sl1d];
        ddNdtByVdcTp[sl0d*(this->alignN)+sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl1d*(this->alignN)+sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl2d*(this->alignN)+sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl0d*(this->alignN)+sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sl1d*(this->alignN)+sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sl2d*(this->alignN)+sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sl0d*(this->alignN)+sl2d] -= Kf3*dCfdC2d;
        ddNdtByVdcTp[sl1d*(this->alignN)+sl2d] -= Kf3*dCfdC2d;
        ddNdtByVdcTp[sl2d*(this->alignN)+sl2d] -= Kf3*dCfdC2d;
        ddNdtByVdcTp[sr0d*(this->alignN)+sl0d] += Kf3*dCfdC0d;
        ddNdtByVdcTp[sr1d*(this->alignN)+sl0d] += Kf3*dCfdC0d;
        ddNdtByVdcTp[sr0d*(this->alignN)+sl1d] += Kf3*dCfdC1d;
        ddNdtByVdcTp[sr1d*(this->alignN)+sl1d] += Kf3*dCfdC1d;
        ddNdtByVdcTp[sr0d*(this->alignN)+sl2d] += Kf3*dCfdC2d;
        ddNdtByVdcTp[sr1d*(this->alignN)+sl2d] += Kf3*dCfdC2d;

        const double Kr3 = this->Kf_[i3 - Ikf[1] + Ikf[9]];
        const double dCrdC0d = C[sr1d];
        const double dCrdC1d = C[sr0d];
        ddNdtByVdcTp[sl0d*(this->alignN)+sr0d] -= (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sl1d*(this->alignN)+sr0d] -= (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sl2d*(this->alignN)+sr0d] -= (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sl0d*(this->alignN)+sr1d] -= (-Kr3*dCrdC1d);
        ddNdtByVdcTp[sl1d*(this->alignN)+sr1d] -= (-Kr3*dCrdC1d);
        ddNdtByVdcTp[sl2d*(this->alignN)+sr1d] -= (-Kr3*dCrdC1d);
        ddNdtByVdcTp[sr0d*(this->alignN)+sr0d] += (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sr1d*(this->alignN)+sr0d] += (-Kr3*dCrdC0d);
        ddNdtByVdcTp[sr0d*(this->alignN)+sr1d] += (-Kr3*dCrdC1d);
        ddNdtByVdcTp[sr1d*(this->alignN)+sr1d] += (-Kr3*dCrdC1d);

        const double dKfdT3 = this->dKfdT_[i3];
        const double dKrdT3 = this->dKfdT_[i3 - Ikf[1] + Ikf[9]];
        const double CF3 = C[sl0d]*C[sl1d]*C[sl2d];
        const double CR3 = C[sr0d]*C[sr1d];
        const double dqdT3 = (dKfdT3*CF3)-(dKrdT3*CR3);
        ddNdtByVdcTp[sl0d*(this->alignN)+this->nSpecies] -= dqdT3;
        ddNdtByVdcTp[sl1d*(this->alignN)+this->nSpecies] -= dqdT3;
        ddNdtByVdcTp[sl2d*(this->alignN)+this->nSpecies] -= dqdT3;
        ddNdtByVdcTp[sr0d*(this->alignN)+this->nSpecies] += dqdT3;
        ddNdtByVdcTp[sr1d*(this->alignN)+this->nSpecies] += dqdT3;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

        if(i3>=this->Ikf[2] && i3<this->Ikf[3])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*CF3 - this->dKfdC_[idx3+this->Itbr[4]]*CR3;
            __m256d Wv3 = _mm256_set1_pd(W3);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv3,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2d*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2d*(this->alignN)+j],sl2v);
                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j],sr0v);
                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j],sr1v);
            }
        }
    }
    for(std::size_t k=endNER32Dup-remainNER32Dup; k<endNER32Dup; k=k+1)
    {
        const unsigned int i0 = this->NER32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->NER32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->NER32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->NER32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->NER32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->NER32DupAllIndex[(k+0)*6+5];

        const double Kf0 = this->Kf_[i0];
        const double dCfdC0a = C[sl1a]*C[sl2a];
        const double dCfdC1a = C[sl0a]*C[sl2a];
        const double dCfdC2a = C[sl0a]*C[sl1a];
        ddNdtByVdcTp[sl0a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl0a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sl1a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sl2a*(this->alignN)+sl2a] -= Kf0*dCfdC2a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+sl2a] += Kf0*dCfdC2a;
        ddNdtByVdcTp[sr1a*(this->alignN)+sl2a] += Kf0*dCfdC2a;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a];
        ddNdtByVdcTp[sl0a*(this->alignN)+sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl1a*(this->alignN)+sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl2a*(this->alignN)+sr0a] -= (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sl0a*(this->alignN)+sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sl1a*(this->alignN)+sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sl2a*(this->alignN)+sr1a] -= (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sr0a*(this->alignN)+sr0a] += (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sr1a*(this->alignN)+sr0a] += (-Kr0*dCrdC0a);
        ddNdtByVdcTp[sr0a*(this->alignN)+sr1a] += (-Kr0*dCrdC1a);
        ddNdtByVdcTp[sr1a*(this->alignN)+sr1a] += (-Kr0*dCrdC1a);

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a]*C[sl1a]*C[sl2a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        ddNdtByVdcTp[sl0a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sl1a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sl2a*(this->alignN)+this->nSpecies] -= dqdT0;
        ddNdtByVdcTp[sr0a*(this->alignN)+this->nSpecies] += dqdT0;
        ddNdtByVdcTp[sr1a*(this->alignN)+this->nSpecies] += dqdT0;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        if(i0>=this->Ikf[2] && i0<this->Ikf[3])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0 - this->dKfdC_[idx0+this->Itbr[4]]*CR0;
            __m256d Wv0 = _mm256_set1_pd(W0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(Wv0,dMdC);
                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);
                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);
                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2a*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2a*(this->alignN)+j],sl2v);
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
