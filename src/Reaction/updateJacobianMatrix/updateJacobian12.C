/*---------------------------------------------------------------------------*\
  Description
      Update the net rates of production (dNdtByV) and the molar-concentration
      based Jacobian matrix (ddNdtByVdcTp) for one-reactant / two-product
      (1-2) reactions, e.g. A -> B + B and A -> B + C.

      The Jacobian matrix is stored row-wise with one column per species
      plus one extra column for the temperature derivative (column index
      nSpecies).

      RR : reversible reaction; the reverse rate constant is obtained from
           the equilibrium constant
      IR : irreversible reaction; the reverse rate constant is zero
      NER: non-equilibrium reaction; the reverse rate constant is computed
           with the Arrhenius form instead of the equilibrium constant

  Functions
      JF12RR : reversible 1-2 reactions, incl. duplicate products (A=B+B)
      JF12IR : irreversible 1-2 reactions, incl. duplicate products (A=B+B)
      JF12NER: non-equilibrium 1-2 reactions

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
void FastChemistry::OptReaction::JF12RR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
) const noexcept
{
    const double PowPByRTSumVki = this->Pow_pByRT_SumVki[1];
    const double invKcLimiter = FastChemistry::invKcLimiter;
    std::size_t endRR12 = this->RR12Size;
    std::size_t remainRR12 = endRR12 % 4;
    for(std::size_t k=0; k<endRR12 - remainRR12; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->RR12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->RR12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->RR12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->RR12AllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0*dCrdC0a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0*dCrdC1a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0*dCrdC1a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a]) - (dBdT[sl0a]);
        const double dKcdTByKc0 = sumVdBdT0 - invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CFa = C[sl0a];
        const double CRa = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CFa) - (dKrdT0*CRa);

        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;


        const double q0 = (Kf0*CFa) - (Kr0*CRa);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*(CFa-invKc0*CRa);
            __m256d Wv0 = _mm256_set1_pd(W0);
            unsigned int remain80 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain80;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv0);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+4]),WdMdC4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
            }
            if(remain80==4)
            {
                unsigned int j = this->AlignSpecies-remain80;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR12AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->RR12AllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->RR12AllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->RR12AllIndex[(k+1)*4+3];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;

        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double dCrdC0b = C[sr1b];
        const double dCrdC1b = C[sr0b];
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1*dCrdC0b;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1*dCrdC0b;
        const double J_sr1b_sr0b_increment = J_sr1b_ptr[sr0b] - Kr1*dCrdC0b;
        const double J_sl0b_sr1b_increment = J_sl0b_ptr[sr1b] + Kr1*dCrdC1b;
        const double J_sr0b_sr1b_increment = J_sr0b_ptr[sr1b] - Kr1*dCrdC1b;
        const double J_sr1b_sr1b_increment = J_sr1b_ptr[sr1b] - Kr1*dCrdC1b;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;
        J_sr1b_ptr[sr0b] = J_sr1b_sr0b_increment;
        J_sl0b_ptr[sr1b] = J_sl0b_sr1b_increment;
        J_sr0b_ptr[sr1b] = J_sr0b_sr1b_increment;
        J_sr1b_ptr[sr1b] = J_sr1b_sr1b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double sumVdBdT1 = (dBdT[sr0b] + dBdT[sr1b]) - (dBdT[sl0b]);
        const double dKcdTByKc1 = sumVdBdT1 - invT;
        const double dKrdT1 = (dKfdT1*invKc1 - (invKc1 < invKcLimiter ? Kr1*dKcdTByKc1 : 0));
        const double CFb = C[sl0b];
        const double CRb = C[sr0b]*C[sr1b];
        const double dqdT1 = (dKfdT1*CFb) - (dKrdT1*CRb);

        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;


        const double q1 = (Kf1*CFb) - (Kr1*CRb);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx]*(CFb-invKc1*CRb);
            __m256d Wv = _mm256_set1_pd(W0);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+4]),WdMdC4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+4]),WdMdC4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+4]),WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sr1b_ptr[j+0],sr1v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
                store256d(&J_sr0b_ptr[j+4],sr0v4);
                store256d(&J_sr1b_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0b_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0b_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1b_ptr[j]),WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sr0b_ptr[j],sr0v);
                store256d(&J_sr1b_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR12AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->RR12AllIndex[(k+2)*4+1];
        const unsigned int sr0c = this->RR12AllIndex[(k+2)*4+2];
        const unsigned int sr1c = this->RR12AllIndex[(k+2)*4+3];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;

        const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*(ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double dCrdC0c = C[sr1c];
        const double dCrdC1c = C[sr0c];
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2*dCrdC0c;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2*dCrdC0c;
        const double J_sr1c_sr0c_increment = J_sr1c_ptr[sr0c] - Kr2*dCrdC0c;
        const double J_sl0c_sr1c_increment = J_sl0c_ptr[sr1c] + Kr2*dCrdC1c;
        const double J_sr0c_sr1c_increment = J_sr0c_ptr[sr1c] - Kr2*dCrdC1c;
        const double J_sr1c_sr1c_increment = J_sr1c_ptr[sr1c] - Kr2*dCrdC1c;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;
        J_sr1c_ptr[sr0c] = J_sr1c_sr0c_increment;
        J_sl0c_ptr[sr1c] = J_sl0c_sr1c_increment;
        J_sr0c_ptr[sr1c] = J_sr0c_sr1c_increment;
        J_sr1c_ptr[sr1c] = J_sr1c_sr1c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double sumVdBdT2 = (dBdT[sr0c] + dBdT[sr1c]) - (dBdT[sl0c]);
        const double dKcdTByKc2 = sumVdBdT2 - invT;
        const double dKrdT2 = (dKfdT2*invKc2 - (invKc2 < invKcLimiter ? Kr2*dKcdTByKc2 : 0));
        const double CFc = C[sl0c];
        const double CRc = C[sr0c]*C[sr1c];
        const double dqdT2 = (dKfdT2*CFc) - (dKrdT2*CRc);

        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;


        const double q2 = (Kf2*CFc) - (Kr2*CRc);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx]*(CFc-invKc2*CRc);
            __m256d Wv = _mm256_set1_pd(W0);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+4]),WdMdC4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+4]),WdMdC4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+4]),WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sr1c_ptr[j+0],sr1v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
                store256d(&J_sr0c_ptr[j+4],sr0v4);
                store256d(&J_sr1c_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0c_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0c_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1c_ptr[j]),WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sr0c_ptr[j],sr0v);
                store256d(&J_sr1c_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR12AllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->RR12AllIndex[(k+3)*4+1];
        const unsigned int sr0d = this->RR12AllIndex[(k+3)*4+2];
        const unsigned int sr1d = this->RR12AllIndex[(k+3)*4+3];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;

        const double invKp3 = (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d])*(ExpNegGbyRT[sl0d]);
        double invKc3 = invKp3*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double dCrdC0d = C[sr1d];
        const double dCrdC1d = C[sr0d];
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3*dCrdC0d;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3*dCrdC0d;
        const double J_sr1d_sr0d_increment = J_sr1d_ptr[sr0d] - Kr3*dCrdC0d;
        const double J_sl0d_sr1d_increment = J_sl0d_ptr[sr1d] + Kr3*dCrdC1d;
        const double J_sr0d_sr1d_increment = J_sr0d_ptr[sr1d] - Kr3*dCrdC1d;
        const double J_sr1d_sr1d_increment = J_sr1d_ptr[sr1d] - Kr3*dCrdC1d;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;
        J_sr1d_ptr[sr0d] = J_sr1d_sr0d_increment;
        J_sl0d_ptr[sr1d] = J_sl0d_sr1d_increment;
        J_sr0d_ptr[sr1d] = J_sr0d_sr1d_increment;
        J_sr1d_ptr[sr1d] = J_sr1d_sr1d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double sumVdBdT3 = (dBdT[sr0d] + dBdT[sr1d]) - (dBdT[sl0d]);
        const double dKcdTByKc3 = sumVdBdT3 - invT;
        const double dKrdT3 = (dKfdT3*invKc3 - (invKc3 < invKcLimiter ? Kr3*dKcdTByKc3 : 0));
        const double CFd = C[sl0d];
        const double CRd = C[sr0d]*C[sr1d];
        const double dqdT3 = (dKfdT3*CFd) - (dKrdT3*CRd);

        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;


        const double q3 = (Kf3*CFd) - (Kr3*CRd);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx]*(CFd-invKc3*CRd);
            __m256d Wv = _mm256_set1_pd(W0);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+4]),WdMdC4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+4]),WdMdC4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+4]),WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sr1d_ptr[j+0],sr1v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
                store256d(&J_sr0d_ptr[j+4],sr0v4);
                store256d(&J_sr1d_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0d_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0d_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1d_ptr[j]),WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sr0d_ptr[j],sr0v);
                store256d(&J_sr1d_ptr[j],sr1v);
            }
        }
    }
    for(std::size_t k=endRR12-remainRR12; k<endRR12; k=k+1)
    {
        const unsigned int i0 = this->RR12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->RR12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->RR12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->RR12AllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0*dCrdC0a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0*dCrdC1a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0*dCrdC1a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a]) - (dBdT[sl0a]);
        const double dKcdTByKc0 = sumVdBdT0 - invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CFa = C[sl0a];
        const double CRa = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CFa) - (dKrdT0*CRa);

        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;


        const double q0 = (Kf0*CFa) - (Kr0*CRa);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx0]*(CFa-invKc0*CRa);
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+4]),WdMdC4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }
    }

    std::size_t endRRDup12 = this->RR12DupSize;
    std::size_t remainRRDup12 = endRRDup12 % 4;
    for(std::size_t k=0; k<endRRDup12-remainRRDup12; k=k+4)
    {

        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->RR12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR12DupAllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0 + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double C_sr0a = C[sr0a];

        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + (Kr0*C_sr0a + Kr0*C_sr0a);
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - 4*Kr0*C_sr0a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr0a]) - (dBdT[sl0a]);
        const double dKcdTByKc0 = sumVdBdT0 - invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CFa = C[sl0a];
        const double CRa = C_sr0a*C_sr0a;

        const double q0 = (Kf0*CFa) - (Kr0*CRa);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const double dqdT0 = (dKfdT0*CFa) - (dKrdT0*CRa);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0 + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*(CFa-invKc0*CRa);
            __m256d Wv0 = _mm256_set1_pd(W0);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv0);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }


        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR12DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->RR12DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->RR12DupAllIndex[(k+1)*3+2];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1 + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;

        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double C_sr0b = C[sr0b];

        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + (Kr1*C_sr0b + Kr1*C_sr0b);
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - 4*Kr1*C_sr0b;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double sumVdBdT1 = (dBdT[sr0b] + dBdT[sr0b]) - (dBdT[sl0b]);
        const double dKcdTByKc1 = sumVdBdT1 - invT;
        const double dKrdT1 = (dKfdT1*invKc1 - (invKc1 < invKcLimiter ? Kr1*dKcdTByKc1 : 0));
        const double CFb = C[sl0b];
        const double CRb = C_sr0b*C_sr0b;

        const double q1 = (Kf1*CFb) - (Kr1*CRb);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        const double dqdT1 = (dKfdT1*CFb) - (dKrdT1*CRb);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1 + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*(CFb-invKc1*CRb);
            __m256d Wv1 = _mm256_set1_pd(W1);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv1);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv1);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+4]),WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+4]),WdMdC4);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sr0b_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv1);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0b_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0b_ptr[j]),WdMdC);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sr0b_ptr[j],sr0v);
            }
        }


        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR12DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->RR12DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->RR12DupAllIndex[(k+2)*3+2];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2 + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;

        const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr0c])*(ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double C_sr0c = C[sr0c];

        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + (Kr2*C_sr0c + Kr2*C_sr0c);
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - 4*Kr2*C_sr0c;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double sumVdBdT2 = (dBdT[sr0c] + dBdT[sr0c]) - (dBdT[sl0c]);
        const double dKcdTByKc2 = sumVdBdT2 - invT;
        const double dKrdT2 = (dKfdT2*invKc2 - (invKc2 < invKcLimiter ? Kr2*dKcdTByKc2 : 0));
        const double CFc = C[sl0c];
        const double CRc = C_sr0c*C_sr0c;

        const double q2 = (Kf2*CFc) - (Kr2*CRc);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2 + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        const double dqdT2 = (dKfdT2*CFc) - (dKrdT2*CRc);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2 + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*(CFc-invKc2*CRc);
            __m256d Wv2 = _mm256_set1_pd(W2);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv2);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv2);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+4]),WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+4]),WdMdC4);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sr0c_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv2);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0c_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0c_ptr[j]),WdMdC);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sr0c_ptr[j],sr0v);
            }
        }


        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR12DupAllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->RR12DupAllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->RR12DupAllIndex[(k+3)*3+2];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3 + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;

        const double invKp3 = (invNegGstdByRT[sr0d]*invNegGstdByRT[sr0d])*(ExpNegGbyRT[sl0d]);
        double invKc3 = invKp3*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double C_sr0d = C[sr0d];

        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + (Kr3*C_sr0d + Kr3*C_sr0d);
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - 4*Kr3*C_sr0d;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double sumVdBdT3 = (dBdT[sr0d] + dBdT[sr0d]) - (dBdT[sl0d]);
        const double dKcdTByKc3 = sumVdBdT3 - invT;
        const double dKrdT3 = (dKfdT3*invKc3 - (invKc3 < invKcLimiter ? Kr3*dKcdTByKc3 : 0));
        const double CFd = C[sl0d];
        const double CRd = C_sr0d*C_sr0d;

        const double q3 = (Kf3*CFd) - (Kr3*CRd);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3 + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;

        const double dqdT3 = (dKfdT3*CFd) - (dKrdT3*CRd);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3 + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*(CFd-invKc3*CRd);
            __m256d Wv3 = _mm256_set1_pd(W3);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv3);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv3);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+4]),WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+4]),WdMdC4);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sr0d_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv3);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0d_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0d_ptr[j]),WdMdC);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sr0d_ptr[j],sr0v);
            }
        }

    }
    for(std::size_t k=endRRDup12-remainRRDup12; k<endRRDup12; k=k+1)
    {

        const unsigned int i0 = this->RR12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR12DupAllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0 + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double C_sr0a = C[sr0a];

        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + (Kr0*C_sr0a + Kr0*C_sr0a);
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - 4*Kr0*C_sr0a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr0a]) - (dBdT[sl0a]);
        const double dKcdTByKc0 = sumVdBdT0 - invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CFa = C[sl0a];
        const double CRa = C_sr0a*C_sr0a;

        const double q0 = (Kf0*CFa) - (Kr0*CRa);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const double dqdT0 = (dKfdT0*CFa) - (dKrdT0*CRa);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0 + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*(CFa-invKc0*CRa);
            __m256d Wv0 = _mm256_set1_pd(W0);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv0);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv0);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }

    }
}

void FastChemistry::OptReaction::JF12IR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
) const noexcept
{

    std::size_t endIR12 = this->IR12Size;
    std::size_t remainIR12 = endIR12 % 4;
    for(std::size_t k=0; k<endIR12-remainIR12; k=k+4)
    {

        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->IR12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->IR12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->IR12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->IR12AllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C[sl0a];
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF0);
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR12AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->IR12AllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->IR12AllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->IR12AllIndex[(k+1)*4+3];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double CF1 = C[sl0b];
        const double dqdT1 = (dKfdT1*CF1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;

        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF1);
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+4]),WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+4]),WdMdC4);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sr0b_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+4]),WdMdC4);
                store256d(&J_sr1b_ptr[j+0],sr1v0);
                store256d(&J_sr1b_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+0]),WdMdC0);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                store256d(&J_sr1b_ptr[j+0],sr1v0);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR12AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->IR12AllIndex[(k+2)*4+1];
        const unsigned int sr0c = this->IR12AllIndex[(k+2)*4+2];
        const unsigned int sr1c = this->IR12AllIndex[(k+2)*4+3];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double CF2 = C[sl0c];
        const double dqdT2 = (dKfdT2*CF2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;

        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF2);
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+4]),WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+4]),WdMdC4);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sr0c_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+4]),WdMdC4);
                store256d(&J_sr1c_ptr[j+0],sr1v0);
                store256d(&J_sr1c_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+0]),WdMdC0);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                store256d(&J_sr1c_ptr[j+0],sr1v0);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR12AllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->IR12AllIndex[(k+3)*4+1];
        const unsigned int sr0d = this->IR12AllIndex[(k+3)*4+2];
        const unsigned int sr1d = this->IR12AllIndex[(k+3)*4+3];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double CF3 = C[sl0d];
        const double dqdT3 = (dKfdT3*CF3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;

        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF3);
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+4]),WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+4]),WdMdC4);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sr0d_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+4]),WdMdC4);
                store256d(&J_sr1d_ptr[j+0],sr1v0);
                store256d(&J_sr1d_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+0]),WdMdC0);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                store256d(&J_sr1d_ptr[j+0],sr1v0);
            }
        }
    }
    for(std::size_t k=endIR12-remainIR12; k<endIR12; k=k+1)
    {

        const unsigned int i0 = this->IR12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->IR12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->IR12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->IR12AllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C[sl0a];
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CF0);
            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
            }
        }
    }

    std::size_t endIR12Dup = this->IR12DupSize;
    std::size_t remainIR12Dup = endIR12Dup % 4;
    for(std::size_t k=0; k<endIR12Dup-remainIR12Dup; k=k+4)
    {

        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->IR12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR12DupAllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];


        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0 + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C[sl0a];
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0 + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;



        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];

            const double W = this->dKfdC_[idx]*(CF0);

            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);

                __m256d sl0v0 = load256d(&J_sl0a_ptr[j+0]);
                __m256d sl0v4 = load256d(&J_sl0a_ptr[j+4]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                sl0v4 = _mm256_sub_pd(sl0v4,WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);

                __m256d sr0v0 = load256d(&J_sr0a_ptr[j+0]);
                __m256d sr0v4 = load256d(&J_sr0a_ptr[j+4]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);

                __m256d sl0v0 = load256d(&J_sl0a_ptr[j+0]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);

                __m256d sr0v0 = load256d(&J_sr0a_ptr[j+0]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR12DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->IR12DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->IR12DupAllIndex[(k+1)*3+2];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];


        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1 + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double CF1 = C[sl0b];
        const double dqdT1 = (dKfdT1*CF1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1 + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;

        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;



        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx = i1 - this->Ikf[2];

            const double W = this->dKfdC_[idx]*(CF1);

            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);

                __m256d sl0v0 = load256d(&J_sl0b_ptr[j+0]);
                __m256d sl0v4 = load256d(&J_sl0b_ptr[j+4]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                sl0v4 = _mm256_sub_pd(sl0v4,WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);

                __m256d sr0v0 = load256d(&J_sr0b_ptr[j+0]);
                __m256d sr0v4 = load256d(&J_sr0b_ptr[j+4]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sr0b_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);

                __m256d sl0v0 = load256d(&J_sl0b_ptr[j+0]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                store256d(&J_sl0b_ptr[j+0],sl0v0);

                __m256d sr0v0 = load256d(&J_sr0b_ptr[j+0]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR12DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->IR12DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->IR12DupAllIndex[(k+2)*3+2];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];


        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2 + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double CF2 = C[sl0c];
        const double dqdT2 = (dKfdT2*CF2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2 + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;

        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2 + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;



        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx = i2 - this->Ikf[2];

            const double W = this->dKfdC_[idx]*(CF2);

            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);

                __m256d sl0v0 = load256d(&J_sl0c_ptr[j+0]);
                __m256d sl0v4 = load256d(&J_sl0c_ptr[j+4]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                sl0v4 = _mm256_sub_pd(sl0v4,WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);

                __m256d sr0v0 = load256d(&J_sr0c_ptr[j+0]);
                __m256d sr0v4 = load256d(&J_sr0c_ptr[j+4]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sr0c_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);

                __m256d sl0v0 = load256d(&J_sl0c_ptr[j+0]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                store256d(&J_sl0c_ptr[j+0],sl0v0);

                __m256d sr0v0 = load256d(&J_sr0c_ptr[j+0]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR12DupAllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->IR12DupAllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->IR12DupAllIndex[(k+3)*3+2];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];


        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3 + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double CF3 = C[sl0d];
        const double dqdT3 = (dKfdT3*CF3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3 + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;

        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3 + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;



        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx = i3 - this->Ikf[2];

            const double W = this->dKfdC_[idx]*(CF3);

            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);

                __m256d sl0v0 = load256d(&J_sl0d_ptr[j+0]);
                __m256d sl0v4 = load256d(&J_sl0d_ptr[j+4]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                sl0v4 = _mm256_sub_pd(sl0v4,WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);

                __m256d sr0v0 = load256d(&J_sr0d_ptr[j+0]);
                __m256d sr0v4 = load256d(&J_sr0d_ptr[j+4]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sr0d_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);

                __m256d sl0v0 = load256d(&J_sl0d_ptr[j+0]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                store256d(&J_sl0d_ptr[j+0],sl0v0);

                __m256d sr0v0 = load256d(&J_sr0d_ptr[j+0]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
            }
        }
    }
    for(std::size_t k=endIR12Dup-remainIR12Dup; k<endIR12Dup; k=k+1)
    {

        const unsigned int i0 = this->IR12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR12DupAllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];


        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0 + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C[sl0a];
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0 + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;



        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];

            const double W = this->dKfdC_[idx]*(CF0);

            __m256d Wv = _mm256_set1_pd(W);
            unsigned int remain8 = this->AlignSpecies%8;
            for(unsigned int j = 0; j < this->AlignSpecies-remain8;j=j+8)
            {
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d dMdC4 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+4]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);

                __m256d sl0v0 = load256d(&J_sl0a_ptr[j+0]);
                __m256d sl0v4 = load256d(&J_sl0a_ptr[j+4]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                sl0v4 = _mm256_sub_pd(sl0v4,WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);

                __m256d sr0v0 = load256d(&J_sr0a_ptr[j+0]);
                __m256d sr0v4 = load256d(&J_sr0a_ptr[j+4]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);

                __m256d sl0v0 = load256d(&J_sl0a_ptr[j+0]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);

                __m256d sr0v0 = load256d(&J_sr0a_ptr[j+0]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
            }
        }
    }
}


void  FastChemistry::OptReaction::JF12NER
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t endNER12 = this->NER12Size;
    std::size_t remainNER12 = endNER12 % 4;

    for(std::size_t k=0; k<endNER12-remainNER12; k=k+4)
    {

        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->NER12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->NER12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->NER12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->NER12AllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];

        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0*dCrdC0a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0*dCrdC1a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0*dCrdC1a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
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
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->NER12AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->NER12AllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->NER12AllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->NER12AllIndex[(k+1)*4+3];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;

        const double Kr1 = this->Kf_[i1 - Ikf[1] + Ikf[9]];

        const double dCrdC0b = C[sr1b];
        const double dCrdC1b = C[sr0b];
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1*dCrdC0b;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1*dCrdC0b;
        const double J_sr1b_sr0b_increment = J_sr1b_ptr[sr0b] - Kr1*dCrdC0b;
        const double J_sl0b_sr1b_increment = J_sl0b_ptr[sr1b] + Kr1*dCrdC1b;
        const double J_sr0b_sr1b_increment = J_sr0b_ptr[sr1b] - Kr1*dCrdC1b;
        const double J_sr1b_sr1b_increment = J_sr1b_ptr[sr1b] - Kr1*dCrdC1b;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;
        J_sr1b_ptr[sr0b] = J_sr1b_sr0b_increment;
        J_sl0b_ptr[sr1b] = J_sl0b_sr1b_increment;
        J_sr0b_ptr[sr1b] = J_sr0b_sr1b_increment;
        J_sr1b_ptr[sr1b] = J_sr1b_sr1b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double dKrdT1 = this->dKfdT_[i1 - Ikf[1] + Ikf[9]];
        const double CF1 = C[sl0b];
        const double CR1 = C[sr0b]*C[sr1b];
        const double dqdT1 = (dKfdT1*CF1)-(dKrdT1*CR1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
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
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+4]),WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+4]),WdMdC4);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sr0b_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+4]),WdMdC4);
                store256d(&J_sr1b_ptr[j+0],sr1v0);
                store256d(&J_sr1b_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0b_ptr[j+0]),WdMdC0);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                store256d(&J_sr1b_ptr[j+0],sr1v0);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->NER12AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->NER12AllIndex[(k+2)*4+1];
        const unsigned int sr0c = this->NER12AllIndex[(k+2)*4+2];
        const unsigned int sr1c = this->NER12AllIndex[(k+2)*4+3];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;

        const double Kr2 = this->Kf_[i2 - Ikf[1] + Ikf[9]];

        const double dCrdC0c = C[sr1c];
        const double dCrdC1c = C[sr0c];
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2*dCrdC0c;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2*dCrdC0c;
        const double J_sr1c_sr0c_increment = J_sr1c_ptr[sr0c] - Kr2*dCrdC0c;
        const double J_sl0c_sr1c_increment = J_sl0c_ptr[sr1c] + Kr2*dCrdC1c;
        const double J_sr0c_sr1c_increment = J_sr0c_ptr[sr1c] - Kr2*dCrdC1c;
        const double J_sr1c_sr1c_increment = J_sr1c_ptr[sr1c] - Kr2*dCrdC1c;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;
        J_sr1c_ptr[sr0c] = J_sr1c_sr0c_increment;
        J_sl0c_ptr[sr1c] = J_sl0c_sr1c_increment;
        J_sr0c_ptr[sr1c] = J_sr0c_sr1c_increment;
        J_sr1c_ptr[sr1c] = J_sr1c_sr1c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double dKrdT2 = this->dKfdT_[i2 - Ikf[1] + Ikf[9]];
        const double CF2 = C[sl0c];
        const double CR2 = C[sr0c]*C[sr1c];
        const double dqdT2 = (dKfdT2*CF2)-(dKrdT2*CR2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
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
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+4]),WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+4]),WdMdC4);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sr0c_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+4]),WdMdC4);
                store256d(&J_sr1c_ptr[j+0],sr1v0);
                store256d(&J_sr1c_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0c_ptr[j+0]),WdMdC0);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                store256d(&J_sr1c_ptr[j+0],sr1v0);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->NER12AllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->NER12AllIndex[(k+3)*4+1];
        const unsigned int sr0d = this->NER12AllIndex[(k+3)*4+2];
        const unsigned int sr1d = this->NER12AllIndex[(k+3)*4+3];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;

        const double Kr3 = this->Kf_[i3 - Ikf[1] + Ikf[9]];

        const double dCrdC0d = C[sr1d];
        const double dCrdC1d = C[sr0d];
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3*dCrdC0d;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3*dCrdC0d;
        const double J_sr1d_sr0d_increment = J_sr1d_ptr[sr0d] - Kr3*dCrdC0d;
        const double J_sl0d_sr1d_increment = J_sl0d_ptr[sr1d] + Kr3*dCrdC1d;
        const double J_sr0d_sr1d_increment = J_sr0d_ptr[sr1d] - Kr3*dCrdC1d;
        const double J_sr1d_sr1d_increment = J_sr1d_ptr[sr1d] - Kr3*dCrdC1d;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;
        J_sr1d_ptr[sr0d] = J_sr1d_sr0d_increment;
        J_sl0d_ptr[sr1d] = J_sl0d_sr1d_increment;
        J_sr0d_ptr[sr1d] = J_sr0d_sr1d_increment;
        J_sr1d_ptr[sr1d] = J_sr1d_sr1d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double dKrdT3 = this->dKfdT_[i3 - Ikf[1] + Ikf[9]];
        const double CF3 = C[sl0d];
        const double CR3 = C[sr0d]*C[sr1d];
        const double dqdT3 = (dKfdT3*CF3)-(dKrdT3*CR3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
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
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+4]),WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+4]),WdMdC4);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sr0d_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+4]),WdMdC4);
                store256d(&J_sr1d_ptr[j+0],sr1v0);
                store256d(&J_sr1d_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0d_ptr[j+0]),WdMdC0);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                store256d(&J_sr1d_ptr[j+0],sr1v0);
            }
        }

    }
    for(std::size_t k=endNER12-remainNER12; k<endNER12; k=k+1)
    {

        const unsigned int i0 = this->NER12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->NER12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->NER12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->NER12AllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];

        const double dCrdC0a = C[sr1a];
        const double dCrdC1a = C[sr0a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0*dCrdC0a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0*dCrdC1a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0*dCrdC1a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a];
        const double CR0 = C[sr0a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
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
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                __m256d sl0v4 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d sl0v0 = _mm256_sub_pd(load256d(&J_sl0a_ptr[j+0]),WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
            }
        }

    }

    std::size_t endNERDup12 = this->NER12DupSize;
    std::size_t remainNERDup12 = endNERDup12 % 4;

    for(std::size_t k=0; k<endNERDup12-remainNERDup12; k=k+4)
    {

        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->NER12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER12DupAllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0 + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];

        const double dCrdC0a = C[sr0a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + (Kr0*dCrdC0a + Kr0*dCrdC0a);
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - 4*Kr0*dCrdC0a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a];
        const double CR0 = C[sr0a]*C[sr0a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0 + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

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
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);

                __m256d sl0v0 = load256d(&J_sl0a_ptr[j+0]);
                __m256d sl0v4 = load256d(&J_sl0a_ptr[j+4]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                sl0v4 = _mm256_sub_pd(sl0v4,WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);

                __m256d sr0v0 = load256d(&J_sr0a_ptr[j+0]);
                __m256d sr0v4 = load256d(&J_sr0a_ptr[j+4]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);

                __m256d sl0v0 = load256d(&J_sl0a_ptr[j+0]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);

                __m256d sr0v0 = load256d(&J_sr0a_ptr[j+0]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->NER12DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->NER12DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->NER12DupAllIndex[(k+1)*3+2];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1 + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;

        const double Kr1 = this->Kf_[i1 - Ikf[1] + Ikf[9]];

        const double dCrdC0b = C[sr0b];
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + (Kr1*dCrdC0b + Kr1*dCrdC0b);
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - 4*Kr1*dCrdC0b;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double dKrdT1 = this->dKfdT_[i1 - Ikf[1] + Ikf[9]];
        const double CF1 = C[sl0b];
        const double CR1 = C[sr0b]*C[sr0b];
        const double dqdT1 = (dKfdT1*CF1)-(dKrdT1*CR1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1 + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

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
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);

                __m256d sl0v0 = load256d(&J_sl0b_ptr[j+0]);
                __m256d sl0v4 = load256d(&J_sl0b_ptr[j+4]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                sl0v4 = _mm256_sub_pd(sl0v4,WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);

                __m256d sr0v0 = load256d(&J_sr0b_ptr[j+0]);
                __m256d sr0v4 = load256d(&J_sr0b_ptr[j+4]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sr0b_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);

                __m256d sl0v0 = load256d(&J_sl0b_ptr[j+0]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                store256d(&J_sl0b_ptr[j+0],sl0v0);

                __m256d sr0v0 = load256d(&J_sr0b_ptr[j+0]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->NER12DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->NER12DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->NER12DupAllIndex[(k+2)*3+2];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2 + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;

        const double Kr2 = this->Kf_[i2 - Ikf[1] + Ikf[9]];

        const double dCrdC0c = C[sr0c];
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + (Kr2*dCrdC0c + Kr2*dCrdC0c);
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - 4*Kr2*dCrdC0c;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double dKrdT2 = this->dKfdT_[i2 - Ikf[1] + Ikf[9]];
        const double CF2 = C[sl0c];
        const double CR2 = C[sr0c]*C[sr0c];
        const double dqdT2 = (dKfdT2*CF2)-(dKrdT2*CR2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2 + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2 + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

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
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);

                __m256d sl0v0 = load256d(&J_sl0c_ptr[j+0]);
                __m256d sl0v4 = load256d(&J_sl0c_ptr[j+4]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                sl0v4 = _mm256_sub_pd(sl0v4,WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);

                __m256d sr0v0 = load256d(&J_sr0c_ptr[j+0]);
                __m256d sr0v4 = load256d(&J_sr0c_ptr[j+4]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sr0c_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);

                __m256d sl0v0 = load256d(&J_sl0c_ptr[j+0]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                store256d(&J_sl0c_ptr[j+0],sl0v0);

                __m256d sr0v0 = load256d(&J_sr0c_ptr[j+0]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->NER12DupAllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->NER12DupAllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->NER12DupAllIndex[(k+3)*3+2];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3 + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;

        const double Kr3 = this->Kf_[i3 - Ikf[1] + Ikf[9]];

        const double dCrdC0d = C[sr0d];
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + (Kr3*dCrdC0d + Kr3*dCrdC0d);
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - 4*Kr3*dCrdC0d;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double dKrdT3 = this->dKfdT_[i3 - Ikf[1] + Ikf[9]];
        const double CF3 = C[sl0d];
        const double CR3 = C[sr0d]*C[sr0d];
        const double dqdT3 = (dKfdT3*CF3)-(dKrdT3*CR3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3 + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3 + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;

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
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);

                __m256d sl0v0 = load256d(&J_sl0d_ptr[j+0]);
                __m256d sl0v4 = load256d(&J_sl0d_ptr[j+4]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                sl0v4 = _mm256_sub_pd(sl0v4,WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);

                __m256d sr0v0 = load256d(&J_sr0d_ptr[j+0]);
                __m256d sr0v4 = load256d(&J_sr0d_ptr[j+4]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sr0d_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);

                __m256d sl0v0 = load256d(&J_sl0d_ptr[j+0]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                store256d(&J_sl0d_ptr[j+0],sl0v0);

                __m256d sr0v0 = load256d(&J_sr0d_ptr[j+0]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
            }
        }

    }
    for(std::size_t k=endNERDup12-remainNERDup12; k<endNERDup12; k=k+1)
    {

        const unsigned int i0 = this->NER12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER12DupAllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0 + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];

        const double dCrdC0a = C[sr0a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + (Kr0*dCrdC0a + Kr0*dCrdC0a);
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - 4*Kr0*dCrdC0a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a];
        const double CR0 = C[sr0a]*C[sr0a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0 + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

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
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);
                __m256d WdMdC4 = _mm256_mul_pd(dMdC4,Wv);

                __m256d sl0v0 = load256d(&J_sl0a_ptr[j+0]);
                __m256d sl0v4 = load256d(&J_sl0a_ptr[j+4]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                sl0v4 = _mm256_sub_pd(sl0v4,WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);

                __m256d sr0v0 = load256d(&J_sr0a_ptr[j+0]);
                __m256d sr0v4 = load256d(&J_sr0a_ptr[j+4]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                sr0v4 = _mm256_add_pd(sr0v4,WdMdC4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC0 = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j+0]);
                __m256d WdMdC0 = _mm256_mul_pd(dMdC0,Wv);

                __m256d sl0v0 = load256d(&J_sl0a_ptr[j+0]);
                sl0v0 = _mm256_sub_pd(sl0v0,WdMdC0);
                store256d(&J_sl0a_ptr[j+0],sl0v0);

                __m256d sr0v0 = load256d(&J_sr0a_ptr[j+0]);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                sr0v0 = _mm256_add_pd(sr0v0,WdMdC0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
            }
        }

    }
}   
