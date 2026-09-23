/*---------------------------------------------------------------------------*\
  Description
      Update the net rates of production (dNdtByV) and the molar-concentration
      based Jacobian matrix (ddNdtByVdcTp) for one-reactant / one-product
      (1-1) reactions, e.g. A -> B and A -> C.

      The Jacobian matrix is stored row-wise with one column per species
      plus one extra column for the temperature derivative (column index
      nSpecies).

      RR : reversible reaction; the reverse rate constant is obtained from
           the equilibrium constant
      IR : irreversible reaction; the reverse rate constant is zero
      NER: non-equilibrium reaction; the reverse rate constant is computed
           with the Arrhenius form instead of the equilibrium constant

  Functions
      JF11RR : reversible 1-1 reactions
      JF11IR : irreversible 1-1 reactions
      JF11NER: non-equilibrium 1-1 reactions

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

void FastChemistry::OptReaction::JF11RR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
) const noexcept
{
    const double invKcLimiter = FastChemistry::invKcLimiter;
    std::size_t endRR11 = this->RR11Size;
    std::size_t remainRR11 = endRR11 % 4;

    for(std::size_t k=0; k<endRR11-remainRR11; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->RR11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR11AllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a]) - (dBdT[sl0a]);
        const double dKcdTByKc0 = sumVdBdT0;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CFa = C[sl0a];
        const double CRa = C[sr0a];
        const double dqdT0 = (dKfdT0*CFa) - (dKrdT0*CRa);

        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CFa) - (Kr0*CRa);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
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
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR11AllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->RR11AllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->RR11AllIndex[(k+1)*3+2];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;

        const double invKp1 = (invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double sumVdBdT1 = (dBdT[sr0b]) - (dBdT[sl0b]);
        const double dKcdTByKc1 = sumVdBdT1;
        const double dKrdT1 = (dKfdT1*invKc1 - (invKc1 < invKcLimiter ? Kr1*dKcdTByKc1 : 0));
        const double CFb = C[sl0b];
        const double CRb = C[sr0b];
        const double dqdT1 = (dKfdT1*CFb) - (dKrdT1*CRb);

        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;

        const double q1 = (Kf1*CFb) - (Kr1*CRb);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+4]),WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
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
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sr0b_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR11AllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->RR11AllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->RR11AllIndex[(k+2)*3+2];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;

        const double invKp2 = (invNegGstdByRT[sr0c])*(ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2;
        invKc2 = std::min(invKc2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double sumVdBdT2 = (dBdT[sr0c]) - (dBdT[sl0c]);
        const double dKcdTByKc2 = sumVdBdT2;
        const double dKrdT2 = (dKfdT2*invKc2 - (invKc2 < invKcLimiter ? Kr2*dKcdTByKc2 : 0));
        const double CFc = C[sl0c];
        const double CRc = C[sr0c];
        const double dqdT2 = (dKfdT2*CFc) - (dKrdT2*CRc);

        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;

        const double q2 = (Kf2*CFc) - (Kr2*CRc);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+4]),WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
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
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sr0c_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR11AllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->RR11AllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->RR11AllIndex[(k+3)*3+2];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;

        const double invKp3 = (invNegGstdByRT[sr0d])*(ExpNegGbyRT[sl0d]);
        double invKc3 = invKp3;
        invKc3 = std::min(invKc3, invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double sumVdBdT3 = (dBdT[sr0d]) - (dBdT[sl0d]);
        const double dKcdTByKc3 = sumVdBdT3;
        const double dKrdT3 = (dKfdT3*invKc3 - (invKc3 < invKcLimiter ? Kr3*dKcdTByKc3 : 0));
        const double CFd = C[sl0d];
        const double CRd = C[sr0d];
        const double dqdT3 = (dKfdT3*CFd) - (dKrdT3*CRd);

        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;

        const double q3 = (Kf3*CFd) - (Kr3*CRd);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;

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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+4]),WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
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
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sr0d_ptr[j],sr0v);
            }
        }

    }

    for(std::size_t k=endRR11-remainRR11; k<endRR11; k=k+1)
    {
        const unsigned int i0 = this->RR11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR11AllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a]) - (dBdT[sl0a]);
        const double dKcdTByKc0 = sumVdBdT0;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CFa = C[sl0a];
        const double CRa = C[sr0a];
        const double dqdT0 = (dKfdT0*CFa) - (dKrdT0*CRa);

        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CFa) - (Kr0*CRa);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CFa-invKc0*CRa);
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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }
    }
}


void FastChemistry::OptReaction::JF11IR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
) const noexcept
{
    std::size_t endIR11 = this->IR11Size;
    std::size_t remainIR11 = endIR11 % 4;

    for(std::size_t k=0; k<endIR11-remainIR11; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->IR11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR11AllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CFa = C[sl0a];
        const double dqdT0 = (dKfdT0*CFa);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CFa);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*(CFa);
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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
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
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR11AllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->IR11AllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->IR11AllIndex[(k+1)*3+2];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double CFb = C[sl0b];
        const double dqdT1 = (dKfdT1*CFb);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;

        const double q1 = (Kf1*CFb);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*(CFb);
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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+4]),WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
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
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sr0b_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR11AllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->IR11AllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->IR11AllIndex[(k+2)*3+2];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double CFc = C[sl0c];
        const double dqdT2 = (dKfdT2*CFc);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;

        const double q2 = (Kf2*CFc);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*(CFc);
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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+4]),WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
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
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sr0c_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR11AllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->IR11AllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->IR11AllIndex[(k+3)*3+2];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double CFd = C[sl0d];
        const double dqdT3 = (dKfdT3*CFd);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;

        const double q3 = (Kf3*CFd);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*(CFd);
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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+4]),WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
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
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sr0d_ptr[j],sr0v);
            }
        }

    }

    for(std::size_t k=endIR11-remainIR11; k<endIR11; k=k+1)
    {
        const unsigned int i0 = this->IR11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR11AllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CFa = C[sl0a];
        const double dqdT0 = (dKfdT0*CFa);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CFa);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*(CFa);
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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }
    }
}


void FastChemistry::OptReaction::JF11NER
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
) const noexcept
{
    std::size_t endNER11 = this->NER11Size;
    std::size_t remainNER11 = endNER11 % 4;

    for(std::size_t k=0; k<endNER11-remainNER11; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->NER11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER11AllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CFa = C[sl0a];
        const double CRa = C[sr0a];
        const double dqdT0 = (dKfdT0*CFa) - (dKrdT0*CRa);

        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CFa) - (Kr0*CRa);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[2] && i0<this->Ikf[3])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CFa - this->dKfdC_[idx0+this->Itbr[4]]*CRa;
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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
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
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->NER11AllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->NER11AllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->NER11AllIndex[(k+1)*3+2];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;

        const double Kr1 = this->Kf_[i1 - Ikf[1] + Ikf[9]];
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double dKrdT1 = this->dKfdT_[i1 - Ikf[1] + Ikf[9]];
        const double CFb = C[sl0b];
        const double CRb = C[sr0b];
        const double dqdT1 = (dKfdT1*CFb) - (dKrdT1*CRb);

        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;

        const double q1 = (Kf1*CFb) - (Kr1*CRb);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        if(i1>=this->Ikf[2] && i1<this->Ikf[3])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*CFb - this->dKfdC_[idx1+this->Itbr[4]]*CRb;
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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0b_ptr[j+4]),WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
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
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sr0b_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->NER11AllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->NER11AllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->NER11AllIndex[(k+2)*3+2];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;

        const double Kr2 = this->Kf_[i2 - Ikf[1] + Ikf[9]];
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double dKrdT2 = this->dKfdT_[i2 - Ikf[1] + Ikf[9]];
        const double CFc = C[sl0c];
        const double CRc = C[sr0c];
        const double dqdT2 = (dKfdT2*CFc) - (dKrdT2*CRc);

        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;

        const double q2 = (Kf2*CFc) - (Kr2*CRc);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        if(i2>=this->Ikf[2] && i2<this->Ikf[3])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*CFc - this->dKfdC_[idx2+this->Itbr[4]]*CRc;
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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0c_ptr[j+4]),WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
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
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sr0c_ptr[j],sr0v);
            }
        }

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->NER11AllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->NER11AllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->NER11AllIndex[(k+3)*3+2];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;

        const double Kr3 = this->Kf_[i3 - Ikf[1] + Ikf[9]];
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double dKrdT3 = this->dKfdT_[i3 - Ikf[1] + Ikf[9]];
        const double CFd = C[sl0d];
        const double CRd = C[sr0d];
        const double dqdT3 = (dKfdT3*CFd) - (dKrdT3*CRd);

        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;

        const double q3 = (Kf3*CFd) - (Kr3*CRd);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;

        if(i3>=this->Ikf[2] && i3<this->Ikf[3])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*CFd - this->dKfdC_[idx3+this->Itbr[4]]*CRd;
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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0d_ptr[j+4]),WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
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
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sr0d_ptr[j],sr0v);
            }
        }

    }

    for(std::size_t k=endNER11-remainNER11; k<endNER11; k=k+1)
    {
        const unsigned int i0 = this->NER11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER11AllIndex[(k+0)*3+2];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CFa = C[sl0a];
        const double CRa = C[sr0a];
        const double dqdT0 = (dKfdT0*CFa) - (dKrdT0*CRa);

        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;

        const double q0 = (Kf0*CFa) - (Kr0*CRa);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        if(i0>=this->Ikf[2] && i0<this->Ikf[3])
        {
            const unsigned int idx = i0 - this->Ikf[2];
            const double W = this->dKfdC_[idx]*CFa - this->dKfdC_[idx+this->Itbr[4]]*CRa;
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
                __m256d sr0v0 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+0]),WdMdC0);
                __m256d sr0v4 = _mm256_add_pd(load256d(&J_sr0a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
            }
        }
    }
}
