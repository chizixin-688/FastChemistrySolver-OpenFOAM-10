/*---------------------------------------------------------------------------*\
  Description
      Computing the molar concentration based jacobian matrix. The function 
      is used for one-three reaction, e.g. A=B+B+B. A=B+C+D
      
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

//=============================================================================//


void  FastChemistry::OptReaction::JF13RR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{
    const double PowPByRTSumVki = this->Pow_pByRT_SumVki[0];
    const double invKcLimiter = FastChemistry::invKcLimiter;
    std::size_t endRR13 = this->RR13Size;
    std::size_t remainRR13 = endRR13%4;
    for(std::size_t k=0; k<endRR13-remainRR13; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->RR13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->RR13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->RR13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->RR13AllIndex[(k+0)*5+4];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];
        double* __restrict__ J_sr2a_ptr = &ddNdtByVdcTp[sr2a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0;
        const double J_sr2a_sl0a_increment = J_sr2a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        J_sr2a_ptr[sl0a] = J_sr2a_sl0a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
            (invNegGstdByRT[sr2a]*ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double dCrdC0a = C[sr1a]*C[sr2a];
        const double dCrdC1a = C[sr0a]*C[sr2a];
        const double dCrdC2a = C[sr0a]*C[sr1a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0*dCrdC0a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr2a_sr0a_increment = J_sr2a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0*dCrdC1a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr2a_sr1a_increment = J_sr2a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sl0a_sr2a_increment = J_sl0a_ptr[sr2a] + Kr0*dCrdC2a;
        const double J_sr0a_sr2a_increment = J_sr0a_ptr[sr2a] - Kr0*dCrdC2a;
        const double J_sr1a_sr2a_increment = J_sr1a_ptr[sr2a] - Kr0*dCrdC2a;
        const double J_sr2a_sr2a_increment = J_sr2a_ptr[sr2a] - Kr0*dCrdC2a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sr2a_ptr[sr0a] = J_sr2a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;
        J_sr2a_ptr[sr1a] = J_sr2a_sr1a_increment;
        J_sl0a_ptr[sr2a] = J_sl0a_sr2a_increment;
        J_sr0a_ptr[sr2a] = J_sr0a_sr2a_increment;
        J_sr1a_ptr[sr2a] = J_sr1a_sr2a_increment;
        J_sr2a_ptr[sr2a] = J_sr2a_sr2a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a] + dBdT[sr2a]) - (dBdT[sl0a]);
        const double dKcdTByKc0 = sumVdBdT0 - invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C[sl0a];
        const double CR0 = C[sr0a]*C[sr1a]*C[sr2a];
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);

        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        const double J_sr2a_T_increment = J_sr2a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;
        J_sr2a_ptr[this->nSpecies] = J_sr2a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*(CF0-invKc0*CR0);
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2a_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sr2a_ptr[j+0],sr2v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
                store256d(&J_sr2a_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
                store256d(&J_sr2a_ptr[j],sr2v);
            }
        }

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->RR13AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR13AllIndex[(k+1)*5+1];
        const unsigned int sr0b = this->RR13AllIndex[(k+1)*5+2];
        const unsigned int sr1b = this->RR13AllIndex[(k+1)*5+3];
        const unsigned int sr2b = this->RR13AllIndex[(k+1)*5+4];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];
        double* __restrict__ J_sr2b_ptr = &ddNdtByVdcTp[sr2b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1;
        const double J_sr2b_sl0b_increment = J_sr2b_ptr[sl0b] + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;
        J_sr2b_ptr[sl0b] = J_sr2b_sl0b_increment;

        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*
            (invNegGstdByRT[sr2b]*ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double dCrdC0b = C[sr1b]*C[sr2b];
        const double dCrdC1b = C[sr0b]*C[sr2b];
        const double dCrdC2b = C[sr0b]*C[sr1b];
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1*dCrdC0b;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1*dCrdC0b;
        const double J_sr1b_sr0b_increment = J_sr1b_ptr[sr0b] - Kr1*dCrdC0b;
        const double J_sr2b_sr0b_increment = J_sr2b_ptr[sr0b] - Kr1*dCrdC0b;
        const double J_sl0b_sr1b_increment = J_sl0b_ptr[sr1b] + Kr1*dCrdC1b;
        const double J_sr0b_sr1b_increment = J_sr0b_ptr[sr1b] - Kr1*dCrdC1b;
        const double J_sr1b_sr1b_increment = J_sr1b_ptr[sr1b] - Kr1*dCrdC1b;
        const double J_sr2b_sr1b_increment = J_sr2b_ptr[sr1b] - Kr1*dCrdC1b;
        const double J_sl0b_sr2b_increment = J_sl0b_ptr[sr2b] + Kr1*dCrdC2b;
        const double J_sr0b_sr2b_increment = J_sr0b_ptr[sr2b] - Kr1*dCrdC2b;
        const double J_sr1b_sr2b_increment = J_sr1b_ptr[sr2b] - Kr1*dCrdC2b;
        const double J_sr2b_sr2b_increment = J_sr2b_ptr[sr2b] - Kr1*dCrdC2b;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;
        J_sr1b_ptr[sr0b] = J_sr1b_sr0b_increment;
        J_sr2b_ptr[sr0b] = J_sr2b_sr0b_increment;
        J_sl0b_ptr[sr1b] = J_sl0b_sr1b_increment;
        J_sr0b_ptr[sr1b] = J_sr0b_sr1b_increment;
        J_sr1b_ptr[sr1b] = J_sr1b_sr1b_increment;
        J_sr2b_ptr[sr1b] = J_sr2b_sr1b_increment;
        J_sl0b_ptr[sr2b] = J_sl0b_sr2b_increment;
        J_sr0b_ptr[sr2b] = J_sr0b_sr2b_increment;
        J_sr1b_ptr[sr2b] = J_sr1b_sr2b_increment;
        J_sr2b_ptr[sr2b] = J_sr2b_sr2b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double sumVdBdT1 = (dBdT[sr0b] + dBdT[sr1b] + dBdT[sr2b]) - (dBdT[sl0b]);
        const double dKcdTByKc1 = sumVdBdT1 - invT;
        const double dKrdT1 = (dKfdT1*invKc1 - (invKc1 < invKcLimiter ? Kr1*dKcdTByKc1 : 0));
        const double CF1 = C[sl0b];
        const double CR1 = C[sr0b]*C[sr1b]*C[sr2b];
        const double dqdT1 = (dKfdT1*CF1) - (dKrdT1*CR1);

        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1;
        const double J_sr2b_T_increment = J_sr2b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;
        J_sr2b_ptr[this->nSpecies] = J_sr2b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*(CF1-invKc1*CR1);
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2b_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2b_ptr[j+4]),WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sr1b_ptr[j+0],sr1v0);
                store256d(&J_sr2b_ptr[j+0],sr2v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
                store256d(&J_sr0b_ptr[j+4],sr0v4);
                store256d(&J_sr1b_ptr[j+4],sr1v4);
                store256d(&J_sr2b_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv1);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0b_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0b_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1b_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2b_ptr[j]),WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sr0b_ptr[j],sr0v);
                store256d(&J_sr1b_ptr[j],sr1v);
                store256d(&J_sr2b_ptr[j],sr2v);
            }
        }

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->RR13AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->RR13AllIndex[(k+2)*5+1];
        const unsigned int sr0c = this->RR13AllIndex[(k+2)*5+2];
        const unsigned int sr1c = this->RR13AllIndex[(k+2)*5+3];
        const unsigned int sr2c = this->RR13AllIndex[(k+2)*5+4];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];
        double* __restrict__ J_sr2c_ptr = &ddNdtByVdcTp[sr2c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2;
        const double J_sr2c_sl0c_increment = J_sr2c_ptr[sl0c] + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;
        J_sr2c_ptr[sl0c] = J_sr2c_sl0c_increment;

        const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*
            (invNegGstdByRT[sr2c]*ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double dCrdC0c = C[sr1c]*C[sr2c];
        const double dCrdC1c = C[sr0c]*C[sr2c];
        const double dCrdC2c = C[sr0c]*C[sr1c];
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2*dCrdC0c;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2*dCrdC0c;
        const double J_sr1c_sr0c_increment = J_sr1c_ptr[sr0c] - Kr2*dCrdC0c;
        const double J_sr2c_sr0c_increment = J_sr2c_ptr[sr0c] - Kr2*dCrdC0c;
        const double J_sl0c_sr1c_increment = J_sl0c_ptr[sr1c] + Kr2*dCrdC1c;
        const double J_sr0c_sr1c_increment = J_sr0c_ptr[sr1c] - Kr2*dCrdC1c;
        const double J_sr1c_sr1c_increment = J_sr1c_ptr[sr1c] - Kr2*dCrdC1c;
        const double J_sr2c_sr1c_increment = J_sr2c_ptr[sr1c] - Kr2*dCrdC1c;
        const double J_sl0c_sr2c_increment = J_sl0c_ptr[sr2c] + Kr2*dCrdC2c;
        const double J_sr0c_sr2c_increment = J_sr0c_ptr[sr2c] - Kr2*dCrdC2c;
        const double J_sr1c_sr2c_increment = J_sr1c_ptr[sr2c] - Kr2*dCrdC2c;
        const double J_sr2c_sr2c_increment = J_sr2c_ptr[sr2c] - Kr2*dCrdC2c;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;
        J_sr1c_ptr[sr0c] = J_sr1c_sr0c_increment;
        J_sr2c_ptr[sr0c] = J_sr2c_sr0c_increment;
        J_sl0c_ptr[sr1c] = J_sl0c_sr1c_increment;
        J_sr0c_ptr[sr1c] = J_sr0c_sr1c_increment;
        J_sr1c_ptr[sr1c] = J_sr1c_sr1c_increment;
        J_sr2c_ptr[sr1c] = J_sr2c_sr1c_increment;
        J_sl0c_ptr[sr2c] = J_sl0c_sr2c_increment;
        J_sr0c_ptr[sr2c] = J_sr0c_sr2c_increment;
        J_sr1c_ptr[sr2c] = J_sr1c_sr2c_increment;
        J_sr2c_ptr[sr2c] = J_sr2c_sr2c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double sumVdBdT2 = (dBdT[sr0c] + dBdT[sr1c] + dBdT[sr2c]) - (dBdT[sl0c]);
        const double dKcdTByKc2 = sumVdBdT2 - invT;
        const double dKrdT2 = (dKfdT2*invKc2 - (invKc2 < invKcLimiter ? Kr2*dKcdTByKc2 : 0));
        const double CF2 = C[sl0c];
        const double CR2 = C[sr0c]*C[sr1c]*C[sr2c];
        const double dqdT2 = (dKfdT2*CF2) - (dKrdT2*CR2);

        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2;
        const double J_sr2c_T_increment = J_sr2c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;
        J_sr2c_ptr[this->nSpecies] = J_sr2c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*(CF2-invKc2*CR2);
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2c_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2c_ptr[j+4]),WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sr1c_ptr[j+0],sr1v0);
                store256d(&J_sr2c_ptr[j+0],sr2v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
                store256d(&J_sr0c_ptr[j+4],sr0v4);
                store256d(&J_sr1c_ptr[j+4],sr1v4);
                store256d(&J_sr2c_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv2);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0c_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0c_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1c_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2c_ptr[j]),WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sr0c_ptr[j],sr0v);
                store256d(&J_sr1c_ptr[j],sr1v);
                store256d(&J_sr2c_ptr[j],sr2v);
            }
        }

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->RR13AllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->RR13AllIndex[(k+3)*5+1];
        const unsigned int sr0d = this->RR13AllIndex[(k+3)*5+2];
        const unsigned int sr1d = this->RR13AllIndex[(k+3)*5+3];
        const unsigned int sr2d = this->RR13AllIndex[(k+3)*5+4];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];
        double* __restrict__ J_sr2d_ptr = &ddNdtByVdcTp[sr2d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3;
        const double J_sr2d_sl0d_increment = J_sr2d_ptr[sl0d] + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;
        J_sr2d_ptr[sl0d] = J_sr2d_sl0d_increment;

        const double invKp3 = (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d])*
            (invNegGstdByRT[sr2d]*ExpNegGbyRT[sl0d]);
        double invKc3 = invKp3*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double dCrdC0d = C[sr1d]*C[sr2d];
        const double dCrdC1d = C[sr0d]*C[sr2d];
        const double dCrdC2d = C[sr0d]*C[sr1d];
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3*dCrdC0d;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3*dCrdC0d;
        const double J_sr1d_sr0d_increment = J_sr1d_ptr[sr0d] - Kr3*dCrdC0d;
        const double J_sr2d_sr0d_increment = J_sr2d_ptr[sr0d] - Kr3*dCrdC0d;
        const double J_sl0d_sr1d_increment = J_sl0d_ptr[sr1d] + Kr3*dCrdC1d;
        const double J_sr0d_sr1d_increment = J_sr0d_ptr[sr1d] - Kr3*dCrdC1d;
        const double J_sr1d_sr1d_increment = J_sr1d_ptr[sr1d] - Kr3*dCrdC1d;
        const double J_sr2d_sr1d_increment = J_sr2d_ptr[sr1d] - Kr3*dCrdC1d;
        const double J_sl0d_sr2d_increment = J_sl0d_ptr[sr2d] + Kr3*dCrdC2d;
        const double J_sr0d_sr2d_increment = J_sr0d_ptr[sr2d] - Kr3*dCrdC2d;
        const double J_sr1d_sr2d_increment = J_sr1d_ptr[sr2d] - Kr3*dCrdC2d;
        const double J_sr2d_sr2d_increment = J_sr2d_ptr[sr2d] - Kr3*dCrdC2d;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;
        J_sr1d_ptr[sr0d] = J_sr1d_sr0d_increment;
        J_sr2d_ptr[sr0d] = J_sr2d_sr0d_increment;
        J_sl0d_ptr[sr1d] = J_sl0d_sr1d_increment;
        J_sr0d_ptr[sr1d] = J_sr0d_sr1d_increment;
        J_sr1d_ptr[sr1d] = J_sr1d_sr1d_increment;
        J_sr2d_ptr[sr1d] = J_sr2d_sr1d_increment;
        J_sl0d_ptr[sr2d] = J_sl0d_sr2d_increment;
        J_sr0d_ptr[sr2d] = J_sr0d_sr2d_increment;
        J_sr1d_ptr[sr2d] = J_sr1d_sr2d_increment;
        J_sr2d_ptr[sr2d] = J_sr2d_sr2d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double sumVdBdT3 = (dBdT[sr0d] + dBdT[sr1d] + dBdT[sr2d]) - (dBdT[sl0d]);
        const double dKcdTByKc3 = sumVdBdT3 - invT;
        const double dKrdT3 = (dKfdT3*invKc3 - (invKc3 < invKcLimiter ? Kr3*dKcdTByKc3 : 0));
        const double CF3 = C[sl0d];
        const double CR3 = C[sr0d]*C[sr1d]*C[sr2d];
        const double dqdT3 = (dKfdT3*CF3) - (dKrdT3*CR3);

        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3;
        const double J_sr2d_T_increment = J_sr2d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;
        J_sr2d_ptr[this->nSpecies] = J_sr2d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        const double omega_sr2d_increment = dNdtByV[sr2d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
        dNdtByV[sr2d] = omega_sr2d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*(CF3-invKc3*CR3);
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2d_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2d_ptr[j+4]),WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sr1d_ptr[j+0],sr1v0);
                store256d(&J_sr2d_ptr[j+0],sr2v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
                store256d(&J_sr0d_ptr[j+4],sr0v4);
                store256d(&J_sr1d_ptr[j+4],sr1v4);
                store256d(&J_sr2d_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv3);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0d_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0d_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1d_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2d_ptr[j]),WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sr0d_ptr[j],sr0v);
                store256d(&J_sr1d_ptr[j],sr1v);
                store256d(&J_sr2d_ptr[j],sr2v);
            }
        }
    }
    for(std::size_t k=endRR13-remainRR13; k<endRR13; k=k+1)
    {
        const unsigned int i0 = this->RR13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->RR13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->RR13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->RR13AllIndex[(k+0)*5+4];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];
        double* __restrict__ J_sr2a_ptr = &ddNdtByVdcTp[sr2a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0;
        const double J_sr2a_sl0a_increment = J_sr2a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        J_sr2a_ptr[sl0a] = J_sr2a_sl0a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
            (invNegGstdByRT[sr2a]*ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double dCrdC0a = C[sr1a]*C[sr2a];
        const double dCrdC1a = C[sr0a]*C[sr2a];
        const double dCrdC2a = C[sr0a]*C[sr1a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0*dCrdC0a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr2a_sr0a_increment = J_sr2a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0*dCrdC1a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr2a_sr1a_increment = J_sr2a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sl0a_sr2a_increment = J_sl0a_ptr[sr2a] + Kr0*dCrdC2a;
        const double J_sr0a_sr2a_increment = J_sr0a_ptr[sr2a] - Kr0*dCrdC2a;
        const double J_sr1a_sr2a_increment = J_sr1a_ptr[sr2a] - Kr0*dCrdC2a;
        const double J_sr2a_sr2a_increment = J_sr2a_ptr[sr2a] - Kr0*dCrdC2a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sr2a_ptr[sr0a] = J_sr2a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;
        J_sr2a_ptr[sr1a] = J_sr2a_sr1a_increment;
        J_sl0a_ptr[sr2a] = J_sl0a_sr2a_increment;
        J_sr0a_ptr[sr2a] = J_sr0a_sr2a_increment;
        J_sr1a_ptr[sr2a] = J_sr1a_sr2a_increment;
        J_sr2a_ptr[sr2a] = J_sr2a_sr2a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a] + dBdT[sr2a]) - (dBdT[sl0a]);
        const double dKcdTByKc0 = sumVdBdT0 - invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C[sl0a];
        const double CR0 = C[sr0a]*C[sr1a]*C[sr2a];
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);

        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        const double J_sr2a_T_increment = J_sr2a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;
        J_sr2a_ptr[this->nSpecies] = J_sr2a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*(CF0-invKc0*CR0);
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2a_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sr2a_ptr[j+0],sr2v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
                store256d(&J_sr2a_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
                store256d(&J_sr2a_ptr[j],sr2v);
            }
        }
    }

    std::size_t endRRDup13 = this->RR13DupSize;
    std::size_t remainRRDup13 = endRRDup13%4;
    for(std::size_t k=0; k<endRRDup13-remainRRDup13; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->RR13DupAllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->RR13DupAllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->RR13DupAllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->RR13DupAllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0 + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
            (invNegGstdByRT[sr1a]*ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double dCrdC0a = C[sr1a]*C[sr1a];
        const double dCrdC1a = 2*C[sr0a]*C[sr1a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0*dCrdC0a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0*dCrdC0a - Kr0*dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0*dCrdC1a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0*dCrdC1a - Kr0*dCrdC1a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a] + dBdT[sr1a]) - (dBdT[sl0a]);
        const double dKcdTByKc0 = sumVdBdT0 - invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C[sl0a];
        const double CR0 = C[sr0a]*C[sr1a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);

        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0 + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*(CF0-invKc0*CR0);
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->RR13DupAllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->RR13DupAllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->RR13DupAllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->RR13DupAllIndex[(k+1)*4+3];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1 + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;

        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*
            (invNegGstdByRT[sr1b]*ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double dCrdC0b = C[sr1b]*C[sr1b];
        const double dCrdC1b = 2*C[sr0b]*C[sr1b];
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1*dCrdC0b;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1*dCrdC0b;
        const double J_sr1b_sr0b_increment = J_sr1b_ptr[sr0b] - Kr1*dCrdC0b - Kr1*dCrdC0b;
        const double J_sl0b_sr1b_increment = J_sl0b_ptr[sr1b] + Kr1*dCrdC1b;
        const double J_sr0b_sr1b_increment = J_sr0b_ptr[sr1b] - Kr1*dCrdC1b;
        const double J_sr1b_sr1b_increment = J_sr1b_ptr[sr1b] - Kr1*dCrdC1b - Kr1*dCrdC1b;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;
        J_sr1b_ptr[sr0b] = J_sr1b_sr0b_increment;
        J_sl0b_ptr[sr1b] = J_sl0b_sr1b_increment;
        J_sr0b_ptr[sr1b] = J_sr0b_sr1b_increment;
        J_sr1b_ptr[sr1b] = J_sr1b_sr1b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double sumVdBdT1 = (dBdT[sr0b] + dBdT[sr1b] + dBdT[sr1b]) - (dBdT[sl0b]);
        const double dKcdTByKc1 = sumVdBdT1 - invT;
        const double dKrdT1 = (dKfdT1*invKc1 - (invKc1 < invKcLimiter ? Kr1*dKcdTByKc1 : 0));
        const double CF1 = C[sl0b];
        const double CR1 = C[sr0b]*C[sr1b]*C[sr1b];
        const double dqdT1 = (dKfdT1*CF1) - (dKrdT1*CR1);

        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1 + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*(CF1-invKc1*CR1);
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv1);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0b_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0b_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1b_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sr0b_ptr[j],sr0v);
                store256d(&J_sr1b_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->RR13DupAllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->RR13DupAllIndex[(k+2)*4+1];
        const unsigned int sr0c = this->RR13DupAllIndex[(k+2)*4+2];
        const unsigned int sr1c = this->RR13DupAllIndex[(k+2)*4+3];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2 + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;

        const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*
            (invNegGstdByRT[sr1c]*ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double dCrdC0c = C[sr1c]*C[sr1c];
        const double dCrdC1c = 2*C[sr0c]*C[sr1c];
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2*dCrdC0c;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2*dCrdC0c;
        const double J_sr1c_sr0c_increment = J_sr1c_ptr[sr0c] - Kr2*dCrdC0c - Kr2*dCrdC0c;
        const double J_sl0c_sr1c_increment = J_sl0c_ptr[sr1c] + Kr2*dCrdC1c;
        const double J_sr0c_sr1c_increment = J_sr0c_ptr[sr1c] - Kr2*dCrdC1c;
        const double J_sr1c_sr1c_increment = J_sr1c_ptr[sr1c] - Kr2*dCrdC1c - Kr2*dCrdC1c;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;
        J_sr1c_ptr[sr0c] = J_sr1c_sr0c_increment;
        J_sl0c_ptr[sr1c] = J_sl0c_sr1c_increment;
        J_sr0c_ptr[sr1c] = J_sr0c_sr1c_increment;
        J_sr1c_ptr[sr1c] = J_sr1c_sr1c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double sumVdBdT2 = (dBdT[sr0c] + dBdT[sr1c] + dBdT[sr1c]) - (dBdT[sl0c]);
        const double dKcdTByKc2 = sumVdBdT2 - invT;
        const double dKrdT2 = (dKfdT2*invKc2 - (invKc2 < invKcLimiter ? Kr2*dKcdTByKc2 : 0));
        const double CF2 = C[sl0c];
        const double CR2 = C[sr0c]*C[sr1c]*C[sr1c];
        const double dqdT2 = (dKfdT2*CF2) - (dKrdT2*CR2);

        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2 + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2 + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*(CF2-invKc2*CR2);
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv2);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0c_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0c_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1c_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sr0c_ptr[j],sr0v);
                store256d(&J_sr1c_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->RR13DupAllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->RR13DupAllIndex[(k+3)*4+1];
        const unsigned int sr0d = this->RR13DupAllIndex[(k+3)*4+2];
        const unsigned int sr1d = this->RR13DupAllIndex[(k+3)*4+3];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3 + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;

        const double invKp3 = (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d])*
            (invNegGstdByRT[sr1d]*ExpNegGbyRT[sl0d]);
        double invKc3 = invKp3*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double dCrdC0d = C[sr1d]*C[sr1d];
        const double dCrdC1d = 2*C[sr0d]*C[sr1d];
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3*dCrdC0d;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3*dCrdC0d;
        const double J_sr1d_sr0d_increment = J_sr1d_ptr[sr0d] - Kr3*dCrdC0d - Kr3*dCrdC0d;
        const double J_sl0d_sr1d_increment = J_sl0d_ptr[sr1d] + Kr3*dCrdC1d;
        const double J_sr0d_sr1d_increment = J_sr0d_ptr[sr1d] - Kr3*dCrdC1d;
        const double J_sr1d_sr1d_increment = J_sr1d_ptr[sr1d] - Kr3*dCrdC1d - Kr3*dCrdC1d;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;
        J_sr1d_ptr[sr0d] = J_sr1d_sr0d_increment;
        J_sl0d_ptr[sr1d] = J_sl0d_sr1d_increment;
        J_sr0d_ptr[sr1d] = J_sr0d_sr1d_increment;
        J_sr1d_ptr[sr1d] = J_sr1d_sr1d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double sumVdBdT3 = (dBdT[sr0d] + dBdT[sr1d] + dBdT[sr1d]) - (dBdT[sl0d]);
        const double dKcdTByKc3 = sumVdBdT3 - invT;
        const double dKrdT3 = (dKfdT3*invKc3 - (invKc3 < invKcLimiter ? Kr3*dKcdTByKc3 : 0));
        const double CF3 = C[sl0d];
        const double CR3 = C[sr0d]*C[sr1d]*C[sr1d];
        const double dqdT3 = (dKfdT3*CF3) - (dKrdT3*CR3);

        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3 + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3 + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*(CF3-invKc3*CR3);
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv3);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0d_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0d_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1d_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sr0d_ptr[j],sr0v);
                store256d(&J_sr1d_ptr[j],sr1v);
            }
        }
    }
    for(std::size_t k=endRRDup13-remainRRDup13; k<endRRDup13; k=k+1)
    {
        const unsigned int i0 = this->RR13DupAllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->RR13DupAllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->RR13DupAllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->RR13DupAllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0 + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;

        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
            (invNegGstdByRT[sr1a]*ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double dCrdC0a = C[sr1a]*C[sr1a];
        const double dCrdC1a = 2*C[sr0a]*C[sr1a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0*dCrdC0a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0*dCrdC0a - Kr0*dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0*dCrdC1a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0*dCrdC1a - Kr0*dCrdC1a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a] + dBdT[sr1a]) - (dBdT[sl0a]);
        const double dKcdTByKc0 = sumVdBdT0 - invT;
        const double dKrdT0 = (dKfdT0*invKc0 - (invKc0 < invKcLimiter ? Kr0*dKcdTByKc0 : 0));
        const double CF0 = C[sl0a];
        const double CR0 = C[sr0a]*C[sr1a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0) - (dKrdT0*CR0);

        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0 + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*(CF0-invKc0*CR0);
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }
    }
}


void  FastChemistry::OptReaction::JF13IR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{
    std::size_t endIR13 = this->IR13Size;
    std::size_t remainIR13 = endIR13%4;
    for(std::size_t k=0; k<endIR13-remainIR13; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->IR13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->IR13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->IR13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->IR13AllIndex[(k+0)*5+4];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];
        double* __restrict__ J_sr2a_ptr = &ddNdtByVdcTp[sr2a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0;
        const double J_sr2a_sl0a_increment = J_sr2a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        J_sr2a_ptr[sl0a] = J_sr2a_sl0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C[sl0a];
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        const double J_sr2a_T_increment = J_sr2a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;
        J_sr2a_ptr[this->nSpecies] = J_sr2a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2a_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sr2a_ptr[j+0],sr2v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
                store256d(&J_sr2a_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
                store256d(&J_sr2a_ptr[j],sr2v);
            }
        }

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->IR13AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR13AllIndex[(k+1)*5+1];
        const unsigned int sr0b = this->IR13AllIndex[(k+1)*5+2];
        const unsigned int sr1b = this->IR13AllIndex[(k+1)*5+3];
        const unsigned int sr2b = this->IR13AllIndex[(k+1)*5+4];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];
        double* __restrict__ J_sr2b_ptr = &ddNdtByVdcTp[sr2b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1;
        const double J_sr2b_sl0b_increment = J_sr2b_ptr[sl0b] + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;
        J_sr2b_ptr[sl0b] = J_sr2b_sl0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double CF1 = C[sl0b];
        const double dqdT1 = (dKfdT1*CF1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1;
        const double J_sr2b_T_increment = J_sr2b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;
        J_sr2b_ptr[this->nSpecies] = J_sr2b_T_increment;

        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*CF1;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2b_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2b_ptr[j+4]),WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sr1b_ptr[j+0],sr1v0);
                store256d(&J_sr2b_ptr[j+0],sr2v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
                store256d(&J_sr0b_ptr[j+4],sr0v4);
                store256d(&J_sr1b_ptr[j+4],sr1v4);
                store256d(&J_sr2b_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv1);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0b_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0b_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1b_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2b_ptr[j]),WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sr0b_ptr[j],sr0v);
                store256d(&J_sr1b_ptr[j],sr1v);
                store256d(&J_sr2b_ptr[j],sr2v);
            }
        }

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->IR13AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->IR13AllIndex[(k+2)*5+1];
        const unsigned int sr0c = this->IR13AllIndex[(k+2)*5+2];
        const unsigned int sr1c = this->IR13AllIndex[(k+2)*5+3];
        const unsigned int sr2c = this->IR13AllIndex[(k+2)*5+4];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];
        double* __restrict__ J_sr2c_ptr = &ddNdtByVdcTp[sr2c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2;
        const double J_sr2c_sl0c_increment = J_sr2c_ptr[sl0c] + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;
        J_sr2c_ptr[sl0c] = J_sr2c_sl0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double CF2 = C[sl0c];
        const double dqdT2 = (dKfdT2*CF2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2;
        const double J_sr2c_T_increment = J_sr2c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;
        J_sr2c_ptr[this->nSpecies] = J_sr2c_T_increment;

        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*CF2;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2c_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2c_ptr[j+4]),WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sr1c_ptr[j+0],sr1v0);
                store256d(&J_sr2c_ptr[j+0],sr2v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
                store256d(&J_sr0c_ptr[j+4],sr0v4);
                store256d(&J_sr1c_ptr[j+4],sr1v4);
                store256d(&J_sr2c_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv2);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0c_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0c_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1c_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2c_ptr[j]),WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sr0c_ptr[j],sr0v);
                store256d(&J_sr1c_ptr[j],sr1v);
                store256d(&J_sr2c_ptr[j],sr2v);
            }
        }

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->IR13AllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->IR13AllIndex[(k+3)*5+1];
        const unsigned int sr0d = this->IR13AllIndex[(k+3)*5+2];
        const unsigned int sr1d = this->IR13AllIndex[(k+3)*5+3];
        const unsigned int sr2d = this->IR13AllIndex[(k+3)*5+4];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];
        double* __restrict__ J_sr2d_ptr = &ddNdtByVdcTp[sr2d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3;
        const double J_sr2d_sl0d_increment = J_sr2d_ptr[sl0d] + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;
        J_sr2d_ptr[sl0d] = J_sr2d_sl0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double CF3 = C[sl0d];
        const double dqdT3 = (dKfdT3*CF3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3;
        const double J_sr2d_T_increment = J_sr2d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;
        J_sr2d_ptr[this->nSpecies] = J_sr2d_T_increment;

        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        const double omega_sr2d_increment = dNdtByV[sr2d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
        dNdtByV[sr2d] = omega_sr2d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*CF3;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2d_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2d_ptr[j+4]),WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sr1d_ptr[j+0],sr1v0);
                store256d(&J_sr2d_ptr[j+0],sr2v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
                store256d(&J_sr0d_ptr[j+4],sr0v4);
                store256d(&J_sr1d_ptr[j+4],sr1v4);
                store256d(&J_sr2d_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv3);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0d_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0d_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1d_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2d_ptr[j]),WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sr0d_ptr[j],sr0v);
                store256d(&J_sr1d_ptr[j],sr1v);
                store256d(&J_sr2d_ptr[j],sr2v);
            }
        }
    }
    for(std::size_t k=endIR13-remainIR13; k<endIR13; k=k+1)
    {
        const unsigned int i0 = this->IR13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->IR13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->IR13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->IR13AllIndex[(k+0)*5+4];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];
        double* __restrict__ J_sr2a_ptr = &ddNdtByVdcTp[sr2a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0;
        const double J_sr2a_sl0a_increment = J_sr2a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        J_sr2a_ptr[sl0a] = J_sr2a_sl0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C[sl0a];
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        const double J_sr2a_T_increment = J_sr2a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;
        J_sr2a_ptr[this->nSpecies] = J_sr2a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2a_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sr2a_ptr[j+0],sr2v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
                store256d(&J_sr2a_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
                store256d(&J_sr2a_ptr[j],sr2v);
            }
        }
    }

    std::size_t endIRDup13 = this->IR13DupSize;
    std::size_t remainIRDup13 = endIRDup13%4;
    for(std::size_t k=0; k<endIRDup13-remainIRDup13; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->IR13DupAllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->IR13DupAllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->IR13DupAllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->IR13DupAllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0 + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C[sl0a];
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0 + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->IR13DupAllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->IR13DupAllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->IR13DupAllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->IR13DupAllIndex[(k+1)*4+3];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1 + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double CF1 = C[sl0b];
        const double dqdT1 = (dKfdT1*CF1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1 + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;

        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        if(i1>=this->Ikf[3] && i1<this->Ikf[6])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*CF1;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv1);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0b_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0b_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1b_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sr0b_ptr[j],sr0v);
                store256d(&J_sr1b_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->IR13DupAllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->IR13DupAllIndex[(k+2)*4+1];
        const unsigned int sr0c = this->IR13DupAllIndex[(k+2)*4+2];
        const unsigned int sr1c = this->IR13DupAllIndex[(k+2)*4+3];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2 + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double CF2 = C[sl0c];
        const double dqdT2 = (dKfdT2*CF2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2 + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;

        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2 + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        if(i2>=this->Ikf[3] && i2<this->Ikf[6])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*CF2;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv2);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0c_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0c_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1c_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sr0c_ptr[j],sr0v);
                store256d(&J_sr1c_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->IR13DupAllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->IR13DupAllIndex[(k+3)*4+1];
        const unsigned int sr0d = this->IR13DupAllIndex[(k+3)*4+2];
        const unsigned int sr1d = this->IR13DupAllIndex[(k+3)*4+3];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3 + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double CF3 = C[sl0d];
        const double dqdT3 = (dKfdT3*CF3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3 + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;

        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3 + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;

        if(i3>=this->Ikf[3] && i3<this->Ikf[6])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*CF3;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv3);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0d_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0d_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1d_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sr0d_ptr[j],sr0v);
                store256d(&J_sr1d_ptr[j],sr1v);
            }
        }
    }
    for(std::size_t k=endIRDup13-remainIRDup13; k<endIRDup13; k=k+1)
    {
        const unsigned int i0 = this->IR13DupAllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->IR13DupAllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->IR13DupAllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->IR13DupAllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0 + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double CF0 = C[sl0a];
        const double dqdT0 = (dKfdT0*CF0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0 + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[3] && i0<this->Ikf[6])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }
    }
}


void  FastChemistry::OptReaction::JF13NER
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{
    std::size_t endNER13 = this->NER13Size;
    std::size_t remainNER13 = endNER13%4;
    for(std::size_t k=0; k<endNER13-remainNER13; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->NER13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->NER13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->NER13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->NER13AllIndex[(k+0)*5+4];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];
        double* __restrict__ J_sr2a_ptr = &ddNdtByVdcTp[sr2a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0;
        const double J_sr2a_sl0a_increment = J_sr2a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        J_sr2a_ptr[sl0a] = J_sr2a_sl0a_increment;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double dCrdC0a = C[sr1a]*C[sr2a];
        const double dCrdC1a = C[sr0a]*C[sr2a];
        const double dCrdC2a = C[sr0a]*C[sr1a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0*dCrdC0a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr2a_sr0a_increment = J_sr2a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0*dCrdC1a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr2a_sr1a_increment = J_sr2a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sl0a_sr2a_increment = J_sl0a_ptr[sr2a] + Kr0*dCrdC2a;
        const double J_sr0a_sr2a_increment = J_sr0a_ptr[sr2a] - Kr0*dCrdC2a;
        const double J_sr1a_sr2a_increment = J_sr1a_ptr[sr2a] - Kr0*dCrdC2a;
        const double J_sr2a_sr2a_increment = J_sr2a_ptr[sr2a] - Kr0*dCrdC2a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sr2a_ptr[sr0a] = J_sr2a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;
        J_sr2a_ptr[sr1a] = J_sr2a_sr1a_increment;
        J_sl0a_ptr[sr2a] = J_sl0a_sr2a_increment;
        J_sr0a_ptr[sr2a] = J_sr0a_sr2a_increment;
        J_sr1a_ptr[sr2a] = J_sr1a_sr2a_increment;
        J_sr2a_ptr[sr2a] = J_sr2a_sr2a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a];
        const double CR0 = C[sr0a]*C[sr1a]*C[sr2a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        const double J_sr2a_T_increment = J_sr2a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;
        J_sr2a_ptr[this->nSpecies] = J_sr2a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        if(i0>=this->Ikf[2] && i0<this->Ikf[3])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0 - this->dKfdC_[idx0+this->Itbr[4]]*CR0;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2a_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sr2a_ptr[j+0],sr2v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
                store256d(&J_sr2a_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
                store256d(&J_sr2a_ptr[j],sr2v);
            }
        }

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->NER13AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER13AllIndex[(k+1)*5+1];
        const unsigned int sr0b = this->NER13AllIndex[(k+1)*5+2];
        const unsigned int sr1b = this->NER13AllIndex[(k+1)*5+3];
        const unsigned int sr2b = this->NER13AllIndex[(k+1)*5+4];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];
        double* __restrict__ J_sr2b_ptr = &ddNdtByVdcTp[sr2b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1;
        const double J_sr2b_sl0b_increment = J_sr2b_ptr[sl0b] + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;
        J_sr2b_ptr[sl0b] = J_sr2b_sl0b_increment;

        const double Kr1 = this->Kf_[i1 - Ikf[1] + Ikf[9]];
        const double dCrdC0b = C[sr1b]*C[sr2b];
        const double dCrdC1b = C[sr0b]*C[sr2b];
        const double dCrdC2b = C[sr0b]*C[sr1b];
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1*dCrdC0b;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1*dCrdC0b;
        const double J_sr1b_sr0b_increment = J_sr1b_ptr[sr0b] - Kr1*dCrdC0b;
        const double J_sr2b_sr0b_increment = J_sr2b_ptr[sr0b] - Kr1*dCrdC0b;
        const double J_sl0b_sr1b_increment = J_sl0b_ptr[sr1b] + Kr1*dCrdC1b;
        const double J_sr0b_sr1b_increment = J_sr0b_ptr[sr1b] - Kr1*dCrdC1b;
        const double J_sr1b_sr1b_increment = J_sr1b_ptr[sr1b] - Kr1*dCrdC1b;
        const double J_sr2b_sr1b_increment = J_sr2b_ptr[sr1b] - Kr1*dCrdC1b;
        const double J_sl0b_sr2b_increment = J_sl0b_ptr[sr2b] + Kr1*dCrdC2b;
        const double J_sr0b_sr2b_increment = J_sr0b_ptr[sr2b] - Kr1*dCrdC2b;
        const double J_sr1b_sr2b_increment = J_sr1b_ptr[sr2b] - Kr1*dCrdC2b;
        const double J_sr2b_sr2b_increment = J_sr2b_ptr[sr2b] - Kr1*dCrdC2b;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;
        J_sr1b_ptr[sr0b] = J_sr1b_sr0b_increment;
        J_sr2b_ptr[sr0b] = J_sr2b_sr0b_increment;
        J_sl0b_ptr[sr1b] = J_sl0b_sr1b_increment;
        J_sr0b_ptr[sr1b] = J_sr0b_sr1b_increment;
        J_sr1b_ptr[sr1b] = J_sr1b_sr1b_increment;
        J_sr2b_ptr[sr1b] = J_sr2b_sr1b_increment;
        J_sl0b_ptr[sr2b] = J_sl0b_sr2b_increment;
        J_sr0b_ptr[sr2b] = J_sr0b_sr2b_increment;
        J_sr1b_ptr[sr2b] = J_sr1b_sr2b_increment;
        J_sr2b_ptr[sr2b] = J_sr2b_sr2b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double dKrdT1 = this->dKfdT_[i1 - Ikf[1] + Ikf[9]];
        const double CF1 = C[sl0b];
        const double CR1 = C[sr0b]*C[sr1b]*C[sr2b];
        const double dqdT1 = (dKfdT1*CF1)-(dKrdT1*CR1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1;
        const double J_sr2b_T_increment = J_sr2b_ptr[this->nSpecies] + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;
        J_sr2b_ptr[this->nSpecies] = J_sr2b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        if(i1>=this->Ikf[2] && i1<this->Ikf[3])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*CF1 - this->dKfdC_[idx1+this->Itbr[4]]*CR1;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2b_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2b_ptr[j+4]),WdMdC4);
                store256d(&J_sl0b_ptr[j+0],sl0v0);
                store256d(&J_sr0b_ptr[j+0],sr0v0);
                store256d(&J_sr1b_ptr[j+0],sr1v0);
                store256d(&J_sr2b_ptr[j+0],sr2v0);
                store256d(&J_sl0b_ptr[j+4],sl0v4);
                store256d(&J_sr0b_ptr[j+4],sr0v4);
                store256d(&J_sr1b_ptr[j+4],sr1v4);
                store256d(&J_sr2b_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv1);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0b_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0b_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1b_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2b_ptr[j]),WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sr0b_ptr[j],sr0v);
                store256d(&J_sr1b_ptr[j],sr1v);
                store256d(&J_sr2b_ptr[j],sr2v);
            }
        }

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->NER13AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER13AllIndex[(k+2)*5+1];
        const unsigned int sr0c = this->NER13AllIndex[(k+2)*5+2];
        const unsigned int sr1c = this->NER13AllIndex[(k+2)*5+3];
        const unsigned int sr2c = this->NER13AllIndex[(k+2)*5+4];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];
        double* __restrict__ J_sr2c_ptr = &ddNdtByVdcTp[sr2c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2;
        const double J_sr2c_sl0c_increment = J_sr2c_ptr[sl0c] + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;
        J_sr2c_ptr[sl0c] = J_sr2c_sl0c_increment;

        const double Kr2 = this->Kf_[i2 - Ikf[1] + Ikf[9]];
        const double dCrdC0c = C[sr1c]*C[sr2c];
        const double dCrdC1c = C[sr0c]*C[sr2c];
        const double dCrdC2c = C[sr0c]*C[sr1c];
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2*dCrdC0c;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2*dCrdC0c;
        const double J_sr1c_sr0c_increment = J_sr1c_ptr[sr0c] - Kr2*dCrdC0c;
        const double J_sr2c_sr0c_increment = J_sr2c_ptr[sr0c] - Kr2*dCrdC0c;
        const double J_sl0c_sr1c_increment = J_sl0c_ptr[sr1c] + Kr2*dCrdC1c;
        const double J_sr0c_sr1c_increment = J_sr0c_ptr[sr1c] - Kr2*dCrdC1c;
        const double J_sr1c_sr1c_increment = J_sr1c_ptr[sr1c] - Kr2*dCrdC1c;
        const double J_sr2c_sr1c_increment = J_sr2c_ptr[sr1c] - Kr2*dCrdC1c;
        const double J_sl0c_sr2c_increment = J_sl0c_ptr[sr2c] + Kr2*dCrdC2c;
        const double J_sr0c_sr2c_increment = J_sr0c_ptr[sr2c] - Kr2*dCrdC2c;
        const double J_sr1c_sr2c_increment = J_sr1c_ptr[sr2c] - Kr2*dCrdC2c;
        const double J_sr2c_sr2c_increment = J_sr2c_ptr[sr2c] - Kr2*dCrdC2c;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;
        J_sr1c_ptr[sr0c] = J_sr1c_sr0c_increment;
        J_sr2c_ptr[sr0c] = J_sr2c_sr0c_increment;
        J_sl0c_ptr[sr1c] = J_sl0c_sr1c_increment;
        J_sr0c_ptr[sr1c] = J_sr0c_sr1c_increment;
        J_sr1c_ptr[sr1c] = J_sr1c_sr1c_increment;
        J_sr2c_ptr[sr1c] = J_sr2c_sr1c_increment;
        J_sl0c_ptr[sr2c] = J_sl0c_sr2c_increment;
        J_sr0c_ptr[sr2c] = J_sr0c_sr2c_increment;
        J_sr1c_ptr[sr2c] = J_sr1c_sr2c_increment;
        J_sr2c_ptr[sr2c] = J_sr2c_sr2c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double dKrdT2 = this->dKfdT_[i2 - Ikf[1] + Ikf[9]];
        const double CF2 = C[sl0c];
        const double CR2 = C[sr0c]*C[sr1c]*C[sr2c];
        const double dqdT2 = (dKfdT2*CF2)-(dKrdT2*CR2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2;
        const double J_sr2c_T_increment = J_sr2c_ptr[this->nSpecies] + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;
        J_sr2c_ptr[this->nSpecies] = J_sr2c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;

        if(i2>=this->Ikf[2] && i2<this->Ikf[3])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*CF2 - this->dKfdC_[idx2+this->Itbr[4]]*CR2;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2c_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2c_ptr[j+4]),WdMdC4);
                store256d(&J_sl0c_ptr[j+0],sl0v0);
                store256d(&J_sr0c_ptr[j+0],sr0v0);
                store256d(&J_sr1c_ptr[j+0],sr1v0);
                store256d(&J_sr2c_ptr[j+0],sr2v0);
                store256d(&J_sl0c_ptr[j+4],sl0v4);
                store256d(&J_sr0c_ptr[j+4],sr0v4);
                store256d(&J_sr1c_ptr[j+4],sr1v4);
                store256d(&J_sr2c_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv2);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0c_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0c_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1c_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2c_ptr[j]),WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sr0c_ptr[j],sr0v);
                store256d(&J_sr1c_ptr[j],sr1v);
                store256d(&J_sr2c_ptr[j],sr2v);
            }
        }

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->NER13AllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->NER13AllIndex[(k+3)*5+1];
        const unsigned int sr0d = this->NER13AllIndex[(k+3)*5+2];
        const unsigned int sr1d = this->NER13AllIndex[(k+3)*5+3];
        const unsigned int sr2d = this->NER13AllIndex[(k+3)*5+4];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];
        double* __restrict__ J_sr2d_ptr = &ddNdtByVdcTp[sr2d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3;
        const double J_sr2d_sl0d_increment = J_sr2d_ptr[sl0d] + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;
        J_sr2d_ptr[sl0d] = J_sr2d_sl0d_increment;

        const double Kr3 = this->Kf_[i3 - Ikf[1] + Ikf[9]];
        const double dCrdC0d = C[sr1d]*C[sr2d];
        const double dCrdC1d = C[sr0d]*C[sr2d];
        const double dCrdC2d = C[sr0d]*C[sr1d];
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3*dCrdC0d;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3*dCrdC0d;
        const double J_sr1d_sr0d_increment = J_sr1d_ptr[sr0d] - Kr3*dCrdC0d;
        const double J_sr2d_sr0d_increment = J_sr2d_ptr[sr0d] - Kr3*dCrdC0d;
        const double J_sl0d_sr1d_increment = J_sl0d_ptr[sr1d] + Kr3*dCrdC1d;
        const double J_sr0d_sr1d_increment = J_sr0d_ptr[sr1d] - Kr3*dCrdC1d;
        const double J_sr1d_sr1d_increment = J_sr1d_ptr[sr1d] - Kr3*dCrdC1d;
        const double J_sr2d_sr1d_increment = J_sr2d_ptr[sr1d] - Kr3*dCrdC1d;
        const double J_sl0d_sr2d_increment = J_sl0d_ptr[sr2d] + Kr3*dCrdC2d;
        const double J_sr0d_sr2d_increment = J_sr0d_ptr[sr2d] - Kr3*dCrdC2d;
        const double J_sr1d_sr2d_increment = J_sr1d_ptr[sr2d] - Kr3*dCrdC2d;
        const double J_sr2d_sr2d_increment = J_sr2d_ptr[sr2d] - Kr3*dCrdC2d;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;
        J_sr1d_ptr[sr0d] = J_sr1d_sr0d_increment;
        J_sr2d_ptr[sr0d] = J_sr2d_sr0d_increment;
        J_sl0d_ptr[sr1d] = J_sl0d_sr1d_increment;
        J_sr0d_ptr[sr1d] = J_sr0d_sr1d_increment;
        J_sr1d_ptr[sr1d] = J_sr1d_sr1d_increment;
        J_sr2d_ptr[sr1d] = J_sr2d_sr1d_increment;
        J_sl0d_ptr[sr2d] = J_sl0d_sr2d_increment;
        J_sr0d_ptr[sr2d] = J_sr0d_sr2d_increment;
        J_sr1d_ptr[sr2d] = J_sr1d_sr2d_increment;
        J_sr2d_ptr[sr2d] = J_sr2d_sr2d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double dKrdT3 = this->dKfdT_[i3 - Ikf[1] + Ikf[9]];
        const double CF3 = C[sl0d];
        const double CR3 = C[sr0d]*C[sr1d]*C[sr2d];
        const double dqdT3 = (dKfdT3*CF3)-(dKrdT3*CR3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3;
        const double J_sr2d_T_increment = J_sr2d_ptr[this->nSpecies] + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;
        J_sr2d_ptr[this->nSpecies] = J_sr2d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        const double omega_sr2d_increment = dNdtByV[sr2d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
        dNdtByV[sr2d] = omega_sr2d_increment;

        if(i3>=this->Ikf[2] && i3<this->Ikf[3])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*CF3 - this->dKfdC_[idx3+this->Itbr[4]]*CR3;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2d_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2d_ptr[j+4]),WdMdC4);
                store256d(&J_sl0d_ptr[j+0],sl0v0);
                store256d(&J_sr0d_ptr[j+0],sr0v0);
                store256d(&J_sr1d_ptr[j+0],sr1v0);
                store256d(&J_sr2d_ptr[j+0],sr2v0);
                store256d(&J_sl0d_ptr[j+4],sl0v4);
                store256d(&J_sr0d_ptr[j+4],sr0v4);
                store256d(&J_sr1d_ptr[j+4],sr1v4);
                store256d(&J_sr2d_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv3);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0d_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0d_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1d_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2d_ptr[j]),WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sr0d_ptr[j],sr0v);
                store256d(&J_sr1d_ptr[j],sr1v);
                store256d(&J_sr2d_ptr[j],sr2v);
            }
        }
    }
    for(std::size_t k=endNER13-remainNER13; k<endNER13; k=k+1)
    {
        const unsigned int i0 = this->NER13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->NER13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->NER13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->NER13AllIndex[(k+0)*5+4];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];
        double* __restrict__ J_sr2a_ptr = &ddNdtByVdcTp[sr2a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0;
        const double J_sr2a_sl0a_increment = J_sr2a_ptr[sl0a] + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;
        J_sr2a_ptr[sl0a] = J_sr2a_sl0a_increment;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double dCrdC0a = C[sr1a]*C[sr2a];
        const double dCrdC1a = C[sr0a]*C[sr2a];
        const double dCrdC2a = C[sr0a]*C[sr1a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0*dCrdC0a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr2a_sr0a_increment = J_sr2a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0*dCrdC1a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr2a_sr1a_increment = J_sr2a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sl0a_sr2a_increment = J_sl0a_ptr[sr2a] + Kr0*dCrdC2a;
        const double J_sr0a_sr2a_increment = J_sr0a_ptr[sr2a] - Kr0*dCrdC2a;
        const double J_sr1a_sr2a_increment = J_sr1a_ptr[sr2a] - Kr0*dCrdC2a;
        const double J_sr2a_sr2a_increment = J_sr2a_ptr[sr2a] - Kr0*dCrdC2a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sr2a_ptr[sr0a] = J_sr2a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;
        J_sr2a_ptr[sr1a] = J_sr2a_sr1a_increment;
        J_sl0a_ptr[sr2a] = J_sl0a_sr2a_increment;
        J_sr0a_ptr[sr2a] = J_sr0a_sr2a_increment;
        J_sr1a_ptr[sr2a] = J_sr1a_sr2a_increment;
        J_sr2a_ptr[sr2a] = J_sr2a_sr2a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a];
        const double CR0 = C[sr0a]*C[sr1a]*C[sr2a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0;
        const double J_sr2a_T_increment = J_sr2a_ptr[this->nSpecies] + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;
        J_sr2a_ptr[this->nSpecies] = J_sr2a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        if(i0>=this->Ikf[2] && i0<this->Ikf[3])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0 - this->dKfdC_[idx0+this->Itbr[4]]*CR0;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                __m256d sr2v0 = _mm256_add_pd(load256d(&J_sr2a_ptr[j+0]),WdMdC0);
                __m256d sr2v4 = _mm256_add_pd(load256d(&J_sr2a_ptr[j+4]),WdMdC4);
                store256d(&J_sl0a_ptr[j+0],sl0v0);
                store256d(&J_sr0a_ptr[j+0],sr0v0);
                store256d(&J_sr1a_ptr[j+0],sr1v0);
                store256d(&J_sr2a_ptr[j+0],sr2v0);
                store256d(&J_sl0a_ptr[j+4],sl0v4);
                store256d(&J_sr0a_ptr[j+4],sr0v4);
                store256d(&J_sr1a_ptr[j+4],sr1v4);
                store256d(&J_sr2a_ptr[j+4],sr2v4);
            }
            if(remain8==4)
            {
                unsigned int j = this->AlignSpecies-remain8;
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx0*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                __m256d sr2v = _mm256_add_pd(load256d(&J_sr2a_ptr[j]),WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
                store256d(&J_sr2a_ptr[j],sr2v);
            }
        }
    }

    std::size_t endNERDup13 = this->NER13DupSize;
    std::size_t remainNERDup13 = endNERDup13%4;
    for(std::size_t k=0; k<endNERDup13-remainNERDup13; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->NER13DupAllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->NER13DupAllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->NER13DupAllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->NER13DupAllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0 + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double dCrdC0a = C[sr1a]*C[sr1a];
        const double dCrdC1a = 2*C[sr0a]*C[sr1a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0*dCrdC0a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0*dCrdC0a - Kr0*dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0*dCrdC1a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0*dCrdC1a - Kr0*dCrdC1a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a];
        const double CR0 = C[sr0a]*C[sr1a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0 + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[2] && i0<this->Ikf[3])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0 - this->dKfdC_[idx0+this->Itbr[4]]*CR0;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->NER13DupAllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->NER13DupAllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->NER13DupAllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->NER13DupAllIndex[(k+1)*4+3];
        double* __restrict__ J_sl0b_ptr = &ddNdtByVdcTp[sl0b*(this->alignN)];
        double* __restrict__ J_sr0b_ptr = &ddNdtByVdcTp[sr0b*(this->alignN)];
        double* __restrict__ J_sr1b_ptr = &ddNdtByVdcTp[sr1b*(this->alignN)];

        const double Kf1 = this->Kf_[i1];
        const double J_sl0b_sl0b_increment = J_sl0b_ptr[sl0b] - Kf1;
        const double J_sr0b_sl0b_increment = J_sr0b_ptr[sl0b] + Kf1;
        const double J_sr1b_sl0b_increment = J_sr1b_ptr[sl0b] + Kf1 + Kf1;
        J_sl0b_ptr[sl0b] = J_sl0b_sl0b_increment;
        J_sr0b_ptr[sl0b] = J_sr0b_sl0b_increment;
        J_sr1b_ptr[sl0b] = J_sr1b_sl0b_increment;

        const double Kr1 = this->Kf_[i1 - Ikf[1] + Ikf[9]];
        const double dCrdC0b = C[sr1b]*C[sr1b];
        const double dCrdC1b = 2*C[sr0b]*C[sr1b];
        const double J_sl0b_sr0b_increment = J_sl0b_ptr[sr0b] + Kr1*dCrdC0b;
        const double J_sr0b_sr0b_increment = J_sr0b_ptr[sr0b] - Kr1*dCrdC0b;
        const double J_sr1b_sr0b_increment = J_sr1b_ptr[sr0b] - Kr1*dCrdC0b - Kr1*dCrdC0b;
        const double J_sl0b_sr1b_increment = J_sl0b_ptr[sr1b] + Kr1*dCrdC1b;
        const double J_sr0b_sr1b_increment = J_sr0b_ptr[sr1b] - Kr1*dCrdC1b;
        const double J_sr1b_sr1b_increment = J_sr1b_ptr[sr1b] - Kr1*dCrdC1b - Kr1*dCrdC1b;
        J_sl0b_ptr[sr0b] = J_sl0b_sr0b_increment;
        J_sr0b_ptr[sr0b] = J_sr0b_sr0b_increment;
        J_sr1b_ptr[sr0b] = J_sr1b_sr0b_increment;
        J_sl0b_ptr[sr1b] = J_sl0b_sr1b_increment;
        J_sr0b_ptr[sr1b] = J_sr0b_sr1b_increment;
        J_sr1b_ptr[sr1b] = J_sr1b_sr1b_increment;

        const double dKfdT1 = this->dKfdT_[i1];
        const double dKrdT1 = this->dKfdT_[i1 - Ikf[1] + Ikf[9]];
        const double CF1 = C[sl0b];
        const double CR1 = C[sr0b]*C[sr1b]*C[sr1b];
        const double dqdT1 = (dKfdT1*CF1)-(dKrdT1*CR1);
        const double J_sl0b_T_increment = J_sl0b_ptr[this->nSpecies] - dqdT1;
        const double J_sr0b_T_increment = J_sr0b_ptr[this->nSpecies] + dqdT1;
        const double J_sr1b_T_increment = J_sr1b_ptr[this->nSpecies] + dqdT1 + dqdT1;
        J_sl0b_ptr[this->nSpecies] = J_sl0b_T_increment;
        J_sr0b_ptr[this->nSpecies] = J_sr0b_T_increment;
        J_sr1b_ptr[this->nSpecies] = J_sr1b_T_increment;

        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        if(i1>=this->Ikf[2] && i1<this->Ikf[3])
        {
            const unsigned int idx1 = i1 - this->Ikf[2];
            const double W1 = this->dKfdC_[idx1]*CF1 - this->dKfdC_[idx1+this->Itbr[4]]*CR1;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1b_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx1*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv1);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0b_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0b_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1b_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0b_ptr[j],sl0v);
                store256d(&J_sr0b_ptr[j],sr0v);
                store256d(&J_sr1b_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->NER13DupAllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->NER13DupAllIndex[(k+2)*4+1];
        const unsigned int sr0c = this->NER13DupAllIndex[(k+2)*4+2];
        const unsigned int sr1c = this->NER13DupAllIndex[(k+2)*4+3];
        double* __restrict__ J_sl0c_ptr = &ddNdtByVdcTp[sl0c*(this->alignN)];
        double* __restrict__ J_sr0c_ptr = &ddNdtByVdcTp[sr0c*(this->alignN)];
        double* __restrict__ J_sr1c_ptr = &ddNdtByVdcTp[sr1c*(this->alignN)];

        const double Kf2 = this->Kf_[i2];
        const double J_sl0c_sl0c_increment = J_sl0c_ptr[sl0c] - Kf2;
        const double J_sr0c_sl0c_increment = J_sr0c_ptr[sl0c] + Kf2;
        const double J_sr1c_sl0c_increment = J_sr1c_ptr[sl0c] + Kf2 + Kf2;
        J_sl0c_ptr[sl0c] = J_sl0c_sl0c_increment;
        J_sr0c_ptr[sl0c] = J_sr0c_sl0c_increment;
        J_sr1c_ptr[sl0c] = J_sr1c_sl0c_increment;

        const double Kr2 = this->Kf_[i2 - Ikf[1] + Ikf[9]];
        const double dCrdC0c = C[sr1c]*C[sr1c];
        const double dCrdC1c = 2*C[sr0c]*C[sr1c];
        const double J_sl0c_sr0c_increment = J_sl0c_ptr[sr0c] + Kr2*dCrdC0c;
        const double J_sr0c_sr0c_increment = J_sr0c_ptr[sr0c] - Kr2*dCrdC0c;
        const double J_sr1c_sr0c_increment = J_sr1c_ptr[sr0c] - Kr2*dCrdC0c - Kr2*dCrdC0c;
        const double J_sl0c_sr1c_increment = J_sl0c_ptr[sr1c] + Kr2*dCrdC1c;
        const double J_sr0c_sr1c_increment = J_sr0c_ptr[sr1c] - Kr2*dCrdC1c;
        const double J_sr1c_sr1c_increment = J_sr1c_ptr[sr1c] - Kr2*dCrdC1c - Kr2*dCrdC1c;
        J_sl0c_ptr[sr0c] = J_sl0c_sr0c_increment;
        J_sr0c_ptr[sr0c] = J_sr0c_sr0c_increment;
        J_sr1c_ptr[sr0c] = J_sr1c_sr0c_increment;
        J_sl0c_ptr[sr1c] = J_sl0c_sr1c_increment;
        J_sr0c_ptr[sr1c] = J_sr0c_sr1c_increment;
        J_sr1c_ptr[sr1c] = J_sr1c_sr1c_increment;

        const double dKfdT2 = this->dKfdT_[i2];
        const double dKrdT2 = this->dKfdT_[i2 - Ikf[1] + Ikf[9]];
        const double CF2 = C[sl0c];
        const double CR2 = C[sr0c]*C[sr1c]*C[sr1c];
        const double dqdT2 = (dKfdT2*CF2)-(dKrdT2*CR2);
        const double J_sl0c_T_increment = J_sl0c_ptr[this->nSpecies] - dqdT2;
        const double J_sr0c_T_increment = J_sr0c_ptr[this->nSpecies] + dqdT2;
        const double J_sr1c_T_increment = J_sr1c_ptr[this->nSpecies] + dqdT2 + dqdT2;
        J_sl0c_ptr[this->nSpecies] = J_sl0c_T_increment;
        J_sr0c_ptr[this->nSpecies] = J_sr0c_T_increment;
        J_sr1c_ptr[this->nSpecies] = J_sr1c_T_increment;

        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2 + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        if(i2>=this->Ikf[2] && i2<this->Ikf[3])
        {
            const unsigned int idx2 = i2 - this->Ikf[2];
            const double W2 = this->dKfdC_[idx2]*CF2 - this->dKfdC_[idx2+this->Itbr[4]]*CR2;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1c_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx2*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv2);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0c_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0c_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1c_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0c_ptr[j],sl0v);
                store256d(&J_sr0c_ptr[j],sr0v);
                store256d(&J_sr1c_ptr[j],sr1v);
            }
        }

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->NER13DupAllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->NER13DupAllIndex[(k+3)*4+1];
        const unsigned int sr0d = this->NER13DupAllIndex[(k+3)*4+2];
        const unsigned int sr1d = this->NER13DupAllIndex[(k+3)*4+3];
        double* __restrict__ J_sl0d_ptr = &ddNdtByVdcTp[sl0d*(this->alignN)];
        double* __restrict__ J_sr0d_ptr = &ddNdtByVdcTp[sr0d*(this->alignN)];
        double* __restrict__ J_sr1d_ptr = &ddNdtByVdcTp[sr1d*(this->alignN)];

        const double Kf3 = this->Kf_[i3];
        const double J_sl0d_sl0d_increment = J_sl0d_ptr[sl0d] - Kf3;
        const double J_sr0d_sl0d_increment = J_sr0d_ptr[sl0d] + Kf3;
        const double J_sr1d_sl0d_increment = J_sr1d_ptr[sl0d] + Kf3 + Kf3;
        J_sl0d_ptr[sl0d] = J_sl0d_sl0d_increment;
        J_sr0d_ptr[sl0d] = J_sr0d_sl0d_increment;
        J_sr1d_ptr[sl0d] = J_sr1d_sl0d_increment;

        const double Kr3 = this->Kf_[i3 - Ikf[1] + Ikf[9]];
        const double dCrdC0d = C[sr1d]*C[sr1d];
        const double dCrdC1d = 2*C[sr0d]*C[sr1d];
        const double J_sl0d_sr0d_increment = J_sl0d_ptr[sr0d] + Kr3*dCrdC0d;
        const double J_sr0d_sr0d_increment = J_sr0d_ptr[sr0d] - Kr3*dCrdC0d;
        const double J_sr1d_sr0d_increment = J_sr1d_ptr[sr0d] - Kr3*dCrdC0d - Kr3*dCrdC0d;
        const double J_sl0d_sr1d_increment = J_sl0d_ptr[sr1d] + Kr3*dCrdC1d;
        const double J_sr0d_sr1d_increment = J_sr0d_ptr[sr1d] - Kr3*dCrdC1d;
        const double J_sr1d_sr1d_increment = J_sr1d_ptr[sr1d] - Kr3*dCrdC1d - Kr3*dCrdC1d;
        J_sl0d_ptr[sr0d] = J_sl0d_sr0d_increment;
        J_sr0d_ptr[sr0d] = J_sr0d_sr0d_increment;
        J_sr1d_ptr[sr0d] = J_sr1d_sr0d_increment;
        J_sl0d_ptr[sr1d] = J_sl0d_sr1d_increment;
        J_sr0d_ptr[sr1d] = J_sr0d_sr1d_increment;
        J_sr1d_ptr[sr1d] = J_sr1d_sr1d_increment;

        const double dKfdT3 = this->dKfdT_[i3];
        const double dKrdT3 = this->dKfdT_[i3 - Ikf[1] + Ikf[9]];
        const double CF3 = C[sl0d];
        const double CR3 = C[sr0d]*C[sr1d]*C[sr1d];
        const double dqdT3 = (dKfdT3*CF3)-(dKrdT3*CR3);
        const double J_sl0d_T_increment = J_sl0d_ptr[this->nSpecies] - dqdT3;
        const double J_sr0d_T_increment = J_sr0d_ptr[this->nSpecies] + dqdT3;
        const double J_sr1d_T_increment = J_sr1d_ptr[this->nSpecies] + dqdT3 + dqdT3;
        J_sl0d_ptr[this->nSpecies] = J_sl0d_T_increment;
        J_sr0d_ptr[this->nSpecies] = J_sr0d_T_increment;
        J_sr1d_ptr[this->nSpecies] = J_sr1d_T_increment;

        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3 + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;

        if(i3>=this->Ikf[2] && i3<this->Ikf[3])
        {
            const unsigned int idx3 = i3 - this->Ikf[2];
            const double W3 = this->dKfdC_[idx3]*CF3 - this->dKfdC_[idx3+this->Itbr[4]]*CR3;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1d_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[idx3*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv3);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0d_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0d_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1d_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0d_ptr[j],sl0v);
                store256d(&J_sr0d_ptr[j],sr0v);
                store256d(&J_sr1d_ptr[j],sr1v);
            }
        }
    }
    for(std::size_t k=endNERDup13-remainNERDup13; k<endNERDup13; k=k+1)
    {
        const unsigned int i0 = this->NER13DupAllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->NER13DupAllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->NER13DupAllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->NER13DupAllIndex[(k+0)*4+3];
        double* __restrict__ J_sl0a_ptr = &ddNdtByVdcTp[sl0a*(this->alignN)];
        double* __restrict__ J_sr0a_ptr = &ddNdtByVdcTp[sr0a*(this->alignN)];
        double* __restrict__ J_sr1a_ptr = &ddNdtByVdcTp[sr1a*(this->alignN)];

        const double Kf0 = this->Kf_[i0];
        const double J_sl0a_sl0a_increment = J_sl0a_ptr[sl0a] - Kf0;
        const double J_sr0a_sl0a_increment = J_sr0a_ptr[sl0a] + Kf0;
        const double J_sr1a_sl0a_increment = J_sr1a_ptr[sl0a] + Kf0 + Kf0;
        J_sl0a_ptr[sl0a] = J_sl0a_sl0a_increment;
        J_sr0a_ptr[sl0a] = J_sr0a_sl0a_increment;
        J_sr1a_ptr[sl0a] = J_sr1a_sl0a_increment;

        const double Kr0 = this->Kf_[i0 - Ikf[1] + Ikf[9]];
        const double dCrdC0a = C[sr1a]*C[sr1a];
        const double dCrdC1a = 2*C[sr0a]*C[sr1a];
        const double J_sl0a_sr0a_increment = J_sl0a_ptr[sr0a] + Kr0*dCrdC0a;
        const double J_sr0a_sr0a_increment = J_sr0a_ptr[sr0a] - Kr0*dCrdC0a;
        const double J_sr1a_sr0a_increment = J_sr1a_ptr[sr0a] - Kr0*dCrdC0a - Kr0*dCrdC0a;
        const double J_sl0a_sr1a_increment = J_sl0a_ptr[sr1a] + Kr0*dCrdC1a;
        const double J_sr0a_sr1a_increment = J_sr0a_ptr[sr1a] - Kr0*dCrdC1a;
        const double J_sr1a_sr1a_increment = J_sr1a_ptr[sr1a] - Kr0*dCrdC1a - Kr0*dCrdC1a;
        J_sl0a_ptr[sr0a] = J_sl0a_sr0a_increment;
        J_sr0a_ptr[sr0a] = J_sr0a_sr0a_increment;
        J_sr1a_ptr[sr0a] = J_sr1a_sr0a_increment;
        J_sl0a_ptr[sr1a] = J_sl0a_sr1a_increment;
        J_sr0a_ptr[sr1a] = J_sr0a_sr1a_increment;
        J_sr1a_ptr[sr1a] = J_sr1a_sr1a_increment;

        const double dKfdT0 = this->dKfdT_[i0];
        const double dKrdT0 = this->dKfdT_[i0 - Ikf[1] + Ikf[9]];
        const double CF0 = C[sl0a];
        const double CR0 = C[sr0a]*C[sr1a]*C[sr1a];
        const double dqdT0 = (dKfdT0*CF0)-(dKrdT0*CR0);
        const double J_sl0a_T_increment = J_sl0a_ptr[this->nSpecies] - dqdT0;
        const double J_sr0a_T_increment = J_sr0a_ptr[this->nSpecies] + dqdT0;
        const double J_sr1a_T_increment = J_sr1a_ptr[this->nSpecies] + dqdT0 + dqdT0;
        J_sl0a_ptr[this->nSpecies] = J_sl0a_T_increment;
        J_sr0a_ptr[this->nSpecies] = J_sr0a_T_increment;
        J_sr1a_ptr[this->nSpecies] = J_sr1a_T_increment;

        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        if(i0>=this->Ikf[2] && i0<this->Ikf[3])
        {
            const unsigned int idx0 = i0 - this->Ikf[2];
            const double W0 = this->dKfdC_[idx0]*CF0 - this->dKfdC_[idx0+this->Itbr[4]]*CR0;
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
                __m256d sr1v0 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+0]),WdMdC0);
                __m256d sr1v4 = _mm256_add_pd(load256d(&J_sr1a_ptr[j+4]),WdMdC4);
                sr1v0 = _mm256_add_pd(sr1v0,WdMdC0);
                sr1v4 = _mm256_add_pd(sr1v4,WdMdC4);
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
                __m256d WdMdC = _mm256_mul_pd(dMdC,Wv0);
                __m256d sl0v = _mm256_sub_pd(load256d(&J_sl0a_ptr[j]),WdMdC);
                __m256d sr0v = _mm256_add_pd(load256d(&J_sr0a_ptr[j]),WdMdC);
                __m256d sr1v = _mm256_add_pd(load256d(&J_sr1a_ptr[j]),WdMdC);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&J_sl0a_ptr[j],sl0v);
                store256d(&J_sr0a_ptr[j],sr0v);
                store256d(&J_sr1a_ptr[j],sr1v);
            }
        }
    }
}
