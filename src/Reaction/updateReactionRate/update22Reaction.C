/*---------------------------------------------------------------------------*\
  Description
      Update the net rates of production (dNdtByV) for two-reactant /
      two-product (2-2) reactions, e.g. A + A -> B + B and A + B -> C + D.

      RR : reversible reaction; the reverse rate constant is obtained from
           the equilibrium constant
      IR : irreversible reaction; the reverse rate constant is zero
      NER: non-equilibrium reaction; the reverse rate constant is computed
           with the Arrhenius form instead of the equilibrium constant

  Functions
      RF22RR : reversible 2-2 reactions, incl. duplicate reactants/products
      RF22IR : irreversible 2-2 reactions, incl. duplicate reactants/products
      RF22NER: non-equilibrium 2-2 reactions

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. Standard C++ library headers
//---------------------------------
#include <algorithm>

//---------------------------------
// 2. FastChemistry headers
//---------------------------------
#include "OptReaction.H"

//=============================================================================//
void
FastChemistry::OptReaction::RF22RR
(
    const double* __restrict__ c,
    double* __restrict__ dNdtByV,
    const double* __restrict__ ExpNegGbyRT
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
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]);
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR22AllIndex[k*5+5];
        const unsigned int sl0b = this->RR22AllIndex[k*5+6];
        const unsigned int sl1b = this->RR22AllIndex[k*5+7];
        const unsigned int sr0b = this->RR22AllIndex[k*5+8];
        const unsigned int sr1b = this->RR22AllIndex[k*5+9];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b])*(invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]);
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR22AllIndex[k*5+10];
        const unsigned int sl0c = this->RR22AllIndex[k*5+11];
        const unsigned int sl1c = this->RR22AllIndex[k*5+12];
        const unsigned int sr0c = this->RR22AllIndex[k*5+13];
        const unsigned int sr1c = this->RR22AllIndex[k*5+14];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c])*(invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c]);
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR22AllIndex[k*5+15];
        const unsigned int sl0d = this->RR22AllIndex[k*5+16];
        const unsigned int sl1d = this->RR22AllIndex[k*5+17];
        const unsigned int sr0d = this->RR22AllIndex[k*5+18];
        const unsigned int sr1d = this->RR22AllIndex[k*5+19];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        double invKc3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d])*(invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d]);
        invKc3 = std::min(invKc3, invKcLimiter);
        Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d]*c[sl1d];
        const double CR3 = c[sr0d]*c[sr1d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
    }
    if(remainRR22==1)
    {
        std::size_t k = endRR22-1;
        const unsigned int i0 = this->RR22AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR22AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR22AllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->RR22AllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->RR22AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]);
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
    }
    else if(remainRR22==2)
    {
        std::size_t k = endRR22-2;
        const unsigned int i0 = this->RR22AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR22AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR22AllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->RR22AllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->RR22AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]);
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        const unsigned int i1 = this->RR22AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR22AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->RR22AllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->RR22AllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->RR22AllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b])*(invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]);
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
    }
    else if(remainRR22==3)
    {
        std::size_t k = endRR22-3;
        const unsigned int i0 = this->RR22AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR22AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR22AllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->RR22AllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->RR22AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]);
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        const unsigned int i1 = this->RR22AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR22AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->RR22AllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->RR22AllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->RR22AllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b])*(invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]);
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        const unsigned int i2 = this->RR22AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->RR22AllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->RR22AllIndex[(k+2)*5+2];
        const unsigned int sr0c = this->RR22AllIndex[(k+2)*5+3];
        const unsigned int sr1c = this->RR22AllIndex[(k+2)*5+4];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c])*(invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c]);
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
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
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]);
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR22DupAllIndex[k*5+5];
        const unsigned int sl0b = this->RR22DupAllIndex[k*5+6];
        const unsigned int sl1b = this->RR22DupAllIndex[k*5+7];
        const unsigned int sr0b = this->RR22DupAllIndex[k*5+8];
        const unsigned int sr1b = this->RR22DupAllIndex[k*5+9];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b])*(invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]);
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR22DupAllIndex[k*5+10];
        const unsigned int sl0c = this->RR22DupAllIndex[k*5+11];
        const unsigned int sl1c = this->RR22DupAllIndex[k*5+12];
        const unsigned int sr0c = this->RR22DupAllIndex[k*5+13];
        const unsigned int sr1c = this->RR22DupAllIndex[k*5+14];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c])*(invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c]);
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR22DupAllIndex[k*5+15];
        const unsigned int sl0d = this->RR22DupAllIndex[k*5+16];
        const unsigned int sl1d = this->RR22DupAllIndex[k*5+17];
        const unsigned int sr0d = this->RR22DupAllIndex[k*5+18];
        const unsigned int sr1d = this->RR22DupAllIndex[k*5+19];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        double invKc3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d])*(invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d]);
        invKc3 = std::min(invKc3, invKcLimiter);
        Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d]*c[sl1d];
        const double CR3 = c[sr0d]*c[sr1d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
    }
    if(remainRRDup22==1)
    {
        std::size_t k = endRRDup22-1;
        const unsigned int i0 = this->RR22DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR22DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR22DupAllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->RR22DupAllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->RR22DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = this->Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]);
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
    }
    else if(remainRRDup22==2)
    {
        std::size_t k = endRRDup22-2;
        const unsigned int i0 = this->RR22DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR22DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR22DupAllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->RR22DupAllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->RR22DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = this->Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]);
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->RR22DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR22DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->RR22DupAllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->RR22DupAllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->RR22DupAllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = this->Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b])*(invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]);
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
    }
    else if(remainRRDup22==3)
    {
        std::size_t k = endRRDup22-3;
        const unsigned int i0 = this->RR22DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR22DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR22DupAllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->RR22DupAllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->RR22DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = this->Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]);
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->RR22DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR22DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->RR22DupAllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->RR22DupAllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->RR22DupAllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = this->Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b])*(invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]);
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        const unsigned int i2 = this->RR22DupAllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->RR22DupAllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->RR22DupAllIndex[(k+2)*5+2];
        const unsigned int sr0c = this->RR22DupAllIndex[(k+2)*5+3];
        const unsigned int sr1c = this->RR22DupAllIndex[(k+2)*5+4];
        double Kr2 = 0;
        const double Kf2 = this->Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c])*(invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c]);
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
    }
    }
}

void
FastChemistry::OptReaction::RF22IR
(
    const double* __restrict__ c,
    double* __restrict__ dNdtByV,
    const double* __restrict__ ExpNegGbyRT
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
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR22AllIndex[k*5+5];
        const unsigned int sl0b = this->IR22AllIndex[k*5+6];
        const unsigned int sl1b = this->IR22AllIndex[k*5+7];
        const unsigned int sr0b = this->IR22AllIndex[k*5+8];
        const unsigned int sr1b = this->IR22AllIndex[k*5+9];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR22AllIndex[k*5+10];
        const unsigned int sl0c = this->IR22AllIndex[k*5+11];
        const unsigned int sl1c = this->IR22AllIndex[k*5+12];
        const unsigned int sr0c = this->IR22AllIndex[k*5+13];
        const unsigned int sr1c = this->IR22AllIndex[k*5+14];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR22AllIndex[k*5+15];
        const unsigned int sl0d = this->IR22AllIndex[k*5+16];
        const unsigned int sl1d = this->IR22AllIndex[k*5+17];
        const unsigned int sr0d = this->IR22AllIndex[k*5+18];
        const unsigned int sr1d = this->IR22AllIndex[k*5+19];
        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d]*c[sl1d];
        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
    }

    if(remainIR22==1)
    {
        std::size_t k = endIR22-1;
        const unsigned int i0 = this->IR22AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR22AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR22AllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->IR22AllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->IR22AllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
    }
    else if(remainIR22==2)
    {
        std::size_t k = endIR22-2;
        const unsigned int i0 = this->IR22AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR22AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR22AllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->IR22AllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->IR22AllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        const unsigned int i1 = this->IR22AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR22AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->IR22AllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->IR22AllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->IR22AllIndex[(k+1)*5+4];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
    }
    else if(remainIR22==3)
    {
        std::size_t k = endIR22-3;
        const unsigned int i0 = this->IR22AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR22AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR22AllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->IR22AllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->IR22AllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        const unsigned int i1 = this->IR22AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR22AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->IR22AllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->IR22AllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->IR22AllIndex[(k+1)*5+4];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        const unsigned int i2 = this->IR22AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->IR22AllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->IR22AllIndex[(k+2)*5+2];
        const unsigned int sr0c = this->IR22AllIndex[(k+2)*5+3];
        const unsigned int sr1c = this->IR22AllIndex[(k+2)*5+4];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
    }

    {
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

        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR22DupAllIndex[k*5+5];
        const unsigned int sl0b = this->IR22DupAllIndex[k*5+6];
        const unsigned int sl1b = this->IR22DupAllIndex[k*5+7];
        const unsigned int sr0b = this->IR22DupAllIndex[k*5+8];
        const unsigned int sr1b = this->IR22DupAllIndex[k*5+9];

        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR22DupAllIndex[k*5+10];
        const unsigned int sl0c = this->IR22DupAllIndex[k*5+11];
        const unsigned int sl1c = this->IR22DupAllIndex[k*5+12];
        const unsigned int sr0c = this->IR22DupAllIndex[k*5+13];
        const unsigned int sr1c = this->IR22DupAllIndex[k*5+14];

        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c];
        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR22DupAllIndex[k*5+15];
        const unsigned int sl0d = this->IR22DupAllIndex[k*5+16];
        const unsigned int sl1d = this->IR22DupAllIndex[k*5+17];
        const unsigned int sr0d = this->IR22DupAllIndex[k*5+18];
        const unsigned int sr1d = this->IR22DupAllIndex[k*5+19];

        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d]*c[sl1d];
        const double q3 = (Kf3*CF3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
    }

    if(remainIR22Dup==1)
    {
        std::size_t k = endIR22Dup-1;
        const unsigned int i0 = this->IR22DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR22DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR22DupAllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->IR22DupAllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->IR22DupAllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
    }
    else if(remainIR22Dup==2)
    {
        std::size_t k = endIR22Dup-2;
        const unsigned int i0 = this->IR22DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR22DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR22DupAllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->IR22DupAllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->IR22DupAllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->IR22DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR22DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->IR22DupAllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->IR22DupAllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->IR22DupAllIndex[(k+1)*5+4];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
    }
    else if(remainIR22Dup==3)
    {
        std::size_t k = endIR22Dup-3;
        const unsigned int i0 = this->IR22DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR22DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR22DupAllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->IR22DupAllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->IR22DupAllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->IR22DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR22DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->IR22DupAllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->IR22DupAllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->IR22DupAllIndex[(k+1)*5+4];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        const unsigned int i2 = this->IR22DupAllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->IR22DupAllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->IR22DupAllIndex[(k+2)*5+2];
        const unsigned int sr0c = this->IR22DupAllIndex[(k+2)*5+3];
        const unsigned int sr1c = this->IR22DupAllIndex[(k+2)*5+4];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c];
        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
    }
    }
}

void
FastChemistry::OptReaction::RF22NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
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
        double Kr0 = 0;
        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->NER22AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER22AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER22AllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->NER22AllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->NER22AllIndex[(k+1)*5+4];
        double Kr1 = 0;
        double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        const unsigned int i2 = this->NER22AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER22AllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->NER22AllIndex[(k+2)*5+2];
        const unsigned int sr0c = this->NER22AllIndex[(k+2)*5+3];
        const unsigned int sr1c = this->NER22AllIndex[(k+2)*5+4];
        double Kr2 = 0;
        double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        const unsigned int i3 = this->NER22AllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->NER22AllIndex[(k+3)*5+1];
        const unsigned int sl1d = this->NER22AllIndex[(k+3)*5+2];
        const unsigned int sr0d = this->NER22AllIndex[(k+3)*5+3];
        const unsigned int sr1d = this->NER22AllIndex[(k+3)*5+4];
        double Kr3 = 0;
        double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d]*c[sl1d];
        const double CR3 = c[sr0d]*c[sr1d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
    }
    if(remainNER22==1)
    {
        std::size_t k = endNER22-1;
        const unsigned int i0 = this->NER22AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER22AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER22AllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->NER22AllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->NER22AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
    }
    else if(remainNER22==2)
    {
        std::size_t k = endNER22-2;
        const unsigned int i0 = this->NER22AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER22AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER22AllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->NER22AllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->NER22AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->NER22AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER22AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER22AllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->NER22AllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->NER22AllIndex[(k+1)*5+4];
        double Kr1 = 0;
        double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
    }
    else if(remainNER22==3)
    {
        std::size_t k = endNER22-3;
        const unsigned int i0 = this->NER22AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER22AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER22AllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->NER22AllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->NER22AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->NER22AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER22AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER22AllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->NER22AllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->NER22AllIndex[(k+1)*5+4];
        double Kr1 = 0;
        double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        const unsigned int i2 = this->NER22AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER22AllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->NER22AllIndex[(k+2)*5+2];
        const unsigned int sr0c = this->NER22AllIndex[(k+2)*5+3];
        const unsigned int sr1c = this->NER22AllIndex[(k+2)*5+4];
        double Kr2 = 0;
        double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
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
        double Kr0 = 0;
        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->NER22DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER22DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER22DupAllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->NER22DupAllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->NER22DupAllIndex[(k+1)*5+4];
        double Kr1 = 0;
        double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        const unsigned int i2 = this->NER22DupAllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER22DupAllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->NER22DupAllIndex[(k+2)*5+2];
        const unsigned int sr0c = this->NER22DupAllIndex[(k+2)*5+3];
        const unsigned int sr1c = this->NER22DupAllIndex[(k+2)*5+4];
        double Kr2 = 0;
        double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        const unsigned int i3 = this->NER22DupAllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->NER22DupAllIndex[(k+3)*5+1];
        const unsigned int sl1d = this->NER22DupAllIndex[(k+3)*5+2];
        const unsigned int sr0d = this->NER22DupAllIndex[(k+3)*5+3];
        const unsigned int sr1d = this->NER22DupAllIndex[(k+3)*5+4];
        double Kr3 = 0;
        double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d]*c[sl1d];
        const double CR3 = c[sr0d]*c[sr1d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
    }
    if(remainNER22Dup==1)
    {
        std::size_t k = endNER22Dup-1;
        const unsigned int i0 = this->NER22DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER22DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER22DupAllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->NER22DupAllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->NER22DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
    }
    else if(remainNER22Dup==2)
    {
        std::size_t k = endNER22Dup-2;
        const unsigned int i0 = this->NER22DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER22DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER22DupAllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->NER22DupAllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->NER22DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->NER22DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER22DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER22DupAllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->NER22DupAllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->NER22DupAllIndex[(k+1)*5+4];
        double Kr1 = 0;
        double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
    }
    else if(remainNER22Dup==3)
    {
        std::size_t k = endNER22Dup-3;
        const unsigned int i0 = this->NER22DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER22DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER22DupAllIndex[(k+0)*5+2];
        const unsigned int sr0a = this->NER22DupAllIndex[(k+0)*5+3];
        const unsigned int sr1a = this->NER22DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->NER22DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER22DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER22DupAllIndex[(k+1)*5+2];
        const unsigned int sr0b = this->NER22DupAllIndex[(k+1)*5+3];
        const unsigned int sr1b = this->NER22DupAllIndex[(k+1)*5+4];
        double Kr1 = 0;
        double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        const unsigned int i2 = this->NER22DupAllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER22DupAllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->NER22DupAllIndex[(k+2)*5+2];
        const unsigned int sr0c = this->NER22DupAllIndex[(k+2)*5+3];
        const unsigned int sr1c = this->NER22DupAllIndex[(k+2)*5+4];
        double Kr2 = 0;
        double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
    }
}
