/*---------------------------------------------------------------------------*\
  Description
      Update the net rates of production (dNdtByV) for three-reactant /
      two-product (3-2) reactions, e.g. A + A + A -> B + B and A + B + C -> D + E.

      RR : reversible reaction; the reverse rate constant is obtained from
           the equilibrium constant
      IR : irreversible reaction; the reverse rate constant is zero
      NER: non-equilibrium reaction; the reverse rate constant is computed
           with the Arrhenius form instead of the equilibrium constant

  Functions
      RF32RR : reversible 3-2 reactions, incl. duplicate species (A+A+A=B+B)
      RF32IR : irreversible 3-2 reactions, incl. duplicate species (A+A+A=B+B)
      RF32NER: non-equilibrium 3-2 reactions

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
FastChemistry::OptReaction::RF32RR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
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
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
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

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->RR32AllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->RR32AllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->RR32AllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->RR32AllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->RR32AllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->RR32AllIndex[(k+1)*6+5];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b];
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

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->RR32AllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->RR32AllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->RR32AllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->RR32AllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->RR32AllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->RR32AllIndex[(k+2)*6+5];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*
            (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c];
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

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->RR32AllIndex[(k+3)*6+0];
        const unsigned int sl0d = this->RR32AllIndex[(k+3)*6+1];
        const unsigned int sl1d = this->RR32AllIndex[(k+3)*6+2];
        const unsigned int sl2d = this->RR32AllIndex[(k+3)*6+3];
        const unsigned int sr0d = this->RR32AllIndex[(k+3)*6+4];
        const unsigned int sr1d = this->RR32AllIndex[(k+3)*6+5];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        double invKc3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]*ExpNegGbyRT[sl2d])*
            (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d])*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d]*c[sr1d];
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
    }
    if(remainRR32==1)
    {
        std::size_t k = endRR32-1;
        const unsigned int i0 = this->RR32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->RR32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->RR32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->RR32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->RR32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->RR32AllIndex[(k+0)*6+5];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
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
    }
    else if(remainRR32==2)
    {
        std::size_t k = endRR32-2;
        const unsigned int i0 = this->RR32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->RR32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->RR32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->RR32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->RR32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->RR32AllIndex[(k+0)*6+5];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
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

        const unsigned int i1 = this->RR32AllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->RR32AllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->RR32AllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->RR32AllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->RR32AllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->RR32AllIndex[(k+1)*6+5];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b];
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
    }
    else if(remainRR32==3)
    {
        std::size_t k = endRR32-3;
        const unsigned int i0 = this->RR32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->RR32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->RR32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->RR32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->RR32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->RR32AllIndex[(k+0)*6+5];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
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

        const unsigned int i1 = this->RR32AllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->RR32AllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->RR32AllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->RR32AllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->RR32AllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->RR32AllIndex[(k+1)*6+5];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b];
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

        const unsigned int i2 = this->RR32AllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->RR32AllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->RR32AllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->RR32AllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->RR32AllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->RR32AllIndex[(k+2)*6+5];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*
            (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c];
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
    }
    {
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
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->RR32DupAllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->RR32DupAllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->RR32DupAllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->RR32DupAllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->RR32DupAllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->RR32DupAllIndex[(k+1)*6+5];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->RR32DupAllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->RR32DupAllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->RR32DupAllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->RR32DupAllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->RR32DupAllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->RR32DupAllIndex[(k+2)*6+5];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*
            (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->RR32DupAllIndex[(k+3)*6+0];
        const unsigned int sl0d = this->RR32DupAllIndex[(k+3)*6+1];
        const unsigned int sl1d = this->RR32DupAllIndex[(k+3)*6+2];
        const unsigned int sl2d = this->RR32DupAllIndex[(k+3)*6+3];
        const unsigned int sr0d = this->RR32DupAllIndex[(k+3)*6+4];
        const unsigned int sr1d = this->RR32DupAllIndex[(k+3)*6+5];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        double invKc3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]*ExpNegGbyRT[sl2d])*
            (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d])*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d]*c[sr1d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
    }
    if(remainRR32Dup==1)
    {
        std::size_t k = endRR32Dup-1;
        const unsigned int i0 = this->RR32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->RR32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->RR32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->RR32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->RR32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->RR32DupAllIndex[(k+0)*6+5];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
    }
    else if(remainRR32Dup==2)
    {
        std::size_t k = endRR32Dup-2;
        const unsigned int i0 = this->RR32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->RR32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->RR32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->RR32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->RR32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->RR32DupAllIndex[(k+0)*6+5];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->RR32DupAllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->RR32DupAllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->RR32DupAllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->RR32DupAllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->RR32DupAllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->RR32DupAllIndex[(k+1)*6+5];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
    }
    else if(remainRR32Dup==3)
    {
        std::size_t k = endRR32Dup-3;
        const unsigned int i0 = this->RR32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->RR32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->RR32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->RR32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->RR32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->RR32DupAllIndex[(k+0)*6+5];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->RR32DupAllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->RR32DupAllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->RR32DupAllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->RR32DupAllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->RR32DupAllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->RR32DupAllIndex[(k+1)*6+5];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        const unsigned int i2 = this->RR32DupAllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->RR32DupAllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->RR32DupAllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->RR32DupAllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->RR32DupAllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->RR32DupAllIndex[(k+2)*6+5];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*
            (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
    }
    }
}


void 
FastChemistry::OptReaction::RF32IR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
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
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
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

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->IR32AllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->IR32AllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->IR32AllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->IR32AllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->IR32AllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->IR32AllIndex[(k+1)*6+5];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
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

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->IR32AllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->IR32AllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->IR32AllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->IR32AllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->IR32AllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->IR32AllIndex[(k+2)*6+5];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
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

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->IR32AllIndex[(k+3)*6+0];
        const unsigned int sl0d = this->IR32AllIndex[(k+3)*6+1];
        const unsigned int sl1d = this->IR32AllIndex[(k+3)*6+2];
        const unsigned int sl2d = this->IR32AllIndex[(k+3)*6+3];
        const unsigned int sr0d = this->IR32AllIndex[(k+3)*6+4];
        const unsigned int sr1d = this->IR32AllIndex[(k+3)*6+5];
        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
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
    }
    if(remainIR32==1)
    {
        std::size_t k = endIR32-1;
        const unsigned int i0 = this->IR32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->IR32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->IR32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->IR32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->IR32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->IR32AllIndex[(k+0)*6+5];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
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
    }
    else if(remainIR32==2)
    {
        std::size_t k = endIR32-2;
        const unsigned int i0 = this->IR32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->IR32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->IR32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->IR32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->IR32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->IR32AllIndex[(k+0)*6+5];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
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

        const unsigned int i1 = this->IR32AllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->IR32AllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->IR32AllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->IR32AllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->IR32AllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->IR32AllIndex[(k+1)*6+5];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
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
    }
    else if(remainIR32==3)
    {
        std::size_t k = endIR32-3;
        const unsigned int i0 = this->IR32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->IR32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->IR32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->IR32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->IR32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->IR32AllIndex[(k+0)*6+5];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
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

        const unsigned int i1 = this->IR32AllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->IR32AllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->IR32AllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->IR32AllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->IR32AllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->IR32AllIndex[(k+1)*6+5];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
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

        const unsigned int i2 = this->IR32AllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->IR32AllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->IR32AllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->IR32AllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->IR32AllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->IR32AllIndex[(k+2)*6+5];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
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
    }
    {
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
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->IR32DupAllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->IR32DupAllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->IR32DupAllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->IR32DupAllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->IR32DupAllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->IR32DupAllIndex[(k+1)*6+5];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->IR32DupAllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->IR32DupAllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->IR32DupAllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->IR32DupAllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->IR32DupAllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->IR32DupAllIndex[(k+2)*6+5];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->IR32DupAllIndex[(k+3)*6+0];
        const unsigned int sl0d = this->IR32DupAllIndex[(k+3)*6+1];
        const unsigned int sl1d = this->IR32DupAllIndex[(k+3)*6+2];
        const unsigned int sl2d = this->IR32DupAllIndex[(k+3)*6+3];
        const unsigned int sr0d = this->IR32DupAllIndex[(k+3)*6+4];
        const unsigned int sr1d = this->IR32DupAllIndex[(k+3)*6+5];
        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double q3 = (Kf3*CF3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
    }
    if(remainIR32Dup==1)
    {
        std::size_t k = endIR32Dup-1;
        const unsigned int i0 = this->IR32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->IR32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->IR32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->IR32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->IR32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->IR32DupAllIndex[(k+0)*6+5];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
    }
    else if(remainIR32Dup==2)
    {
        std::size_t k = endIR32Dup-2;
        const unsigned int i0 = this->IR32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->IR32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->IR32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->IR32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->IR32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->IR32DupAllIndex[(k+0)*6+5];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->IR32DupAllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->IR32DupAllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->IR32DupAllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->IR32DupAllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->IR32DupAllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->IR32DupAllIndex[(k+1)*6+5];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
    }
    else if(remainIR32Dup==3)
    {
        std::size_t k = endIR32Dup-3;
        const unsigned int i0 = this->IR32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->IR32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->IR32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->IR32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->IR32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->IR32DupAllIndex[(k+0)*6+5];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->IR32DupAllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->IR32DupAllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->IR32DupAllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->IR32DupAllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->IR32DupAllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->IR32DupAllIndex[(k+1)*6+5];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        const unsigned int i2 = this->IR32DupAllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->IR32DupAllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->IR32DupAllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->IR32DupAllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->IR32DupAllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->IR32DupAllIndex[(k+2)*6+5];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
    }
    }
}


void 
FastChemistry::OptReaction::RF32NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
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
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
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

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->NER32AllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->NER32AllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->NER32AllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->NER32AllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->NER32AllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->NER32AllIndex[(k+1)*6+5];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b];
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

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->NER32AllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->NER32AllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->NER32AllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->NER32AllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->NER32AllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->NER32AllIndex[(k+2)*6+5];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c];
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

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->NER32AllIndex[(k+3)*6+0];
        const unsigned int sl0d = this->NER32AllIndex[(k+3)*6+1];
        const unsigned int sl1d = this->NER32AllIndex[(k+3)*6+2];
        const unsigned int sl2d = this->NER32AllIndex[(k+3)*6+3];
        const unsigned int sr0d = this->NER32AllIndex[(k+3)*6+4];
        const unsigned int sr1d = this->NER32AllIndex[(k+3)*6+5];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d]*c[sr1d];
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
    }
    if(remainNER32==1)
    {
        std::size_t k = endNER32-1;
        const unsigned int i0 = this->NER32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->NER32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->NER32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->NER32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->NER32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->NER32AllIndex[(k+0)*6+5];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
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
    }
    else if(remainNER32==2)
    {
        std::size_t k = endNER32-2;
        const unsigned int i0 = this->NER32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->NER32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->NER32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->NER32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->NER32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->NER32AllIndex[(k+0)*6+5];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
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

        const unsigned int i1 = this->NER32AllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->NER32AllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->NER32AllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->NER32AllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->NER32AllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->NER32AllIndex[(k+1)*6+5];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b];
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
    }
    else if(remainNER32==3)
    {
        std::size_t k = endNER32-3;
        const unsigned int i0 = this->NER32AllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->NER32AllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->NER32AllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->NER32AllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->NER32AllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->NER32AllIndex[(k+0)*6+5];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
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

        const unsigned int i1 = this->NER32AllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->NER32AllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->NER32AllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->NER32AllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->NER32AllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->NER32AllIndex[(k+1)*6+5];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b];
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

        const unsigned int i2 = this->NER32AllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->NER32AllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->NER32AllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->NER32AllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->NER32AllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->NER32AllIndex[(k+2)*6+5];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c];
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
    }
    {
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
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->NER32DupAllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->NER32DupAllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->NER32DupAllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->NER32DupAllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->NER32DupAllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->NER32DupAllIndex[(k+1)*6+5];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->NER32DupAllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->NER32DupAllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->NER32DupAllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->NER32DupAllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->NER32DupAllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->NER32DupAllIndex[(k+2)*6+5];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->NER32DupAllIndex[(k+3)*6+0];
        const unsigned int sl0d = this->NER32DupAllIndex[(k+3)*6+1];
        const unsigned int sl1d = this->NER32DupAllIndex[(k+3)*6+2];
        const unsigned int sl2d = this->NER32DupAllIndex[(k+3)*6+3];
        const unsigned int sr0d = this->NER32DupAllIndex[(k+3)*6+4];
        const unsigned int sr1d = this->NER32DupAllIndex[(k+3)*6+5];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d]*c[sr1d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
    }
    if(remainNER32Dup==1)
    {
        std::size_t k = endNER32Dup-1;
        const unsigned int i0 = this->NER32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->NER32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->NER32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->NER32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->NER32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->NER32DupAllIndex[(k+0)*6+5];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
    }
    else if(remainNER32Dup==2)
    {
        std::size_t k = endNER32Dup-2;
        const unsigned int i0 = this->NER32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->NER32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->NER32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->NER32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->NER32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->NER32DupAllIndex[(k+0)*6+5];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->NER32DupAllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->NER32DupAllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->NER32DupAllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->NER32DupAllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->NER32DupAllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->NER32DupAllIndex[(k+1)*6+5];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
    }
    else if(remainNER32Dup==3)
    {
        std::size_t k = endNER32Dup-3;
        const unsigned int i0 = this->NER32DupAllIndex[(k+0)*6+0];
        const unsigned int sl0a = this->NER32DupAllIndex[(k+0)*6+1];
        const unsigned int sl1a = this->NER32DupAllIndex[(k+0)*6+2];
        const unsigned int sl2a = this->NER32DupAllIndex[(k+0)*6+3];
        const unsigned int sr0a = this->NER32DupAllIndex[(k+0)*6+4];
        const unsigned int sr1a = this->NER32DupAllIndex[(k+0)*6+5];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

        const unsigned int i1 = this->NER32DupAllIndex[(k+1)*6+0];
        const unsigned int sl0b = this->NER32DupAllIndex[(k+1)*6+1];
        const unsigned int sl1b = this->NER32DupAllIndex[(k+1)*6+2];
        const unsigned int sl2b = this->NER32DupAllIndex[(k+1)*6+3];
        const unsigned int sr0b = this->NER32DupAllIndex[(k+1)*6+4];
        const unsigned int sr1b = this->NER32DupAllIndex[(k+1)*6+5];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

        const unsigned int i2 = this->NER32DupAllIndex[(k+2)*6+0];
        const unsigned int sl0c = this->NER32DupAllIndex[(k+2)*6+1];
        const unsigned int sl1c = this->NER32DupAllIndex[(k+2)*6+2];
        const unsigned int sl2c = this->NER32DupAllIndex[(k+2)*6+3];
        const unsigned int sr0c = this->NER32DupAllIndex[(k+2)*6+4];
        const unsigned int sr1c = this->NER32DupAllIndex[(k+2)*6+5];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
    }
    }
}
