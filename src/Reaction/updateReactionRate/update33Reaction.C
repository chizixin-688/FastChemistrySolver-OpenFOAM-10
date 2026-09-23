/*---------------------------------------------------------------------------*\
  Description
      Update the net rates of production (dNdtByV) for three-reactant /
      three-product (3-3) reactions, e.g. A + A + A -> B + B + B and A + B + C -> D + E + F.

      RR : reversible reaction; the reverse rate constant is obtained from
           the equilibrium constant
      IR : irreversible reaction; the reverse rate constant is zero
      NER: non-equilibrium reaction; the reverse rate constant is computed
           with the Arrhenius form instead of the equilibrium constant

  Functions
      RF33RR : reversible 3-3 reactions, incl. duplicate species
      RF33IR : irreversible 3-3 reactions, incl. duplicate species
      RF33NER: non-equilibrium 3-3 reactions

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
FastChemistry::OptReaction::RF33RR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    const double PowPByRTSumVki = this->Pow_pByRT_SumVki[2];
    const double invKcLimiter = FastChemistry::invKcLimiter;
    std::size_t endRR33 = this->RR33Size;
    std::size_t remainRR33 = endRR33%4;
    for(std::size_t k=0; k<endRR33-remainRR33; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->RR33AllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->RR33AllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->RR33AllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->RR33AllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->RR33AllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->RR33AllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->RR33AllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]*invNegGstdByRT[sr2a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->RR33AllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->RR33AllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->RR33AllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->RR33AllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->RR33AllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->RR33AllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->RR33AllIndex[(k+1)*7+6];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]*invNegGstdByRT[sr2b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->RR33AllIndex[(k+2)*7+0];
        const unsigned int sl0c = this->RR33AllIndex[(k+2)*7+1];
        const unsigned int sl1c = this->RR33AllIndex[(k+2)*7+2];
        const unsigned int sl2c = this->RR33AllIndex[(k+2)*7+3];
        const unsigned int sr0c = this->RR33AllIndex[(k+2)*7+4];
        const unsigned int sr1c = this->RR33AllIndex[(k+2)*7+5];
        const unsigned int sr2c = this->RR33AllIndex[(k+2)*7+6];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*
            (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c]*invNegGstdByRT[sr2c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->RR33AllIndex[(k+3)*7+0];
        const unsigned int sl0d = this->RR33AllIndex[(k+3)*7+1];
        const unsigned int sl1d = this->RR33AllIndex[(k+3)*7+2];
        const unsigned int sl2d = this->RR33AllIndex[(k+3)*7+3];
        const unsigned int sr0d = this->RR33AllIndex[(k+3)*7+4];
        const unsigned int sr1d = this->RR33AllIndex[(k+3)*7+5];
        const unsigned int sr2d = this->RR33AllIndex[(k+3)*7+6];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        double invKc3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]*ExpNegGbyRT[sl2d])*
            (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d]*invNegGstdByRT[sr2d])*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d]*c[sr1d]*c[sr2d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sl2d_increment = dNdtByV[sl2d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        const double omega_sr2d_increment = dNdtByV[sr2d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sl2d] = omega_sl2d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
        dNdtByV[sr2d] = omega_sr2d_increment;
    }
    if(remainRR33==1)
    {
        std::size_t k = endRR33-1;
        const unsigned int i0 = this->RR33AllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->RR33AllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->RR33AllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->RR33AllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->RR33AllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->RR33AllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->RR33AllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]*invNegGstdByRT[sr2a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;
    }
    else if(remainRR33==2)
    {
        std::size_t k = endRR33-2;
        const unsigned int i0 = this->RR33AllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->RR33AllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->RR33AllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->RR33AllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->RR33AllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->RR33AllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->RR33AllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]*invNegGstdByRT[sr2a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->RR33AllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->RR33AllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->RR33AllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->RR33AllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->RR33AllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->RR33AllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->RR33AllIndex[(k+1)*7+6];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]*invNegGstdByRT[sr2b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;
    }
    else if(remainRR33==3)
    {
        std::size_t k = endRR33-3;
        const unsigned int i0 = this->RR33AllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->RR33AllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->RR33AllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->RR33AllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->RR33AllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->RR33AllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->RR33AllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]*invNegGstdByRT[sr2a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->RR33AllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->RR33AllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->RR33AllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->RR33AllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->RR33AllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->RR33AllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->RR33AllIndex[(k+1)*7+6];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]*invNegGstdByRT[sr2b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        const unsigned int i2 = this->RR33AllIndex[(k+2)*7+0];
        const unsigned int sl0c = this->RR33AllIndex[(k+2)*7+1];
        const unsigned int sl1c = this->RR33AllIndex[(k+2)*7+2];
        const unsigned int sl2c = this->RR33AllIndex[(k+2)*7+3];
        const unsigned int sr0c = this->RR33AllIndex[(k+2)*7+4];
        const unsigned int sr1c = this->RR33AllIndex[(k+2)*7+5];
        const unsigned int sr2c = this->RR33AllIndex[(k+2)*7+6];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*
            (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c]*invNegGstdByRT[sr2c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;
    }
    {
    std::size_t endRR33Dup = this->RR33DupSize;
    std::size_t remainRR33Dup = endRR33Dup%4;
    for(std::size_t k=0; k<endRR33Dup-remainRR33Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->RR33DupAllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->RR33DupAllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->RR33DupAllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->RR33DupAllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->RR33DupAllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->RR33DupAllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->RR33DupAllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]*invNegGstdByRT[sr2a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->RR33DupAllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->RR33DupAllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->RR33DupAllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->RR33DupAllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->RR33DupAllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->RR33DupAllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->RR33DupAllIndex[(k+1)*7+6];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]*invNegGstdByRT[sr2b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->RR33DupAllIndex[(k+2)*7+0];
        const unsigned int sl0c = this->RR33DupAllIndex[(k+2)*7+1];
        const unsigned int sl1c = this->RR33DupAllIndex[(k+2)*7+2];
        const unsigned int sl2c = this->RR33DupAllIndex[(k+2)*7+3];
        const unsigned int sr0c = this->RR33DupAllIndex[(k+2)*7+4];
        const unsigned int sr1c = this->RR33DupAllIndex[(k+2)*7+5];
        const unsigned int sr2c = this->RR33DupAllIndex[(k+2)*7+6];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*
            (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c]*invNegGstdByRT[sr2c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->RR33DupAllIndex[(k+3)*7+0];
        const unsigned int sl0d = this->RR33DupAllIndex[(k+3)*7+1];
        const unsigned int sl1d = this->RR33DupAllIndex[(k+3)*7+2];
        const unsigned int sl2d = this->RR33DupAllIndex[(k+3)*7+3];
        const unsigned int sr0d = this->RR33DupAllIndex[(k+3)*7+4];
        const unsigned int sr1d = this->RR33DupAllIndex[(k+3)*7+5];
        const unsigned int sr2d = this->RR33DupAllIndex[(k+3)*7+6];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        double invKc3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]*ExpNegGbyRT[sl2d])*
            (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d]*invNegGstdByRT[sr2d])*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d]*c[sr1d]*c[sr2d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
        dNdtByV[sr2d] = dNdtByV[sr2d] + q3;
    }
    if(remainRR33Dup==1)
    {
        std::size_t k = endRR33Dup-1;
        const unsigned int i0 = this->RR33DupAllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->RR33DupAllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->RR33DupAllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->RR33DupAllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->RR33DupAllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->RR33DupAllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->RR33DupAllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]*invNegGstdByRT[sr2a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;
    }
    else if(remainRR33Dup==2)
    {
        std::size_t k = endRR33Dup-2;
        const unsigned int i0 = this->RR33DupAllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->RR33DupAllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->RR33DupAllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->RR33DupAllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->RR33DupAllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->RR33DupAllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->RR33DupAllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]*invNegGstdByRT[sr2a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        const unsigned int i1 = this->RR33DupAllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->RR33DupAllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->RR33DupAllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->RR33DupAllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->RR33DupAllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->RR33DupAllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->RR33DupAllIndex[(k+1)*7+6];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]*invNegGstdByRT[sr2b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;
    }
    else if(remainRR33Dup==3)
    {
        std::size_t k = endRR33Dup-3;
        const unsigned int i0 = this->RR33DupAllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->RR33DupAllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->RR33DupAllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->RR33DupAllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->RR33DupAllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->RR33DupAllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->RR33DupAllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]*invNegGstdByRT[sr2a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        const unsigned int i1 = this->RR33DupAllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->RR33DupAllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->RR33DupAllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->RR33DupAllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->RR33DupAllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->RR33DupAllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->RR33DupAllIndex[(k+1)*7+6];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]*invNegGstdByRT[sr2b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;

        const unsigned int i2 = this->RR33DupAllIndex[(k+2)*7+0];
        const unsigned int sl0c = this->RR33DupAllIndex[(k+2)*7+1];
        const unsigned int sl1c = this->RR33DupAllIndex[(k+2)*7+2];
        const unsigned int sl2c = this->RR33DupAllIndex[(k+2)*7+3];
        const unsigned int sr0c = this->RR33DupAllIndex[(k+2)*7+4];
        const unsigned int sr1c = this->RR33DupAllIndex[(k+2)*7+5];
        const unsigned int sr2c = this->RR33DupAllIndex[(k+2)*7+6];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*
            (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c]*invNegGstdByRT[sr2c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;
    }
    }
}


void 
FastChemistry::OptReaction::RF33IR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t endIR33 = this->IR33Size;
    std::size_t remainIR33 = endIR33%4;
    for(std::size_t k=0; k<endIR33-remainIR33; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->IR33AllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->IR33AllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->IR33AllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->IR33AllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->IR33AllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->IR33AllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->IR33AllIndex[(k+0)*7+6];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->IR33AllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->IR33AllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->IR33AllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->IR33AllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->IR33AllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->IR33AllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->IR33AllIndex[(k+1)*7+6];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->IR33AllIndex[(k+2)*7+0];
        const unsigned int sl0c = this->IR33AllIndex[(k+2)*7+1];
        const unsigned int sl1c = this->IR33AllIndex[(k+2)*7+2];
        const unsigned int sl2c = this->IR33AllIndex[(k+2)*7+3];
        const unsigned int sr0c = this->IR33AllIndex[(k+2)*7+4];
        const unsigned int sr1c = this->IR33AllIndex[(k+2)*7+5];
        const unsigned int sr2c = this->IR33AllIndex[(k+2)*7+6];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->IR33AllIndex[(k+3)*7+0];
        const unsigned int sl0d = this->IR33AllIndex[(k+3)*7+1];
        const unsigned int sl1d = this->IR33AllIndex[(k+3)*7+2];
        const unsigned int sl2d = this->IR33AllIndex[(k+3)*7+3];
        const unsigned int sr0d = this->IR33AllIndex[(k+3)*7+4];
        const unsigned int sr1d = this->IR33AllIndex[(k+3)*7+5];
        const unsigned int sr2d = this->IR33AllIndex[(k+3)*7+6];
        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sl2d_increment = dNdtByV[sl2d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        const double omega_sr2d_increment = dNdtByV[sr2d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sl2d] = omega_sl2d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
        dNdtByV[sr2d] = omega_sr2d_increment;
    }
    if(remainIR33==1)
    {
        std::size_t k = endIR33-1;
        const unsigned int i0 = this->IR33AllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->IR33AllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->IR33AllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->IR33AllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->IR33AllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->IR33AllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->IR33AllIndex[(k+0)*7+6];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;
    }
    else if(remainIR33==2)
    {
        std::size_t k = endIR33-2;
        const unsigned int i0 = this->IR33AllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->IR33AllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->IR33AllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->IR33AllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->IR33AllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->IR33AllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->IR33AllIndex[(k+0)*7+6];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->IR33AllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->IR33AllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->IR33AllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->IR33AllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->IR33AllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->IR33AllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->IR33AllIndex[(k+1)*7+6];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;
    }
    else if(remainIR33==3)
    {
        std::size_t k = endIR33-3;
        const unsigned int i0 = this->IR33AllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->IR33AllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->IR33AllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->IR33AllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->IR33AllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->IR33AllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->IR33AllIndex[(k+0)*7+6];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->IR33AllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->IR33AllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->IR33AllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->IR33AllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->IR33AllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->IR33AllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->IR33AllIndex[(k+1)*7+6];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        const unsigned int i2 = this->IR33AllIndex[(k+2)*7+0];
        const unsigned int sl0c = this->IR33AllIndex[(k+2)*7+1];
        const unsigned int sl1c = this->IR33AllIndex[(k+2)*7+2];
        const unsigned int sl2c = this->IR33AllIndex[(k+2)*7+3];
        const unsigned int sr0c = this->IR33AllIndex[(k+2)*7+4];
        const unsigned int sr1c = this->IR33AllIndex[(k+2)*7+5];
        const unsigned int sr2c = this->IR33AllIndex[(k+2)*7+6];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;
    }
    {
    std::size_t endIR33Dup = this->IR33DupSize;
    std::size_t remainIR33Dup = endIR33Dup%4;
    for(std::size_t k=0; k<endIR33Dup-remainIR33Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->IR33DupAllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->IR33DupAllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->IR33DupAllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->IR33DupAllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->IR33DupAllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->IR33DupAllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->IR33DupAllIndex[(k+0)*7+6];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->IR33DupAllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->IR33DupAllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->IR33DupAllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->IR33DupAllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->IR33DupAllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->IR33DupAllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->IR33DupAllIndex[(k+1)*7+6];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->IR33DupAllIndex[(k+2)*7+0];
        const unsigned int sl0c = this->IR33DupAllIndex[(k+2)*7+1];
        const unsigned int sl1c = this->IR33DupAllIndex[(k+2)*7+2];
        const unsigned int sl2c = this->IR33DupAllIndex[(k+2)*7+3];
        const unsigned int sr0c = this->IR33DupAllIndex[(k+2)*7+4];
        const unsigned int sr1c = this->IR33DupAllIndex[(k+2)*7+5];
        const unsigned int sr2c = this->IR33DupAllIndex[(k+2)*7+6];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->IR33DupAllIndex[(k+3)*7+0];
        const unsigned int sl0d = this->IR33DupAllIndex[(k+3)*7+1];
        const unsigned int sl1d = this->IR33DupAllIndex[(k+3)*7+2];
        const unsigned int sl2d = this->IR33DupAllIndex[(k+3)*7+3];
        const unsigned int sr0d = this->IR33DupAllIndex[(k+3)*7+4];
        const unsigned int sr1d = this->IR33DupAllIndex[(k+3)*7+5];
        const unsigned int sr2d = this->IR33DupAllIndex[(k+3)*7+6];
        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double q3 = (Kf3*CF3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
        dNdtByV[sr2d] = dNdtByV[sr2d] + q3;
    }
    if(remainIR33Dup==1)
    {
        std::size_t k = endIR33Dup-1;
        const unsigned int i0 = this->IR33DupAllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->IR33DupAllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->IR33DupAllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->IR33DupAllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->IR33DupAllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->IR33DupAllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->IR33DupAllIndex[(k+0)*7+6];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;
    }
    else if(remainIR33Dup==2)
    {
        std::size_t k = endIR33Dup-2;
        const unsigned int i0 = this->IR33DupAllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->IR33DupAllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->IR33DupAllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->IR33DupAllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->IR33DupAllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->IR33DupAllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->IR33DupAllIndex[(k+0)*7+6];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        const unsigned int i1 = this->IR33DupAllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->IR33DupAllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->IR33DupAllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->IR33DupAllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->IR33DupAllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->IR33DupAllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->IR33DupAllIndex[(k+1)*7+6];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;
    }
    else if(remainIR33Dup==3)
    {
        std::size_t k = endIR33Dup-3;
        const unsigned int i0 = this->IR33DupAllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->IR33DupAllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->IR33DupAllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->IR33DupAllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->IR33DupAllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->IR33DupAllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->IR33DupAllIndex[(k+0)*7+6];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        const unsigned int i1 = this->IR33DupAllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->IR33DupAllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->IR33DupAllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->IR33DupAllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->IR33DupAllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->IR33DupAllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->IR33DupAllIndex[(k+1)*7+6];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;

        const unsigned int i2 = this->IR33DupAllIndex[(k+2)*7+0];
        const unsigned int sl0c = this->IR33DupAllIndex[(k+2)*7+1];
        const unsigned int sl1c = this->IR33DupAllIndex[(k+2)*7+2];
        const unsigned int sl2c = this->IR33DupAllIndex[(k+2)*7+3];
        const unsigned int sr0c = this->IR33DupAllIndex[(k+2)*7+4];
        const unsigned int sr1c = this->IR33DupAllIndex[(k+2)*7+5];
        const unsigned int sr2c = this->IR33DupAllIndex[(k+2)*7+6];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;
    }
    }
}


void 
FastChemistry::OptReaction::RF33NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t endNER33 = this->NER33Size;
    std::size_t remainNER33 = endNER33%4;
    for(std::size_t k=0; k<endNER33-remainNER33; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->NER33AllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->NER33AllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->NER33AllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->NER33AllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->NER33AllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->NER33AllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->NER33AllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->NER33AllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->NER33AllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->NER33AllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->NER33AllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->NER33AllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->NER33AllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->NER33AllIndex[(k+1)*7+6];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->NER33AllIndex[(k+2)*7+0];
        const unsigned int sl0c = this->NER33AllIndex[(k+2)*7+1];
        const unsigned int sl1c = this->NER33AllIndex[(k+2)*7+2];
        const unsigned int sl2c = this->NER33AllIndex[(k+2)*7+3];
        const unsigned int sr0c = this->NER33AllIndex[(k+2)*7+4];
        const unsigned int sr1c = this->NER33AllIndex[(k+2)*7+5];
        const unsigned int sr2c = this->NER33AllIndex[(k+2)*7+6];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->NER33AllIndex[(k+3)*7+0];
        const unsigned int sl0d = this->NER33AllIndex[(k+3)*7+1];
        const unsigned int sl1d = this->NER33AllIndex[(k+3)*7+2];
        const unsigned int sl2d = this->NER33AllIndex[(k+3)*7+3];
        const unsigned int sr0d = this->NER33AllIndex[(k+3)*7+4];
        const unsigned int sr1d = this->NER33AllIndex[(k+3)*7+5];
        const unsigned int sr2d = this->NER33AllIndex[(k+3)*7+6];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d]*c[sr1d]*c[sr2d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sl2d_increment = dNdtByV[sl2d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        const double omega_sr2d_increment = dNdtByV[sr2d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sl2d] = omega_sl2d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
        dNdtByV[sr2d] = omega_sr2d_increment;
    }
    if(remainNER33==1)
    {
        std::size_t k = endNER33-1;
        const unsigned int i0 = this->NER33AllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->NER33AllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->NER33AllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->NER33AllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->NER33AllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->NER33AllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->NER33AllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;
    }
    else if(remainNER33==2)
    {
        std::size_t k = endNER33-2;
        const unsigned int i0 = this->NER33AllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->NER33AllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->NER33AllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->NER33AllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->NER33AllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->NER33AllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->NER33AllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->NER33AllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->NER33AllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->NER33AllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->NER33AllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->NER33AllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->NER33AllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->NER33AllIndex[(k+1)*7+6];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;
    }
    else if(remainNER33==3)
    {
        std::size_t k = endNER33-3;
        const unsigned int i0 = this->NER33AllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->NER33AllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->NER33AllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->NER33AllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->NER33AllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->NER33AllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->NER33AllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->NER33AllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->NER33AllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->NER33AllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->NER33AllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->NER33AllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->NER33AllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->NER33AllIndex[(k+1)*7+6];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        const unsigned int i2 = this->NER33AllIndex[(k+2)*7+0];
        const unsigned int sl0c = this->NER33AllIndex[(k+2)*7+1];
        const unsigned int sl1c = this->NER33AllIndex[(k+2)*7+2];
        const unsigned int sl2c = this->NER33AllIndex[(k+2)*7+3];
        const unsigned int sr0c = this->NER33AllIndex[(k+2)*7+4];
        const unsigned int sr1c = this->NER33AllIndex[(k+2)*7+5];
        const unsigned int sr2c = this->NER33AllIndex[(k+2)*7+6];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;
    }
    {
    std::size_t endNER33Dup = this->NER33DupSize;
    std::size_t remainNER33Dup = endNER33Dup%4;
    for(std::size_t k=0; k<endNER33Dup-remainNER33Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->NER33DupAllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->NER33DupAllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->NER33DupAllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->NER33DupAllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->NER33DupAllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->NER33DupAllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->NER33DupAllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->NER33DupAllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->NER33DupAllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->NER33DupAllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->NER33DupAllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->NER33DupAllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->NER33DupAllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->NER33DupAllIndex[(k+1)*7+6];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->NER33DupAllIndex[(k+2)*7+0];
        const unsigned int sl0c = this->NER33DupAllIndex[(k+2)*7+1];
        const unsigned int sl1c = this->NER33DupAllIndex[(k+2)*7+2];
        const unsigned int sl2c = this->NER33DupAllIndex[(k+2)*7+3];
        const unsigned int sr0c = this->NER33DupAllIndex[(k+2)*7+4];
        const unsigned int sr1c = this->NER33DupAllIndex[(k+2)*7+5];
        const unsigned int sr2c = this->NER33DupAllIndex[(k+2)*7+6];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->NER33DupAllIndex[(k+3)*7+0];
        const unsigned int sl0d = this->NER33DupAllIndex[(k+3)*7+1];
        const unsigned int sl1d = this->NER33DupAllIndex[(k+3)*7+2];
        const unsigned int sl2d = this->NER33DupAllIndex[(k+3)*7+3];
        const unsigned int sr0d = this->NER33DupAllIndex[(k+3)*7+4];
        const unsigned int sr1d = this->NER33DupAllIndex[(k+3)*7+5];
        const unsigned int sr2d = this->NER33DupAllIndex[(k+3)*7+6];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d]*c[sr1d]*c[sr2d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
        dNdtByV[sr2d] = dNdtByV[sr2d] + q3;
    }
    if(remainNER33Dup==1)
    {
        std::size_t k = endNER33Dup-1;
        const unsigned int i0 = this->NER33DupAllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->NER33DupAllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->NER33DupAllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->NER33DupAllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->NER33DupAllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->NER33DupAllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->NER33DupAllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;
    }
    else if(remainNER33Dup==2)
    {
        std::size_t k = endNER33Dup-2;
        const unsigned int i0 = this->NER33DupAllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->NER33DupAllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->NER33DupAllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->NER33DupAllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->NER33DupAllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->NER33DupAllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->NER33DupAllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        const unsigned int i1 = this->NER33DupAllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->NER33DupAllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->NER33DupAllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->NER33DupAllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->NER33DupAllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->NER33DupAllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->NER33DupAllIndex[(k+1)*7+6];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;
    }
    else if(remainNER33Dup==3)
    {
        std::size_t k = endNER33Dup-3;
        const unsigned int i0 = this->NER33DupAllIndex[(k+0)*7+0];
        const unsigned int sl0a = this->NER33DupAllIndex[(k+0)*7+1];
        const unsigned int sl1a = this->NER33DupAllIndex[(k+0)*7+2];
        const unsigned int sl2a = this->NER33DupAllIndex[(k+0)*7+3];
        const unsigned int sr0a = this->NER33DupAllIndex[(k+0)*7+4];
        const unsigned int sr1a = this->NER33DupAllIndex[(k+0)*7+5];
        const unsigned int sr2a = this->NER33DupAllIndex[(k+0)*7+6];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        const unsigned int i1 = this->NER33DupAllIndex[(k+1)*7+0];
        const unsigned int sl0b = this->NER33DupAllIndex[(k+1)*7+1];
        const unsigned int sl1b = this->NER33DupAllIndex[(k+1)*7+2];
        const unsigned int sl2b = this->NER33DupAllIndex[(k+1)*7+3];
        const unsigned int sr0b = this->NER33DupAllIndex[(k+1)*7+4];
        const unsigned int sr1b = this->NER33DupAllIndex[(k+1)*7+5];
        const unsigned int sr2b = this->NER33DupAllIndex[(k+1)*7+6];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;

        const unsigned int i2 = this->NER33DupAllIndex[(k+2)*7+0];
        const unsigned int sl0c = this->NER33DupAllIndex[(k+2)*7+1];
        const unsigned int sl1c = this->NER33DupAllIndex[(k+2)*7+2];
        const unsigned int sl2c = this->NER33DupAllIndex[(k+2)*7+3];
        const unsigned int sr0c = this->NER33DupAllIndex[(k+2)*7+4];
        const unsigned int sr1c = this->NER33DupAllIndex[(k+2)*7+5];
        const unsigned int sr2c = this->NER33DupAllIndex[(k+2)*7+6];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;
    }
    }
}
