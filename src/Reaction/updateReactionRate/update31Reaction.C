/*---------------------------------------------------------------------------*\
  Description
      Update the net rates of production (dNdtByV) for three-reactant /
      one-product (3-1) reactions, e.g. A + A + A -> B and A + B + C -> D.

      RR : reversible reaction; the reverse rate constant is obtained from
           the equilibrium constant
      IR : irreversible reaction; the reverse rate constant is zero
      NER: non-equilibrium reaction; the reverse rate constant is computed
           with the Arrhenius form instead of the equilibrium constant

  Functions
      RF31RR : reversible 3-1 reactions, incl. duplicate reactants (A+A+A=B)
      RF31IR : irreversible 3-1 reactions, incl. duplicate reactants (A+A+A=B)
      RF31NER: non-equilibrium 3-1 reactions

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
FastChemistry::OptReaction::RF31RR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    const double PowPByRTSumVki = this->Pow_pByRT_SumVki[4];
    const double invKcLimiter = FastChemistry::invKcLimiter;
    std::size_t endRR31 = this->RR31Size;
    std::size_t remainRR31 = endRR31%4;
    for(std::size_t k=0; k<endRR31-remainRR31; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->RR31AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR31AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR31AllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->RR31AllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->RR31AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->RR31AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR31AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->RR31AllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->RR31AllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->RR31AllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->RR31AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->RR31AllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->RR31AllIndex[(k+2)*5+2];
        const unsigned int sl2c = this->RR31AllIndex[(k+2)*5+3];
        const unsigned int sr0c = this->RR31AllIndex[(k+2)*5+4];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*
            (invNegGstdByRT[sr0c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->RR31AllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->RR31AllIndex[(k+3)*5+1];
        const unsigned int sl1d = this->RR31AllIndex[(k+3)*5+2];
        const unsigned int sl2d = this->RR31AllIndex[(k+3)*5+3];
        const unsigned int sr0d = this->RR31AllIndex[(k+3)*5+4];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        double invKc3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]*ExpNegGbyRT[sl2d])*
            (invNegGstdByRT[sr0d])*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sl2d_increment = dNdtByV[sl2d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sl2d] = omega_sl2d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
    }
    if(remainRR31==1)
    {
        std::size_t k = endRR31-1;
        const unsigned int i0 = this->RR31AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR31AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR31AllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->RR31AllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->RR31AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }
    else if(remainRR31==2)
    {
        std::size_t k = endRR31-2;
        const unsigned int i0 = this->RR31AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR31AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR31AllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->RR31AllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->RR31AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->RR31AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR31AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->RR31AllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->RR31AllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->RR31AllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
    }
    else if(remainRR31==3)
    {
        std::size_t k = endRR31-3;
        const unsigned int i0 = this->RR31AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR31AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR31AllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->RR31AllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->RR31AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->RR31AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR31AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->RR31AllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->RR31AllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->RR31AllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        const unsigned int i2 = this->RR31AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->RR31AllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->RR31AllIndex[(k+2)*5+2];
        const unsigned int sl2c = this->RR31AllIndex[(k+2)*5+3];
        const unsigned int sr0c = this->RR31AllIndex[(k+2)*5+4];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*
            (invNegGstdByRT[sr0c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }
    {
    std::size_t endRR31Dup = this->RR31DupSize;
    std::size_t remainRR31Dup = endRR31Dup%4;
    for(std::size_t k=0; k<endRR31Dup-remainRR31Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->RR31DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR31DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR31DupAllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->RR31DupAllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->RR31DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->RR31DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR31DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->RR31DupAllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->RR31DupAllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->RR31DupAllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->RR31DupAllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->RR31DupAllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->RR31DupAllIndex[(k+2)*5+2];
        const unsigned int sl2c = this->RR31DupAllIndex[(k+2)*5+3];
        const unsigned int sr0c = this->RR31DupAllIndex[(k+2)*5+4];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*
            (invNegGstdByRT[sr0c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->RR31DupAllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->RR31DupAllIndex[(k+3)*5+1];
        const unsigned int sl1d = this->RR31DupAllIndex[(k+3)*5+2];
        const unsigned int sl2d = this->RR31DupAllIndex[(k+3)*5+3];
        const unsigned int sr0d = this->RR31DupAllIndex[(k+3)*5+4];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        double invKc3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]*ExpNegGbyRT[sl2d])*
            (invNegGstdByRT[sr0d])*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
    }
    if(remainRR31Dup==1)
    {
        std::size_t k = endRR31Dup-1;
        const unsigned int i0 = this->RR31DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR31DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR31DupAllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->RR31DupAllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->RR31DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
    }
    else if(remainRR31Dup==2)
    {
        std::size_t k = endRR31Dup-2;
        const unsigned int i0 = this->RR31DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR31DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR31DupAllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->RR31DupAllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->RR31DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        const unsigned int i1 = this->RR31DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR31DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->RR31DupAllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->RR31DupAllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->RR31DupAllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
    }
    else if(remainRR31Dup==3)
    {
        std::size_t k = endRR31Dup-3;
        const unsigned int i0 = this->RR31DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR31DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->RR31DupAllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->RR31DupAllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->RR31DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a])*
            (invNegGstdByRT[sr0a])*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        const unsigned int i1 = this->RR31DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR31DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->RR31DupAllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->RR31DupAllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->RR31DupAllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b])*
            (invNegGstdByRT[sr0b])*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;

        const unsigned int i2 = this->RR31DupAllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->RR31DupAllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->RR31DupAllIndex[(k+2)*5+2];
        const unsigned int sl2c = this->RR31DupAllIndex[(k+2)*5+3];
        const unsigned int sr0c = this->RR31DupAllIndex[(k+2)*5+4];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c])*
            (invNegGstdByRT[sr0c])*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
    }
    }
}


void 
FastChemistry::OptReaction::RF31IR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t endIR31 = this->IR31Size;
    std::size_t remainIR31 = endIR31%4;
    for(std::size_t k=0; k<endIR31-remainIR31; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->IR31AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR31AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR31AllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->IR31AllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->IR31AllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->IR31AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR31AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->IR31AllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->IR31AllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->IR31AllIndex[(k+1)*5+4];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->IR31AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->IR31AllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->IR31AllIndex[(k+2)*5+2];
        const unsigned int sl2c = this->IR31AllIndex[(k+2)*5+3];
        const unsigned int sr0c = this->IR31AllIndex[(k+2)*5+4];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->IR31AllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->IR31AllIndex[(k+3)*5+1];
        const unsigned int sl1d = this->IR31AllIndex[(k+3)*5+2];
        const unsigned int sl2d = this->IR31AllIndex[(k+3)*5+3];
        const unsigned int sr0d = this->IR31AllIndex[(k+3)*5+4];
        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sl2d_increment = dNdtByV[sl2d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sl2d] = omega_sl2d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
    }
    if(remainIR31==1)
    {
        std::size_t k = endIR31-1;
        const unsigned int i0 = this->IR31AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR31AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR31AllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->IR31AllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->IR31AllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }
    else if(remainIR31==2)
    {
        std::size_t k = endIR31-2;
        const unsigned int i0 = this->IR31AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR31AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR31AllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->IR31AllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->IR31AllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->IR31AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR31AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->IR31AllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->IR31AllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->IR31AllIndex[(k+1)*5+4];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
    }
    else if(remainIR31==3)
    {
        std::size_t k = endIR31-3;
        const unsigned int i0 = this->IR31AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR31AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR31AllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->IR31AllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->IR31AllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->IR31AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR31AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->IR31AllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->IR31AllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->IR31AllIndex[(k+1)*5+4];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        const unsigned int i2 = this->IR31AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->IR31AllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->IR31AllIndex[(k+2)*5+2];
        const unsigned int sl2c = this->IR31AllIndex[(k+2)*5+3];
        const unsigned int sr0c = this->IR31AllIndex[(k+2)*5+4];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }
    {
    std::size_t endIR31Dup = this->IR31DupSize;
    std::size_t remainIR31Dup = endIR31Dup%4;
    for(std::size_t k=0; k<endIR31Dup-remainIR31Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->IR31DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR31DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR31DupAllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->IR31DupAllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->IR31DupAllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->IR31DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR31DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->IR31DupAllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->IR31DupAllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->IR31DupAllIndex[(k+1)*5+4];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->IR31DupAllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->IR31DupAllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->IR31DupAllIndex[(k+2)*5+2];
        const unsigned int sl2c = this->IR31DupAllIndex[(k+2)*5+3];
        const unsigned int sr0c = this->IR31DupAllIndex[(k+2)*5+4];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->IR31DupAllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->IR31DupAllIndex[(k+3)*5+1];
        const unsigned int sl1d = this->IR31DupAllIndex[(k+3)*5+2];
        const unsigned int sl2d = this->IR31DupAllIndex[(k+3)*5+3];
        const unsigned int sr0d = this->IR31DupAllIndex[(k+3)*5+4];
        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double q3 = (Kf3*CF3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
    }
    if(remainIR31Dup==1)
    {
        std::size_t k = endIR31Dup-1;
        const unsigned int i0 = this->IR31DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR31DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR31DupAllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->IR31DupAllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->IR31DupAllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
    }
    else if(remainIR31Dup==2)
    {
        std::size_t k = endIR31Dup-2;
        const unsigned int i0 = this->IR31DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR31DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR31DupAllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->IR31DupAllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->IR31DupAllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        const unsigned int i1 = this->IR31DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR31DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->IR31DupAllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->IR31DupAllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->IR31DupAllIndex[(k+1)*5+4];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
    }
    else if(remainIR31Dup==3)
    {
        std::size_t k = endIR31Dup-3;
        const unsigned int i0 = this->IR31DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR31DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->IR31DupAllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->IR31DupAllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->IR31DupAllIndex[(k+0)*5+4];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        const unsigned int i1 = this->IR31DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR31DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->IR31DupAllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->IR31DupAllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->IR31DupAllIndex[(k+1)*5+4];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;

        const unsigned int i2 = this->IR31DupAllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->IR31DupAllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->IR31DupAllIndex[(k+2)*5+2];
        const unsigned int sl2c = this->IR31DupAllIndex[(k+2)*5+3];
        const unsigned int sr0c = this->IR31DupAllIndex[(k+2)*5+4];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
    }
    }
}


void 
FastChemistry::OptReaction::RF31NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t endNER31 = this->NER31Size;
    std::size_t remainNER31 = endNER31%4;
    for(std::size_t k=0; k<endNER31-remainNER31; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->NER31AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER31AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER31AllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->NER31AllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->NER31AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->NER31AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER31AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER31AllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->NER31AllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->NER31AllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->NER31AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER31AllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->NER31AllIndex[(k+2)*5+2];
        const unsigned int sl2c = this->NER31AllIndex[(k+2)*5+3];
        const unsigned int sr0c = this->NER31AllIndex[(k+2)*5+4];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->NER31AllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->NER31AllIndex[(k+3)*5+1];
        const unsigned int sl1d = this->NER31AllIndex[(k+3)*5+2];
        const unsigned int sl2d = this->NER31AllIndex[(k+3)*5+3];
        const unsigned int sr0d = this->NER31AllIndex[(k+3)*5+4];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sl2d_increment = dNdtByV[sl2d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sl2d] = omega_sl2d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
    }
    if(remainNER31==1)
    {
        std::size_t k = endNER31-1;
        const unsigned int i0 = this->NER31AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER31AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER31AllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->NER31AllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->NER31AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }
    else if(remainNER31==2)
    {
        std::size_t k = endNER31-2;
        const unsigned int i0 = this->NER31AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER31AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER31AllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->NER31AllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->NER31AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->NER31AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER31AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER31AllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->NER31AllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->NER31AllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
    }
    else if(remainNER31==3)
    {
        std::size_t k = endNER31-3;
        const unsigned int i0 = this->NER31AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER31AllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER31AllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->NER31AllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->NER31AllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sl2a_increment = dNdtByV[sl2a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sl2a] = omega_sl2a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->NER31AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER31AllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER31AllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->NER31AllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->NER31AllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sl2b_increment = dNdtByV[sl2b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sl2b] = omega_sl2b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        const unsigned int i2 = this->NER31AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER31AllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->NER31AllIndex[(k+2)*5+2];
        const unsigned int sl2c = this->NER31AllIndex[(k+2)*5+3];
        const unsigned int sr0c = this->NER31AllIndex[(k+2)*5+4];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sl2c_increment = dNdtByV[sl2c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sl2c] = omega_sl2c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }
    {
    std::size_t endNER31Dup = this->NER31DupSize;
    std::size_t remainNER31Dup = endNER31Dup%4;
    for(std::size_t k=0; k<endNER31Dup-remainNER31Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========

        const unsigned int i0 = this->NER31DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER31DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER31DupAllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->NER31DupAllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->NER31DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        // ========== Reaction block k+1 (b) ==========

        const unsigned int i1 = this->NER31DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER31DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER31DupAllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->NER31DupAllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->NER31DupAllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;

        // ========== Reaction block k+2 (c) ==========

        const unsigned int i2 = this->NER31DupAllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER31DupAllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->NER31DupAllIndex[(k+2)*5+2];
        const unsigned int sl2c = this->NER31DupAllIndex[(k+2)*5+3];
        const unsigned int sr0c = this->NER31DupAllIndex[(k+2)*5+4];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;

        // ========== Reaction block k+3 (d) ==========

        const unsigned int i3 = this->NER31DupAllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->NER31DupAllIndex[(k+3)*5+1];
        const unsigned int sl1d = this->NER31DupAllIndex[(k+3)*5+2];
        const unsigned int sl2d = this->NER31DupAllIndex[(k+3)*5+3];
        const unsigned int sr0d = this->NER31DupAllIndex[(k+3)*5+4];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
    }
    if(remainNER31Dup==1)
    {
        std::size_t k = endNER31Dup-1;
        const unsigned int i0 = this->NER31DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER31DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER31DupAllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->NER31DupAllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->NER31DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
    }
    else if(remainNER31Dup==2)
    {
        std::size_t k = endNER31Dup-2;
        const unsigned int i0 = this->NER31DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER31DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER31DupAllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->NER31DupAllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->NER31DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        const unsigned int i1 = this->NER31DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER31DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER31DupAllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->NER31DupAllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->NER31DupAllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
    }
    else if(remainNER31Dup==3)
    {
        std::size_t k = endNER31Dup-3;
        const unsigned int i0 = this->NER31DupAllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER31DupAllIndex[(k+0)*5+1];
        const unsigned int sl1a = this->NER31DupAllIndex[(k+0)*5+2];
        const unsigned int sl2a = this->NER31DupAllIndex[(k+0)*5+3];
        const unsigned int sr0a = this->NER31DupAllIndex[(k+0)*5+4];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        const unsigned int i1 = this->NER31DupAllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER31DupAllIndex[(k+1)*5+1];
        const unsigned int sl1b = this->NER31DupAllIndex[(k+1)*5+2];
        const unsigned int sl2b = this->NER31DupAllIndex[(k+1)*5+3];
        const unsigned int sr0b = this->NER31DupAllIndex[(k+1)*5+4];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;

        const unsigned int i2 = this->NER31DupAllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER31DupAllIndex[(k+2)*5+1];
        const unsigned int sl1c = this->NER31DupAllIndex[(k+2)*5+2];
        const unsigned int sl2c = this->NER31DupAllIndex[(k+2)*5+3];
        const unsigned int sr0c = this->NER31DupAllIndex[(k+2)*5+4];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
    }
    }
}
