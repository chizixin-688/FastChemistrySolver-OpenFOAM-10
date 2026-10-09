/*---------------------------------------------------------------------------*\
  Description
      Update the net rates of production (dNdtByV) for two-reactant /
      one-product (2-1) reactions, e.g. A + A -> B and A + B -> C.

      RR : reversible reaction; the reverse rate constant is obtained from
           the equilibrium constant
      IR : irreversible reaction; the reverse rate constant is zero
      NER: non-equilibrium reaction; the reverse rate constant is computed
           with the Arrhenius form instead of the equilibrium constant

  Functions
      RF21RR : reversible 2-1 reactions, incl. duplicate reactants (A+A=B)
      RF21IR : irreversible 2-1 reactions, incl. duplicate reactants (A+A=B)
      RF21NER: non-equilibrium 2-1 reactions

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
FastChemistry::OptReaction::RF21RR
(
    const double* __restrict__ c,
    double* __restrict__ dNdtByV,
    const double* __restrict__ ExpNegGbyRT
) const noexcept
{
    const double PowPByRT_SumVki = this->Pow_pByRT_SumVki[3];
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
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        const double invKp0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]);
        double invKc0 = invKp0*PowPByRT_SumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR21AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->RR21AllIndex[(k+1)*4+1];
        const unsigned int sl1b = this->RR21AllIndex[(k+1)*4+2];
        const unsigned int sr0b = this->RR21AllIndex[(k+1)*4+3];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        const double invKp1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b])*(invNegGstdByRT[sr0b]);
        double invKc1 = invKp1*PowPByRT_SumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR21AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->RR21AllIndex[(k+2)*4+1];
        const unsigned int sl1c = this->RR21AllIndex[(k+2)*4+2];
        const unsigned int sr0c = this->RR21AllIndex[(k+2)*4+3];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        const double invKp2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c])*(invNegGstdByRT[sr0c]);
        double invKc2 = invKp2*PowPByRT_SumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR21AllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->RR21AllIndex[(k+3)*4+1];
        const unsigned int sl1d = this->RR21AllIndex[(k+3)*4+2];
        const unsigned int sr0d = this->RR21AllIndex[(k+3)*4+3];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        const double invKp3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d])*(invNegGstdByRT[sr0d]);
        double invKc3 = invKp3*PowPByRT_SumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d]*c[sl1d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
    }
    if(remainRR21==1)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->RR21AllIndex[(endRR21-1)*4+0];
        const unsigned int sl0a = this->RR21AllIndex[(endRR21-1)*4+1];
        const unsigned int sl1a = this->RR21AllIndex[(endRR21-1)*4+2];
        const unsigned int sr0a = this->RR21AllIndex[(endRR21-1)*4+3];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        const double invKp0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]);
        double invKc0 = invKp0*PowPByRT_SumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }

    else if(remainRR21==2)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->RR21AllIndex[(endRR21-2)*4+0];
        const unsigned int sl0a = this->RR21AllIndex[(endRR21-2)*4+1];
        const unsigned int sl1a = this->RR21AllIndex[(endRR21-2)*4+2];
        const unsigned int sr0a = this->RR21AllIndex[(endRR21-2)*4+3];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        const double invKp0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]);
        double invKc0 = invKp0*PowPByRT_SumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Remainder block 2 (b) ==========
        const unsigned int i1 = this->RR21AllIndex[(endRR21-1)*4+0];
        const unsigned int sl0b = this->RR21AllIndex[(endRR21-1)*4+1];
        const unsigned int sl1b = this->RR21AllIndex[(endRR21-1)*4+2];
        const unsigned int sr0b = this->RR21AllIndex[(endRR21-1)*4+3];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        const double invKp1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b])*(invNegGstdByRT[sr0b]);
        double invKc1 = invKp1*PowPByRT_SumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
    }

    else if(remainRR21==3)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->RR21AllIndex[(endRR21-3)*4+0];
        const unsigned int sl0a = this->RR21AllIndex[(endRR21-3)*4+1];
        const unsigned int sl1a = this->RR21AllIndex[(endRR21-3)*4+2];
        const unsigned int sr0a = this->RR21AllIndex[(endRR21-3)*4+3];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        const double invKp0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]);
        double invKc0 = invKp0*PowPByRT_SumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Remainder block 2 (b) ==========
        const unsigned int i1 = this->RR21AllIndex[(endRR21-2)*4+0];
        const unsigned int sl0b = this->RR21AllIndex[(endRR21-2)*4+1];
        const unsigned int sl1b = this->RR21AllIndex[(endRR21-2)*4+2];
        const unsigned int sr0b = this->RR21AllIndex[(endRR21-2)*4+3];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        const double invKp1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b])*(invNegGstdByRT[sr0b]);
        double invKc1 = invKp1*PowPByRT_SumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Remainder block 3 (c) ==========
        const unsigned int i2 = this->RR21AllIndex[(endRR21-1)*4+0];
        const unsigned int sl0c = this->RR21AllIndex[(endRR21-1)*4+1];
        const unsigned int sl1c = this->RR21AllIndex[(endRR21-1)*4+2];
        const unsigned int sr0c = this->RR21AllIndex[(endRR21-1)*4+3];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        const double invKp2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c])*(invNegGstdByRT[sr0c]);
        double invKc2 = invKp2*PowPByRT_SumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }

    std::size_t endRR21Dup = this->RR21DupSize;
    std::size_t remainRR21Dup = endRR21Dup % 4;

    for(std::size_t k=0; k<endRR21Dup-remainRR21Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->RR21DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR21DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR21DupAllIndex[(k+0)*3+2];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        const double invKp0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl0a])*(invNegGstdByRT[sr0a]);
        double invKc0 = invKp0*PowPByRT_SumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0 - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR21DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->RR21DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->RR21DupAllIndex[(k+1)*3+2];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        const double invKp1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl0b])*(invNegGstdByRT[sr0b]);
        double invKc1 = invKp1*PowPByRT_SumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1 - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR21DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->RR21DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->RR21DupAllIndex[(k+2)*3+2];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        const double invKp2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl0c])*(invNegGstdByRT[sr0c]);
        double invKc2 = invKp2*PowPByRT_SumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl0c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2 - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR21DupAllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->RR21DupAllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->RR21DupAllIndex[(k+3)*3+2];
        double Kr3 = 0;
        const double Kf3 = Kf_[i3];
        const double invKp3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl0d])*(invNegGstdByRT[sr0d]);
        double invKc3 = invKp3*PowPByRT_SumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d]*c[sl0d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3 - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
    }
    if(remainRR21Dup==1)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->RR21DupAllIndex[(endRR21Dup-1)*3+0];
        const unsigned int sl0a = this->RR21DupAllIndex[(endRR21Dup-1)*3+1];
        const unsigned int sr0a = this->RR21DupAllIndex[(endRR21Dup-1)*3+2];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        const double invKp0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl0a])*(invNegGstdByRT[sr0a]);
        double invKc0 = invKp0*PowPByRT_SumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0 - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }

    else if(remainRR21Dup==2)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->RR21DupAllIndex[(endRR21Dup-2)*3+0];
        const unsigned int sl0a = this->RR21DupAllIndex[(endRR21Dup-2)*3+1];
        const unsigned int sr0a = this->RR21DupAllIndex[(endRR21Dup-2)*3+2];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        const double invKp0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl0a])*(invNegGstdByRT[sr0a]);
        double invKc0 = invKp0*PowPByRT_SumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0 - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Remainder block 2 (b) ==========
        const unsigned int i1 = this->RR21DupAllIndex[(endRR21Dup-1)*3+0];
        const unsigned int sl0b = this->RR21DupAllIndex[(endRR21Dup-1)*3+1];
        const unsigned int sr0b = this->RR21DupAllIndex[(endRR21Dup-1)*3+2];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        const double invKp1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl0b])*(invNegGstdByRT[sr0b]);
        double invKc1 = invKp1*PowPByRT_SumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1 - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
    }

    else if(remainRR21Dup==3)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->RR21DupAllIndex[(endRR21Dup-3)*3+0];
        const unsigned int sl0a = this->RR21DupAllIndex[(endRR21Dup-3)*3+1];
        const unsigned int sr0a = this->RR21DupAllIndex[(endRR21Dup-3)*3+2];
        double Kr0 = 0;
        const double Kf0 = Kf_[i0];
        const double invKp0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl0a])*(invNegGstdByRT[sr0a]);
        double invKc0 = invKp0*PowPByRT_SumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a]*c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0 - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Remainder block 2 (b) ==========
        const unsigned int i1 = this->RR21DupAllIndex[(endRR21Dup-2)*3+0];
        const unsigned int sl0b = this->RR21DupAllIndex[(endRR21Dup-2)*3+1];
        const unsigned int sr0b = this->RR21DupAllIndex[(endRR21Dup-2)*3+2];
        double Kr1 = 0;
        const double Kf1 = Kf_[i1];
        const double invKp1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl0b])*(invNegGstdByRT[sr0b]);
        double invKc1 = invKp1*PowPByRT_SumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b]*c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1 - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Remainder block 3 (c) ==========
        const unsigned int i2 = this->RR21DupAllIndex[(endRR21Dup-1)*3+0];
        const unsigned int sl0c = this->RR21DupAllIndex[(endRR21Dup-1)*3+1];
        const unsigned int sr0c = this->RR21DupAllIndex[(endRR21Dup-1)*3+2];
        double Kr2 = 0;
        const double Kf2 = Kf_[i2];
        const double invKp2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl0c])*(invNegGstdByRT[sr0c]);
        double invKc2 = invKp2*PowPByRT_SumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c]*c[sl0c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2 - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }
}


void
FastChemistry::OptReaction::RF21IR
(
    const double* __restrict__ c,
    double* __restrict__ dNdtByV,
    const double* __restrict__ ExpNegGbyRT
) const noexcept
{
    std::size_t endIR21 = this->IR21Size;
    std::size_t remainIR21 = endIR21 % 4;

    for(std::size_t k=0; k<endIR21-remainIR21; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->IR21AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->IR21AllIndex[(k+0)*4+1];
        const unsigned int sl1a = this->IR21AllIndex[(k+0)*4+2];
        const unsigned int sr0a = this->IR21AllIndex[(k+0)*4+3];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR21AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->IR21AllIndex[(k+1)*4+1];
        const unsigned int sl1b = this->IR21AllIndex[(k+1)*4+2];
        const unsigned int sr0b = this->IR21AllIndex[(k+1)*4+3];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR21AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->IR21AllIndex[(k+2)*4+1];
        const unsigned int sl1c = this->IR21AllIndex[(k+2)*4+2];
        const unsigned int sr0c = this->IR21AllIndex[(k+2)*4+3];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR21AllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->IR21AllIndex[(k+3)*4+1];
        const unsigned int sl1d = this->IR21AllIndex[(k+3)*4+2];
        const unsigned int sr0d = this->IR21AllIndex[(k+3)*4+3];
        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d]*c[sl1d];
        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
    }

    if(remainIR21==1)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->IR21AllIndex[(endIR21-1)*4+0];
        const unsigned int sl0a = this->IR21AllIndex[(endIR21-1)*4+1];
        const unsigned int sl1a = this->IR21AllIndex[(endIR21-1)*4+2];
        const unsigned int sr0a = this->IR21AllIndex[(endIR21-1)*4+3];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }

    else if(remainIR21==2)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->IR21AllIndex[(endIR21-2)*4+0];
        const unsigned int sl0a = this->IR21AllIndex[(endIR21-2)*4+1];
        const unsigned int sl1a = this->IR21AllIndex[(endIR21-2)*4+2];
        const unsigned int sr0a = this->IR21AllIndex[(endIR21-2)*4+3];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Remainder block 2 (b) ==========
        const unsigned int i1 = this->IR21AllIndex[(endIR21-1)*4+0];
        const unsigned int sl0b = this->IR21AllIndex[(endIR21-1)*4+1];
        const unsigned int sl1b = this->IR21AllIndex[(endIR21-1)*4+2];
        const unsigned int sr0b = this->IR21AllIndex[(endIR21-1)*4+3];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
    }

    else if(remainIR21==3)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->IR21AllIndex[(endIR21-3)*4+0];
        const unsigned int sl0a = this->IR21AllIndex[(endIR21-3)*4+1];
        const unsigned int sl1a = this->IR21AllIndex[(endIR21-3)*4+2];
        const unsigned int sr0a = this->IR21AllIndex[(endIR21-3)*4+3];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Remainder block 2 (b) ==========
        const unsigned int i1 = this->IR21AllIndex[(endIR21-2)*4+0];
        const unsigned int sl0b = this->IR21AllIndex[(endIR21-2)*4+1];
        const unsigned int sl1b = this->IR21AllIndex[(endIR21-2)*4+2];
        const unsigned int sr0b = this->IR21AllIndex[(endIR21-2)*4+3];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Remainder block 3 (c) ==========
        const unsigned int i2 = this->IR21AllIndex[(endIR21-1)*4+0];
        const unsigned int sl0c = this->IR21AllIndex[(endIR21-1)*4+1];
        const unsigned int sl1c = this->IR21AllIndex[(endIR21-1)*4+2];
        const unsigned int sr0c = this->IR21AllIndex[(endIR21-1)*4+3];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl1c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }

    
    std::size_t endDup = this->IR21DupSize;
    std::size_t remainDup = endDup % 4;

    for(std::size_t k=0; k<endDup-remainDup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->IR21DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR21DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR21DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0 - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR21DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->IR21DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->IR21DupAllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl0b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1 - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR21DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->IR21DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->IR21DupAllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl0c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2 - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR21DupAllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->IR21DupAllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->IR21DupAllIndex[(k+3)*3+2];
        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d]*c[sl0d];
        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3 - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
    }

    if(remainDup==1)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->IR21DupAllIndex[(endDup-1)*3+0];
        const unsigned int sl0a = this->IR21DupAllIndex[(endDup-1)*3+1];
        const unsigned int sr0a = this->IR21DupAllIndex[(endDup-1)*3+2];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0 - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }

    else if(remainDup==2)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->IR21DupAllIndex[(endDup-2)*3+0];
        const unsigned int sl0a = this->IR21DupAllIndex[(endDup-2)*3+1];
        const unsigned int sr0a = this->IR21DupAllIndex[(endDup-2)*3+2];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0 - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Remainder block 2 (b) ==========
        const unsigned int i1 = this->IR21DupAllIndex[(endDup-1)*3+0];
        const unsigned int sl0b = this->IR21DupAllIndex[(endDup-1)*3+1];
        const unsigned int sr0b = this->IR21DupAllIndex[(endDup-1)*3+2];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl0b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1 - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
    }

    else if(remainDup==3)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->IR21DupAllIndex[(endDup-3)*3+0];
        const unsigned int sl0a = this->IR21DupAllIndex[(endDup-3)*3+1];
        const unsigned int sr0a = this->IR21DupAllIndex[(endDup-3)*3+2];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a]*c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0 - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Remainder block 2 (b) ==========
        const unsigned int i1 = this->IR21DupAllIndex[(endDup-2)*3+0];
        const unsigned int sl0b = this->IR21DupAllIndex[(endDup-2)*3+1];
        const unsigned int sr0b = this->IR21DupAllIndex[(endDup-2)*3+2];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b]*c[sl0b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1 - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Remainder block 3 (c) ==========
        const unsigned int i2 = this->IR21DupAllIndex[(endDup-1)*3+0];
        const unsigned int sl0c = this->IR21DupAllIndex[(endDup-1)*3+1];
        const unsigned int sr0c = this->IR21DupAllIndex[(endDup-1)*3+2];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c]*c[sl0c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2 - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }
}


void
FastChemistry::OptReaction::RF21NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
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
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->NER21AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->NER21AllIndex[(k+1)*4+1];
        const unsigned int sl1b = this->NER21AllIndex[(k+1)*4+2];
        const unsigned int sr0b = this->NER21AllIndex[(k+1)*4+3];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->NER21AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->NER21AllIndex[(k+2)*4+1];
        const unsigned int sl1c = this->NER21AllIndex[(k+2)*4+2];
        const unsigned int sr0c = this->NER21AllIndex[(k+2)*4+3];
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->NER21AllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->NER21AllIndex[(k+3)*4+1];
        const unsigned int sl1d = this->NER21AllIndex[(k+3)*4+2];
        const unsigned int sr0d = this->NER21AllIndex[(k+3)*4+3];
        const double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        const double Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d]*c[sl1d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sl1d_increment = dNdtByV[sl1d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sl1d] = omega_sl1d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
    }
    if(remainNER21==1)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->NER21AllIndex[(endNER21-1)*4+0];
        const unsigned int sl0a = this->NER21AllIndex[(endNER21-1)*4+1];
        const unsigned int sl1a = this->NER21AllIndex[(endNER21-1)*4+2];
        const unsigned int sr0a = this->NER21AllIndex[(endNER21-1)*4+3];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }

    else if(remainNER21==2)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->NER21AllIndex[(endNER21-2)*4+0];
        const unsigned int sl0a = this->NER21AllIndex[(endNER21-2)*4+1];
        const unsigned int sl1a = this->NER21AllIndex[(endNER21-2)*4+2];
        const unsigned int sr0a = this->NER21AllIndex[(endNER21-2)*4+3];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Remainder block 2 (b) ==========
        const unsigned int i1 = this->NER21AllIndex[(endNER21-1)*4+0];
        const unsigned int sl0b = this->NER21AllIndex[(endNER21-1)*4+1];
        const unsigned int sl1b = this->NER21AllIndex[(endNER21-1)*4+2];
        const unsigned int sr0b = this->NER21AllIndex[(endNER21-1)*4+3];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
    }

    else if(remainNER21==3)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->NER21AllIndex[(endNER21-3)*4+0];
        const unsigned int sl0a = this->NER21AllIndex[(endNER21-3)*4+1];
        const unsigned int sl1a = this->NER21AllIndex[(endNER21-3)*4+2];
        const unsigned int sr0a = this->NER21AllIndex[(endNER21-3)*4+3];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sl1a_increment = dNdtByV[sl1a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sl1a] = omega_sl1a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Remainder block 2 (b) ==========
        const unsigned int i1 = this->NER21AllIndex[(endNER21-2)*4+0];
        const unsigned int sl0b = this->NER21AllIndex[(endNER21-2)*4+1];
        const unsigned int sl1b = this->NER21AllIndex[(endNER21-2)*4+2];
        const unsigned int sr0b = this->NER21AllIndex[(endNER21-2)*4+3];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sl1b_increment = dNdtByV[sl1b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sl1b] = omega_sl1b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Remainder block 3 (c) ==========
        const unsigned int i2 = this->NER21AllIndex[(endNER21-1)*4+0];
        const unsigned int sl0c = this->NER21AllIndex[(endNER21-1)*4+1];
        const unsigned int sl1c = this->NER21AllIndex[(endNER21-1)*4+2];
        const unsigned int sr0c = this->NER21AllIndex[(endNER21-1)*4+3];
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sl1c_increment = dNdtByV[sl1c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sl1c] = omega_sl1c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }

    std::size_t endNER21Dup = this->NER21DupAllIndex.size()/3;
    std::size_t remainNER21Dup = endNER21Dup % 4;

    for(std::size_t k=0; k<endNER21Dup-remainNER21Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->NER21DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER21DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER21DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0 - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->NER21DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->NER21DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->NER21DupAllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1 - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->NER21DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->NER21DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->NER21DupAllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl0c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2 - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->NER21DupAllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->NER21DupAllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->NER21DupAllIndex[(k+3)*3+2];
        const double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        const double Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d]*c[sl0d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3 - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
    }
    if(remainNER21Dup==1)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->NER21DupAllIndex[(endNER21Dup-1)*3+0];
        const unsigned int sl0a = this->NER21DupAllIndex[(endNER21Dup-1)*3+1];
        const unsigned int sr0a = this->NER21DupAllIndex[(endNER21Dup-1)*3+2];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0 - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }

    else if(remainNER21Dup==2)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->NER21DupAllIndex[(endNER21Dup-2)*3+0];
        const unsigned int sl0a = this->NER21DupAllIndex[(endNER21Dup-2)*3+1];
        const unsigned int sr0a = this->NER21DupAllIndex[(endNER21Dup-2)*3+2];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0 - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Remainder block 2 (b) ==========
        const unsigned int i1 = this->NER21DupAllIndex[(endNER21Dup-1)*3+0];
        const unsigned int sl0b = this->NER21DupAllIndex[(endNER21Dup-1)*3+1];
        const unsigned int sr0b = this->NER21DupAllIndex[(endNER21Dup-1)*3+2];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1 - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
    }

    else if(remainNER21Dup==3)
    {
        // ========== Remainder block 1 (a) ==========
        const unsigned int i0 = this->NER21DupAllIndex[(endNER21Dup-3)*3+0];
        const unsigned int sl0a = this->NER21DupAllIndex[(endNER21Dup-3)*3+1];
        const unsigned int sr0a = this->NER21DupAllIndex[(endNER21Dup-3)*3+2];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a]*c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0 - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Remainder block 2 (b) ==========
        const unsigned int i1 = this->NER21DupAllIndex[(endNER21Dup-2)*3+0];
        const unsigned int sl0b = this->NER21DupAllIndex[(endNER21Dup-2)*3+1];
        const unsigned int sr0b = this->NER21DupAllIndex[(endNER21Dup-2)*3+2];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b]*c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1 - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Remainder block 3 (c) ==========
        const unsigned int i2 = this->NER21DupAllIndex[(endNER21Dup-1)*3+0];
        const unsigned int sl0c = this->NER21DupAllIndex[(endNER21Dup-1)*3+1];
        const unsigned int sr0c = this->NER21DupAllIndex[(endNER21Dup-1)*3+2];
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c]*c[sl0c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2 - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }
}
