/*---------------------------------------------------------------------------*\
  Description
      Update the net rates of production (dNdtByV) for one-reactant /
      two-product (1-2) reactions, e.g. A -> B + B and A -> B + C.

      RR : reversible reaction; the reverse rate constant is obtained from
           the equilibrium constant
      IR : irreversible reaction; the reverse rate constant is zero
      NER: non-equilibrium reaction; the reverse rate constant is computed
           with the Arrhenius form instead of the equilibrium constant

  Functions
      RF12RR : reversible 1-2 reactions, incl. duplicate products (A=B+B)
      RF12IR : irreversible 1-2 reactions, incl. duplicate products (A=B+B)
      RF12NER: non-equilibrium 1-2 reactions

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
FastChemistry::OptReaction::RF12RR
(
    const double* __restrict__ c,
    double* __restrict__ dNdtByV,
    const double* __restrict__ ExpNegGbyRT
) const noexcept
{
    const double PowPByRTSumVki = this->Pow_pByRT_SumVki[1];
    const double invKcLimiter = FastChemistry::invKcLimiter;
    std::size_t endRR12 = this->RR12Size;
    std::size_t remainRR12 = endRR12 % 4;
    for(std::size_t k=0; k<endRR12-remainRR12; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->RR12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->RR12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->RR12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->RR12AllIndex[(k+0)*4+3];

        const double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR12AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->RR12AllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->RR12AllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->RR12AllIndex[(k+1)*4+3];
        const double Kf1 = Kf_[i1];
        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR12AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->RR12AllIndex[(k+2)*4+1];
        const unsigned int sr0c = this->RR12AllIndex[(k+2)*4+2];
        const unsigned int sr1c = this->RR12AllIndex[(k+2)*4+3];
        const double Kf2 = Kf_[i2];
        const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*(ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR12AllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->RR12AllIndex[(k+3)*4+1];
        const unsigned int sr0d = this->RR12AllIndex[(k+3)*4+2];
        const unsigned int sr1d = this->RR12AllIndex[(k+3)*4+3];
        const double Kf3 = Kf_[i3];
        const double invKp3 = (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d])*(ExpNegGbyRT[sl0d]);
        double invKc3 = invKp3*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d];
        const double CR3 = c[sr0d]*c[sr1d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
    }
    if(remainRR12==1)
    {
        std::size_t k = endRR12-1;
        const unsigned int i0 = this->RR12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->RR12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->RR12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->RR12AllIndex[(k+0)*4+3];
        const double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
    }
    else if(remainRR12==2)
    {
        std::size_t k = endRR12-2;
        const unsigned int i0 = this->RR12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->RR12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->RR12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->RR12AllIndex[(k+0)*4+3];
        const double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        const unsigned int i1 = this->RR12AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->RR12AllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->RR12AllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->RR12AllIndex[(k+1)*4+3];
        const double Kf1 = Kf_[i1];
        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
    }
    else if(remainRR12==3)
    {
        std::size_t k = endRR12-3;
        const unsigned int i0 = this->RR12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->RR12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->RR12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->RR12AllIndex[(k+0)*4+3];
        const double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        const unsigned int i1 = this->RR12AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->RR12AllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->RR12AllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->RR12AllIndex[(k+1)*4+3];
        const double Kf1 = Kf_[i1];
        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        const unsigned int i2 = this->RR12AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->RR12AllIndex[(k+2)*4+1];
        const unsigned int sr0c = this->RR12AllIndex[(k+2)*4+2];
        const unsigned int sr1c = this->RR12AllIndex[(k+2)*4+3];
        const double Kf2 = Kf_[i2];
        const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*(ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
    }

    std::size_t endRRDup12 = this->RR12DupSize;
    std::size_t remainRRDup12 = endRRDup12 % 4;
    for(std::size_t k=0; k<endRRDup12-remainRRDup12; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->RR12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR12DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR12DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->RR12DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->RR12DupAllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR12DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->RR12DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->RR12DupAllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr0c])*(ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2 + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR12DupAllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->RR12DupAllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->RR12DupAllIndex[(k+3)*3+2];
        const double Kf3 = Kf_[i3];
        const double invKp3 = (invNegGstdByRT[sr0d]*invNegGstdByRT[sr0d])*(ExpNegGbyRT[sl0d]);
        double invKc3 = invKp3*PowPByRTSumVki;
        invKc3 = std::min(invKc3, invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d];
        const double CR3 = c[sr0d]*c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3 + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
    }
    if(remainRRDup12==1)
    {
        std::size_t k = endRRDup12-1;
        const unsigned int i0 = this->RR12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR12DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }
    else if(remainRRDup12==2)
    {
        std::size_t k = endRRDup12-2;
        const unsigned int i0 = this->RR12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR12DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->RR12DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->RR12DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->RR12DupAllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
    }
    else if(remainRRDup12==3)
    {
        std::size_t k = endRRDup12-3;
        const unsigned int i0 = this->RR12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR12DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->RR12DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->RR12DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->RR12DupAllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1, invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        const unsigned int i2 = this->RR12DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->RR12DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->RR12DupAllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr0c])*(ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2*PowPByRTSumVki;
        invKc2 = std::min(invKc2, invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2 + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }
}

void
FastChemistry::OptReaction::RF12IR
(
    const double* __restrict__ c,
    double* __restrict__ dNdtByV,
    const double* __restrict__ ExpNegGbyRT
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
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR12AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->IR12AllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->IR12AllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->IR12AllIndex[(k+1)*4+3];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR12AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->IR12AllIndex[(k+2)*4+1];
        const unsigned int sr0c = this->IR12AllIndex[(k+2)*4+2];
        const unsigned int sr1c = this->IR12AllIndex[(k+2)*4+3];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR12AllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->IR12AllIndex[(k+3)*4+1];
        const unsigned int sr0d = this->IR12AllIndex[(k+3)*4+2];
        const unsigned int sr1d = this->IR12AllIndex[(k+3)*4+3];
        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d];
        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
    }
    if(remainIR12==1)
    {
        std::size_t k=endIR12-1;

        const unsigned int i0 = this->IR12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->IR12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->IR12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->IR12AllIndex[(k+0)*4+3];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

    }
    else if(remainIR12==2)
    {
        std::size_t k=endIR12-2;

        const unsigned int i0 = this->IR12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->IR12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->IR12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->IR12AllIndex[(k+0)*4+3];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        const unsigned int i1 = this->IR12AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->IR12AllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->IR12AllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->IR12AllIndex[(k+1)*4+3];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

    }
    else if(remainIR12==3)
    {
        std::size_t k=endIR12-3;

        const unsigned int i0 = this->IR12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->IR12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->IR12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->IR12AllIndex[(k+0)*4+3];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        const unsigned int i1 = this->IR12AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->IR12AllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->IR12AllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->IR12AllIndex[(k+1)*4+3];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        const unsigned int i2 = this->IR12AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->IR12AllIndex[(k+2)*4+1];
        const unsigned int sr0c = this->IR12AllIndex[(k+2)*4+2];
        const unsigned int sr1c = this->IR12AllIndex[(k+2)*4+3];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
    }

    std::size_t endIR12Dup = this->IR12DupSize;
    std::size_t remainIR12Dup = endIR12Dup % 4;
    for(std::size_t k=0; k<endIR12Dup-remainIR12Dup; k=k+4)
    {

        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->IR12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR12DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR12DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->IR12DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->IR12DupAllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR12DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->IR12DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->IR12DupAllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2 + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR12DupAllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->IR12DupAllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->IR12DupAllIndex[(k+3)*3+2];
        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d];
        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3 + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;

    }
    if(remainIR12Dup==1)
    {
        std::size_t k=endIR12Dup-1;

        const unsigned int i0 = this->IR12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR12DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }
    else if(remainIR12Dup==2)
    {
        std::size_t k=endIR12Dup-2;

        const unsigned int i0 = this->IR12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR12DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->IR12DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->IR12DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->IR12DupAllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
    }
    else if(remainIR12Dup==3)
    {
        std::size_t k=endIR12Dup-3;

        const unsigned int i0 = this->IR12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR12DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->IR12DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->IR12DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->IR12DupAllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        const unsigned int i2 = this->IR12DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->IR12DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->IR12DupAllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2 + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }
}

void
FastChemistry::OptReaction::RF12NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
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
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->NER12AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->NER12AllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->NER12AllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->NER12AllIndex[(k+1)*4+3];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->NER12AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->NER12AllIndex[(k+2)*4+1];
        const unsigned int sr0c = this->NER12AllIndex[(k+2)*4+2];
        const unsigned int sr1c = this->NER12AllIndex[(k+2)*4+3];
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->NER12AllIndex[(k+3)*4+0];
        const unsigned int sl0d = this->NER12AllIndex[(k+3)*4+1];
        const unsigned int sr0d = this->NER12AllIndex[(k+3)*4+2];
        const unsigned int sr1d = this->NER12AllIndex[(k+3)*4+3];
        const double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        const double Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d];
        const double CR3 = c[sr0d]*c[sr1d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
    }
    if(remainNER12==1)
    {
        std::size_t k = endNER12-1;

        const unsigned int i0 = this->NER12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->NER12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->NER12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->NER12AllIndex[(k+0)*4+3];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

    }
    else if(remainNER12==2)
    {
        std::size_t k = endNER12-2;

        const unsigned int i0 = this->NER12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->NER12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->NER12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->NER12AllIndex[(k+0)*4+3];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        const unsigned int i1 = this->NER12AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->NER12AllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->NER12AllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->NER12AllIndex[(k+1)*4+3];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
    }
    else if(remainNER12==3)
    {
        std::size_t k = endNER12-3;

        const unsigned int i0 = this->NER12AllIndex[(k+0)*4+0];
        const unsigned int sl0a = this->NER12AllIndex[(k+0)*4+1];
        const unsigned int sr0a = this->NER12AllIndex[(k+0)*4+2];
        const unsigned int sr1a = this->NER12AllIndex[(k+0)*4+3];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;

        const unsigned int i1 = this->NER12AllIndex[(k+1)*4+0];
        const unsigned int sl0b = this->NER12AllIndex[(k+1)*4+1];
        const unsigned int sr0b = this->NER12AllIndex[(k+1)*4+2];
        const unsigned int sr1b = this->NER12AllIndex[(k+1)*4+3];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;

        const unsigned int i2 = this->NER12AllIndex[(k+2)*4+0];
        const unsigned int sl0c = this->NER12AllIndex[(k+2)*4+1];
        const unsigned int sr0c = this->NER12AllIndex[(k+2)*4+2];
        const unsigned int sr1c = this->NER12AllIndex[(k+2)*4+3];
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr1c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
    }

    std::size_t endNERDup12 = this->NER12DupSize;
    std::size_t remainNER12Dup = endNERDup12 % 4;
    for(std::size_t k=0; k<endNERDup12-remainNER12Dup; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->NER12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER12DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->NER12DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->NER12DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->NER12DupAllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->NER12DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->NER12DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->NER12DupAllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2 + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->NER12DupAllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->NER12DupAllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->NER12DupAllIndex[(k+3)*3+2];
        const double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        const double Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d];
        const double CR3 = c[sr0d]*c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3 + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
    }
    if(remainNER12Dup==1)
    {
        std::size_t k = endNERDup12-1;

        const unsigned int i0 = this->NER12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER12DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }
    else if(remainNER12Dup==2)
    {
        std::size_t k = endNERDup12-2;

        const unsigned int i0 = this->NER12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER12DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->NER12DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->NER12DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->NER12DupAllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
    }
    else if(remainNER12Dup==3)
    {
        std::size_t k = endNERDup12-3;

        const unsigned int i0 = this->NER12DupAllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER12DupAllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER12DupAllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0 + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->NER12DupAllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->NER12DupAllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->NER12DupAllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1 + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        const unsigned int i2 = this->NER12DupAllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->NER12DupAllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->NER12DupAllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2 + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }
}
