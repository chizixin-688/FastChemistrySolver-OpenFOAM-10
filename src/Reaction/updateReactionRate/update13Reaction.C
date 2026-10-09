/*---------------------------------------------------------------------------*\
  Description
      Update the net rates of production (dNdtByV) for one-reactant /
      three-product (1-3) reactions, e.g. A -> B + B + C and A -> B + C + D.

      RR : reversible reaction; the reverse rate constant is obtained from
           the equilibrium constant
      IR : irreversible reaction; the reverse rate constant is zero
      NER: non-equilibrium reaction; the reverse rate constant is computed
           with the Arrhenius form instead of the equilibrium constant

  Functions
      RF13RR : reversible 1-3 reactions, incl. duplicate products (A=B+B+B)
      RF13IR : irreversible 1-3 reactions, incl. duplicate products (A=B+B+B)
      RF13NER: non-equilibrium 1-3 reactions

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
FastChemistry::OptReaction::RF13RR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    const double PowPByRTSumVki = this->Pow_pByRT_SumVki[0];
    const double invKcLimiter = FastChemistry::invKcLimiter;
    std::size_t endRR13 = this->RR13AllIndex.size()/5;
    std::size_t remainRR13 = endRR13%4;
    for(std::size_t k=0; k<endRR13-remainRR13; k=k+4)
    {
        const unsigned int i0 = this->RR13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->RR13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->RR13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->RR13AllIndex[(k+0)*5+4];
        double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
            (invNegGstdByRT[sr2a]*ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0,invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->RR13AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR13AllIndex[(k+1)*5+1];
        const unsigned int sr0b = this->RR13AllIndex[(k+1)*5+2];
        const unsigned int sr1b = this->RR13AllIndex[(k+1)*5+3];
        const unsigned int sr2b = this->RR13AllIndex[(k+1)*5+4];
        double Kf1 = Kf_[i1];
        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*
            (invNegGstdByRT[sr2b]*ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1,invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        const unsigned int i2 = this->RR13AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->RR13AllIndex[(k+2)*5+1];
        const unsigned int sr0c = this->RR13AllIndex[(k+2)*5+2];
        const unsigned int sr1c = this->RR13AllIndex[(k+2)*5+3];
        const unsigned int sr2c = this->RR13AllIndex[(k+2)*5+4];
        double Kf2 = Kf_[i2];
        const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*
            (invNegGstdByRT[sr2c]*ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2*PowPByRTSumVki;
        invKc2 = std::min(invKc2,invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;

        const unsigned int i3 = this->RR13AllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->RR13AllIndex[(k+3)*5+1];
        const unsigned int sr0d = this->RR13AllIndex[(k+3)*5+2];
        const unsigned int sr1d = this->RR13AllIndex[(k+3)*5+3];
        const unsigned int sr2d = this->RR13AllIndex[(k+3)*5+4];
        double Kf3 = Kf_[i3];
        const double invKp3 = (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d])*
            (invNegGstdByRT[sr2d]*ExpNegGbyRT[sl0d]);
        double invKc3 = invKp3*PowPByRTSumVki;
        invKc3 = std::min(invKc3,invKcLimiter);
        const double Kr3 = Kf3*invKc3;
        const double CF3 = c[sl0d];
        const double CR3 = c[sr0d]*c[sr1d]*c[sr2d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        const double omega_sr2d_increment = dNdtByV[sr2d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
        dNdtByV[sr2d] = omega_sr2d_increment;

    }
    if(remainRR13==1)
    {
        std::size_t k = endRR13-1;

        const unsigned int i0 = this->RR13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->RR13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->RR13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->RR13AllIndex[(k+0)*5+4];
        double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
            (invNegGstdByRT[sr2a]*ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0,invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

    }
    else if(remainRR13==2)
    {
        std::size_t k = endRR13-2;

        const unsigned int i0 = this->RR13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->RR13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->RR13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->RR13AllIndex[(k+0)*5+4];
        double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
            (invNegGstdByRT[sr2a]*ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0,invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->RR13AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR13AllIndex[(k+1)*5+1];
        const unsigned int sr0b = this->RR13AllIndex[(k+1)*5+2];
        const unsigned int sr1b = this->RR13AllIndex[(k+1)*5+3];
        const unsigned int sr2b = this->RR13AllIndex[(k+1)*5+4];
        double Kf1 = Kf_[i1];
        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*
            (invNegGstdByRT[sr2b]*ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1,invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

    }
    else if(remainRR13==3)
    {
        std::size_t k = endRR13-3;

        const unsigned int i0 = this->RR13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->RR13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->RR13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->RR13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->RR13AllIndex[(k+0)*5+4];
        double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
            (invNegGstdByRT[sr2a]*ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0*PowPByRTSumVki;
        invKc0 = std::min(invKc0,invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->RR13AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->RR13AllIndex[(k+1)*5+1];
        const unsigned int sr0b = this->RR13AllIndex[(k+1)*5+2];
        const unsigned int sr1b = this->RR13AllIndex[(k+1)*5+3];
        const unsigned int sr2b = this->RR13AllIndex[(k+1)*5+4];
        double Kf1 = Kf_[i1];
        const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b])*
            (invNegGstdByRT[sr2b]*ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1*PowPByRTSumVki;
        invKc1 = std::min(invKc1,invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        const unsigned int i2 = this->RR13AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->RR13AllIndex[(k+2)*5+1];
        const unsigned int sr0c = this->RR13AllIndex[(k+2)*5+2];
        const unsigned int sr1c = this->RR13AllIndex[(k+2)*5+3];
        const unsigned int sr2c = this->RR13AllIndex[(k+2)*5+4];
        double Kf2 = Kf_[i2];
        const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c])*
            (invNegGstdByRT[sr2c]*ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2*PowPByRTSumVki;
        invKc2 = std::min(invKc2,invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;

    }

    // ========== RF13RR: duplicate-product reactions (e.g. A = B + B + C) ==========
    {
        const std::size_t endRRDup13 = this->RR13DupSize;
        const double invKcLimiterDup = FastChemistry::invKcLimiter;
        for(std::size_t k=0; k<endRRDup13; ++k)
        {
            const unsigned int i0 = this->RR13DupAllIndex[(k+0)*4+0];
            const unsigned int sl0a = this->RR13DupAllIndex[(k+0)*4+1];
            const unsigned int sr0a = this->RR13DupAllIndex[(k+0)*4+2];
            const unsigned int sr1a = this->RR13DupAllIndex[(k+0)*4+3];
            const double Kf0 = Kf_[i0];
            // sr1a is the duplicated product (appears twice)
            const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a])*
                (invNegGstdByRT[sr1a]*ExpNegGbyRT[sl0a]);
            double invKc0 = invKp0*PowPByRTSumVki;
            invKc0 = std::min(invKc0, invKcLimiterDup);
            const double Kr0 = Kf0*invKc0;
            const double CF0 = c[sl0a];
            const double CR0 = c[sr0a]*c[sr1a]*c[sr1a];
            const double q0 = (Kf0*CF0) - (Kr0*CR0);
            dNdtByV[sl0a] -= q0;
            dNdtByV[sr0a] += q0;
            dNdtByV[sr1a] += q0 + q0;
        }
    }

}


void 
FastChemistry::OptReaction::RF13IR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t endIR13 = this->IR13AllIndex.size()/5;
    std::size_t remainIR13 = endIR13%4;
    for(std::size_t k=0; k<endIR13-remainIR13; k=k+4)
    {
        const unsigned int i0 = this->IR13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->IR13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->IR13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->IR13AllIndex[(k+0)*5+4];
        double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->IR13AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR13AllIndex[(k+1)*5+1];
        const unsigned int sr0b = this->IR13AllIndex[(k+1)*5+2];
        const unsigned int sr1b = this->IR13AllIndex[(k+1)*5+3];
        const unsigned int sr2b = this->IR13AllIndex[(k+1)*5+4];
        double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        const unsigned int i2 = this->IR13AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->IR13AllIndex[(k+2)*5+1];
        const unsigned int sr0c = this->IR13AllIndex[(k+2)*5+2];
        const unsigned int sr1c = this->IR13AllIndex[(k+2)*5+3];
        const unsigned int sr2c = this->IR13AllIndex[(k+2)*5+4];
        double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;

        const unsigned int i3 = this->IR13AllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->IR13AllIndex[(k+3)*5+1];
        const unsigned int sr0d = this->IR13AllIndex[(k+3)*5+2];
        const unsigned int sr1d = this->IR13AllIndex[(k+3)*5+3];
        const unsigned int sr2d = this->IR13AllIndex[(k+3)*5+4];
        double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d];
        const double q3 = (Kf3*CF3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        const double omega_sr2d_increment = dNdtByV[sr2d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
        dNdtByV[sr2d] = omega_sr2d_increment;

    }
    if(remainIR13==1)
    {
        std::size_t k = endIR13-1;

        const unsigned int i0 = this->IR13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->IR13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->IR13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->IR13AllIndex[(k+0)*5+4];
        double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

    }
    else if(remainIR13==2)
    {
        std::size_t k = endIR13-2;

        const unsigned int i0 = this->IR13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->IR13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->IR13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->IR13AllIndex[(k+0)*5+4];
        double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->IR13AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR13AllIndex[(k+1)*5+1];
        const unsigned int sr0b = this->IR13AllIndex[(k+1)*5+2];
        const unsigned int sr1b = this->IR13AllIndex[(k+1)*5+3];
        const unsigned int sr2b = this->IR13AllIndex[(k+1)*5+4];
        double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

    }
    else if(remainIR13==3)
    {
        std::size_t k = endIR13-3;

        const unsigned int i0 = this->IR13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->IR13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->IR13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->IR13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->IR13AllIndex[(k+0)*5+4];
        double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->IR13AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->IR13AllIndex[(k+1)*5+1];
        const unsigned int sr0b = this->IR13AllIndex[(k+1)*5+2];
        const unsigned int sr1b = this->IR13AllIndex[(k+1)*5+3];
        const unsigned int sr2b = this->IR13AllIndex[(k+1)*5+4];
        double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b];
        const double q1 = (Kf1*CF1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        const unsigned int i2 = this->IR13AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->IR13AllIndex[(k+2)*5+1];
        const unsigned int sr0c = this->IR13AllIndex[(k+2)*5+2];
        const unsigned int sr1c = this->IR13AllIndex[(k+2)*5+3];
        const unsigned int sr2c = this->IR13AllIndex[(k+2)*5+4];
        double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c];
        const double q2 = (Kf2*CF2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;

    }

    // ========== RF13IR: duplicate-product reactions (e.g. A = B + B + C) ==========
    {
        const std::size_t endIRDup13 = this->IR13DupSize;
        for(std::size_t k=0; k<endIRDup13; ++k)
        {
            const unsigned int i0 = this->IR13DupAllIndex[(k+0)*4+0];
            const unsigned int sl0a = this->IR13DupAllIndex[(k+0)*4+1];
            const unsigned int sr0a = this->IR13DupAllIndex[(k+0)*4+2];
            const unsigned int sr1a = this->IR13DupAllIndex[(k+0)*4+3];
            const double Kf0 = Kf_[i0];
            const double CF0 = c[sl0a];
            const double q0 = (Kf0*CF0);
            dNdtByV[sl0a] -= q0;
            dNdtByV[sr0a] += q0;
            dNdtByV[sr1a] += q0 + q0;
        }
    }

}

void 
FastChemistry::OptReaction::RF13NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t endNER13 = this->NER13AllIndex.size()/5;
    std::size_t remainNER13 = endNER13%4;
    for(std::size_t k=0; k<endNER13-remainNER13; k=k+4)
    {
        const unsigned int i0 = this->NER13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->NER13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->NER13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->NER13AllIndex[(k+0)*5+4];
        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->NER13AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER13AllIndex[(k+1)*5+1];
        const unsigned int sr0b = this->NER13AllIndex[(k+1)*5+2];
        const unsigned int sr1b = this->NER13AllIndex[(k+1)*5+3];
        const unsigned int sr2b = this->NER13AllIndex[(k+1)*5+4];
        double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        const unsigned int i2 = this->NER13AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER13AllIndex[(k+2)*5+1];
        const unsigned int sr0c = this->NER13AllIndex[(k+2)*5+2];
        const unsigned int sr1c = this->NER13AllIndex[(k+2)*5+3];
        const unsigned int sr2c = this->NER13AllIndex[(k+2)*5+4];
        double Kf2 = Kf_[i2]; 
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;

        const unsigned int i3 = this->NER13AllIndex[(k+3)*5+0];
        const unsigned int sl0d = this->NER13AllIndex[(k+3)*5+1];
        const unsigned int sr0d = this->NER13AllIndex[(k+3)*5+2];
        const unsigned int sr1d = this->NER13AllIndex[(k+3)*5+3];
        const unsigned int sr2d = this->NER13AllIndex[(k+3)*5+4];
        double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        const double Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d];
        const double CR3 = c[sr0d]*c[sr1d]*c[sr2d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        const double omega_sr1d_increment = dNdtByV[sr1d] + q3;
        const double omega_sr2d_increment = dNdtByV[sr2d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
        dNdtByV[sr1d] = omega_sr1d_increment;
        dNdtByV[sr2d] = omega_sr2d_increment;

    }
    if(remainNER13==1)
    {
        std::size_t k = endNER13-1;

        const unsigned int i0 = this->NER13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->NER13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->NER13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->NER13AllIndex[(k+0)*5+4];
        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

    }
    else if(remainNER13==2)
    {
        std::size_t k = endNER13-2;

        const unsigned int i0 = this->NER13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->NER13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->NER13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->NER13AllIndex[(k+0)*5+4];
        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->NER13AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER13AllIndex[(k+1)*5+1];
        const unsigned int sr0b = this->NER13AllIndex[(k+1)*5+2];
        const unsigned int sr1b = this->NER13AllIndex[(k+1)*5+3];
        const unsigned int sr2b = this->NER13AllIndex[(k+1)*5+4];
        double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

    }
    else if(remainNER13==3)
    {
        std::size_t k = endNER13-3;

        const unsigned int i0 = this->NER13AllIndex[(k+0)*5+0];
        const unsigned int sl0a = this->NER13AllIndex[(k+0)*5+1];
        const unsigned int sr0a = this->NER13AllIndex[(k+0)*5+2];
        const unsigned int sr1a = this->NER13AllIndex[(k+0)*5+3];
        const unsigned int sr2a = this->NER13AllIndex[(k+0)*5+4];
        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        const double omega_sr1a_increment = dNdtByV[sr1a] + q0;
        const double omega_sr2a_increment = dNdtByV[sr2a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
        dNdtByV[sr1a] = omega_sr1a_increment;
        dNdtByV[sr2a] = omega_sr2a_increment;

        const unsigned int i1 = this->NER13AllIndex[(k+1)*5+0];
        const unsigned int sl0b = this->NER13AllIndex[(k+1)*5+1];
        const unsigned int sr0b = this->NER13AllIndex[(k+1)*5+2];
        const unsigned int sr1b = this->NER13AllIndex[(k+1)*5+3];
        const unsigned int sr2b = this->NER13AllIndex[(k+1)*5+4];
        double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        const double omega_sr1b_increment = dNdtByV[sr1b] + q1;
        const double omega_sr2b_increment = dNdtByV[sr2b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;
        dNdtByV[sr1b] = omega_sr1b_increment;
        dNdtByV[sr2b] = omega_sr2b_increment;

        const unsigned int i2 = this->NER13AllIndex[(k+2)*5+0];
        const unsigned int sl0c = this->NER13AllIndex[(k+2)*5+1];
        const unsigned int sr0c = this->NER13AllIndex[(k+2)*5+2];
        const unsigned int sr1c = this->NER13AllIndex[(k+2)*5+3];
        const unsigned int sr2c = this->NER13AllIndex[(k+2)*5+4];
        double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        const double omega_sr1c_increment = dNdtByV[sr1c] + q2;
        const double omega_sr2c_increment = dNdtByV[sr2c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
        dNdtByV[sr1c] = omega_sr1c_increment;
        dNdtByV[sr2c] = omega_sr2c_increment;

    }

    // ========== RF13NER: duplicate-product reactions (e.g. A = B + B + C) ==========
    {
        const std::size_t endNERDup13 = this->NER13DupSize;
        for(std::size_t k=0; k<endNERDup13; ++k)
        {
            const unsigned int i0 = this->NER13DupAllIndex[(k+0)*4+0];
            const unsigned int sl0a = this->NER13DupAllIndex[(k+0)*4+1];
            const unsigned int sr0a = this->NER13DupAllIndex[(k+0)*4+2];
            const unsigned int sr1a = this->NER13DupAllIndex[(k+0)*4+3];
            const double Kf0 = Kf_[i0];
            unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
            const double Kr0 = this->Kf_[l0];
            const double CF0 = c[sl0a];
            const double CR0 = c[sr0a]*c[sr1a]*c[sr1a];
            const double q0 = (Kf0*CF0) - (Kr0*CR0);
            dNdtByV[sl0a] -= q0;
            dNdtByV[sr0a] += q0;
            dNdtByV[sr1a] += q0 + q0;
        }
    }
}
