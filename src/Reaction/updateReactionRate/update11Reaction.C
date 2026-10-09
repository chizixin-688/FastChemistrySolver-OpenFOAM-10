/*---------------------------------------------------------------------------*\
  Description
      Update the net rates of production (dNdtByV) for one-reactant /
      one-product (1-1) reactions, e.g. A -> B.

      RR : reversible reaction; the reverse rate constant is obtained from
           the equilibrium constant
      IR : irreversible reaction; the reverse rate constant is zero
      NER: non-equilibrium reaction; the reverse rate constant is computed
           with the Arrhenius form instead of the equilibrium constant

  Functions
      RF11RR : reversible 1-1 reactions
      RF11IR : irreversible 1-1 reactions
      RF11NER: non-equilibrium 1-1 reactions

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
FastChemistry::OptReaction::RF11RR
(
    const double* __restrict__ c,
    double* __restrict__ dNdtByV,
    const double* __restrict__ ExpNegGbyRT
) const noexcept
{
    std::size_t endRR11 = this->RR11Size;
    std::size_t remainRR11 = endRR11 % 4;
    const double invKcLimiter = FastChemistry::invKcLimiter;
    for(std::size_t k=0; k<endRR11-remainRR11; k=k+4)
    {
        // ========== Reaction block k+0 (a) ==========
        const unsigned int i0 = this->RR11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR11AllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->RR11AllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->RR11AllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->RR11AllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double invKp1 = (invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1;
        invKc1 = std::min(invKc1, FastChemistry::invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->RR11AllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->RR11AllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->RR11AllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        const double invKp2 = (invNegGstdByRT[sr0c])*(ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2;
        invKc2 = std::min(invKc2, FastChemistry::invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->RR11AllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->RR11AllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->RR11AllIndex[(k+3)*3+2];
        const double Kf3 = Kf_[i3];
        const double Kp3 = (ExpNegGbyRT[sr0d])*(invNegGstdByRT[sl0d]);
        double Kc3 = Kp3;
        Kc3 = std::max(Kc3, KcLimiter);
        const double Kr3 = Kf3/Kc3;
        const double CF3 = c[sl0d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        const double omega_sl0d_increment = dNdtByV[sl0d] - q3;
        const double omega_sr0d_increment = dNdtByV[sr0d] + q3;
        dNdtByV[sl0d] = omega_sl0d_increment;
        dNdtByV[sr0d] = omega_sr0d_increment;
    }
    if(remainRR11==1)
    {
        std::size_t k = endRR11-1;
        const unsigned int i0 = this->RR11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR11AllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;
    }
    else if(remainRR11==2)
    {
        std::size_t k = endRR11-2;
        const unsigned int i0 = this->RR11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR11AllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->RR11AllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->RR11AllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->RR11AllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double invKp1 = (invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1;
        invKc1 = std::min(invKc1, FastChemistry::invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

    }
    else if(remainRR11==3)
    {
        std::size_t k = endRR11-3;
        const unsigned int i0 = this->RR11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->RR11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->RR11AllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0;
        invKc0 = std::min(invKc0, invKcLimiter);
        const double Kr0 = Kf0*invKc0;
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        const double omega_sl0a_increment = dNdtByV[sl0a] - q0;
        const double omega_sr0a_increment = dNdtByV[sr0a] + q0;
        dNdtByV[sl0a] = omega_sl0a_increment;
        dNdtByV[sr0a] = omega_sr0a_increment;

        const unsigned int i1 = this->RR11AllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->RR11AllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->RR11AllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double invKp1 = (invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1;
        invKc1 = std::min(invKc1, FastChemistry::invKcLimiter);
        const double Kr1 = Kf1*invKc1;
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        const double omega_sl0b_increment = dNdtByV[sl0b] - q1;
        const double omega_sr0b_increment = dNdtByV[sr0b] + q1;
        dNdtByV[sl0b] = omega_sl0b_increment;
        dNdtByV[sr0b] = omega_sr0b_increment;

        const unsigned int i2 = this->RR11AllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->RR11AllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->RR11AllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        const double invKp2 = (invNegGstdByRT[sr0c])*(ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2;
        invKc2 = std::min(invKc2, FastChemistry::invKcLimiter);
        const double Kr2 = Kf2*invKc2;
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        const double omega_sl0c_increment = dNdtByV[sl0c] - q2;
        const double omega_sr0c_increment = dNdtByV[sr0c] + q2;
        dNdtByV[sl0c] = omega_sl0c_increment;
        dNdtByV[sr0c] = omega_sr0c_increment;
    }
}

void
FastChemistry::OptReaction::RF11IR
(
    const double* __restrict__ c,
    double* __restrict__ dNdtByV,
    const double* __restrict__ ExpNegGbyRT
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
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->IR11AllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->IR11AllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->IR11AllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;

        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->IR11AllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->IR11AllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->IR11AllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c];
        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->IR11AllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->IR11AllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->IR11AllIndex[(k+3)*3+2];
        const double Kf3 = Kf_[i3];
        const double CF3 = c[sl0d];
        const double q3 = (Kf3*CF3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
    }
    if(remainIR11==1)
    {
        std::size_t k = endIR11-1;
        const unsigned int i0 = this->IR11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR11AllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
    }
    else if(remainIR11==2)
    {
        std::size_t k = endIR11-2;
        const unsigned int i0 = this->IR11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR11AllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        const unsigned int i1 = this->IR11AllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->IR11AllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->IR11AllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
    }
    else if(remainIR11==3)
    {
        std::size_t k = endIR11-3;

        const unsigned int i0 = this->IR11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->IR11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->IR11AllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        const unsigned int i1 = this->IR11AllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->IR11AllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->IR11AllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        const double CF1 = c[sl0b];
        const double q1 = (Kf1*CF1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;

        const unsigned int i2 = this->IR11AllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->IR11AllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->IR11AllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        const double CF2 = c[sl0c];
        const double q2 = (Kf2*CF2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
    }
}

void
FastChemistry::OptReaction::RF11NER
(
    const double* __restrict__ c,
    double* __restrict__ dNdtByV,
    const double* __restrict__ ExpNegGbyRT
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
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        // ========== Reaction block k+1 (b) ==========
        const unsigned int i1 = this->NER11AllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->NER11AllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->NER11AllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        // ========== Reaction block k+2 (c) ==========
        const unsigned int i2 = this->NER11AllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->NER11AllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->NER11AllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;

        // ========== Reaction block k+3 (d) ==========
        const unsigned int i3 = this->NER11AllIndex[(k+3)*3+0];
        const unsigned int sl0d = this->NER11AllIndex[(k+3)*3+1];
        const unsigned int sr0d = this->NER11AllIndex[(k+3)*3+2];
        const double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        const double Kr3 = this->Kf_[l3];
        const double CF3 = c[sl0d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);
        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
    }
    if(remainNER11==1)
    {
        std::size_t k = endNER11-1;

        const unsigned int i0 = this->NER11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER11AllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
    }
    else if(remainNER11==2)
    {
        std::size_t k = endNER11-2;

        const unsigned int i0 = this->NER11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER11AllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        const unsigned int i1 = this->NER11AllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->NER11AllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->NER11AllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
    }
    else if(remainNER11==3)
    {
        std::size_t k = endNER11-3;

        const unsigned int i0 = this->NER11AllIndex[(k+0)*3+0];
        const unsigned int sl0a = this->NER11AllIndex[(k+0)*3+1];
        const unsigned int sr0a = this->NER11AllIndex[(k+0)*3+2];
        const double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);
        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        const unsigned int i1 = this->NER11AllIndex[(k+1)*3+0];
        const unsigned int sl0b = this->NER11AllIndex[(k+1)*3+1];
        const unsigned int sr0b = this->NER11AllIndex[(k+1)*3+2];
        const double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);
        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;

        const unsigned int i2 = this->NER11AllIndex[(k+2)*3+0];
        const unsigned int sl0c = this->NER11AllIndex[(k+2)*3+1];
        const unsigned int sr0c = this->NER11AllIndex[(k+2)*3+2];
        const double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);
        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
    }
}