/*---------------------------------------------------------------------------*\
  Description
      Computing the molar concentration reaction rate. The function 
      is used for two-one reaction, e.g. A+A=B. A+B=C
      
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

/*void 
FastChemistry::OptReaction::update21Reaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reaction21index.size();
    //std::size_t remain = end%4;
    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    for(std::size_t i=0; i<end; i++)
    {
        const unsigned int KfIndex = this->reaction21index[i];
        const unsigned int sl0 = lhsSpeciesIndex1D21[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D21[lhsIndex+1];

        const unsigned int sr0 = rhsSpeciesIndex1D21[rhsIndex+0];


        double Kr = 0;
        double Kf = Kf_[KfIndex];
        if(this->isIrreversible[KfIndex]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
            double Kc = Kp*this->Pow_pByRT_SumVki[1];
            Kc = std::max(Kc,KcLimiter);
            Kr = Kf/Kc;         
        }
        else if(this->isIrreversible[KfIndex]==2)
        {
            unsigned int l = KfIndex - this->Ikf[1] + this->Ikf[9];
            Kr = this->Kf_[l];
        }

        const double CF = c[sl0]*c[sl1];
        const double CR = c[sr0];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q; 
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;

        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 1;
    }
}*/

void 
FastChemistry::OptReaction::RF21RR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reversibleReaction21index.size();
    std::size_t remain = end%4;
    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->reversibleReaction21index[k+0];
        const unsigned int i1 = this->reversibleReaction21index[k+1];
        const unsigned int i2 = this->reversibleReaction21index[k+2];
        const unsigned int i3 = this->reversibleReaction21index[k+3];

        const unsigned int sl0a = lhsSpeciesIndex1D21RR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D21RR[lhsIndex+1];
        const unsigned int sr0a = rhsSpeciesIndex1D21RR[rhsIndex+0];

        double Kr0 = 0;
        double Kf0 = Kf_[i0];

        
        //const double Kp0 = (ExpNegGbyRT[sr0a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
        //double Kc0 = Kp0*this->Pow_pByRT_SumVki[1];
        //Kc0 = std::max(Kc0,KcLimiter);
        //Kr0 = Kf0/Kc0;
        const double invKp0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]);
        double invKc0 = invKp0*this->Pow_pByRT_SumVki[3];
        invKc0 = std::min(invKc0,FastChemistry::invKcLimiter);
        Kr0 = Kf0*invKc0;

        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;



        const unsigned int sl0b = lhsSpeciesIndex1D21RR[lhsIndex+2];
        const unsigned int sl1b = lhsSpeciesIndex1D21RR[lhsIndex+3];
        const unsigned int sr0b = rhsSpeciesIndex1D21RR[rhsIndex+1];

        double Kr1 = 0;
        double Kf1 = Kf_[i1];

        
        //const double Kp1 = (ExpNegGbyRT[sr0b])/(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]);
        //double Kc1 = Kp1*this->Pow_pByRT_SumVki[1];
        //Kc1 = std::max(Kc1,KcLimiter);
        //Kr1 = Kf1/Kc1;         
        const double invKp1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b])*(invNegGstdByRT[sr0b]);
        double invKc1 = invKp1*this->Pow_pByRT_SumVki[3];
        invKc1 = std::min(invKc1,FastChemistry::invKcLimiter);
        Kr1 = Kf1*invKc1;

        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;



        const unsigned int sl0c = lhsSpeciesIndex1D21RR[lhsIndex+4];
        const unsigned int sl1c = lhsSpeciesIndex1D21RR[lhsIndex+5];
        const unsigned int sr0c = rhsSpeciesIndex1D21RR[rhsIndex+2];

        double Kr2 = 0;
        double Kf2 = Kf_[i2];

        
        //const double Kp2 = (ExpNegGbyRT[sr0c])/(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]);
        //double Kc2 = Kp2*this->Pow_pByRT_SumVki[1];
        //Kc2 = std::max(Kc2,KcLimiter);
        //Kr2 = Kf2/Kc2;
        const double invKp2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c])*(invNegGstdByRT[sr0c]);
        double invKc2 = invKp2*this->Pow_pByRT_SumVki[3];
        invKc2 = std::min(invKc2,FastChemistry::invKcLimiter);
        Kr2 = Kf2*invKc2;
        

        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;



        const unsigned int sl0d = lhsSpeciesIndex1D21RR[lhsIndex+6];
        const unsigned int sl1d = lhsSpeciesIndex1D21RR[lhsIndex+7];
        const unsigned int sr0d = rhsSpeciesIndex1D21RR[rhsIndex+3];

        double Kr3 = 0;
        double Kf3 = Kf_[i3];

        
        //const double Kp3 = (ExpNegGbyRT[sr0d])/(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]);
        //double Kc3 = Kp3*this->Pow_pByRT_SumVki[1];
        //Kc3 = std::max(Kc3,KcLimiter);
        //Kr3 = Kf3/Kc3;
        const double invKp3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d])*(invNegGstdByRT[sr0d]);
        double invKc3 = invKp3*this->Pow_pByRT_SumVki[3];
        invKc3 = std::min(invKc3,FastChemistry::invKcLimiter);
        Kr3 = Kf3*invKc3;

        const double CF3 = c[sl0d]*c[sl1d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;

        lhsIndex = lhsIndex + 8;
        rhsIndex = rhsIndex + 4;
    }
    for(std::size_t k=end-remain; k<end; k=k+1)
    {
        const unsigned int i0 = this->reversibleReaction21index[k];
        const unsigned int sl0 = lhsSpeciesIndex1D21RR[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D21RR[lhsIndex+1];
        const unsigned int sr0 = rhsSpeciesIndex1D21RR[rhsIndex+0];


        double Kr = 0;
        double Kf = Kf_[i0];

        
        //const double Kp0 = (ExpNegGbyRT[sr0a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
        //double Kc0 = Kp0*this->Pow_pByRT_SumVki[1];
        //Kc0 = std::max(Kc0,KcLimiter);
        //Kr0 = Kf0/Kc0;         
        const double invKp = (ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1])*(invNegGstdByRT[sr0]);
        double invKc = invKp*this->Pow_pByRT_SumVki[3];
        invKc = std::min(invKc,FastChemistry::invKcLimiter);
        Kr = Kf*invKc;

        const double CF = c[sl0]*c[sl1];
        const double CR = c[sr0];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q; 
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;

        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 1;
    }
}


void 
FastChemistry::OptReaction::RF21IR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = irreversibleReaction21index.size();
    std::size_t remain = end%4;
    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->irreversibleReaction21index[k+0];
        const unsigned int i1 = this->irreversibleReaction21index[k+1];
        const unsigned int i2 = this->irreversibleReaction21index[k+2];
        const unsigned int i3 = this->irreversibleReaction21index[k+3];


        const unsigned int sl0a = lhsSpeciesIndex1D21IR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D21IR[lhsIndex+1];

        const unsigned int sr0a = rhsSpeciesIndex1D21IR[rhsIndex+0];

        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a]*c[sl1a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;



        const unsigned int sl0b = lhsSpeciesIndex1D21IR[lhsIndex+2];
        const unsigned int sl1b = lhsSpeciesIndex1D21IR[lhsIndex+3];

        const unsigned int sr0b = rhsSpeciesIndex1D21IR[rhsIndex+1];

        double Kf1 = Kf_[i1];

        const double CF1 = c[sl0b]*c[sl1b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;


        const unsigned int sl0c = lhsSpeciesIndex1D21IR[lhsIndex+4];
        const unsigned int sl1c = lhsSpeciesIndex1D21IR[lhsIndex+5];

        const unsigned int sr0c = rhsSpeciesIndex1D21IR[rhsIndex+2];

        double Kf2 = Kf_[i2];

        const double CF2 = c[sl0c]*c[sl1c];

        const double q2 = (Kf2*CF2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;



        const unsigned int sl0d = lhsSpeciesIndex1D21IR[lhsIndex+6];
        const unsigned int sl1d = lhsSpeciesIndex1D21IR[lhsIndex+7];

        const unsigned int sr0d = rhsSpeciesIndex1D21IR[rhsIndex+3];

        double Kf3 = Kf_[i3];

        const double CF3 = c[sl0d]*c[sl1d];

        const double q3 = (Kf3*CF3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;

        lhsIndex = lhsIndex + 8;
        rhsIndex = rhsIndex + 4;
    }
    for(std::size_t k=end-remain; k<end; k=k+1)
    {
        const unsigned int i0 = this->irreversibleReaction21index[k+0];


        const unsigned int sl0a = lhsSpeciesIndex1D21IR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D21IR[lhsIndex+1];

        const unsigned int sr0a = rhsSpeciesIndex1D21IR[rhsIndex+0];

        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a]*c[sl1a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 1;
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
    std::size_t end = nonEquilibriumReaction21index.size();
    std::size_t remain = end%4;
    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->nonEquilibriumReaction21index[k+0];
        const unsigned int i1 = this->nonEquilibriumReaction21index[k+1];
        const unsigned int i2 = this->nonEquilibriumReaction21index[k+2];
        const unsigned int i3 = this->nonEquilibriumReaction21index[k+3];


        const unsigned int sl0a = lhsSpeciesIndex1D21NER[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D21NER[lhsIndex+1];

        const unsigned int sr0a = rhsSpeciesIndex1D21NER[rhsIndex+0];

        double Kr0 = 0;
        double Kf0 = Kf_[i0];

        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];

        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;



        const unsigned int sl0b = lhsSpeciesIndex1D21NER[lhsIndex+2];
        const unsigned int sl1b = lhsSpeciesIndex1D21NER[lhsIndex+3];

        const unsigned int sr0b = rhsSpeciesIndex1D21NER[rhsIndex+1];

        double Kr1 = 0;
        double Kf1 = Kf_[i1];

        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        Kr1 = this->Kf_[l1];

        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;



        const unsigned int sl0c = lhsSpeciesIndex1D21NER[lhsIndex+4];
        const unsigned int sl1c = lhsSpeciesIndex1D21NER[lhsIndex+5];

        const unsigned int sr0c = rhsSpeciesIndex1D21NER[rhsIndex+2];

        double Kr2 = 0;
        double Kf2 = Kf_[i2];

        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        Kr2 = this->Kf_[l2];

        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;


        const unsigned int sl0d = lhsSpeciesIndex1D21NER[lhsIndex+6];
        const unsigned int sl1d = lhsSpeciesIndex1D21NER[lhsIndex+7];

        const unsigned int sr0d = rhsSpeciesIndex1D21NER[rhsIndex+3];

        double Kr3 = 0;
        double Kf3 = Kf_[i3];

        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        Kr3 = this->Kf_[l3];

        const double CF3 = c[sl0d]*c[sl1d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;

        lhsIndex = lhsIndex + 8;
        rhsIndex = rhsIndex + 4;
    }
    for(std::size_t k=end-remain; k<end; k=k+1)
    {
        const unsigned int i0 = this->nonEquilibriumReaction21index[k+0];


        const unsigned int sl0a = lhsSpeciesIndex1D21NER[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D21NER[lhsIndex+1];

        const unsigned int sr0a = rhsSpeciesIndex1D21NER[rhsIndex+0];

        double Kr0 = 0;
        double Kf0 = Kf_[i0];

        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        Kr0 = this->Kf_[l0];

        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 1;
    }
}