#include "OptReaction.H"

void 
OptReaction::update31Reaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reaction31index.size();
    std::size_t remain = end%4;

    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    
    for(std::size_t i=0; i<end; i++)
    {
        const unsigned int KfIndex = this->reaction31index[i];

        const unsigned int sl0 = lhsSpeciesIndex1D31[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D31[lhsIndex+1];
        const unsigned int sl2 = lhsSpeciesIndex1D31[lhsIndex+2];

        const unsigned int sr0 = rhsSpeciesIndex1D31[rhsIndex+0];


        double Kr = 0;
        double Kf = Kf_[KfIndex];
        if(this->isIrreversible[KfIndex]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
            double Kc = Kp*this->Pow_pByRT_SumVki[0];
            Kc = std::max(Kc,KcLimiter);
            Kr = Kf/Kc;         
        }
        else if(this->isIrreversible[KfIndex]==2)
        {
            unsigned int l = KfIndex - this->Ikf[1] + this->Ikf[9];
            Kr = this->Kf_[l];
        }

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;


        lhsIndex = lhsIndex + 3;
        rhsIndex = rhsIndex + 1;

    }
}