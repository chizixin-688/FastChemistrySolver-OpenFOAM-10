#include "OptReaction.H"

void 
OptReaction::update33Reaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reaction33index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    std::size_t lhsIndex = 0;
    for(std::size_t i=0; i<end; i++)
    {
        const unsigned int KfIndex = this->reaction33index[i];

        const unsigned int sl0 = lhsSpeciesIndex1D33[rhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D33[rhsIndex+1];
        const unsigned int sl2 = lhsSpeciesIndex1D33[rhsIndex+2];

        const unsigned int sr0 = rhsSpeciesIndex1D33[lhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D33[lhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D33[lhsIndex+2];

        double Kr = 0;
        double Kf = Kf_[KfIndex];
        if(this->isIrreversible[KfIndex]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
            double Kc = Kp;
            Kc = std::max(Kc,KcLimiter);
            Kr = Kf/Kc;         
        }
        else if(this->isIrreversible[KfIndex]==2)
        {
            unsigned int l = KfIndex - this->Ikf[1] + this->Ikf[9];
            Kr = this->Kf_[l];
        }

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        rhsIndex = rhsIndex + 3;
        lhsIndex = lhsIndex + 3;
    }
}