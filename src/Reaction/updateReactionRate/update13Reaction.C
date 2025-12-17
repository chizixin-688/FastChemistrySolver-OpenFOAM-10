#include "OptReaction.H"

void 
OptReaction::update13Reaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reaction13index.size();
    //std::size_t remain = end%4;

    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;

    for(std::size_t i=0; i<end; i++)
    {
        const unsigned int KfIndex = this->reaction13index[i];

        const unsigned int sl0 = lhsSpeciesIndex1D13[lhsIndex+0];

        const unsigned int sr0 = rhsSpeciesIndex1D13[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D13[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D13[rhsIndex+2];


        double Kr = 0;
        double Kf = Kf_[KfIndex];
        if(this->isIrreversible[KfIndex]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]);
            double Kc = Kp*this->Pow_pByRT_SumVki[4];
            Kc = std::max(Kc,KcLimiter);
            Kr = Kf/Kc;         
        }
        else if(this->isIrreversible[KfIndex]==2)
        {
            unsigned int l = KfIndex - this->Ikf[1] + this->Ikf[9];
            Kr = this->Kf_[l];
        }

        const double CF = c[sl0];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 3;

    }
}