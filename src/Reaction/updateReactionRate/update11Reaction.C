#include "OptReaction.H"

void 
OptReaction::update11Reaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reaction11index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    for(std::size_t i=0; i<end; i++)
    {
        const unsigned int KfIndex = this->reaction11index[i];

        const unsigned int sl0 = lhsSpeciesIndex1D11[rhsIndex+0];
        
        const unsigned int sr0 = rhsSpeciesIndex1D11[rhsIndex+0];


        double Kr = 0;
        double Kf = Kf_[KfIndex];
        if(this->isIrreversible[KfIndex]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0])/(ExpNegGbyRT[sl0]);
            double Kc = Kp;
            Kc = std::max(Kc,KcLimiter);
            Kr = Kf/Kc;         
        }
        else if(this->isIrreversible[KfIndex]==2)
        {
            unsigned int l = KfIndex - this->Ikf[1] + this->Ikf[9];
            Kr = this->Kf_[l];
        }

        const double CF = c[sl0];
        const double CR = c[sr0];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q; 
        dNdtByV[sr0] = dNdtByV[sr0] + q;


        rhsIndex = rhsIndex+1;
    }
}