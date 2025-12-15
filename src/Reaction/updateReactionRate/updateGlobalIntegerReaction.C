#include "OptReaction.H"

void 
OptReaction::updateGlobalIntegerReaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reactionGIindex.size();
    for(std::size_t k=0; k<end; k++)
    {

        int i = this->reactionGIindex[k];

        int sumVki = 0;
        double Kp = 1.0;
        double CR = 1.0;
        double CF = 1.0;
        const double Kf = this->Kf_[i];
        double Kc_ = 1.0;
        double Kr = 0;
        if(this->isIrreversible[i]==0)
        {
            for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->lhsSpeciesIndex[i][j];
                Kp = Kp / ExpNegGbyRT[si];
                sumVki = sumVki - 1;
                CF = CF * c[si];
            }
            
            for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->rhsSpeciesIndex[i][j];
                Kp = Kp * ExpNegGbyRT[si];
                sumVki = sumVki + 1;
                CR = CR * c[si];
            }

            Kc_ = Kp*this->Pow_pByRT_SumVki_I[sumVki];
            Kc_ = std::max(Kc_,KcLimiter);
            Kr = Kf/Kc_;
        }
        else if(this->isIrreversible[i]==2)
        {
            unsigned int J = i - this->Ikf[1] + this->Ikf[9];

            for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->lhsSpeciesIndex[i][j];
                CF = CF * c[si];
            }
            
            for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->rhsSpeciesIndex[i][j];
                CR = CR * c[si];
            }
            Kr = this->Kf_[J];
        }
        else
        {
            for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->lhsSpeciesIndex[i][j];
                CF = CF * c[si];
            }
            
            for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->rhsSpeciesIndex[i][j];
                CR = CR * c[si];
            }
        }

        const double q = Kf*CF-Kr*CR;

        for(unsigned int j = 0; j < lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = lhsSpeciesIndex[i][j];
            dNdtByV[si] = dNdtByV[si] - q;
        }
        
        for(unsigned int j = 0; j < rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = rhsSpeciesIndex[i][j];
            dNdtByV[si] = dNdtByV[si] + q;
        }
    }
}