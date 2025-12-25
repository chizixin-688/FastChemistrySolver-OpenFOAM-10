/*---------------------------------------------------------------------------*\
  Description
      Computing the molar concentration reaction rate. The function 
      is used for general reaction with non-integer stoichiometric number,e.g.
      4.5A+5.5B=6.5C

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"

//=============================================================================//

void 
FastChemistry::OptReaction::RFGNI
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reactionGNIindex.size();
    for(std::size_t k=0; k<end; k++)
    {
        int i = this->reactionGNIindex[k];

        double CR = 1.0;
        double CF = 1.0;
        const double Kf = this->Kf_[i];
        double Kc_ = 1.0;
        double Kr = 0;
        if(this->isIrreversible[i]==0) // Reversible reaction
        {
            double Kp = 0;
            double sumVki = 0;
            for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->lhsSpeciesIndex[i][j];
                const double sl = this->lhsStoichCoeff[i][j];
                const double el = this->lhsReactionOrder[i][j];
                Kp += sl*this->negGstdByRT[si];
                sumVki = sumVki - sl;
                CF = CF * (c[si] >= FastChemistry::ConcentrationLimiter || el >= 1 ? std::pow(c[si], el) : 0.0);
            }
            
            for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->rhsSpeciesIndex[i][j];
                const double sr = this->rhsStoichCoeff[i][j];
                const double er = this->rhsReactionOrder[i][j];
                Kp -= sr*this->negGstdByRT[si];
                sumVki = sumVki + sr;
                CR = CR * (c[si] >= FastChemistry::ConcentrationLimiter || er >= 1 ? std::pow(c[si], er) : 0.0);            
            }
            __m256d KpPow = _mm256_setr_pd(Kp,sumVki*(this->logP-this->logRuT),1,1);
            KpPow = vec256_expd(KpPow);

            //Kp = std::exp(Kp);
            //double r = std::pow(this->Pstd/(this->Ru*this->T),sumVki);
            //Kc_ = Kp*r;
            Kp = this->get_elem0(KpPow);
            double r = this->get_elem1(KpPow);
            Kc_ = Kp*r;
            
            Kc_ = std::max(Kc_,KcLimiter);
            Kr = Kf/Kc_;
        }
        else if(this->isIrreversible[i]==2)
        {
            unsigned int J = i - this->Ikf[1] + this->Ikf[9];
            Kr = this->Kf_[J];

            for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->lhsSpeciesIndex[i][j];
                const double el = this->lhsReactionOrder[i][j];
                CF = CF * (c[si] >= FastChemistry::ConcentrationLimiter || el >= 1 ? std::pow(c[si], el) : 0.0);
            }
            
            for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->rhsSpeciesIndex[i][j];
                const double er = this->rhsReactionOrder[i][j];
                CR = CR * (c[si] >= FastChemistry::ConcentrationLimiter || er >= 1 ? std::pow(c[si], er) : 0.0);            
            }
        }
        else
        {
            for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->lhsSpeciesIndex[i][j];
                const double el = this->lhsReactionOrder[i][j];
                CF = CF * (c[si] >= FastChemistry::ConcentrationLimiter || el >= 1 ? std::pow(c[si], el) : 0.0);
            }
            
            for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->rhsSpeciesIndex[i][j];
                const double er = this->rhsReactionOrder[i][j];
                CR = CR * (c[si] >= FastChemistry::ConcentrationLimiter || er >= 1 ? std::pow(c[si], er) : 0.0);
            }
        }

        const double q = Kf*CF-Kr*CR;

        for(unsigned int j = 0; j < lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = lhsSpeciesIndex[i][j];
            const double sl = lhsStoichCoeff[i][j];
            dNdtByV[si] = dNdtByV[si] - sl*q;
        }
        
        for(unsigned int j = 0; j < rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = rhsSpeciesIndex[i][j];
            const double sl = rhsStoichCoeff[i][j];        
            dNdtByV[si] = dNdtByV[si] + sl*q;
        }
    }
}