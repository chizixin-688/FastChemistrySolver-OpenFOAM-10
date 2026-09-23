/*---------------------------------------------------------------------------*\
  Description
      Global Jacobian update of the per-volume molar production rate for a
      reaction with integer stoichiometric coefficients (e.g. 4A+5B=6C).

      For every reaction handled here the routine
        - assembles the forward/reverse rate coefficients (equilibrium
          constant built from the species Gibbs factors invNegGstdByRT /
          tmp_Exp, or stored reverse rates for the paired-reaction case),
        - forms the net rate        q    = Kf*CF - Kr*CR
          and its temperature part  dqdT = dKfdT*CF - dKrdT*CR
        - scatters q / dqdT into the species rates dNdtByV and into the
          temperature column ddNdtByVdcTp[*][nSpecies] of the Jacobian, and
        - scatters the concentration derivatives of Kf*CF and Kr*CR into the
          species-concentration columns of ddNdtByVdcTp.

      CF and CR are the products of the reactant / product concentrations
      (integer stoichiometry is encoded by repeating the species in the
      lists).  Third-body reactions finally sweep the whole species block in
      4-wide SIMD columns.

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"

//=============================================================================//

void  FastChemistry::OptReaction::JFGI
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
) const noexcept
{

    std::size_t end = reactionGIindex.size();

    // Sweep over the reactions that require the global Jacobian update.
    for (std::size_t K=0; K<end; K++)
    {
        // Rate data of the current reaction.  CF and CR below accumulate
        // the concentration products of the reactants / products.
        const unsigned int i = this->reactionGIindex[K];
        const double Kf = this->Kf_[i];
        const double dKfdT = this->dKfdT_[i];
        double Kr = 0;
        double dKrdT = 0;
        double CF = 1.0;
        double CR = 1.0;
        double invKc = 0;
        // Reversible reaction (code 0): assemble Kc from the species Gibbs
        // factors, then Kr = Kf/Kc and its temperature derivative.
        if(this->isIrreversible[i]==0)
        {
            int sumVki = 0;
            double Kp = 1.0;
            double sumVdBdT = 0.0;
            for (unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->lhsSpeciesIndex[i][j];
                Kp = Kp * this->invNegGstdByRT[si];
                sumVki = sumVki - 1;
                sumVdBdT = sumVdBdT - dBdT[si];
                CF = CF * C[si];
            }

            for (unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->rhsSpeciesIndex[i][j];
                Kp = Kp * this->tmp_Exp[si];
                sumVki = sumVki + 1;
                sumVdBdT = sumVdBdT + dBdT[si];
                CR = CR * C[si];
            }

            const double Kc = std::max(Kp*this->Pow_pByRT_SumVki_I[sumVki],FastChemistry::KcLimiter);
            invKc = 1.0/Kc;

            const double dKcdTByKc = sumVdBdT - sumVki*invT;
            Kr = Kf*invKc;
            dKrdT = dKfdT*invKc - (Kc > FastChemistry::KcLimiter ? Kr*dKcdTByKc : 0);
        }
        // Code 2: the reverse rate is the stored forward rate of the
        // paired reverse reaction; only the concentration products remain.
        else if(this->isIrreversible[i]==2)
        {
            Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
            dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
            invKc = Kr/Kf;
            for (unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->rhsSpeciesIndex[i][j];
                CR = CR * C[si];
            }
            for (unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
            {
                const unsigned int sl = this->lhsSpeciesIndex[i][j];
                CF = CF * C[sl];
            }
        }
        else
        {
            for (unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
            {
                const unsigned int si = this->rhsSpeciesIndex[i][j];
                CR = CR * C[si];
            }
            for (unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
            {
                const unsigned int sl = this->lhsSpeciesIndex[i][j];
                CF = CF * C[sl];
            }
        }
        // Net molar production rate of this reaction and its temperature
        // derivative (sign conventions: reactants consume, products form).
        const double q = Kf*CF - Kr*CR;
        // Scatter q into the species rates dNdtByV and dqdT into the
        // temperature column ddNdtByVdcTp[*][nSpecies] of the Jacobian.
        // Integer stoichiometry is encoded by repeated list entries, so a
        // species with coefficient nu appears nu times.
        const double dqdT = (dKfdT*CF-dKrdT*CR);

        for (unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->lhsSpeciesIndex[i][j];
            ddNdtByVdcTp[si*(this->alignN) + (this->nSpecies)] -= dqdT;
            dNdtByV[si] = dNdtByV[si] - q;
        }

        for (unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[i][j];
            ddNdtByVdcTp[si*(this->alignN) + (this->nSpecies)] += dqdT;
            dNdtByV[si] = dNdtByV[si] + q;
        }

        for (unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int sj = this->lhsSpeciesIndex[i][j];

            // Derivative of the forward flux Kf*CF with respect to the
            // reactant concentration C_j (product rule over the reactants;
            // the factor of the differentiated entry j drops out).
            double dCfdCj = 1.0;
            for (unsigned int k = 0; k < this->lhsSpeciesIndex[i].size();k++)
            {
                const unsigned int sk = this->lhsSpeciesIndex[i][k];

                dCfdCj = dCfdCj*
                (
                    (j==k)? 1:C[sk]
                );
            }

            const double KfdCfdCj = Kf*dCfdCj;

            for (unsigned int k = 0; k < this->lhsSpeciesIndex[i].size();k++)
            {
                const unsigned int sk = this->lhsSpeciesIndex[i][k];
                ddNdtByVdcTp[sk*(this->alignN)+sj] -= KfdCfdCj;
            }

            for (unsigned int k = 0; k < rhsSpeciesIndex[i].size();k++)
            {
                const unsigned int sk = rhsSpeciesIndex[i][k];
                ddNdtByVdcTp[sk*(this->alignN)+sj] += KfdCfdCj;
            }
        }
        // Derivative of the reverse flux Kr*CR with respect to the product
        // concentration C_j (reversible reactions and code 2 only).
        if(this->isIrreversible[i]==0||this->isIrreversible[i]==2)
        {
            for (unsigned int j = 0; j < rhsSpeciesIndex[i].size();j++)
            {
                const unsigned int sj = rhsSpeciesIndex[i][j];

                // Product rule over the products; the factor of the
                // differentiated entry j drops out.
                double dCrdCj = 1.0;
                for (unsigned int k = 0; k < rhsSpeciesIndex[i].size();k++)
                {
                    const unsigned int sk = rhsSpeciesIndex[i][k];

                    dCrdCj = dCrdCj*
                    (
                        (j==k)? 1:C[sk]
                    );
                }

                const double negKrdCrdCj = -Kr*dCrdCj;

                for (unsigned int k = 0; k < this->lhsSpeciesIndex[i].size();k++)
                {
                    const unsigned int sk = this->lhsSpeciesIndex[i][k];
                    ddNdtByVdcTp[sk*(this->alignN)+sj] -= negKrdCrdCj;
                }

                for (unsigned int k = 0; k < rhsSpeciesIndex[i].size();k++)
                {
                    const unsigned int sk = rhsSpeciesIndex[i][k];
                    ddNdtByVdcTp[sk*(this->alignN)+sj] += negKrdCrdCj;
                }
            }
        }

        // Third-body (M) reactions: the rate coefficients depend on the
        // collision-partner concentrations.  Accumulate the M-derivative
        // dWdM = dKrdC*CR - dKfdC*CF and sweep the whole species block in
        // 4-wide SIMD columns (lhs rows get -WdMdC, rhs rows +WdMdC).
        if(i>=this->Ikf[2] && i <this->Ikf[6])
        {
            const unsigned int kk = i - this->Ikf[2];

            __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[kk]);
            __m256d CF_ = _mm256_set1_pd(CF);
            __m256d CR_ = _mm256_set1_pd(CR);
            double dKrdC = 0;
            dKrdC = (this->dKfdC_[kk]*invKc);
            __m256d dKrdC_ = _mm256_set1_pd(dKrdC);

            for (unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[kk*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                for (unsigned int k = 0; k < this->lhsSpeciesIndex[i].size();k++)
                {
                    const unsigned int sk = this->lhsSpeciesIndex[i][k];
                    __m256d skv = load256d(&ddNdtByVdcTp[sk*(this->alignN)+j]);
                    skv = _mm256_sub_pd(skv,WdMdC);
                    store256d(&ddNdtByVdcTp[sk*(this->alignN)+j],skv);
                }

                for (unsigned int k = 0; k < rhsSpeciesIndex[i].size();k++)
                {

                    const unsigned int sk = rhsSpeciesIndex[i][k];
                    __m256d skv = load256d(&ddNdtByVdcTp[sk*(this->alignN)+j]);
                    skv = _mm256_add_pd(skv,WdMdC);
                    store256d(&ddNdtByVdcTp[sk*(this->alignN)+j],skv);
                }
            }
        }
    }
}
