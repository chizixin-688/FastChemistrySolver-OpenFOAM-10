/*---------------------------------------------------------------------------*\
  Description
      Global Jacobian update of the per-volume molar production rate for a
      reaction with non-integer stoichiometric coefficients (e.g.
      4.5A+5.5B=6.5C).

      For every reaction handled here the routine
        - assembles the forward/reverse rate coefficients (equilibrium
          constant from the exponent-weighted Gibbs factors negGstdByRT),
        - forms the net rate        q    = Kf*CF - Kr*CR
          and its temperature part  dwdT = dKfdT*CF - dKrdT*CR
        - scatters q / dwdT into the species rates dNdtByV and into the
          temperature column ddNdtByVdcTp[*][nSpecies] of the Jacobian, and
        - scatters the concentration derivatives of Kf*CF and Kr*CR into the
          species-concentration columns of ddNdtByVdcTp.

      Unlike the integer variant, each concentration enters with the
      reaction order given by lhsReactionOrder / rhsReactionOrder and is
      guarded by FastChemistry::ConcentrationLimiter.  Third-body reactions
      finally sweep the whole species block in 4-wide SIMD columns.

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"

//=============================================================================//

void  FastChemistry::OptReaction::JFGNI
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
) const noexcept
{

    std::size_t end = reactionGNIindex.size();

    // Sweep over the reactions that require the global Jacobian update.
    for (std::size_t K=0; K<end; K++)
    {
        // Rate data of the current reaction.  CF and CR below accumulate
        // the concentration products of the reactants / products.
        const unsigned int iii = this->reactionGNIindex[K];
        const double Kf = this->Kf_[iii];
        const double dKfdT = this->dKfdT_[iii];
        double Kr = 0;
        double dKrdT = 0;
        double CF = 1.0;
        double CR = 1.0;
        double invKc = 0;
        // Reversible reaction (code 0): assemble Kc from the species Gibbs
        // factors, then Kr = Kf/Kc and its temperature derivative.
        if(this->isIrreversible[iii]==0)
        {
            double sumVki = 0;
            double Kp = 1.0;
            double sumVdBdT = 0.0;
            for (unsigned int j = 0; j < this->lhsSpeciesIndex[iii].size();j++)
            {
                const unsigned int si = this->lhsSpeciesIndex[iii][j];
                const double sl = this->lhsStoichCoeff[iii][j];
                const double el = this->lhsReactionOrder[iii][j];
                Kp += sl*this->negGstdByRT[si];
                sumVdBdT = sumVdBdT - sl*dBdT[si];
                sumVki = sumVki - sl;
                CF = CF * (C[si] >= FastChemistry::ConcentrationLimiter || el >= 1 ? std::pow((C[si]), el) : 0.0);
            }

            for (unsigned int j = 0; j < this->rhsSpeciesIndex[iii].size();j++)
            {
                const unsigned int si = this->rhsSpeciesIndex[iii][j];
                const double sr = this->rhsStoichCoeff[iii][j];
                const double er = this->rhsReactionOrder[iii][j];
                Kp -= sr*this->negGstdByRT[si];
                sumVdBdT = sumVdBdT + sr*dBdT[si];
                sumVki = sumVki + sr;
                CR = CR * (C[si] >= FastChemistry::ConcentrationLimiter || er >= 1.0 ? std::pow((C[si]), er) : 0.0);
            }
            Kp = std::exp(Kp);
            double Kc = Kp*std::pow(this->Pstd/(this->Ru*this->T),sumVki);
            Kc = std::max(Kc,FastChemistry::KcLimiter);
            invKc = 1.0/Kc;

            const double dKcdTByKc = sumVdBdT - sumVki*invT;

            Kr = Kf*invKc;
            dKrdT = dKfdT*invKc - (Kc > FastChemistry::KcLimiter ? Kr*dKcdTByKc : 0.0);
        }
        // Code 2: the reverse rate is the stored forward rate of the
        // paired reverse reaction; only the concentration products remain.
        else if(this->isIrreversible[iii]==2)
        {
            Kr = this->Kf_[iii - Ikf[1] + Ikf[9]];
            dKrdT = this->dKfdT_[iii - Ikf[1] + Ikf[9]];
            invKc = Kr/Kf;
            for (unsigned int j = 0; j < this->rhsSpeciesIndex[iii].size();j++)
            {
                const unsigned int si = this->rhsSpeciesIndex[iii][j];
                const double er = lhsReactionOrder[iii][j];
                CR = CR * (C[si] >= FastChemistry::ConcentrationLimiter || er >= 1 ? std::pow((C[si]), er) : 0.0);
            }
            for (unsigned int j = 0; j < lhsSpeciesIndex[iii].size();j++)
            {
                const unsigned int si = lhsSpeciesIndex[iii][j];
                const double el = lhsReactionOrder[iii][j];
                CF = CF * (C[si] >= FastChemistry::ConcentrationLimiter || el >= 1 ? std::pow((C[si]), el) : 0.0);
            }
        }
        else
        {
            for (unsigned int j = 0; j < this->rhsSpeciesIndex[iii].size();j++)
            {
                const unsigned int si = this->rhsSpeciesIndex[iii][j];
                const double er = lhsReactionOrder[iii][j];
                CR = CR * (C[si] >= FastChemistry::ConcentrationLimiter || er >= 1 ? std::pow((C[si]), er) : 0.0);
            }
            for (unsigned int j = 0; j < lhsSpeciesIndex[iii].size();j++)
            {
                const unsigned int si = lhsSpeciesIndex[iii][j];
                const double el = lhsReactionOrder[iii][j];
                CF = CF * (C[si] >= FastChemistry::ConcentrationLimiter || el >= 1 ? std::pow((C[si]), el) : 0.0);
            }
        }

        // Net molar production rate of this reaction and its temperature
        // derivative (sign conventions: reactants consume, products form).
        const double q = Kf*CF - Kr*CR;
        // Scatter q into the species rates dNdtByV and dwdT into the
        // temperature column ddNdtByVdcTp[*][nSpecies] of the Jacobian,
        // each weighted by the (possibly non-integer) stoichiometric
        // coefficients sl / sr.
        const double dwdT = (dKfdT*CF-dKrdT*CR);

        for (unsigned int j = 0; j < this->lhsSpeciesIndex[iii].size();j++)
        {
            const unsigned int si = this->lhsSpeciesIndex[iii][j];
            const double sl = this->lhsStoichCoeff[iii][j];
            ddNdtByVdcTp[si*(this->alignN) + (this->nSpecies)] -= sl*dwdT;
            dNdtByV[si] = dNdtByV[si] - sl*q;
        }

        for (unsigned int j = 0; j < this->rhsSpeciesIndex[iii].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[iii][j];
            const double sr = this->rhsStoichCoeff[iii][j];
            ddNdtByVdcTp[si*(this->alignN) + (this->nSpecies)] += sr*dwdT;
            dNdtByV[si] = dNdtByV[si] + sr*q;
        }

        for (unsigned int j = 0; j < lhsSpeciesIndex[iii].size();j++)
        {
            const unsigned int sj = lhsSpeciesIndex[iii][j];

            // Derivative of the forward flux Kf*CF with respect to the
            // reactant concentration C_j.  Each factor is differentiated with
            // its reaction order and guarded by the concentration limiter.
            double dCfdCj = 1.0;
            for (unsigned int i = 0; i < lhsSpeciesIndex[iii].size();i++)
            {
                const unsigned int si = lhsSpeciesIndex[iii][i];
                const double el = lhsReactionOrder[iii][i];

                if (i == j)
                {
                    dCfdCj *=
                    C[si] >= FastChemistry::ConcentrationLimiter || el >= 1
                        ? el*std::pow((C[si]), el - 1.0)
                        : 0;
                }
                else
                {
                    dCfdCj *=
                        C[si] >= FastChemistry::ConcentrationLimiter || el >= 1
                        ? std::pow((C[si]), el)
                        : 0;
                }
            }

            const double KfdCfdCj = Kf*dCfdCj;

            for (unsigned int i = 0; i < lhsSpeciesIndex[iii].size();i++)
            {
                const unsigned int si = lhsSpeciesIndex[iii][i];
                const double sl = lhsStoichCoeff[iii][i];
                ddNdtByVdcTp[si*(this->alignN)+sj] -= sl*KfdCfdCj;
            }

            for (unsigned int i = 0; i < rhsSpeciesIndex[iii].size();i++)
            {
                const unsigned int si = rhsSpeciesIndex[iii][i];
                const double sr = rhsStoichCoeff[iii][i];
                ddNdtByVdcTp[si*(this->alignN)+sj] += sr*KfdCfdCj;
            }
        }

        // Derivative of the reverse flux Kr*CR with respect to the product
        // concentration C_j (reversible reactions and code 2 only).
        if(this->isIrreversible[iii]==0||this->isIrreversible[iii]==2)
        {
            for (unsigned int j = 0; j < rhsSpeciesIndex[iii].size();j++)
            {
                const unsigned int sj = rhsSpeciesIndex[iii][j];

                // Product rule over the products; each factor is differentiated
                // with its reaction order and guarded by the concentration limiter.
                double dCrdCj = 1.0;
                for (unsigned int i = 0; i < rhsSpeciesIndex[iii].size();i++)
                {

                    const unsigned int si = rhsSpeciesIndex[iii][i];
                    const double er = rhsReactionOrder[iii][i];
                    dCrdCj = dCrdCj*((j==i)? 1:C[si]);

                    if (i == j)
                    {
                        dCrdCj *=
                            C[si] >= FastChemistry::ConcentrationLimiter || er >= 1.0
                        ? er*std::pow((C[si]), er - 1.0)
                        : 0;
                    }
                    else
                    {
                        dCrdCj *=
                            C[si] >= FastChemistry::ConcentrationLimiter || er >= 1.0
                        ? std::pow((C[si]), er)
                        : 0;
                    }
                }

                const double negKrdCrdCj = -Kr*dCrdCj;

                for (unsigned int i = 0; i < lhsSpeciesIndex[iii].size();i++)
                {
                    const unsigned int si = lhsSpeciesIndex[iii][i];
                    const double sl = lhsReactionOrder[iii][j];
                    ddNdtByVdcTp[si*(this->alignN)+sj] -= sl*negKrdCrdCj;
                }

                for (unsigned int i = 0; i < rhsSpeciesIndex[iii].size();i++)
                {
                    const unsigned int si = rhsSpeciesIndex[iii][i];
                    const double sr = rhsReactionOrder[iii][j];
                    ddNdtByVdcTp[si*(this->alignN)+sj] += sr*negKrdCrdCj;
                }
            }
        }

        // Third-body (M) reactions: the rate coefficients depend on the
        // collision-partner concentrations.  Accumulate the M-derivative
        // dWdM = dKrdC*CR - dKfdC*CF and sweep the whole species block in
        // 4-wide SIMD columns (lhs rows get -WdMdC, rhs rows +WdMdC).
        if(iii>=this->Ikf[2] && iii <this->Ikf[6])
        {
            const unsigned int kk = iii - this->Ikf[2];

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

                for (unsigned int i = 0; i < lhsSpeciesIndex[iii].size();i++)
                {
                    const unsigned int si = lhsSpeciesIndex[iii][i];
                    const double sl = lhsStoichCoeff[iii][i];
                    __m256d Jac = load256d(&ddNdtByVdcTp[si*(this->alignN)+j]);
                    __m256d slv = _mm256_set1_pd(sl);

                    Jac = _mm256_sub_pd(Jac,_mm256_mul_pd(WdMdC,slv));
                    store256d(&ddNdtByVdcTp[si*(this->alignN)+j],Jac);
                }

                for (unsigned int i = 0; i < rhsSpeciesIndex[iii].size();i++)
                {
                    const unsigned int si = rhsSpeciesIndex[iii][i];
                    const double sr = rhsStoichCoeff[iii][i];
                    __m256d Jac = load256d(&ddNdtByVdcTp[si*(this->alignN)+j]);
                    __m256d srv = _mm256_set1_pd(sr);

                    Jac = _mm256_add_pd(Jac,_mm256_mul_pd(WdMdC,srv));
                    store256d(&ddNdtByVdcTp[si*(this->alignN)+j],Jac);
                }
            }
        }

    }
}
