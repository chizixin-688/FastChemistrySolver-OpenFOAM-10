/*---------------------------------------------------------------------------*\
  Description
      Mass-fraction-based Jacobian: temperature row of d(dT/dt)/d(Y,T).

      The species Jacobian block d(dY_j/dt)/dY_i and the temperature column
      d(dY_j/dt)/dT are assumed to be already assembled in Jac.  This routine
      fills the last row (row Ns, Ns = number of species):

          Jac[Ns][i] = -[ sum_j Ha_j*Jac[j][i] + Cp_i*dTdt ] / CpM
          Jac[Ns][Ns] = -[ sum_i (Cp_i*dY_i/dt + Ha_i*d(dY_i/dt)/dT)
                            + dCpMdT*dTdt ] / CpM

      with  dTdt = -sum_j Ha_j*dY_j/dt / CpM,  CpM = Cp[Ns] the mixture heat
      capacity and dCpMdT = dCpdT[Ns] its temperature derivative.  dTdt is
      also stored in dPhidt[Ns].

      This variant (suffix 0) is specialised for nSpecie() % 8 == 0 and
      therefore has no remainder tail.

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "FastChemistryModel.H"

//---------------------------------
// 2. SIMD / AVX2 headers
//---------------------------------
#include <immintrin.h>

//=============================================================================//

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddTdtdYT_Vec88_0
(
    const double* __restrict__ Cp,
    const double* __restrict__ dCpdT,
    const double* __restrict__ Ha,
    double* __restrict__ dPhidt,
    double* __restrict__ Jac
) const noexcept
{
    // Mixture total heat capacity (stored at index nSpecie) and its derivative.
    const double& CpM    = Cp[this->nSpecie()];
    const double& dCpMdT = dCpdT[this->nSpecie()];
    const double invCpM  = 1.0/CpM;

    // ---- dT/dt = -sum_i Ha_i*dY_i/dt / CpM ----
    __m256d dTdtv = _mm256_setzero_pd();
    for (int i=0; i<this->nSpecie(); i=i+8)
    {
        dTdtv = _mm256_fmadd_pd(-load256d(&Ha[i+0]),load256d(&dPhidt[i+0]),dTdtv);
        dTdtv = _mm256_fmadd_pd(-load256d(&Ha[i+4]),load256d(&dPhidt[i+4]),dTdtv);
    }
    const double dTdt = hsum4(dTdtv)*invCpM;
    dPhidt[this->nSpecie()] = dTdt;

    // ---- temperature-row accumulator and diagonal entry (row/col index nSpecie) ----
    double& ddTdtdT = Jac[this->nSpecie()*(alignN)+this->nSpecie()];
    ddTdtdT = 0;

    const __m256d invCpMv  = _mm256_set1_pd(invCpM);
    const __m256d negdTdtv = _mm256_set1_pd(-dTdt);
    __m256d ddTdtdTv = _mm256_setzero_pd();   // d(dTdt)/dT numerator accumulator
    for (int i=0; i<this->nSpecie(); i=i+8)
    {
        __m256d ddTdtdYi03v = _mm256_setzero_pd();  // row nSpecie, columns i..i+3
        __m256d ddTdtdYi47v = _mm256_setzero_pd();  // row nSpecie, columns i+4..i+7
        for (int j=0; j<this->nSpecie(); j=j+8)
        {
            // Accumulate -sum_j Ha_j*Jac[j][i] over species rows j..j+7.
            for (int r=0; r<8; r++)
            {
                const __m256d negHajv = _mm256_set1_pd(-Ha[j+r]);
                ddTdtdYi03v = _mm256_fmadd_pd(load256d(&Jac[(j+r)*(alignN)+(i+0)]),negHajv,ddTdtdYi03v);
                ddTdtdYi47v = _mm256_fmadd_pd(load256d(&Jac[(j+r)*(alignN)+(i+4)]),negHajv,ddTdtdYi47v);
            }
        }
        // Jac[Ns][i] = -( acc[i] + Cp_i*dTdt ) / CpM   for columns i..i+3
        const __m256d Cpi03v = load256d(&Cp[i+0]);
        ddTdtdYi03v = _mm256_mul_pd(_mm256_fmadd_pd(Cpi03v,negdTdtv,ddTdtdYi03v),invCpMv);
        store256d(&Jac[this->nSpecie()*(alignN)+i+0],ddTdtdYi03v);

        // ... and for columns i+4..i+7
        const __m256d Cpi47v = load256d(&Cp[i+4]);
        ddTdtdYi47v = _mm256_mul_pd(_mm256_fmadd_pd(Cpi47v,negdTdtv,ddTdtdYi47v),invCpMv);
        store256d(&Jac[this->nSpecie()*(alignN)+i+4],ddTdtdYi47v);
        // d(dTdt)/dT numerator terms for species i..i+7:
        //   sum_i [ Cp_i*dY_i/dt + Ha_i*d(dY_i/dt)/dT ]
        const __m256d dYi03dtv = load256d(&dPhidt[i+0]);
        const __m256d dYi47dtv = load256d(&dPhidt[i+4]);

        const __m256d ddYi03dtdTv = _mm256_setr_pd
        (
            Jac[(i+0)*(alignN)+this->nSpecie()],
            Jac[(i+1)*(alignN)+this->nSpecie()],
            Jac[(i+2)*(alignN)+this->nSpecie()],
            Jac[(i+3)*(alignN)+this->nSpecie()]
        );
        ddTdtdTv = _mm256_fmadd_pd(dYi03dtv,Cpi03v,ddTdtdTv);
        ddTdtdTv = _mm256_fmadd_pd(ddYi03dtdTv,load256d(&Ha[i+0]),ddTdtdTv);

        const __m256d ddYi47dtdTv = _mm256_setr_pd
        (
            Jac[(i+4)*(alignN)+this->nSpecie()],
            Jac[(i+5)*(alignN)+this->nSpecie()],
            Jac[(i+6)*(alignN)+this->nSpecie()],
            Jac[(i+7)*(alignN)+this->nSpecie()]
        );
        ddTdtdTv = _mm256_fmadd_pd(dYi47dtv,Cpi47v,ddTdtdTv);
        ddTdtdTv = _mm256_fmadd_pd(ddYi47dtdTv,load256d(&Ha[i+4]),ddTdtdTv);
    }

    // d(dTdt)/dT = -( numerator + dCpMdT*dTdt ) / CpM
    ddTdtdT = ddTdtdT - hsum4(ddTdtdTv);
    ddTdtdT -= dTdt*dCpMdT;
    ddTdtdT *= invCpM;
}