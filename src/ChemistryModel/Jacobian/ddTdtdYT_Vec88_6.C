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

      This variant (suffix 6) is specialised for nSpecie() % 8 == 6: the
      8-wide SIMD blocks cover Ns-6 species and the remaining six columns
      and rows ("tail") are handled by a 4-wide block plus a 2-lane block
      (low 128-bit lanes).

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
void Foam::FastChemistryModel<UnusedThermo>::ddTdtdYT_Vec88_6
(
    const double* __restrict__ Cp,
    const double* __restrict__ dCpdT,
    const double* __restrict__ Ha,
    double* __restrict__ dPhidt,
    double* __restrict__ Jac
) const noexcept
{
    const double& CpM = Cp[this->nSpecie()];
    const double& dCpMdT = dCpdT[this->nSpecie()];
    const double invCpM = 1.0/CpM;
    // ---- dT/dt = -sum_i Ha_i*dY_i/dt / CpM ----
    __m256d dTdtv = _mm256_setzero_pd();
    for (int i=0; i<this->nSpecie()-6; i=i+8)
    {
        dTdtv = _mm256_fmadd_pd(-load256d(&Ha[i+0]),load256d(&dPhidt[i+0]),dTdtv);
        dTdtv = _mm256_fmadd_pd(-load256d(&Ha[i+4]),load256d(&dPhidt[i+4]),dTdtv);
    }
    {
        int i = this->nSpecie()-6;
        __m256d negHai04v = -load256d(&Ha[i+0]);
        __m256d dPhi04dtv = load256d(&dPhidt[i+0]);
        dTdtv = _mm256_fmadd_pd(negHai04v,dPhi04dtv,dTdtv);

        __m256d negHai45v = -_mm256_insertf128_pd (_mm256_setzero_pd(), _mm_loadu_pd(&Ha[i+4]), 0);
        __m256d dPhi45dtv = _mm256_insertf128_pd (_mm256_setzero_pd (), _mm_loadu_pd(&dPhidt[i+4]), 0);

        dTdtv = _mm256_fmadd_pd(negHai45v,dPhi45dtv,dTdtv);
    }
    double dTdt = hsum4(dTdtv);

    dTdt *= invCpM;

    __m256d invCpMv  = _mm256_set1_pd(invCpM);   // 1/CpM broadcast
    __m256d negdTdtv = _mm256_set1_pd(-dTdt);    // -dTdt broadcast

    dPhidt[this->nSpecie()] = dTdt;
    double& ddTdtdT = Jac[this->nSpecie() *(alignN)+ this->nSpecie()];
    ddTdtdT = 0;
    __m256d ddTdtdTv = _mm256_setzero_pd();

    for (int i=0; i<this->nSpecie()-6; i=i+8)
    {
        __m256d ddTdtdYi03v = _mm256_setzero_pd();
        __m256d ddTdtdYi47v = _mm256_setzero_pd();
        for (int j=0; j<this->nSpecie()-6; j=j+8)
        {
            __m256d negHaj0v = _mm256_set1_pd(-Ha[j+0]);
            __m256d ddYj0dtdYi03v = load256d(&Jac[(j+0)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj0dtdYi03v,negHaj0v,ddTdtdYi03v);
            __m256d ddYj0dtdYi47v = load256d(&Jac[(j+0)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj0dtdYi47v,negHaj0v,ddTdtdYi47v);

            __m256d negHaj1v = _mm256_set1_pd(-Ha[j+1]);
            __m256d ddYj1dtdYi03v = load256d(&Jac[(j+1)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj1dtdYi03v,negHaj1v,ddTdtdYi03v);
            __m256d ddYj1dtdYi47v = load256d(&Jac[(j+1)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj1dtdYi47v,negHaj1v,ddTdtdYi47v);

            __m256d negHaj2v = _mm256_set1_pd(-Ha[j+2]);
            __m256d ddYj2dtdYi03v = load256d(&Jac[(j+2)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj2dtdYi03v,negHaj2v,ddTdtdYi03v);
            __m256d ddYj2dtdYi47v = load256d(&Jac[(j+2)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj2dtdYi47v,negHaj2v,ddTdtdYi47v);

            __m256d negHaj3v = _mm256_set1_pd(-Ha[j+3]);
            __m256d ddYj3dtdYi03v = load256d(&Jac[(j+3)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj3dtdYi03v,negHaj3v,ddTdtdYi03v);
            __m256d ddYj3dtdYi47v = load256d(&Jac[(j+3)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj3dtdYi47v,negHaj3v,ddTdtdYi47v);

            __m256d negHaj4v = _mm256_set1_pd(-Ha[j+4]);
            __m256d ddYj4dtdYi03v = load256d(&Jac[(j+4)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj4dtdYi03v,negHaj4v,ddTdtdYi03v);
            __m256d ddYj4dtdYi47v = load256d(&Jac[(j+4)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj4dtdYi47v,negHaj4v,ddTdtdYi47v);

            __m256d negHaj5v = _mm256_set1_pd(-Ha[j+5]);
            __m256d ddYj5dtdYi03v = load256d(&Jac[(j+5)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj5dtdYi03v,negHaj5v,ddTdtdYi03v);
            __m256d ddYj5dtdYi47v = load256d(&Jac[(j+5)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj5dtdYi47v,negHaj5v,ddTdtdYi47v);

            __m256d negHaj6v = _mm256_set1_pd(-Ha[j+6]);
            __m256d ddYj6dtdYi03v = load256d(&Jac[(j+6)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj6dtdYi03v,negHaj6v,ddTdtdYi03v);
            __m256d ddYj6dtdYi47v = load256d(&Jac[(j+6)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj6dtdYi47v,negHaj6v,ddTdtdYi47v);

            __m256d negHaj7v = _mm256_set1_pd(-Ha[j+7]);
            __m256d ddYj7dtdYi03v = load256d(&Jac[(j+7)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj7dtdYi03v,negHaj7v,ddTdtdYi03v);
            __m256d ddYj7dtdYi47v = load256d(&Jac[(j+7)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj7dtdYi47v,negHaj7v,ddTdtdYi47v);
        }
        {
            int j = this->nSpecie()-6;

            __m256d negHaj0v = _mm256_set1_pd(-Ha[j+0]);
            __m256d ddYj0dtdYi03v = load256d(&Jac[(j+0)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj0dtdYi03v,negHaj0v,ddTdtdYi03v);
            __m256d ddYj0dtdYi47v = load256d(&Jac[(j+0)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj0dtdYi47v,negHaj0v,ddTdtdYi47v);

            __m256d negHaj1v = _mm256_set1_pd(-Ha[j+1]);
            __m256d ddYj1dtdYi03v = load256d(&Jac[(j+1)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj1dtdYi03v,negHaj1v,ddTdtdYi03v);
            __m256d ddYj1dtdYi47v = load256d(&Jac[(j+1)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj1dtdYi47v,negHaj1v,ddTdtdYi47v);

            __m256d negHaj2v = _mm256_set1_pd(-Ha[j+2]);
            __m256d ddYj2dtdYi03v = load256d(&Jac[(j+2)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj2dtdYi03v,negHaj2v,ddTdtdYi03v);
            __m256d ddYj2dtdYi47v = load256d(&Jac[(j+2)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj2dtdYi47v,negHaj2v,ddTdtdYi47v);

            __m256d negHaj3v = _mm256_set1_pd(-Ha[j+3]);
            __m256d ddYj3dtdYi03v = load256d(&Jac[(j+3)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj3dtdYi03v,negHaj3v,ddTdtdYi03v);
            __m256d ddYj3dtdYi47v = load256d(&Jac[(j+3)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj3dtdYi47v,negHaj3v,ddTdtdYi47v);

            __m256d negHaj4v = _mm256_set1_pd(-Ha[j+4]);
            __m256d ddYj4dtdYi03v = load256d(&Jac[(j+4) *(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj4dtdYi03v,negHaj4v,ddTdtdYi03v);
            __m256d ddYj4dtdYi47v = load256d(&Jac[(j+4) *(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj4dtdYi47v,negHaj4v,ddTdtdYi47v);

            __m256d negHaj5v = _mm256_set1_pd(-Ha[j+5]);
            __m256d ddYj5dtdYi03v = load256d(&Jac[(j+5)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj5dtdYi03v,negHaj5v,ddTdtdYi03v);
            __m256d ddYj5dtdYi47v = load256d(&Jac[(j+5)*(alignN)+ (i+4)]);
            ddTdtdYi47v = _mm256_fmadd_pd(ddYj5dtdYi47v,negHaj5v,ddTdtdYi47v);
        }

        ddTdtdYi03v = _mm256_fmadd_pd(load256d(&Cp[i+0]),negdTdtv,ddTdtdYi03v);
        ddTdtdYi03v =_mm256_mul_pd(ddTdtdYi03v,invCpMv);
        store256d(&Jac[this->nSpecie() *(alignN)+ i+0],ddTdtdYi03v);
        __m256d dYi03dtv = load256d(&dPhidt[i+0]);
        const double ddYi0dtdT = Jac[(i+0) *(alignN)+ this->nSpecie()];
        const double ddYi1dtdT = Jac[(i+1) *(alignN)+ this->nSpecie()];
        const double ddYi2dtdT = Jac[(i+2) *(alignN)+ this->nSpecie()];
        const double ddYi3dtdT = Jac[(i+3) *(alignN)+ this->nSpecie()];
        __m256d ddYi03dtdTv = _mm256_setr_pd(ddYi0dtdT,ddYi1dtdT,ddYi2dtdT,ddYi3dtdT);
        ddTdtdTv = _mm256_fmadd_pd(dYi03dtv,load256d(&Cp[i+0]),ddTdtdTv);
        ddTdtdTv = _mm256_fmadd_pd(ddYi03dtdTv,load256d(&Ha[i+0]),ddTdtdTv);


        ddTdtdYi47v = _mm256_fmadd_pd(load256d(&Cp[i+4]),negdTdtv,ddTdtdYi47v);
        ddTdtdYi47v =_mm256_mul_pd(ddTdtdYi47v,invCpMv);
        store256d(&Jac[this->nSpecie() *(alignN)+ i+4],ddTdtdYi47v);
        __m256d dYi47dtv = load256d(&dPhidt[i+4]);
        const double ddYi4dtdT = Jac[(i+4) *(alignN)+ this->nSpecie()];
        const double ddYi5dtdT = Jac[(i+5) *(alignN)+ this->nSpecie()];
        const double ddYi6dtdT = Jac[(i+6) *(alignN)+ this->nSpecie()];
        const double ddYi7dtdT = Jac[(i+7) *(alignN)+ this->nSpecie()];
        __m256d ddYi47dtdTv = _mm256_setr_pd(ddYi4dtdT,ddYi5dtdT,ddYi6dtdT,ddYi7dtdT);
        ddTdtdTv = _mm256_fmadd_pd(dYi47dtv,load256d(&Cp[i+4]),ddTdtdTv);
        ddTdtdTv = _mm256_fmadd_pd(ddYi47dtdTv,load256d(&Ha[i+4]),ddTdtdTv);
    }
    // ---- tail: last six columns (species Ns-6 .. Ns-1) ----
    {
        int i = this->nSpecie()-6;

        __m256d ddTdtdYi03v = _mm256_setzero_pd();    // columns i..i+3
        __m256d ddTdtdYi45v = _mm256_setzero_pd();    // columns i+4..i+5 (low 2 lanes)
        for (int j=0; j<this->nSpecie()-6; j=j+8)
        {
            __m256d negHaj0v = _mm256_set1_pd(-Ha[j+0]);
            __m256d ddYj0dtdYi03v = load256d(&Jac[(j+0)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj0dtdYi03v,negHaj0v,ddTdtdYi03v);
            __m256d ddYj0dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+0)*(alignN)+(i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj0dtdYi45v,negHaj0v,ddTdtdYi45v);

            __m256d negHaj1v = _mm256_set1_pd(-Ha[j+1]);
            __m256d ddYj1dtdYi03v = load256d(&Jac[(j+1)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj1dtdYi03v,negHaj1v,ddTdtdYi03v);
            __m256d ddYj1dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+1)*(alignN)+(i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj1dtdYi45v,negHaj1v,ddTdtdYi45v);

            __m256d negHaj2v = _mm256_set1_pd(-Ha[j+2]);
            __m256d ddYj2dtdYi03v = load256d(&Jac[(j+2)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj2dtdYi03v,negHaj2v,ddTdtdYi03v);
            __m256d ddYj2dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+2)*(alignN)+(i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj2dtdYi45v,negHaj2v,ddTdtdYi45v);

            __m256d negHaj3v = _mm256_set1_pd(-Ha[j+3]);
            __m256d ddYj3dtdYi03v = load256d(&Jac[(j+3)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj3dtdYi03v,negHaj3v,ddTdtdYi03v);
            __m256d ddYj3dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+3)*(alignN)+(i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj3dtdYi45v,negHaj3v,ddTdtdYi45v);

            __m256d negHaj4v = _mm256_set1_pd(-Ha[j+4]);
            __m256d ddYj4dtdYi03v = load256d(&Jac[(j+4)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj4dtdYi03v,negHaj4v,ddTdtdYi03v);
            __m256d ddYj4dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+4)*(alignN)+(i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj4dtdYi45v,negHaj4v,ddTdtdYi45v);

            __m256d negHaj5v = _mm256_set1_pd(-Ha[j+5]);
            __m256d ddYj5dtdYi03v = load256d(&Jac[(j+5)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj5dtdYi03v,negHaj5v,ddTdtdYi03v);
            __m256d ddYj5dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+5)*(alignN)+ (i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj5dtdYi45v,negHaj5v,ddTdtdYi45v);

            __m256d negHaj6v = _mm256_set1_pd(-Ha[j+6]);
            __m256d ddYj6dtdYi03v = load256d(&Jac[(j+6)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj6dtdYi03v,negHaj6v,ddTdtdYi03v);
            __m256d ddYj6dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+6)*(alignN)+ (i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj6dtdYi45v,negHaj6v,ddTdtdYi45v);

            __m256d negHaj7v = _mm256_set1_pd(-Ha[j+7]);
            __m256d ddYj7dtdYi03v = load256d(&Jac[(j+7)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj7dtdYi03v,negHaj7v,ddTdtdYi03v);
            __m256d ddYj7dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+7)*(alignN)+ (i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj7dtdYi45v,negHaj7v,ddTdtdYi45v);
        }
        {
            int j = this->nSpecie()-6;

            __m256d negHaj0v = _mm256_set1_pd(-Ha[j+0]);
            __m256d ddYj0dtdYi03v = load256d(&Jac[(j+0)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj0dtdYi03v,negHaj0v,ddTdtdYi03v);
            __m256d ddYj0dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+0)*(alignN)+(i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj0dtdYi45v,negHaj0v,ddTdtdYi45v);

            __m256d negHaj1v = _mm256_set1_pd(-Ha[j+1]);
            __m256d ddYj1dtdYi03v = load256d(&Jac[(j+1)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj1dtdYi03v,negHaj1v,ddTdtdYi03v);
            __m256d ddYj1dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+1)*(alignN)+(i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj1dtdYi45v,negHaj1v,ddTdtdYi45v);

            __m256d negHaj2v = _mm256_set1_pd(-Ha[j+2]);
            __m256d ddYj2dtdYi03v = load256d(&Jac[(j+2)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj2dtdYi03v,negHaj2v,ddTdtdYi03v);
            __m256d ddYj2dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+2)*(alignN)+(i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj2dtdYi45v,negHaj2v,ddTdtdYi45v);

            __m256d negHaj3v = _mm256_set1_pd(-Ha[j+3]);
            __m256d ddYj3dtdYi03v = load256d(&Jac[(j+3)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj3dtdYi03v,negHaj3v,ddTdtdYi03v);
            __m256d ddYj3dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+3)*(alignN)+(i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj3dtdYi45v,negHaj3v,ddTdtdYi45v);

            __m256d negHaj4v = _mm256_set1_pd(-Ha[j+4]);
            __m256d ddYj4dtdYi03v = load256d(&Jac[(j+4)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj4dtdYi03v,negHaj4v,ddTdtdYi03v);
            __m256d ddYj4dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+4)*(alignN)+(i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj4dtdYi45v,negHaj4v,ddTdtdYi45v);

            __m256d negHaj5v = _mm256_set1_pd(-Ha[j+5]);
            __m256d ddYj5dtdYi03v = load256d(&Jac[(j+5)*(alignN)+ (i+0)]);
            ddTdtdYi03v = _mm256_fmadd_pd(ddYj5dtdYi03v,negHaj5v,ddTdtdYi03v);
            __m256d ddYj5dtdYi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Jac[(j+5)*(alignN)+ (i+4)]), 0);
            ddTdtdYi45v = _mm256_fmadd_pd(ddYj5dtdYi45v,negHaj5v,ddTdtdYi45v);
        }
        __m256d Cpi03v = load256d(&Cp[i+0]);
        ddTdtdYi03v = _mm256_fmadd_pd(Cpi03v,negdTdtv,ddTdtdYi03v);
        ddTdtdYi03v =_mm256_mul_pd(ddTdtdYi03v,invCpMv);
        store256d(&Jac[this->nSpecie() *(alignN)+ i+0],ddTdtdYi03v);
        __m256d dYi03dtv = load256d(&dPhidt[i+0]);
        const double ddYi0dtdT = Jac[(i+0) *(alignN)+ this->nSpecie()];
        const double ddYi1dtdT = Jac[(i+1) *(alignN)+ this->nSpecie()];
        const double ddYi2dtdT = Jac[(i+2) *(alignN)+ this->nSpecie()];
        const double ddYi3dtdT = Jac[(i+3) *(alignN)+ this->nSpecie()];
        __m256d ddYi03dtdTv = _mm256_setr_pd(ddYi0dtdT,ddYi1dtdT,ddYi2dtdT,ddYi3dtdT);
        ddTdtdTv = _mm256_fmadd_pd(dYi03dtv,Cpi03v,ddTdtdTv);
        ddTdtdTv = _mm256_fmadd_pd(ddYi03dtdTv,load256d(&Ha[i+0]),ddTdtdTv);

        __m256d Cpi45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Cp[i+4]), 0);
        ddTdtdYi45v = _mm256_fmadd_pd(Cpi45v,negdTdtv,ddTdtdYi45v);
        ddTdtdYi45v =_mm256_mul_pd(ddTdtdYi45v,invCpMv);
        store256d(&Jac[this->nSpecie() *(alignN)+ i+4],ddTdtdYi45v);
        __m256d dYi45dtv = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&dPhidt[i+4]), 0);
        const double ddYi4dtdT = Jac[(i+4) *(alignN)+ this->nSpecie()];
        const double ddYi5dtdT = Jac[(i+5) *(alignN)+ this->nSpecie()];
        __m256d ddYi45dtdTv = _mm256_setr_pd(ddYi4dtdT,ddYi5dtdT,0,0);
        ddTdtdTv = _mm256_fmadd_pd(dYi45dtv,Cpi45v,ddTdtdTv);
        __m256d Hai45v = _mm256_insertf128_pd(_mm256_setzero_pd (), _mm_loadu_pd(&Ha[i+4]), 0);
        ddTdtdTv = _mm256_fmadd_pd(ddYi45dtdTv,Hai45v,ddTdtdTv);
    }
    // d(dTdt)/dT = -( numerator + dCpMdT*dTdt ) / CpM
    ddTdtdT = ddTdtdT - hsum4(ddTdtdTv);
    ddTdtdT -= dTdt*dCpMdT;
    ddTdtdT *= invCpM;
}