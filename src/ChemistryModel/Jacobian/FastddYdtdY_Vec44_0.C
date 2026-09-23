/*---------------------------------------------------------------------------*\
  Description
      Fast (approximate) species block of the mass-fraction based Jacobian:
      Jac[i][j] = d(dY_i/dt)/dY_j  for i,j in [0, Ns).

      The full conversion dc_k/dY_j is approximated by its diagonal part
      (the off-diagonal elements of the dcdY matrix are zero):

          dc_k/dY_j ~ delta_kj * rhoM*invW[j]

      which yields

          Jac[i][j] = WiByrhoM[i]*ddNdtByVdcT[i][j]*rhoM*invW[j]
                      + rhoMByRhoi[j]*(WiByrhoM[i]*dPhidt[i])

      where rhoM is the reference density.  Phi and dcdY are unused by this
      approximation (kept for interface consistency with the exact kernels).

This variant (suffix 0) is specialised for nSpecie() % 4 == 0, so no remainder tail is required.

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

template<class ThermoType>
void Foam::FastChemistryModel<ThermoType>::FastddYdtdY_Vec44_0
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ rhoMByRhoi,
    const double* __restrict__ WiByrhoM,
    const double* __restrict__ dPhidt,
    const double* __restrict__ invW,
    const double* __restrict__ Phi,
    double* __restrict__ dcdY,
    double* __restrict__ Jac,
    double rhoM
) const noexcept
{
    __m256d rhoMv = _mm256_set1_pd(rhoM);
    for(int i=0; i<this->nSpecie(); i=i+4)
    {
        __m256d Wi03ByrhoMv = load256d(&WiByrhoM[i+0]);
        __m256d dPhi03dtv = load256d(&dPhidt[i+0]);
        dPhi03dtv = _mm256_mul_pd(dPhi03dtv,Wi03ByrhoMv);

        __m256d Wi0ByrhoMv = _mm256_permute4x64_pd(Wi03ByrhoMv, 0x00);
        __m256d Wi1ByrhoMv = _mm256_permute4x64_pd(Wi03ByrhoMv, 0x55);
        __m256d Wi2ByrhoMv = _mm256_permute4x64_pd(Wi03ByrhoMv, 0xAA);
        __m256d Wi3ByrhoMv = _mm256_permute4x64_pd(Wi03ByrhoMv, 0xFF);

        __m256d dPhi0dtv = _mm256_permute4x64_pd(dPhi03dtv, 0x00);
        __m256d dPhi1dtv = _mm256_permute4x64_pd(dPhi03dtv, 0x55);
        __m256d dPhi2dtv = _mm256_permute4x64_pd(dPhi03dtv, 0xAA);
        __m256d dPhi3dtv = _mm256_permute4x64_pd(dPhi03dtv, 0xFF);

        for (int j=0; j<this->nSpecie(); j=j+4)
        {
            __m256d invWv = load256d(&invW[j+0]);
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,invWv);
            __m256d rhoMvj_ = load256d(&rhoMByRhoi[j+0]);

            __m256d Jci0j03v = load256d(&ddNdtByVdcT[(i+0)*(alignN)+j+0]);
            __m256d ddNi0dtByVdYj = _mm256_mul_pd(Jci0j03v,dCjdYj);
            __m256d r0 = _mm256_mul_pd(rhoMvj_,dPhi0dtv);
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,r0);
            store256d(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);

            __m256d Jci1j03v = load256d(&ddNdtByVdcT[(i+1)*(alignN)+j+0]);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(Jci1j03v,dCjdYj);
            __m256d r1 = _mm256_mul_pd(rhoMvj_,dPhi1dtv);
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,r1);
            store256d(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);

            __m256d Jci2j03v = load256d(&ddNdtByVdcT[(i+2)*(alignN)+j+0]);
            __m256d ddNi2dtByVdYj = _mm256_mul_pd(Jci2j03v,dCjdYj);
            __m256d r2 = _mm256_mul_pd(rhoMvj_,dPhi2dtv);
            __m256d ddYi2dtdYj = _mm256_fmadd_pd(Wi2ByrhoMv,ddNi2dtByVdYj,r2);
            store256d(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYj);

            __m256d Jci3j03v = load256d(&ddNdtByVdcT[(i+3)*(alignN)+j+0]);
            __m256d ddNi3dtByVdYj = _mm256_mul_pd(Jci3j03v,dCjdYj);
            __m256d r3 = _mm256_mul_pd(rhoMvj_,dPhi3dtv);
            __m256d ddYi3dtdYj = _mm256_fmadd_pd(Wi3ByrhoMv,ddNi3dtByVdYj,r3);
            store256d(&Jac[(i+3)*(alignN) + j+0],ddYi3dtdYj);
        }
    }
}