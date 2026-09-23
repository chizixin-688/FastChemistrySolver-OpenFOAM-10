/*---------------------------------------------------------------------------*\
  Description
      Temperature column of the mass-fraction based Jacobian (exact path):
      Jac[i][Ns] = d(dY_i/dt)/dT  for i in [0, Ns).

      The explicit temperature dependence ddNdtByVdcT[i][Ns] is combined
      with the implicit one through the concentrations (dc_k/dT = -alpha*c_k
      chain-rule term), giving

          dY_i/dt      = WiByrhoM[i]*dPhidt[i]
          Jac[i][Ns]   = WiByrhoM[i]*( ddNdtByVdcT[i][Ns]
                          - alphavM*sum_j ddNdtByVdcT[i][j]*c[j] )
                          + alphavM*dY_i/dt

      dPhidt is also overwritten in place with dY_i/dt.  The species block
      and the temperature row are assembled by the sibling kernels; this
      routine only fills column Ns.

This variant (suffix 0) is specialised for nSpecie() % 8 == 0, so no remainder tail is required.

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
void Foam::FastChemistryModel<UnusedThermo>::ddYdtdT_Vec88_0
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ WiByrhoM,
    const double* __restrict__ c,
    double* __restrict__ dPhidt,
    double* __restrict__ Jac,
    double alphavM
) const noexcept
{
    __m256d alphavMv = _mm256_set1_pd(alphavM);
    for (int i=0; i<this->nSpecie(); i=i+8)
    {

        __m256d Wi03ByrhoMv = load256d(&WiByrhoM[i]);
        __m256d dPhi03dtv = load256d(&dPhidt[i]);
        dPhi03dtv = _mm256_mul_pd(dPhi03dtv,Wi03ByrhoMv);
        store256d(&dPhidt[i],dPhi03dtv);

        __m256d Wi47ByrhoMv = load256d(&WiByrhoM[i+4]);
        __m256d dPhi47dtv = load256d(&dPhidt[i+4]);
        dPhi47dtv = _mm256_mul_pd(dPhi47dtv,Wi47ByrhoMv);
        store256d(&dPhidt[i+4],dPhi47dtv);

        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*(this->alignN)];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*(this->alignN)];
        const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*(this->alignN)];
        const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i+3)*(this->alignN)];
        const double* __restrict__ JcRowi4 = &ddNdtByVdcT[(i+4)*(this->alignN)];
        const double* __restrict__ JcRowi5 = &ddNdtByVdcT[(i+5)*(this->alignN)];
        const double* __restrict__ JcRowi6 = &ddNdtByVdcT[(i+6)*(this->alignN)];
        const double* __restrict__ JcRowi7 = &ddNdtByVdcT[(i+7)*(this->alignN)];

        double ddNi0dtByVdT = JcRowi0[this->nSpecie()];
        double ddNi1dtByVdT = JcRowi1[this->nSpecie()];
        double ddNi2dtByVdT = JcRowi2[this->nSpecie()];
        double ddNi3dtByVdT = JcRowi3[this->nSpecie()];
        double ddNi4dtByVdT = JcRowi4[this->nSpecie()];
        double ddNi5dtByVdT = JcRowi5[this->nSpecie()];
        double ddNi6dtByVdT = JcRowi6[this->nSpecie()];
        double ddNi7dtByVdT = JcRowi7[this->nSpecie()];
        __m256d sumi0 = _mm256_setzero_pd();
        __m256d sumi1 = _mm256_setzero_pd();
        __m256d sumi2 = _mm256_setzero_pd();
        __m256d sumi3 = _mm256_setzero_pd();
        __m256d sumi4 = _mm256_setzero_pd();
        __m256d sumi5 = _mm256_setzero_pd();
        __m256d sumi6 = _mm256_setzero_pd();
        __m256d sumi7 = _mm256_setzero_pd();
        for (int j=0; j<this->nSpecie(); j=j+8)
        {
            __m256d Cj03v = load256d(&c[j+0]);
            __m256d Cj47v = load256d(&c[j+4]);

            __m256d Jci0j03v = load256d(&JcRowi0[j]);
            sumi0 =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j03v,Cj03v),alphavMv,sumi0);
            __m256d Jci0j47v = load256d(&JcRowi0[j+4]);
            sumi0 =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j47v,Cj47v),alphavMv,sumi0);

            __m256d Jci1j03v = load256d(&JcRowi1[j]);
            sumi1 =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j03v,Cj03v),alphavMv,sumi1);
            __m256d Jci1j47v = load256d(&JcRowi1[j+4]);
            sumi1 =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j47v,Cj47v),alphavMv,sumi1);

            __m256d Jci2j03v = load256d(&JcRowi2[j]);
            sumi2 =_mm256_fmadd_pd(_mm256_mul_pd(Jci2j03v,Cj03v),alphavMv,sumi2);
            __m256d Jci2j47v = load256d(&JcRowi2[j+4]);
            sumi2 =_mm256_fmadd_pd(_mm256_mul_pd(Jci2j47v,Cj47v),alphavMv,sumi2);


            __m256d Jci3j03v = load256d(&JcRowi3[j]);
            sumi3 =_mm256_fmadd_pd(_mm256_mul_pd(Jci3j03v,Cj03v),alphavMv,sumi3);
            __m256d Jci3j47v = load256d(&JcRowi3[j+4]);
            sumi3 =_mm256_fmadd_pd(_mm256_mul_pd(Jci3j47v,Cj47v),alphavMv,sumi3);

            __m256d Jci4j03v = load256d(&JcRowi4[j]);
            sumi4 =_mm256_fmadd_pd(_mm256_mul_pd(Jci4j03v,Cj03v),alphavMv,sumi4);
            __m256d Jci4j47v = load256d(&JcRowi4[j+4]);
            sumi4 =_mm256_fmadd_pd(_mm256_mul_pd(Jci4j47v,Cj47v),alphavMv,sumi4);

            __m256d Jci5j03v = load256d(&JcRowi5[j]);
            sumi5 =_mm256_fmadd_pd(_mm256_mul_pd(Jci5j03v,Cj03v),alphavMv,sumi5);
            __m256d Jci5j47v = load256d(&JcRowi5[j+4]);
            sumi5 =_mm256_fmadd_pd(_mm256_mul_pd(Jci5j47v,Cj47v),alphavMv,sumi5);

            __m256d Jci6j03v = load256d(&JcRowi6[j]);
            sumi6 =_mm256_fmadd_pd(_mm256_mul_pd(Jci6j03v,Cj03v),alphavMv,sumi6);
            __m256d Jci6j47v = load256d(&JcRowi6[j+4]);
            sumi6 =_mm256_fmadd_pd(_mm256_mul_pd(Jci6j47v,Cj47v),alphavMv,sumi6);

            __m256d Jci7j03v = load256d(&JcRowi7[j]);
            sumi7 =_mm256_fmadd_pd(_mm256_mul_pd(Jci7j03v,Cj03v),alphavMv,sumi7);
            __m256d Jci7j47v = load256d(&JcRowi7[j+4]);
            sumi7 =_mm256_fmadd_pd(_mm256_mul_pd(Jci7j47v,Cj47v),alphavMv,sumi7);
        }


        ddNi0dtByVdT = ddNi0dtByVdT - (hsum4(sumi0));
        ddNi1dtByVdT = ddNi1dtByVdT - (hsum4(sumi1));
        ddNi2dtByVdT = ddNi2dtByVdT - (hsum4(sumi2));
        ddNi3dtByVdT = ddNi3dtByVdT - (hsum4(sumi3));

        ddNi4dtByVdT = ddNi4dtByVdT - (hsum4(sumi4));
        ddNi5dtByVdT = ddNi5dtByVdT - (hsum4(sumi5));
        ddNi6dtByVdT = ddNi6dtByVdT - (hsum4(sumi6));
        ddNi7dtByVdT = ddNi7dtByVdT - (hsum4(sumi7));

        __m256d dYi03dtv = load256d(&dPhidt[i]);
        __m256d ddNi03dtByVdTv = _mm256_setr_pd(ddNi0dtByVdT,ddNi1dtByVdT,ddNi2dtByVdT,ddNi3dtByVdT);
        __m256d result0 = _mm256_fmadd_pd(Wi03ByrhoMv,ddNi03dtByVdTv,_mm256_mul_pd(alphavMv,dYi03dtv));
        Jac[(i+0)*(alignN) + this->nSpecie()] = get0(result0);
        Jac[(i+1)*(alignN) + this->nSpecie()] = get1(result0);
        Jac[(i+2)*(alignN) + this->nSpecie()] = get2(result0);
        Jac[(i+3)*(alignN) + this->nSpecie()] = get3(result0);



        __m256d dYi47dtv = load256d(&dPhidt[i+4]);
        __m256d ddNi47dtByVdTv = _mm256_setr_pd(ddNi4dtByVdT,ddNi5dtByVdT,ddNi6dtByVdT,ddNi7dtByVdT);
        __m256d result1 = _mm256_fmadd_pd(Wi47ByrhoMv,ddNi47dtByVdTv,_mm256_mul_pd(alphavMv,dYi47dtv));
        Jac[(i+4)*(alignN) + this->nSpecie()] = get0(result1);
        Jac[(i+5)*(alignN) + this->nSpecie()] = get1(result1);
        Jac[(i+6)*(alignN) + this->nSpecie()] = get2(result1);
        Jac[(i+7)*(alignN) + this->nSpecie()] = get3(result1);

    }
}