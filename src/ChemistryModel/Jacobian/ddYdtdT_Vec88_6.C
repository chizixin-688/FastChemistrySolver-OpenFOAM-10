/*---------------------------------------------------------------------------*\
  Description
      Computing mass fraction based jacobian matrix using molar concentration 
      based matrix
      Computing partial derivatives of dYdt w.r.t temperature
      the tail is 6.  

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
void Foam::FastChemistryModel<UnusedThermo>::ddYdtdT_Vec88_6
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
    for (int i=0; i<this->nSpecie()-6; i=i+8)
    {

        __m256d Wi03ByrhoMv = load256d(&WiByrhoM[i+0]);
        __m256d dPhi03dt = load256d(&dPhidt[i+0]);   
        dPhi03dt = _mm256_mul_pd(dPhi03dt,Wi03ByrhoMv);   
        store256d(&dPhidt[i+0],dPhi03dt);

        __m256d Wi47ByrhoMv = load256d(&WiByrhoM[i+4]);
        __m256d dPhi47dt = load256d(&dPhidt[i+4]);   
        dPhi47dt = _mm256_mul_pd(dPhi47dt,Wi47ByrhoMv);   
        store256d(&dPhidt[i+4],dPhi47dt);

        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*(alignN)];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*(alignN)];
        const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*(alignN)];
        const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i+3)*(alignN)];
        const double* __restrict__ JcRowi4 = &ddNdtByVdcT[(i+4)*(alignN)];
        const double* __restrict__ JcRowi5 = &ddNdtByVdcT[(i+5)*(alignN)];
        const double* __restrict__ JcRowi6 = &ddNdtByVdcT[(i+6)*(alignN)];
        const double* __restrict__ JcRowi7 = &ddNdtByVdcT[(i+7)*(alignN)];

        double ddNi0dtByVdT = JcRowi0[this->nSpecie()];
        double ddNi1dtByVdT = JcRowi1[this->nSpecie()];
        double ddNi2dtByVdT = JcRowi2[this->nSpecie()];
        double ddNi3dtByVdT = JcRowi3[this->nSpecie()];
        double ddNi4dtByVdT = JcRowi4[this->nSpecie()];
        double ddNi5dtByVdT = JcRowi5[this->nSpecie()];
        double ddNi6dtByVdT = JcRowi6[this->nSpecie()];
        double ddNi7dtByVdT = JcRowi7[this->nSpecie()];
        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        __m256d sumi2v = _mm256_setzero_pd();
        __m256d sumi3v = _mm256_setzero_pd();
        __m256d sumi4v = _mm256_setzero_pd();
        __m256d sumi5v = _mm256_setzero_pd();
        __m256d sumi6v = _mm256_setzero_pd();
        __m256d sumi7v = _mm256_setzero_pd();
        for (int j=0; j<this->nSpecie()-6; j=j+8)
        {
            __m256d Cj03v = load256d(&c[j+0]);
            __m256d Cj47v = load256d(&c[j+4]);

            __m256d Jci0j03v = load256d(&JcRowi0[j+0]);
            sumi0v =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j03v,Cj03v),alphavMv,sumi0v);
            __m256d Jci0j47v = load256d(&JcRowi0[j+4]);
            sumi0v =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j47v,Cj47v),alphavMv,sumi0v);

            __m256d Jci1j03v = load256d(&JcRowi1[j+0]);
            sumi1v =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j03v,Cj03v),alphavMv,sumi1v);
            __m256d Jci1j47v = load256d(&JcRowi1[j+4]);
            sumi1v =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j47v,Cj47v),alphavMv,sumi1v);

            __m256d Jci2j03v = load256d(&JcRowi2[j+0]);
            sumi2v =_mm256_fmadd_pd(_mm256_mul_pd(Jci2j03v,Cj03v),alphavMv,sumi2v);
            __m256d Jci2j47v = load256d(&JcRowi2[j+4]);
            sumi2v =_mm256_fmadd_pd(_mm256_mul_pd(Jci2j47v,Cj47v),alphavMv,sumi2v);

            __m256d Jci3j03v = load256d(&JcRowi3[j+0]);
            sumi3v =_mm256_fmadd_pd(_mm256_mul_pd(Jci3j03v,Cj03v),alphavMv,sumi3v);
            __m256d Jci3j47v = load256d(&JcRowi3[j+4]);
            sumi3v =_mm256_fmadd_pd(_mm256_mul_pd(Jci3j47v,Cj47v),alphavMv,sumi3v);

            __m256d Jci4j03v = load256d(&JcRowi4[j+0]);
            sumi4v =_mm256_fmadd_pd(_mm256_mul_pd(Jci4j03v,Cj03v),alphavMv,sumi4v);
            __m256d Jci4j47v = load256d(&JcRowi4[j+4]);
            sumi4v =_mm256_fmadd_pd(_mm256_mul_pd(Jci4j47v,Cj47v),alphavMv,sumi4v);

            __m256d Jci5j03v = load256d(&JcRowi5[j+0]);
            sumi5v =_mm256_fmadd_pd(_mm256_mul_pd(Jci5j03v,Cj03v),alphavMv,sumi5v);
            __m256d Jci5j47v = load256d(&JcRowi5[j+4]);
            sumi5v =_mm256_fmadd_pd(_mm256_mul_pd(Jci5j47v,Cj47v),alphavMv,sumi5v);

            __m256d Jci6j03v = load256d(&JcRowi6[j+0]);
            sumi6v =_mm256_fmadd_pd(_mm256_mul_pd(Jci6j03v,Cj03v),alphavMv,sumi6v);
            __m256d Jci6j47v = load256d(&JcRowi6[j+4]);
            sumi6v =_mm256_fmadd_pd(_mm256_mul_pd(Jci6j47v,Cj47v),alphavMv,sumi6v);

            __m256d Jci7j03v = load256d(&JcRowi7[j+0]);
            sumi7v =_mm256_fmadd_pd(_mm256_mul_pd(Jci7j03v,Cj03v),alphavMv,sumi7v);
            __m256d Jci7j47v = load256d(&JcRowi7[j+4]);
            sumi7v =_mm256_fmadd_pd(_mm256_mul_pd(Jci7j47v,Cj47v),alphavMv,sumi7v);
        }
        {
            int j = this->nSpecie()-6;
            __m256d Cj03v = load256d(&c[j+0]);
            __m256d Cj45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&c[j+4]));

            __m256d Jci0j03v = load256d(&JcRowi0[j+0]);
            sumi0v =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j03v,Cj03v),alphavMv,sumi0v);
            __m256d Jci0j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi0[j+4]));
            sumi0v =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j45v,Cj45v),alphavMv,sumi0v);

            __m256d Jci1j03v = load256d(&JcRowi1[j+0]);
            sumi1v =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j03v,Cj03v),alphavMv,sumi1v);
            __m256d Jci1j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi1[j+4]));
            sumi1v =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j45v,Cj45v),alphavMv,sumi1v);

            __m256d Jci2j03v = load256d(&JcRowi2[j+0]);
            sumi2v =_mm256_fmadd_pd(_mm256_mul_pd(Jci2j03v,Cj03v),alphavMv,sumi2v);
            __m256d Jci2j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi2[j+4]));
            sumi2v =_mm256_fmadd_pd(_mm256_mul_pd(Jci2j45v,Cj45v),alphavMv,sumi2v);

            __m256d Jci3j03v = load256d(&JcRowi3[j+0]);
            sumi3v =_mm256_fmadd_pd(_mm256_mul_pd(Jci3j03v,Cj03v),alphavMv,sumi3v);
            __m256d Jci3j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi3[j+4]));
            sumi3v =_mm256_fmadd_pd(_mm256_mul_pd(Jci3j45v,Cj45v),alphavMv,sumi3v);

            __m256d Jci4j03v = load256d(&JcRowi4[j+0]);
            sumi4v =_mm256_fmadd_pd(_mm256_mul_pd(Jci4j03v,Cj03v),alphavMv,sumi4v);
            __m256d Jci4j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi4[j+4]));
            sumi4v =_mm256_fmadd_pd(_mm256_mul_pd(Jci4j45v,Cj45v),alphavMv,sumi4v);

            __m256d Jci5j03v = load256d(&JcRowi5[j+0]);
            sumi5v =_mm256_fmadd_pd(_mm256_mul_pd(Jci5j03v,Cj03v),alphavMv,sumi5v);
            __m256d Jci5j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi5[j+4]));
            sumi5v =_mm256_fmadd_pd(_mm256_mul_pd(Jci5j45v,Cj45v),alphavMv,sumi5v);

            __m256d Jci6j03v = load256d(&JcRowi6[j+0]);
            sumi6v =_mm256_fmadd_pd(_mm256_mul_pd(Jci6j03v,Cj03v),alphavMv,sumi6v);
            __m256d Jci6j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi6[j+4]));
            sumi6v =_mm256_fmadd_pd(_mm256_mul_pd(Jci6j45v,Cj45v),alphavMv,sumi6v);

            __m256d Jci7j03v = load256d(&JcRowi7[j+0]);
            sumi7v =_mm256_fmadd_pd(_mm256_mul_pd(Jci7j03v,Cj03v),alphavMv,sumi7v);
            __m256d Jci7j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi7[j+4]));
            sumi7v =_mm256_fmadd_pd(_mm256_mul_pd(Jci7j45v,Cj45v),alphavMv,sumi7v);
        }

        ddNi0dtByVdT = ddNi0dtByVdT - (hsum4(sumi0v));
        ddNi1dtByVdT = ddNi1dtByVdT - (hsum4(sumi1v));
        ddNi2dtByVdT = ddNi2dtByVdT - (hsum4(sumi2v));
        ddNi3dtByVdT = ddNi3dtByVdT - (hsum4(sumi3v));
        ddNi4dtByVdT = ddNi4dtByVdT - (hsum4(sumi4v));
        ddNi5dtByVdT = ddNi5dtByVdT - (hsum4(sumi5v));
        ddNi6dtByVdT = ddNi6dtByVdT - (hsum4(sumi6v));
        ddNi7dtByVdT = ddNi7dtByVdT - (hsum4(sumi7v));

        __m256d dYi03dtv = load256d(&dPhidt[i+0]);
        __m256d ddNi03dtByVdTv = _mm256_setr_pd(ddNi0dtByVdT,ddNi1dtByVdT,ddNi2dtByVdT,ddNi3dtByVdT);
        __m256d result0 = _mm256_fmadd_pd(Wi03ByrhoMv,ddNi03dtByVdTv,_mm256_mul_pd(alphavMv,dYi03dtv));
        Jac[(i+0)*(alignN) + this->nSpecie()] = get_elem0(result0);
        Jac[(i+1)*(alignN) + this->nSpecie()] = get_elem1(result0);
        Jac[(i+2)*(alignN) + this->nSpecie()] = get_elem2(result0);
        Jac[(i+3)*(alignN) + this->nSpecie()] = get_elem3(result0);

        __m256d dYi47dtv = load256d(&dPhidt[i+4]);
        __m256d ddNi47dtByVdTv = _mm256_setr_pd(ddNi4dtByVdT,ddNi5dtByVdT,ddNi6dtByVdT,ddNi7dtByVdT);
        __m256d result1 = _mm256_fmadd_pd(Wi47ByrhoMv,ddNi47dtByVdTv,_mm256_mul_pd(alphavMv,dYi47dtv));

        Jac[(i+4)*(alignN) + this->nSpecie()] = get_elem0(result1);
        Jac[(i+5)*(alignN) + this->nSpecie()] = get_elem1(result1);
        Jac[(i+6)*(alignN) + this->nSpecie()] = get_elem2(result1);
        Jac[(i+7)*(alignN) + this->nSpecie()] = get_elem3(result1);
    }

    {
        int i = this->nSpecie()-6;
        __m256d Wi03ByrhoMv = load256d(&WiByrhoM[i+0]);
        __m256d dPhi03dt = load256d(&dPhidt[i+0]);   
        dPhi03dt = _mm256_mul_pd(dPhi03dt,Wi03ByrhoMv);   
        store256d(&dPhidt[i+0],dPhi03dt);

        __m256d Wi45ByrhoMv = _mm256_zextpd128_pd256(_mm_loadu_pd(&WiByrhoM[i+4]));
        __m256d dPhi45dt = _mm256_zextpd128_pd256(_mm_loadu_pd(&dPhidt[i+4]));   
        dPhi45dt = _mm256_mul_pd(dPhi45dt,Wi45ByrhoMv);   
        store256d(&dPhidt[i+4],dPhi45dt);

        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*(alignN)];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*(alignN)];
        const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*(alignN)];
        const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i+3)*(alignN)];
        const double* __restrict__ JcRowi4 = &ddNdtByVdcT[(i+4)*(alignN)];
        const double* __restrict__ JcRowi5 = &ddNdtByVdcT[(i+5)*(alignN)];

        double ddNi0dtByVdT = JcRowi0[this->nSpecie()];
        double ddNi1dtByVdT = JcRowi1[this->nSpecie()];
        double ddNi2dtByVdT = JcRowi2[this->nSpecie()];
        double ddNi3dtByVdT = JcRowi3[this->nSpecie()];
        double ddNi4dtByVdT = JcRowi4[this->nSpecie()];
        double ddNi5dtByVdT = JcRowi5[this->nSpecie()];

        __m256d sumi0v = _mm256_setzero_pd();
        __m256d sumi1v = _mm256_setzero_pd();
        __m256d sumi2v = _mm256_setzero_pd();
        __m256d sumi3v = _mm256_setzero_pd();
        __m256d sumi4v = _mm256_setzero_pd();
        __m256d sumi5v = _mm256_setzero_pd();

        for (int j=0; j<this->nSpecie()-6; j=j+8)
        {
            __m256d Cj03v = load256d(&c[j+0]);
            __m256d Cj47v = load256d(&c[j+4]);

            __m256d Jci0j03v = load256d(&JcRowi0[j+0]);
            sumi0v =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j03v,Cj03v),alphavMv,sumi0v);
            __m256d Jci0j47v = load256d(&JcRowi0[j+4]);
            sumi0v =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j47v,Cj47v),alphavMv,sumi0v);

            __m256d Jci1j03v = load256d(&JcRowi1[j+0]);
            sumi1v =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j03v,Cj03v),alphavMv,sumi1v);
            __m256d Jci1j47v = load256d(&JcRowi1[j+4]);
            sumi1v =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j47v,Cj47v),alphavMv,sumi1v);

            __m256d Jci2j03v = load256d(&JcRowi2[j+0]);
            sumi2v =_mm256_fmadd_pd(_mm256_mul_pd(Jci2j03v,Cj03v),alphavMv,sumi2v);
            __m256d Jci2j47v = load256d(&JcRowi2[j+4]);
            sumi2v =_mm256_fmadd_pd(_mm256_mul_pd(Jci2j47v,Cj47v),alphavMv,sumi2v);

            __m256d Jci3j03v = load256d(&JcRowi3[j+0]);
            sumi3v =_mm256_fmadd_pd(_mm256_mul_pd(Jci3j03v,Cj03v),alphavMv,sumi3v);
            __m256d Jci3j47v = load256d(&JcRowi3[j+4]);
            sumi3v =_mm256_fmadd_pd(_mm256_mul_pd(Jci3j47v,Cj47v),alphavMv,sumi3v);

            __m256d Jci4j03v = load256d(&JcRowi4[j+0]);
            sumi4v =_mm256_fmadd_pd(_mm256_mul_pd(Jci4j03v,Cj03v),alphavMv,sumi4v);
            __m256d Jci4j47v = load256d(&JcRowi4[j+4]);
            sumi4v =_mm256_fmadd_pd(_mm256_mul_pd(Jci4j47v,Cj47v),alphavMv,sumi4v);

            __m256d Jci5j03v = load256d(&JcRowi5[j+0]);
            sumi5v =_mm256_fmadd_pd(_mm256_mul_pd(Jci5j03v,Cj03v),alphavMv,sumi5v);
            __m256d Jci5j47v = load256d(&JcRowi5[j+4]);
            sumi5v =_mm256_fmadd_pd(_mm256_mul_pd(Jci5j47v,Cj47v),alphavMv,sumi5v);
        }
        {
            int j = this->nSpecie()-6;
            __m256d Cj03v = load256d(&c[j+0]);
            __m256d Cj45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&c[j+4]));

            __m256d Jci0j03v = load256d(&JcRowi0[j+0]);
            sumi0v =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j03v,Cj03v),alphavMv,sumi0v);
            __m256d Jci0j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi0[j+4]));
            sumi0v =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j45v,Cj45v),alphavMv,sumi0v);

            __m256d Jci1j03v = load256d(&JcRowi1[j+0]);
            sumi1v =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j03v,Cj03v),alphavMv,sumi1v);
            __m256d Jci1j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi1[j+4]));
            sumi1v =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j45v,Cj45v),alphavMv,sumi1v);

            __m256d Jci2j03v = load256d(&JcRowi2[j+0]);
            sumi2v =_mm256_fmadd_pd(_mm256_mul_pd(Jci2j03v,Cj03v),alphavMv,sumi2v);
            __m256d Jci2j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi2[j+4]));
            sumi2v =_mm256_fmadd_pd(_mm256_mul_pd(Jci2j45v,Cj45v),alphavMv,sumi2v);

            __m256d Jci3j03v = load256d(&JcRowi3[j+0]);
            sumi3v =_mm256_fmadd_pd(_mm256_mul_pd(Jci3j03v,Cj03v),alphavMv,sumi3v);
            __m256d Jci3j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi3[j+4]));
            sumi3v =_mm256_fmadd_pd(_mm256_mul_pd(Jci3j45v,Cj45v),alphavMv,sumi3v);

            __m256d Jci4j03v = load256d(&JcRowi4[j+0]);
            sumi4v =_mm256_fmadd_pd(_mm256_mul_pd(Jci4j03v,Cj03v),alphavMv,sumi4v);
            __m256d Jci4j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi4[j+4]));
            sumi4v =_mm256_fmadd_pd(_mm256_mul_pd(Jci4j45v,Cj45v),alphavMv,sumi4v);

            __m256d Jci5j03v = load256d(&JcRowi5[j+0]);
            sumi5v =_mm256_fmadd_pd(_mm256_mul_pd(Jci5j03v,Cj03v),alphavMv,sumi5v);
            __m256d Jci5j45v = _mm256_zextpd128_pd256(_mm_loadu_pd(&JcRowi5[j+4]));
            sumi5v =_mm256_fmadd_pd(_mm256_mul_pd(Jci5j45v,Cj45v),alphavMv,sumi5v);
        }

        ddNi0dtByVdT = ddNi0dtByVdT - (hsum4(sumi0v));
        ddNi1dtByVdT = ddNi1dtByVdT - (hsum4(sumi1v));
        ddNi2dtByVdT = ddNi2dtByVdT - (hsum4(sumi2v));
        ddNi3dtByVdT = ddNi3dtByVdT - (hsum4(sumi3v));
        ddNi4dtByVdT = ddNi4dtByVdT - (hsum4(sumi4v));
        ddNi5dtByVdT = ddNi5dtByVdT - (hsum4(sumi5v));


        __m256d dYi03dtv = load256d(&dPhidt[i+0]);
        __m256d ddNi03dtByVdTv = _mm256_setr_pd(ddNi0dtByVdT,ddNi1dtByVdT,ddNi2dtByVdT,ddNi3dtByVdT);
        __m256d result0 = _mm256_fmadd_pd(Wi03ByrhoMv,ddNi03dtByVdTv,_mm256_mul_pd(alphavMv,dYi03dtv));
        Jac[(i+0)*(alignN) + this->nSpecie()] = get_elem0(result0);
        Jac[(i+1)*(alignN) + this->nSpecie()] = get_elem1(result0);
        Jac[(i+2)*(alignN) + this->nSpecie()] = get_elem2(result0);
        Jac[(i+3)*(alignN) + this->nSpecie()] = get_elem3(result0);

        __m256d dYi45dtv = _mm256_zextpd128_pd256(_mm_loadu_pd(&dPhidt[i+4]));
        __m256d ddNi45dtByVdTv = _mm256_setr_pd(ddNi4dtByVdT,ddNi5dtByVdT,0,0);
        __m256d result1 = _mm256_fmadd_pd(Wi45ByrhoMv,ddNi45dtByVdTv,_mm256_mul_pd(alphavMv,dYi45dtv));

        Jac[(i+4)*(alignN) + this->nSpecie()] = get_elem0(result1);
        Jac[(i+5)*(alignN) + this->nSpecie()] = get_elem1(result1);
    }
}