#include "FastChemistryModel.H"
#include <immintrin.h>  

template<class ThermoType>
void Foam::FastChemistryModel<ThermoType>::ddYdtdT_Vec88_0
(
    const double* __restrict__ ddNdtByVdcT,
    const double* __restrict__ WiByrhoM,
    const double* __restrict__ c,
    double* __restrict__ dPhidt,
    double* __restrict__ Jac,
    double alphavM
) const noexcept
{
    __m256d alphavMv = _mm256_set1_pd(alphavM);//1
    for (int i=0; i<this->nSpecie(); i=i+8)
    {

        __m256d Wi03ByrhoMv = load256d(&WiByrhoM[i]);//2
        __m256d dPhi03dtv = load256d(&dPhidt[i]);   //3
        dPhi03dtv = _mm256_mul_pd(dPhi03dtv,Wi03ByrhoMv);   //3
        store256d(&dPhidt[i],dPhi03dtv);            
        //__m256d dYi03dtv = dPhi03dtv;                       //3

        __m256d Wi47ByrhoMv = load256d(&WiByrhoM[i+4]);//4
        __m256d dPhi47dtv = load256d(&dPhidt[i+4]);   //5
        dPhi47dtv = _mm256_mul_pd(dPhi47dtv,Wi47ByrhoMv);   //5
        store256d(&dPhidt[i+4],dPhi47dtv);            
        //__m256d dYi47dtv = dPhi47dtv;                       //5        

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
        __m256d sumi0 = _mm256_setzero_pd();             //6
        __m256d sumi1 = _mm256_setzero_pd();             //7
        __m256d sumi2 = _mm256_setzero_pd();             //8
        __m256d sumi3 = _mm256_setzero_pd();             //9
        __m256d sumi4 = _mm256_setzero_pd();             //10
        __m256d sumi5 = _mm256_setzero_pd();             //11
        __m256d sumi6 = _mm256_setzero_pd();             //12
        __m256d sumi7 = _mm256_setzero_pd();             //13
        for (int j=0; j<this->nSpecie(); j=j+8)
        {
            __m256d Cj03v = load256d(&c[j+0]);//14
            __m256d Cj47v = load256d(&c[j+4]);//15

            __m256d Jci0j03v = load256d(&JcRowi0[j]);//16
            sumi0 =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j03v,Cj03v),alphavMv,sumi0);//14
            __m256d Jci0j47v = load256d(&JcRowi0[j+4]);//15
            sumi0 =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j47v,Cj47v),alphavMv,sumi0);//14

            __m256d Jci1j03v = load256d(&JcRowi1[j]);//15
            sumi1 =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j03v,Cj03v),alphavMv,sumi1);
            __m256d Jci1j47v = load256d(&JcRowi1[j+4]);//15
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
}