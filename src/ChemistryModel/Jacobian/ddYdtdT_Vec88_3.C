#include "FastChemistryModel.H"
#include <immintrin.h>  

template<class ThermoType>
void Foam::FastChemistryModel<ThermoType>::ddYdtdT_Vec88_3
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
    for (int i=0; i<this->nSpecie()-3; i=i+8)
    {
        
            __m256d Wi03ByrhoMv = load256d(&WiByrhoM[i+0]);
            __m256d dPhi03dtv = load256d(&dPhidt[i+0]);   
            dPhi03dtv = _mm256_mul_pd(dPhi03dtv,Wi03ByrhoMv);   
            store256d(&dPhidt[i+0],dPhi03dtv);

            __m256d Wi47ByrhoMv = load256d(&WiByrhoM[i+4]);
            __m256d dPhi47dtv = load256d(&dPhidt[i+4]);   
            dPhi47dtv = _mm256_mul_pd(dPhi47dtv,Wi47ByrhoMv);   
            store256d(&dPhidt[i+4],dPhi47dtv);
        


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
        for (int j=0; j<this->nSpecie()-3; j=j+8)
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
        ddNi0dtByVdT = ddNi0dtByVdT - (hsum4(sumi0v));
        ddNi1dtByVdT = ddNi1dtByVdT - (hsum4(sumi1v));
        ddNi2dtByVdT = ddNi2dtByVdT - (hsum4(sumi2v));
        ddNi3dtByVdT = ddNi3dtByVdT - (hsum4(sumi3v));
        ddNi4dtByVdT = ddNi4dtByVdT - (hsum4(sumi4v));
        ddNi5dtByVdT = ddNi5dtByVdT - (hsum4(sumi5v));
        ddNi6dtByVdT = ddNi6dtByVdT - (hsum4(sumi6v));
        ddNi7dtByVdT = ddNi7dtByVdT - (hsum4(sumi7v));
        {
            int j = this->nSpecie()-3;
            double cj0 = c[j+0];
            ddNi0dtByVdT -= JcRowi0[j+0]*cj0*alphavM;
            ddNi1dtByVdT -= JcRowi1[j+0]*cj0*alphavM;
            ddNi2dtByVdT -= JcRowi2[j+0]*cj0*alphavM;
            ddNi3dtByVdT -= JcRowi3[j+0]*cj0*alphavM;
            ddNi4dtByVdT -= JcRowi4[j+0]*cj0*alphavM;
            ddNi5dtByVdT -= JcRowi5[j+0]*cj0*alphavM;
            ddNi6dtByVdT -= JcRowi6[j+0]*cj0*alphavM;
            ddNi7dtByVdT -= JcRowi7[j+0]*cj0*alphavM;

            double cj1 = c[j+1];
            ddNi0dtByVdT -= JcRowi0[j+1]*cj1*alphavM;
            ddNi1dtByVdT -= JcRowi1[j+1]*cj1*alphavM;
            ddNi2dtByVdT -= JcRowi2[j+1]*cj1*alphavM;
            ddNi3dtByVdT -= JcRowi3[j+1]*cj1*alphavM;
            ddNi4dtByVdT -= JcRowi4[j+1]*cj1*alphavM;
            ddNi5dtByVdT -= JcRowi5[j+1]*cj1*alphavM;
            ddNi6dtByVdT -= JcRowi6[j+1]*cj1*alphavM;
            ddNi7dtByVdT -= JcRowi7[j+1]*cj1*alphavM;

            double cj2 = c[j+2];
            ddNi0dtByVdT -= JcRowi0[j+2]*cj2*alphavM;
            ddNi1dtByVdT -= JcRowi1[j+2]*cj2*alphavM;
            ddNi2dtByVdT -= JcRowi2[j+2]*cj2*alphavM;
            ddNi3dtByVdT -= JcRowi3[j+2]*cj2*alphavM;
            ddNi4dtByVdT -= JcRowi4[j+2]*cj2*alphavM;
            ddNi5dtByVdT -= JcRowi5[j+2]*cj2*alphavM;
            ddNi6dtByVdT -= JcRowi6[j+2]*cj2*alphavM;
            ddNi7dtByVdT -= JcRowi7[j+2]*cj2*alphavM;
        }
        __m256d dY03dtv = load256d(&dPhidt[i+0]);
        __m256d ddNi03dtByVdTv = _mm256_setr_pd(ddNi0dtByVdT,ddNi1dtByVdT,ddNi2dtByVdT,ddNi3dtByVdT);
        __m256d result0 = _mm256_fmadd_pd(Wi03ByrhoMv,ddNi03dtByVdTv,_mm256_mul_pd(alphavMv,dY03dtv));
        Jac[(i+0)*(alignN) + this->nSpecie()] = get_elem0(result0);
        Jac[(i+1)*(alignN) + this->nSpecie()] = get_elem1(result0);
        Jac[(i+2)*(alignN) + this->nSpecie()] = get_elem2(result0);
        Jac[(i+3)*(alignN) + this->nSpecie()] = get_elem3(result0);

        __m256d dY47dtv = load256d(&dPhidt[i+4]);
        __m256d ddNi47dtByVdTv = _mm256_setr_pd(ddNi4dtByVdT,ddNi5dtByVdT,ddNi6dtByVdT,ddNi7dtByVdT);
        __m256d result1 = _mm256_fmadd_pd(Wi47ByrhoMv,ddNi47dtByVdTv,_mm256_mul_pd(alphavMv,dY47dtv));
        Jac[(i+4)*(alignN) + this->nSpecie()] = get_elem0(result1);
        Jac[(i+5)*(alignN) + this->nSpecie()] = get_elem1(result1);
        Jac[(i+6)*(alignN) + this->nSpecie()] = get_elem2(result1);
        Jac[(i+7)*(alignN) + this->nSpecie()] = get_elem3(result1);
      
    }
    {
        int i = this->nSpecie()-3;

        const double Wi0ByrhoM_ = WiByrhoM[i+0];
        const double Wi1ByrhoM_ = WiByrhoM[i+1];
        const double Wi2ByrhoM_ = WiByrhoM[i+2];

        dPhidt[i+0] = dPhidt[i+0]*Wi0ByrhoM_;
        dPhidt[i+1] = dPhidt[i+1]*Wi1ByrhoM_;
        dPhidt[i+2] = dPhidt[i+2]*Wi2ByrhoM_;

        double dYi0dt = dPhidt[i+0];
        double dYi1dt = dPhidt[i+1];
        double dYi2dt = dPhidt[i+2];

        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*(alignN)];
        const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*(alignN)];
        const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*(alignN)];

        double ddNi0dtByVdT = JcRowi0[this->nSpecie()];
        double ddNi1dtByVdT = JcRowi1[this->nSpecie()];
        double ddNi2dtByVdT = JcRowi2[this->nSpecie()];
        __m256d ddNi0dtByVdTv = _mm256_setzero_pd();
        __m256d ddNi1dtByVdTv = _mm256_setzero_pd();
        __m256d ddNi2dtByVdTv = _mm256_setzero_pd();
        for (int j=0; j<this->nSpecie()-3; j=j+8)
        {
            __m256d Cj03v = load256d(&c[j+0]);
            __m256d Jci0j03v = load256d(&JcRowi0[j+0]);
            ddNi0dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(Jci0j03v,Cj03v),alphavMv,ddNi0dtByVdTv);
            __m256d Jci1j03v = load256d(&JcRowi1[j+0]);
            ddNi1dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(Jci1j03v,Cj03v),alphavMv,ddNi1dtByVdTv);
            __m256d Jci2j03v = load256d(&JcRowi2[j+0]);
            ddNi2dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(Jci2j03v,Cj03v),alphavMv,ddNi2dtByVdTv);  
            
            __m256d Cj47v = load256d(&c[j+4]);
            __m256d Jci0j47v = load256d(&JcRowi0[j+4]);
            ddNi0dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(Jci0j47v,Cj47v),alphavMv,ddNi0dtByVdTv);
            __m256d Jci1j47v = load256d(&JcRowi1[j+4]);
            ddNi1dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(Jci1j47v,Cj47v),alphavMv,ddNi1dtByVdTv);
            __m256d Jci2j47v = load256d(&JcRowi2[j+4]);
            ddNi2dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(Jci2j47v,Cj47v),alphavMv,ddNi2dtByVdTv); 
        }
        {
            int j = this->nSpecie()-3;

            ddNi0dtByVdT -= JcRowi0[j+0]*c[j+0]*alphavM;
            ddNi1dtByVdT -= JcRowi1[j+0]*c[j+0]*alphavM;
            ddNi2dtByVdT -= JcRowi2[j+0]*c[j+0]*alphavM;

            ddNi0dtByVdT -= JcRowi0[j+1]*c[j+1]*alphavM;
            ddNi1dtByVdT -= JcRowi1[j+1]*c[j+1]*alphavM;
            ddNi2dtByVdT -= JcRowi2[j+1]*c[j+1]*alphavM; 

            ddNi0dtByVdT -= JcRowi0[j+2]*c[j+2]*alphavM;
            ddNi1dtByVdT -= JcRowi1[j+2]*c[j+2]*alphavM;
            ddNi2dtByVdT -= JcRowi2[j+2]*c[j+2]*alphavM;
        }
        ddNi0dtByVdT -= hsum4(ddNi0dtByVdTv);
        ddNi1dtByVdT -= hsum4(ddNi1dtByVdTv);
        ddNi2dtByVdT -= hsum4(ddNi2dtByVdTv);
        Jac[(i+0)*(alignN) + this->nSpecie()] = Wi0ByrhoM_*ddNi0dtByVdT + alphavM*dYi0dt;
        Jac[(i+1)*(alignN) + this->nSpecie()] = Wi1ByrhoM_*ddNi1dtByVdT + alphavM*dYi1dt;
        Jac[(i+2)*(alignN) + this->nSpecie()] = Wi2ByrhoM_*ddNi2dtByVdT + alphavM*dYi2dt;
    }
}