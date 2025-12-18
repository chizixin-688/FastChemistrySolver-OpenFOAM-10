#include "FastChemistryModel.H"
#include <immintrin.h>  

template<class UnusedThermo>
void Foam::FastChemistryModel<UnusedThermo>::ddYdtdT_Vec88_1
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
    for (int i=0; i<this->nSpecie()-1; i=i+8)
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
        __m256d sumi0v = _mm256_setzero_pd();             //6
        __m256d sumi1v = _mm256_setzero_pd();             //7
        __m256d sumi2v = _mm256_setzero_pd();             //8
        __m256d sumi3v = _mm256_setzero_pd();             //9
        __m256d sumi4v = _mm256_setzero_pd();             //10
        __m256d sumi5v = _mm256_setzero_pd();             //11
        __m256d sumi6v = _mm256_setzero_pd();             //12
        __m256d sumi7v = _mm256_setzero_pd();             //13
        for (int j=0; j<this->nSpecie()-1; j=j+8)
        {
            __m256d Cj03v = load256d(&c[j+0]);//14
            __m256d Cj47v = load256d(&c[j+4]);//15

            __m256d Jci0j03v = load256d(&JcRowi0[j]);//16
            sumi0v =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j03v,Cj03v),alphavMv,sumi0v);//14
            __m256d Jci0j47v = load256d(&JcRowi0[j+4]);//15
            sumi0v =_mm256_fmadd_pd(_mm256_mul_pd(Jci0j47v,Cj47v),alphavMv,sumi0v);//14

            __m256d Jci1j03v = load256d(&JcRowi1[j]);//15
            sumi1v =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j03v,Cj03v),alphavMv,sumi1v);
            __m256d Jci1j47v = load256d(&JcRowi1[j+4]);//15
            sumi1v =_mm256_fmadd_pd(_mm256_mul_pd(Jci1j47v,Cj47v),alphavMv,sumi1v);

            __m256d Jci2j03v = load256d(&JcRowi2[j]);
            sumi2v =_mm256_fmadd_pd(_mm256_mul_pd(Jci2j03v,Cj03v),alphavMv,sumi2v);
            __m256d Jci2j47v = load256d(&JcRowi2[j+4]);
            sumi2v =_mm256_fmadd_pd(_mm256_mul_pd(Jci2j47v,Cj47v),alphavMv,sumi2v);


            __m256d Jci3j03v = load256d(&JcRowi3[j]);
            sumi3v =_mm256_fmadd_pd(_mm256_mul_pd(Jci3j03v,Cj03v),alphavMv,sumi3v);
            __m256d Jci3j47v = load256d(&JcRowi3[j+4]);
            sumi3v =_mm256_fmadd_pd(_mm256_mul_pd(Jci3j47v,Cj47v),alphavMv,sumi3v);

            __m256d Jci4j03v = load256d(&JcRowi4[j]);
            sumi4v =_mm256_fmadd_pd(_mm256_mul_pd(Jci4j03v,Cj03v),alphavMv,sumi4v);
            __m256d Jci4j47v = load256d(&JcRowi4[j+4]);
            sumi4v =_mm256_fmadd_pd(_mm256_mul_pd(Jci4j47v,Cj47v),alphavMv,sumi4v);

            __m256d Jci5j03v = load256d(&JcRowi5[j]);
            sumi5v =_mm256_fmadd_pd(_mm256_mul_pd(Jci5j03v,Cj03v),alphavMv,sumi5v);
            __m256d Jci5j47v = load256d(&JcRowi5[j+4]);
            sumi5v =_mm256_fmadd_pd(_mm256_mul_pd(Jci5j47v,Cj47v),alphavMv,sumi5v);

            __m256d Jci6j03v = load256d(&JcRowi6[j]);
            sumi6v =_mm256_fmadd_pd(_mm256_mul_pd(Jci6j03v,Cj03v),alphavMv,sumi6v);
            __m256d Jci6j47v = load256d(&JcRowi6[j+4]);
            sumi6v =_mm256_fmadd_pd(_mm256_mul_pd(Jci6j47v,Cj47v),alphavMv,sumi6v);

            __m256d Jci7j03v = load256d(&JcRowi7[j]);
            sumi7v =_mm256_fmadd_pd(_mm256_mul_pd(Jci7j03v,Cj03v),alphavMv,sumi7v);
            __m256d Jci7j47v = load256d(&JcRowi7[j+4]);
            sumi7v =_mm256_fmadd_pd(_mm256_mul_pd(Jci7j47v,Cj47v),alphavMv,sumi7v);
        }
        {
            int j = this->nSpecie()-1;
            double Cj0 = c[j+0];

            ddNi0dtByVdT = ddNi0dtByVdT - JcRowi0[j+0]*Cj0*alphavM;
            ddNi1dtByVdT = ddNi1dtByVdT - JcRowi1[j+0]*Cj0*alphavM;
            ddNi2dtByVdT = ddNi2dtByVdT - JcRowi2[j+0]*Cj0*alphavM;
            ddNi3dtByVdT = ddNi3dtByVdT - JcRowi3[j+0]*Cj0*alphavM;
            ddNi4dtByVdT = ddNi4dtByVdT - JcRowi4[j+0]*Cj0*alphavM;
            ddNi5dtByVdT = ddNi5dtByVdT - JcRowi5[j+0]*Cj0*alphavM;
            ddNi6dtByVdT = ddNi6dtByVdT - JcRowi6[j+0]*Cj0*alphavM;
            ddNi7dtByVdT = ddNi7dtByVdT - JcRowi7[j+0]*Cj0*alphavM;
        }

        ddNi0dtByVdT = ddNi0dtByVdT - (hsum4(sumi0v));
        ddNi1dtByVdT = ddNi1dtByVdT - (hsum4(sumi1v));
        ddNi2dtByVdT = ddNi2dtByVdT - (hsum4(sumi2v));
        ddNi3dtByVdT = ddNi3dtByVdT - (hsum4(sumi3v));
        ddNi4dtByVdT = ddNi4dtByVdT - (hsum4(sumi4v));
        ddNi5dtByVdT = ddNi5dtByVdT - (hsum4(sumi5v));
        ddNi6dtByVdT = ddNi6dtByVdT - (hsum4(sumi6v));
        ddNi7dtByVdT = ddNi7dtByVdT - (hsum4(sumi7v));

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

    //Do the last specie
    {
        int i = this->nSpecie()-1;
        const double Wi0ByrhoM_ = WiByrhoM[i];
        dPhidt[i] = dPhidt[i]*Wi0ByrhoM_;
        double dYi0dt = dPhidt[i];
        const double* __restrict__ JcRowi0 = &ddNdtByVdcT[i*(alignN)];
        double ddNi0dtByVdT = JcRowi0[this->nSpecie()];
        __m256d ddNi0dtByVdTv = _mm256_setzero_pd();
        for (int j=0; j<this->nSpecie()-1; j=j+8)
        {
            __m256d Jci0j03v = load256d(&JcRowi0[j+0]);
            __m256d Cj03v = load256d(&c[j+0]);
            ddNi0dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(Jci0j03v,Cj03v),alphavMv,ddNi0dtByVdTv);
            __m256d Jci0j47v = load256d(&JcRowi0[j+4]);
            __m256d Cj47v = load256d(&c[j+4]);
            ddNi0dtByVdTv = _mm256_fmadd_pd(_mm256_mul_pd(Jci0j47v,Cj47v),alphavMv,ddNi0dtByVdTv);
        }
        {
            int j = this->nSpecie()-1;
            ddNi0dtByVdT = ddNi0dtByVdT - JcRowi0[j]*c[j]*alphavM - hsum4(ddNi0dtByVdTv);
        }
        
        Jac[i*(this->alignN) + this->nSpecie()] = Wi0ByrhoM_*ddNi0dtByVdT + alphavM*dYi0dt;
    }
}