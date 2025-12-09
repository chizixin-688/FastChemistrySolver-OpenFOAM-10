#include "FastChemistryModel.H"
#include <immintrin.h>  

template<class ThermoType>
void Foam::FastChemistryModel<ThermoType>::FastddYdtdY_Vec44_3
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
    __m256d zerov = _mm256_setzero_pd();
    __m256d rhoMv = _mm256_set1_pd(rhoM);
    for(int i=0; i<this->nSpecie()-3; i=i+4)
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

        for (int j=0; j<this->nSpecie()-3; j=j+4)
        {
            __m256d invWv = load256d(&invW[j+0]);//10
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,invWv);//10
            __m256d rhoMvj_ = load256d(&rhoMByRhoi[j+0]);//11

            __m256d Jci0j03v = load256d(&ddNdtByVdcT[(i+0)*(alignN)+j+0]);//12
            __m256d ddNi0dtByVdYj = _mm256_mul_pd(Jci0j03v,dCjdYj);//12
            __m256d r0 = _mm256_mul_pd(rhoMvj_,dPhi0dtv);//13
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,r0);//12
            store256d(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);//11

            __m256d Jci1j03v = load256d(&ddNdtByVdcT[(i+1)*(alignN)+j+0]);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(Jci1j03v,dCjdYj);
            __m256d r1 = _mm256_mul_pd(rhoMvj_,dPhi1dtv);//13
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,r1);
            store256d(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);

            __m256d Jci2j03v = load256d(&ddNdtByVdcT[(i+2)*(alignN)+j+0]);
            __m256d ddNi2dtByVdYj = _mm256_mul_pd(Jci2j03v,dCjdYj);
            __m256d r2 = _mm256_mul_pd(rhoMvj_,dPhi2dtv);//13
            __m256d ddYi2dtdYj = _mm256_fmadd_pd(Wi2ByrhoMv,ddNi2dtByVdYj,r2);
            store256d(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYj);

            __m256d Jci3j03v = load256d(&ddNdtByVdcT[(i+3)*(alignN)+j+0]);
            __m256d ddNi3dtByVdYj = _mm256_mul_pd(Jci3j03v,dCjdYj);//12
            __m256d r3 = _mm256_mul_pd(rhoMvj_,dPhi3dtv);//13
            __m256d ddYi3dtdYj = _mm256_fmadd_pd(Wi3ByrhoMv,ddNi3dtByVdYj,r3);
            store256d(&Jac[(i+3)*(alignN) + j+0],ddYi3dtdYj);
        }
        {
            int j = this->nSpecie()-3;

            __m256d invWv = _mm256_blend_pd(load256d(&invW[j+0]),zerov,0b1000);//10
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,invWv);//10
            __m256d rhoMvj_ = _mm256_blend_pd(load256d(&rhoMByRhoi[j+0]),zerov,0b1000);//11

            __m256d Jci0j03v = _mm256_blend_pd(load256d(&ddNdtByVdcT[(i+0)*(alignN)+j+0]),zerov,0b1000);//12
            __m256d ddNi0dtByVdYj = _mm256_mul_pd(Jci0j03v,dCjdYj);//12
            __m256d r0 = _mm256_mul_pd(rhoMvj_,dPhi0dtv);//13
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,r0);//12
            store256d(&Jac[(i+0)*(alignN)+j+0],ddYi0dtdYj);//11

            __m256d Jci1j03v = _mm256_blend_pd(load256d(&ddNdtByVdcT[(i+1)*(alignN)+j+0]),zerov,0b1000);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(Jci1j03v,dCjdYj);
            __m256d r1 = _mm256_mul_pd(rhoMvj_,dPhi1dtv);//13
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,r1);
            store256d(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);

            __m256d Jci2j03v = _mm256_blend_pd(load256d(&ddNdtByVdcT[(i+2)*(alignN)+j+0]),zerov,0b1000);
            __m256d ddNi2dtByVdYj = _mm256_mul_pd(Jci2j03v,dCjdYj);
            __m256d r2 = _mm256_mul_pd(rhoMvj_,dPhi2dtv);//13
            __m256d ddYi2dtdYj = _mm256_fmadd_pd(Wi2ByrhoMv,ddNi2dtByVdYj,r2);
            store256d(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYj);

            __m256d Jci3j03v = _mm256_blend_pd(load256d(&ddNdtByVdcT[(i+3)*(alignN)+j+0]),zerov,0b1000);
            __m256d ddNi3dtByVdYj = _mm256_mul_pd(Jci3j03v,dCjdYj);//12
            __m256d r3 = _mm256_mul_pd(rhoMvj_,dPhi3dtv);//13
            __m256d ddYi3dtdYj = _mm256_fmadd_pd(Wi3ByrhoMv,ddNi3dtByVdYj,r3);
            store256d(&Jac[(i+3)*(alignN) + j+0],ddYi3dtdYj);
        }
    }
    {
        int i = this->nSpecie()-3;

        __m256d Wi03ByrhoMv = _mm256_blend_pd(load256d(&WiByrhoM[i+0]),zerov,0b1000);
        __m256d dPhi03dtv = _mm256_blend_pd(load256d(&dPhidt[i+0]),zerov,0b1000);
        dPhi03dtv = _mm256_mul_pd(dPhi03dtv,Wi03ByrhoMv);

        __m256d Wi0ByrhoMv = _mm256_permute4x64_pd(Wi03ByrhoMv, 0x00); 
        __m256d Wi1ByrhoMv = _mm256_permute4x64_pd(Wi03ByrhoMv, 0x55); 
        __m256d Wi2ByrhoMv = _mm256_permute4x64_pd(Wi03ByrhoMv, 0xAA); 

        __m256d dPhi0dtv = _mm256_permute4x64_pd(dPhi03dtv, 0x00); 
        __m256d dPhi1dtv = _mm256_permute4x64_pd(dPhi03dtv, 0x55); 
        __m256d dPhi2dtv = _mm256_permute4x64_pd(dPhi03dtv, 0xAA); 

        for (int j=0; j<this->nSpecie()-3; j=j+4)
        {
            __m256d invWv = load256d(&invW[j+0]);//10
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,invWv);//10
            __m256d rhoMvj_ = load256d(&rhoMByRhoi[j+0]);//11

            __m256d Jci0j03v = load256d(&ddNdtByVdcT[(i+0)*(alignN)+j+0]);//12
            __m256d ddNi0dtByVdYj = _mm256_mul_pd(Jci0j03v,dCjdYj);//12
            __m256d r0 = _mm256_mul_pd(rhoMvj_,dPhi0dtv);//13
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,r0);//12
            store256d(&Jac[(i+0)*(alignN) + j+0],ddYi0dtdYj);//11

            __m256d Jci1j03v = load256d(&ddNdtByVdcT[(i+1)*(alignN)+j+0]);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(Jci1j03v,dCjdYj);
            __m256d r1 = _mm256_mul_pd(rhoMvj_,dPhi1dtv);//13
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,r1);
            store256d(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);

            __m256d Jci2j03v = load256d(&ddNdtByVdcT[(i+2)*(alignN)+j+0]);
            __m256d ddNi2dtByVdYj = _mm256_mul_pd(Jci2j03v,dCjdYj);
            __m256d r2 = _mm256_mul_pd(rhoMvj_,dPhi2dtv);//13
            __m256d ddYi2dtdYj = _mm256_fmadd_pd(Wi2ByrhoMv,ddNi2dtByVdYj,r2);
            store256d(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYj);
        }
        {
            int j = this->nSpecie()-3;

            __m256d invWv = _mm256_blend_pd(load256d(&invW[j+0]),zerov,0b1000);//10
            __m256d dCjdYj = _mm256_mul_pd(rhoMv,invWv);//10
            __m256d rhoMvj_ = _mm256_blend_pd(load256d(&rhoMByRhoi[j+0]),zerov,0b1000);//11

            __m256d Jci0j03v = _mm256_blend_pd(load256d(&ddNdtByVdcT[(i+0)*(alignN)+j+0]),zerov,0b1000);//12
            __m256d ddNi0dtByVdYj = _mm256_mul_pd(Jci0j03v,dCjdYj);//12
            __m256d r0 = _mm256_mul_pd(rhoMvj_,dPhi0dtv);//13
            __m256d ddYi0dtdYj = _mm256_fmadd_pd(Wi0ByrhoMv,ddNi0dtByVdYj,r0);//12
            store256d(&Jac[(i+0)*(alignN)+j+0],ddYi0dtdYj);//11

            __m256d Jci1j03v = _mm256_blend_pd(load256d(&ddNdtByVdcT[(i+1)*(alignN)+j+0]),zerov,0b1000);
            __m256d ddNi1dtByVdYj = _mm256_mul_pd(Jci1j03v,dCjdYj);
            __m256d r1 = _mm256_mul_pd(rhoMvj_,dPhi1dtv);//13
            __m256d ddYi1dtdYj = _mm256_fmadd_pd(Wi1ByrhoMv,ddNi1dtByVdYj,r1);
            store256d(&Jac[(i+1)*(alignN) + j+0],ddYi1dtdYj);

            __m256d Jci2j03v = _mm256_blend_pd(load256d(&ddNdtByVdcT[(i+2)*(alignN)+j+0]),zerov,0b1000);
            __m256d ddNi2dtByVdYj = _mm256_mul_pd(Jci2j03v,dCjdYj);
            __m256d r2 = _mm256_mul_pd(rhoMvj_,dPhi2dtv);//13
            __m256d ddYi2dtdYj = _mm256_fmadd_pd(Wi2ByrhoMv,ddNi2dtByVdYj,r2);
            store256d(&Jac[(i+2)*(alignN) + j+0],ddYi2dtdYj);
        }
    }
}