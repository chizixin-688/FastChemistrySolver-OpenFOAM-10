/*---------------------------------------------------------------------------*\
  Description
      Computing mass fraction based jacobian matrix using molar concentration 
      based matrix
      Computing partial derivatives of dYdt w.r.t mass fraction
      the tail is 0.  

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
void Foam::FastChemistryModel<UnusedThermo>::ddYdtdY_Vec844_0
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
    for(int j = 0; j < this->nSpecie(); j=j+8)
    {
        __m256d rhoMbyRhoJ03 = load256d(&rhoMByRhoi[j+0]);
        __m256d rhoMbyRhoJ47 = load256d(&rhoMByRhoi[j+4]);
        for(int i = 0; i < this->nSpecie(); i=i+8)
        {
            __m256d rhoMByWi03YTv = _mm256_mul_pd(load256d(&invW[i+0]), load256d(&Phi[i+0]));
            __m256d rhoMByWi47YTv = _mm256_mul_pd(load256d(&invW[i+4]), load256d(&Phi[i+4]));
            rhoMByWi03YTv = _mm256_mul_pd(-rhoMv,rhoMByWi03YTv);
            rhoMByWi47YTv = _mm256_mul_pd(-rhoMv,rhoMByWi47YTv);

            __m256d rhoMByWi0YTv = _mm256_permute4x64_pd(rhoMByWi03YTv, 0x00);
            __m256d rhoMByWi1YTv = _mm256_permute4x64_pd(rhoMByWi03YTv, 0x55);
            __m256d rhoMByWi2YTv = _mm256_permute4x64_pd(rhoMByWi03YTv, 0xAA);
            __m256d rhoMByWi3YTv = _mm256_permute4x64_pd(rhoMByWi03YTv, 0xFF);
            store256d(&dcdY[i*8+0] ,_mm256_mul_pd(rhoMByWi0YTv,rhoMbyRhoJ03));
            store256d(&dcdY[i*8+4] ,_mm256_mul_pd(rhoMByWi0YTv,rhoMbyRhoJ47));
            store256d(&dcdY[i*8+8] ,_mm256_mul_pd(rhoMByWi1YTv,rhoMbyRhoJ03));
            store256d(&dcdY[i*8+12],_mm256_mul_pd(rhoMByWi1YTv,rhoMbyRhoJ47));
            store256d(&dcdY[i*8+16],_mm256_mul_pd(rhoMByWi2YTv,rhoMbyRhoJ03));
            store256d(&dcdY[i*8+20],_mm256_mul_pd(rhoMByWi2YTv,rhoMbyRhoJ47));
            store256d(&dcdY[i*8+24],_mm256_mul_pd(rhoMByWi3YTv,rhoMbyRhoJ03));
            store256d(&dcdY[i*8+28],_mm256_mul_pd(rhoMByWi3YTv,rhoMbyRhoJ47));

            __m256d rhoMByWi4YTv = _mm256_permute4x64_pd(rhoMByWi47YTv, 0x00);
            __m256d rhoMByWi5YTv = _mm256_permute4x64_pd(rhoMByWi47YTv, 0x55);
            __m256d rhoMByWi6YTv = _mm256_permute4x64_pd(rhoMByWi47YTv, 0xAA);
            __m256d rhoMByWi7YTv = _mm256_permute4x64_pd(rhoMByWi47YTv, 0xFF);
            store256d(&dcdY[i*8+32],_mm256_mul_pd(rhoMByWi4YTv,rhoMbyRhoJ03));
            store256d(&dcdY[i*8+36],_mm256_mul_pd(rhoMByWi4YTv,rhoMbyRhoJ47));
            store256d(&dcdY[i*8+40],_mm256_mul_pd(rhoMByWi5YTv,rhoMbyRhoJ03));
            store256d(&dcdY[i*8+44],_mm256_mul_pd(rhoMByWi5YTv,rhoMbyRhoJ47));
            store256d(&dcdY[i*8+48],_mm256_mul_pd(rhoMByWi6YTv,rhoMbyRhoJ03));
            store256d(&dcdY[i*8+52],_mm256_mul_pd(rhoMByWi6YTv,rhoMbyRhoJ47));
            store256d(&dcdY[i*8+56],_mm256_mul_pd(rhoMByWi7YTv,rhoMbyRhoJ03));
            store256d(&dcdY[i*8+60],_mm256_mul_pd(rhoMByWi7YTv,rhoMbyRhoJ47));            
        }


        dcdY[j*8+0] += rhoM*invW[j+0];
        dcdY[j*8+9] += rhoM*invW[j+1];
        dcdY[j*8+18] += rhoM*invW[j+2];
        dcdY[j*8+27] += rhoM*invW[j+3];
        dcdY[j*8+36] += rhoM*invW[j+4];
        dcdY[j*8+45] += rhoM*invW[j+5];
        dcdY[j*8+54] += rhoM*invW[j+6];
        dcdY[j*8+63] += rhoM*invW[j+7];

        for(int i=0; i<this->nSpecie(); i=i+4)
        {
            const double Wi0ByrhoM_ = WiByrhoM[i+0];
            const double Wi1ByrhoM_ = WiByrhoM[i+1];
            const double Wi2ByrhoM_ = WiByrhoM[i+2];
            const double Wi3ByrhoM_ = WiByrhoM[i+3];
            const double dYi0dt = dPhidt[i+0]*Wi0ByrhoM_;
            const double dYi1dt = dPhidt[i+1]*Wi1ByrhoM_;
            const double dYi2dt = dPhidt[i+2]*Wi2ByrhoM_;
            const double dYi3dt = dPhidt[i+3]*Wi3ByrhoM_;
            __m256d ddNi0dtByVdYj03 = _mm256_setzero_pd();
            __m256d ddNi1dtByVdYj03 = _mm256_setzero_pd();
            __m256d ddNi2dtByVdYj03 = _mm256_setzero_pd();
            __m256d ddNi3dtByVdYj03 = _mm256_setzero_pd();
            __m256d ddNi0dtByVdYj47 = _mm256_setzero_pd();
            __m256d ddNi1dtByVdYj47 = _mm256_setzero_pd();
            __m256d ddNi2dtByVdYj47 = _mm256_setzero_pd();
            __m256d ddNi3dtByVdYj47 = _mm256_setzero_pd();            
            const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*this->alignN];
            const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*this->alignN];
            const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*this->alignN];
            const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i+3)*this->alignN];
            for (int k=0; k<this->nSpecie(); k=k+4)
            {
                __m256d ddNi0dtByVdck0v = _mm256_set1_pd(JcRowi0[k+0]);
                __m256d ddNi1dtByVdck0v = _mm256_set1_pd(JcRowi1[k+0]);
                __m256d ddNi2dtByVdck0v = _mm256_set1_pd(JcRowi2[k+0]);
                __m256d ddNi3dtByVdck0v = _mm256_set1_pd(JcRowi3[k+0]);

                __m256d dCk0dYj03v = load256d(&dcdY[k*8]);
                ddNi0dtByVdYj03 = _mm256_fmadd_pd(ddNi0dtByVdck0v,dCk0dYj03v,ddNi0dtByVdYj03);
                ddNi1dtByVdYj03 = _mm256_fmadd_pd(ddNi1dtByVdck0v,dCk0dYj03v,ddNi1dtByVdYj03);
                ddNi2dtByVdYj03 = _mm256_fmadd_pd(ddNi2dtByVdck0v,dCk0dYj03v,ddNi2dtByVdYj03);
                ddNi3dtByVdYj03 = _mm256_fmadd_pd(ddNi3dtByVdck0v,dCk0dYj03v,ddNi3dtByVdYj03);
                __m256d dCk0dYj47v = load256d(&dcdY[k*8+4]);
                ddNi0dtByVdYj47 = _mm256_fmadd_pd(ddNi0dtByVdck0v,dCk0dYj47v,ddNi0dtByVdYj47);
                ddNi1dtByVdYj47 = _mm256_fmadd_pd(ddNi1dtByVdck0v,dCk0dYj47v,ddNi1dtByVdYj47);
                ddNi2dtByVdYj47 = _mm256_fmadd_pd(ddNi2dtByVdck0v,dCk0dYj47v,ddNi2dtByVdYj47);
                ddNi3dtByVdYj47 = _mm256_fmadd_pd(ddNi3dtByVdck0v,dCk0dYj47v,ddNi3dtByVdYj47);
                
                __m256d ddNi0dtByVdck1v = _mm256_set1_pd(JcRowi0[k+1]);
                __m256d ddNi1dtByVdck1v = _mm256_set1_pd(JcRowi1[k+1]);
                __m256d ddNi2dtByVdck1v = _mm256_set1_pd(JcRowi2[k+1]);
                __m256d ddNi3dtByVdck1v = _mm256_set1_pd(JcRowi3[k+1]);
                __m256d dCk1dYj03v = load256d(&dcdY[k*8+8]);     
                ddNi0dtByVdYj03 = _mm256_fmadd_pd(ddNi0dtByVdck1v,dCk1dYj03v,ddNi0dtByVdYj03);
                ddNi1dtByVdYj03 = _mm256_fmadd_pd(ddNi1dtByVdck1v,dCk1dYj03v,ddNi1dtByVdYj03);
                ddNi2dtByVdYj03 = _mm256_fmadd_pd(ddNi2dtByVdck1v,dCk1dYj03v,ddNi2dtByVdYj03);
                ddNi3dtByVdYj03 = _mm256_fmadd_pd(ddNi3dtByVdck1v,dCk1dYj03v,ddNi3dtByVdYj03);
                __m256d dCk1dYj47v = load256d(&dcdY[k*8+12]);     
                ddNi0dtByVdYj47 = _mm256_fmadd_pd(ddNi0dtByVdck1v,dCk1dYj47v,ddNi0dtByVdYj47);
                ddNi1dtByVdYj47 = _mm256_fmadd_pd(ddNi1dtByVdck1v,dCk1dYj47v,ddNi1dtByVdYj47);
                ddNi2dtByVdYj47 = _mm256_fmadd_pd(ddNi2dtByVdck1v,dCk1dYj47v,ddNi2dtByVdYj47);
                ddNi3dtByVdYj47 = _mm256_fmadd_pd(ddNi3dtByVdck1v,dCk1dYj47v,ddNi3dtByVdYj47);

                __m256d ddNi0dtByVdck2v = _mm256_set1_pd(JcRowi0[k+2]);
                __m256d ddNi1dtByVdck2v = _mm256_set1_pd(JcRowi1[k+2]);
                __m256d ddNi2dtByVdck2v = _mm256_set1_pd(JcRowi2[k+2]);
                __m256d ddNi3dtByVdck2v = _mm256_set1_pd(JcRowi3[k+2]);
                __m256d dCk2dYj03v = load256d(&dcdY[k*8+16]);
                ddNi0dtByVdYj03 = _mm256_fmadd_pd(ddNi0dtByVdck2v,dCk2dYj03v,ddNi0dtByVdYj03);
                ddNi1dtByVdYj03 = _mm256_fmadd_pd(ddNi1dtByVdck2v,dCk2dYj03v,ddNi1dtByVdYj03);
                ddNi2dtByVdYj03 = _mm256_fmadd_pd(ddNi2dtByVdck2v,dCk2dYj03v,ddNi2dtByVdYj03);
                ddNi3dtByVdYj03 = _mm256_fmadd_pd(ddNi3dtByVdck2v,dCk2dYj03v,ddNi3dtByVdYj03);
                __m256d dCk2dYj47v = load256d(&dcdY[k*8+20]);
                ddNi0dtByVdYj47 = _mm256_fmadd_pd(ddNi0dtByVdck2v,dCk2dYj47v,ddNi0dtByVdYj47);
                ddNi1dtByVdYj47 = _mm256_fmadd_pd(ddNi1dtByVdck2v,dCk2dYj47v,ddNi1dtByVdYj47);
                ddNi2dtByVdYj47 = _mm256_fmadd_pd(ddNi2dtByVdck2v,dCk2dYj47v,ddNi2dtByVdYj47);
                ddNi3dtByVdYj47 = _mm256_fmadd_pd(ddNi3dtByVdck2v,dCk2dYj47v,ddNi3dtByVdYj47);

                __m256d ddNi0dtByVdck3v = _mm256_set1_pd(JcRowi0[k+3]);
                __m256d ddNi1dtByVdck3v = _mm256_set1_pd(JcRowi1[k+3]);
                __m256d ddNi2dtByVdck3v = _mm256_set1_pd(JcRowi2[k+3]);
                __m256d ddNi3dtByVdck3v = _mm256_set1_pd(JcRowi3[k+3]);  
                __m256d dCk3dYj03v = load256d(&dcdY[k*8+24]);              
                ddNi0dtByVdYj03 = _mm256_fmadd_pd(ddNi0dtByVdck3v,dCk3dYj03v,ddNi0dtByVdYj03);
                ddNi1dtByVdYj03 = _mm256_fmadd_pd(ddNi1dtByVdck3v,dCk3dYj03v,ddNi1dtByVdYj03);
                ddNi2dtByVdYj03 = _mm256_fmadd_pd(ddNi2dtByVdck3v,dCk3dYj03v,ddNi2dtByVdYj03);
                ddNi3dtByVdYj03 = _mm256_fmadd_pd(ddNi3dtByVdck3v,dCk3dYj03v,ddNi3dtByVdYj03);
                __m256d dCk3dYj47v = load256d(&dcdY[k*8+28]);              
                ddNi0dtByVdYj47 = _mm256_fmadd_pd(ddNi0dtByVdck3v,dCk3dYj47v,ddNi0dtByVdYj47);
                ddNi1dtByVdYj47 = _mm256_fmadd_pd(ddNi1dtByVdck3v,dCk3dYj47v,ddNi1dtByVdYj47);
                ddNi2dtByVdYj47 = _mm256_fmadd_pd(ddNi2dtByVdck3v,dCk3dYj47v,ddNi2dtByVdYj47);
                ddNi3dtByVdYj47 = _mm256_fmadd_pd(ddNi3dtByVdck3v,dCk3dYj47v,ddNi3dtByVdYj47);                
            }
             __m256d WiByrhoM_0 = _mm256_set1_pd(Wi0ByrhoM_);
            __m256d dYidtv0 = _mm256_set1_pd(dYi0dt);
            __m256d ddYi0dtdYj03v = _mm256_mul_pd(rhoMbyRhoJ03,dYidtv0);
            ddYi0dtdYj03v = _mm256_fmadd_pd(WiByrhoM_0,ddNi0dtByVdYj03,ddYi0dtdYj03v);
            __m256d ddYi0dtdYj47v = _mm256_mul_pd(rhoMbyRhoJ47,dYidtv0);
            ddYi0dtdYj47v = _mm256_fmadd_pd(WiByrhoM_0,ddNi0dtByVdYj47,ddYi0dtdYj47v);            


            __m256d WiByrhoM_1 = _mm256_set1_pd(Wi1ByrhoM_);
            __m256d dYidtv1 = _mm256_set1_pd(dYi1dt);
            __m256d ddYi1dtdYj03v = _mm256_mul_pd(rhoMbyRhoJ03,dYidtv1);
            ddYi1dtdYj03v = _mm256_fmadd_pd(WiByrhoM_1,ddNi1dtByVdYj03,ddYi1dtdYj03v);
            __m256d ddYi1dtdYj47v = _mm256_mul_pd(rhoMbyRhoJ47,dYidtv1);
            ddYi1dtdYj47v = _mm256_fmadd_pd(WiByrhoM_1,ddNi1dtByVdYj47,ddYi1dtdYj47v);


            __m256d WiByrhoM_2 = _mm256_set1_pd(Wi2ByrhoM_);
            __m256d dYidtv2 = _mm256_set1_pd(dYi2dt);
            __m256d ddYi2dtdYj03v = _mm256_mul_pd(rhoMbyRhoJ03,dYidtv2);
            ddYi2dtdYj03v = _mm256_fmadd_pd(WiByrhoM_2,ddNi2dtByVdYj03,ddYi2dtdYj03v);
            __m256d ddYi2dtdYj47v = _mm256_mul_pd(rhoMbyRhoJ47,dYidtv2);
            ddYi2dtdYj47v = _mm256_fmadd_pd(WiByrhoM_2,ddNi2dtByVdYj47,ddYi2dtdYj47v);

            __m256d WiByrhoM_3 = _mm256_set1_pd(Wi3ByrhoM_);
            __m256d dYidtv3 = _mm256_set1_pd(dYi3dt);
            __m256d ddYi3dtdYj03v = _mm256_mul_pd(rhoMbyRhoJ03,dYidtv3);
            ddYi3dtdYj03v = _mm256_fmadd_pd(WiByrhoM_3,ddNi3dtByVdYj03,ddYi3dtdYj03v);
            __m256d ddYi3dtdYj47v = _mm256_mul_pd(rhoMbyRhoJ47,dYidtv3);
            ddYi3dtdYj47v = _mm256_fmadd_pd(WiByrhoM_3,ddNi3dtByVdYj47,ddYi3dtdYj47v);

            store256d(&Jac[(i+0)*(this->alignN)+ j+0],ddYi0dtdYj03v);
            store256d(&Jac[(i+1)*(this->alignN)+ j+0],ddYi1dtdYj03v);
            store256d(&Jac[(i+2)*(this->alignN)+ j+0],ddYi2dtdYj03v);
            store256d(&Jac[(i+3)*(this->alignN)+ j+0],ddYi3dtdYj03v);
            store256d(&Jac[(i+0)*(this->alignN)+ j+4],ddYi0dtdYj47v);
            store256d(&Jac[(i+1)*(this->alignN)+ j+4],ddYi1dtdYj47v);
            store256d(&Jac[(i+2)*(this->alignN)+ j+4],ddYi2dtdYj47v);
            store256d(&Jac[(i+3)*(this->alignN)+ j+4],ddYi3dtdYj47v);            
        }     
    }     
}