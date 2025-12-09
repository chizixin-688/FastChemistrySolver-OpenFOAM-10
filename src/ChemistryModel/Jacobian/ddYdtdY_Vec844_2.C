#include "FastChemistryModel.H"
#include <immintrin.h>  

template<class ThermoType>
void Foam::FastChemistryModel<ThermoType>::ddYdtdY_Vec844_2
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
    for(int j = 0; j < this->nSpecie()-2; j = j + 8)
    {
        __m256d rhoMbyRhoJ03 = load256d(&rhoMByRhoi[j+0]);
        __m256d rhoMbyRhoJ47 = load256d(&rhoMByRhoi[j+4]);
        for(int i = 0; i < this->nSpecie()-2; i=i+8)
        {
            __m256d rhoMByWi03YTv = _mm256_mul_pd(load256d(&invW[i+0]), load256d(&Phi[i+0]));
            __m256d rhoMByWi47YTv = _mm256_mul_pd(load256d(&invW[i+4]), load256d(&Phi[i+4]));
            rhoMByWi03YTv = _mm256_mul_pd(-rhoMv,rhoMByWi03YTv);
            rhoMByWi47YTv = _mm256_mul_pd(-rhoMv,rhoMByWi47YTv);

            __m256d rhoMByWi0YTv = _mm256_permute4x64_pd(rhoMByWi03YTv, 0x00);
            __m256d rhoMByWi1YTv = _mm256_permute4x64_pd(rhoMByWi03YTv, 0x55);
            __m256d rhoMByWi2YTv = _mm256_permute4x64_pd(rhoMByWi03YTv, 0xAA);
            __m256d rhoMByWi3YTv = _mm256_permute4x64_pd(rhoMByWi03YTv, 0xFF);
            store256d(&dcdY[i*8+0],_mm256_mul_pd(rhoMByWi0YTv,rhoMbyRhoJ03));
            store256d(&dcdY[i*8+4],_mm256_mul_pd(rhoMByWi0YTv,rhoMbyRhoJ47));
            store256d(&dcdY[i*8+8],_mm256_mul_pd(rhoMByWi1YTv,rhoMbyRhoJ03));
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
        // Do the last remain 2 species
        {

            int i = this->nSpecie()-2;

            const double rhoMByWiYTi4 = -rhoM*invW[i+0]*Phi[i+0];
            __m256d r1 = _mm256_mul_pd(_mm256_set1_pd(rhoMByWiYTi4),rhoMbyRhoJ03);
            __m256d r2 = _mm256_mul_pd(_mm256_set1_pd(rhoMByWiYTi4),rhoMbyRhoJ47);
            store256d(&dcdY[i*8+0],r1);
            store256d(&dcdY[i*8+4],r2);
            const double rhoMByWiYTi5 = -rhoM*invW[i+1]*Phi[i+1];
            __m256d r3 = _mm256_mul_pd(_mm256_set1_pd(rhoMByWiYTi5),rhoMbyRhoJ03);
            __m256d r4 = _mm256_mul_pd(_mm256_set1_pd(rhoMByWiYTi5),rhoMbyRhoJ47);
            store256d(&dcdY[i*8+8],r3);
            store256d(&dcdY[i*8+12],r4);            
        }


        dcdY[j*8+0] += rhoM*invW[j+0];
        dcdY[j*8+9] += rhoM*invW[j+1];
        dcdY[j*8+18] += rhoM*invW[j+2];
        dcdY[j*8+27] += rhoM*invW[j+3];
        dcdY[j*8+36] += rhoM*invW[j+4];
        dcdY[j*8+45] += rhoM*invW[j+5];
        dcdY[j*8+54] += rhoM*invW[j+6];
        dcdY[j*8+63] += rhoM*invW[j+7];

        for(int i=0; i<this->nSpecie()-2; i=i+4)
        {
            const double Wi0ByrhoM_ = WiByrhoM[i+0];
            const double Wi1ByrhoM_ = WiByrhoM[i+1];
            const double Wi2ByrhoM_ = WiByrhoM[i+2];
            const double Wi3ByrhoM_ = WiByrhoM[i+3];
            const double dYi0dt = dPhidt[i+0]*Wi0ByrhoM_;
            const double dYi1dt = dPhidt[i+1]*Wi1ByrhoM_;
            const double dYi2dt = dPhidt[i+2]*Wi2ByrhoM_;
            const double dYi3dt = dPhidt[i+3]*Wi3ByrhoM_;
            __m256d ddNi0dtByVdYj03v = _mm256_setzero_pd();
            __m256d ddNi1dtByVdYj03v = _mm256_setzero_pd();
            __m256d ddNi2dtByVdYj03v = _mm256_setzero_pd();
            __m256d ddNi3dtByVdYj03v = _mm256_setzero_pd();
            __m256d ddNi0dtByVdYj47v = _mm256_setzero_pd();
            __m256d ddNi1dtByVdYj47v = _mm256_setzero_pd();
            __m256d ddNi2dtByVdYj47v = _mm256_setzero_pd();
            __m256d ddNi3dtByVdYj47v = _mm256_setzero_pd();
            const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*this->alignN];
            const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*this->alignN];
            const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*this->alignN];
            const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i+3)*this->alignN];
            for (int k=0; k<this->nSpecie()-2; k=k+4)
            {
                __m256d ddNi0dtByVdck0v = _mm256_set1_pd(JcRowi0[k+0]);
                __m256d ddNi1dtByVdck0v = _mm256_set1_pd(JcRowi1[k+0]);
                __m256d ddNi2dtByVdck0v = _mm256_set1_pd(JcRowi2[k+0]);
                __m256d ddNi3dtByVdck0v = _mm256_set1_pd(JcRowi3[k+0]);

                __m256d dCk0dYj03v = load256d(&dcdY[k*8]);
                ddNi0dtByVdYj03v = _mm256_fmadd_pd(ddNi0dtByVdck0v,dCk0dYj03v,ddNi0dtByVdYj03v);
                ddNi1dtByVdYj03v = _mm256_fmadd_pd(ddNi1dtByVdck0v,dCk0dYj03v,ddNi1dtByVdYj03v);
                ddNi2dtByVdYj03v = _mm256_fmadd_pd(ddNi2dtByVdck0v,dCk0dYj03v,ddNi2dtByVdYj03v);
                ddNi3dtByVdYj03v = _mm256_fmadd_pd(ddNi3dtByVdck0v,dCk0dYj03v,ddNi3dtByVdYj03v);
                __m256d dCk0dYj47v = load256d(&dcdY[k*8+4]);
                ddNi0dtByVdYj47v = _mm256_fmadd_pd(ddNi0dtByVdck0v,dCk0dYj47v,ddNi0dtByVdYj47v);
                ddNi1dtByVdYj47v = _mm256_fmadd_pd(ddNi1dtByVdck0v,dCk0dYj47v,ddNi1dtByVdYj47v);
                ddNi2dtByVdYj47v = _mm256_fmadd_pd(ddNi2dtByVdck0v,dCk0dYj47v,ddNi2dtByVdYj47v);
                ddNi3dtByVdYj47v = _mm256_fmadd_pd(ddNi3dtByVdck0v,dCk0dYj47v,ddNi3dtByVdYj47v);

                
                __m256d ddNi0dtByVdck1v = _mm256_set1_pd(JcRowi0[k+1]);
                __m256d ddNi1dtByVdck1v = _mm256_set1_pd(JcRowi1[k+1]);
                __m256d ddNi2dtByVdck1v = _mm256_set1_pd(JcRowi2[k+1]);
                __m256d ddNi3dtByVdck1v = _mm256_set1_pd(JcRowi3[k+1]);
                __m256d dCk1dYj03v = load256d(&dcdY[k*8+8]);     
                ddNi0dtByVdYj03v = _mm256_fmadd_pd(ddNi0dtByVdck1v,dCk1dYj03v,ddNi0dtByVdYj03v);
                ddNi1dtByVdYj03v = _mm256_fmadd_pd(ddNi1dtByVdck1v,dCk1dYj03v,ddNi1dtByVdYj03v);
                ddNi2dtByVdYj03v = _mm256_fmadd_pd(ddNi2dtByVdck1v,dCk1dYj03v,ddNi2dtByVdYj03v);
                ddNi3dtByVdYj03v = _mm256_fmadd_pd(ddNi3dtByVdck1v,dCk1dYj03v,ddNi3dtByVdYj03v);
                __m256d dCk1dYj47v = load256d(&dcdY[k*8+12]);     
                ddNi0dtByVdYj47v = _mm256_fmadd_pd(ddNi0dtByVdck1v,dCk1dYj47v,ddNi0dtByVdYj47v);
                ddNi1dtByVdYj47v = _mm256_fmadd_pd(ddNi1dtByVdck1v,dCk1dYj47v,ddNi1dtByVdYj47v);
                ddNi2dtByVdYj47v = _mm256_fmadd_pd(ddNi2dtByVdck1v,dCk1dYj47v,ddNi2dtByVdYj47v);
                ddNi3dtByVdYj47v = _mm256_fmadd_pd(ddNi3dtByVdck1v,dCk1dYj47v,ddNi3dtByVdYj47v);

                
                __m256d ddNi0dtByVdck2v = _mm256_set1_pd(JcRowi0[k+2]);
                __m256d ddNi1dtByVdck2v = _mm256_set1_pd(JcRowi1[k+2]);
                __m256d ddNi2dtByVdck2v = _mm256_set1_pd(JcRowi2[k+2]);
                __m256d ddNi3dtByVdck2v = _mm256_set1_pd(JcRowi3[k+2]);
                __m256d dCk2dYj03v = load256d(&dcdY[k*8+16]);
                ddNi0dtByVdYj03v = _mm256_fmadd_pd(ddNi0dtByVdck2v,dCk2dYj03v,ddNi0dtByVdYj03v);
                ddNi1dtByVdYj03v = _mm256_fmadd_pd(ddNi1dtByVdck2v,dCk2dYj03v,ddNi1dtByVdYj03v);
                ddNi2dtByVdYj03v = _mm256_fmadd_pd(ddNi2dtByVdck2v,dCk2dYj03v,ddNi2dtByVdYj03v);
                ddNi3dtByVdYj03v = _mm256_fmadd_pd(ddNi3dtByVdck2v,dCk2dYj03v,ddNi3dtByVdYj03v);
                __m256d dCk2dYj47v = load256d(&dcdY[k*8+20]);
                ddNi0dtByVdYj47v = _mm256_fmadd_pd(ddNi0dtByVdck2v,dCk2dYj47v,ddNi0dtByVdYj47v);
                ddNi1dtByVdYj47v = _mm256_fmadd_pd(ddNi1dtByVdck2v,dCk2dYj47v,ddNi1dtByVdYj47v);
                ddNi2dtByVdYj47v = _mm256_fmadd_pd(ddNi2dtByVdck2v,dCk2dYj47v,ddNi2dtByVdYj47v);
                ddNi3dtByVdYj47v = _mm256_fmadd_pd(ddNi3dtByVdck2v,dCk2dYj47v,ddNi3dtByVdYj47v);

                
                __m256d ddNi0dtByVdck3v = _mm256_set1_pd(JcRowi0[k+3]);
                __m256d ddNi1dtByVdck3v = _mm256_set1_pd(JcRowi1[k+3]);
                __m256d ddNi2dtByVdck3v = _mm256_set1_pd(JcRowi2[k+3]);
                __m256d ddNi3dtByVdck3v = _mm256_set1_pd(JcRowi3[k+3]);  
                __m256d dCk3dYj03v = load256d(&dcdY[k*8+24]);              
                ddNi0dtByVdYj03v = _mm256_fmadd_pd(ddNi0dtByVdck3v,dCk3dYj03v,ddNi0dtByVdYj03v);
                ddNi1dtByVdYj03v = _mm256_fmadd_pd(ddNi1dtByVdck3v,dCk3dYj03v,ddNi1dtByVdYj03v);
                ddNi2dtByVdYj03v = _mm256_fmadd_pd(ddNi2dtByVdck3v,dCk3dYj03v,ddNi2dtByVdYj03v);
                ddNi3dtByVdYj03v = _mm256_fmadd_pd(ddNi3dtByVdck3v,dCk3dYj03v,ddNi3dtByVdYj03v);
                __m256d dCk3dYj47v = load256d(&dcdY[k*8+28]);              
                ddNi0dtByVdYj47v = _mm256_fmadd_pd(ddNi0dtByVdck3v,dCk3dYj47v,ddNi0dtByVdYj47v);
                ddNi1dtByVdYj47v = _mm256_fmadd_pd(ddNi1dtByVdck3v,dCk3dYj47v,ddNi1dtByVdYj47v);
                ddNi2dtByVdYj47v = _mm256_fmadd_pd(ddNi2dtByVdck3v,dCk3dYj47v,ddNi2dtByVdYj47v);
                ddNi3dtByVdYj47v = _mm256_fmadd_pd(ddNi3dtByVdck3v,dCk3dYj47v,ddNi3dtByVdYj47v);                
            }
            // do the remain 2 species
            {
                int k = nSpecie()-2; 
                __m256d ddNi0dtByVdck0v = _mm256_set1_pd(JcRowi0[k+0]);
                __m256d ddNi1dtByVdck0v = _mm256_set1_pd(JcRowi1[k+0]);
                __m256d ddNi2dtByVdck0v = _mm256_set1_pd(JcRowi2[k+0]);
                __m256d ddNi3dtByVdck0v = _mm256_set1_pd(JcRowi3[k+0]);

                __m256d dCk0dYj03v = load256d(&dcdY[k*8]);
                ddNi0dtByVdYj03v = _mm256_fmadd_pd(ddNi0dtByVdck0v,dCk0dYj03v,ddNi0dtByVdYj03v);
                ddNi1dtByVdYj03v = _mm256_fmadd_pd(ddNi1dtByVdck0v,dCk0dYj03v,ddNi1dtByVdYj03v);
                ddNi2dtByVdYj03v = _mm256_fmadd_pd(ddNi2dtByVdck0v,dCk0dYj03v,ddNi2dtByVdYj03v);
                ddNi3dtByVdYj03v = _mm256_fmadd_pd(ddNi3dtByVdck0v,dCk0dYj03v,ddNi3dtByVdYj03v);
                __m256d dCk0dYj47v = load256d(&dcdY[k*8+4]);
                ddNi0dtByVdYj47v = _mm256_fmadd_pd(ddNi0dtByVdck0v,dCk0dYj47v,ddNi0dtByVdYj47v);
                ddNi1dtByVdYj47v = _mm256_fmadd_pd(ddNi1dtByVdck0v,dCk0dYj47v,ddNi1dtByVdYj47v);
                ddNi2dtByVdYj47v = _mm256_fmadd_pd(ddNi2dtByVdck0v,dCk0dYj47v,ddNi2dtByVdYj47v);
                ddNi3dtByVdYj47v = _mm256_fmadd_pd(ddNi3dtByVdck0v,dCk0dYj47v,ddNi3dtByVdYj47v);

                
                __m256d ddNi0dtByVdck1v = _mm256_set1_pd(JcRowi0[k+1]);
                __m256d ddNi1dtByVdck1v = _mm256_set1_pd(JcRowi1[k+1]);
                __m256d ddNi2dtByVdck1v = _mm256_set1_pd(JcRowi2[k+1]);
                __m256d ddNi3dtByVdck1v = _mm256_set1_pd(JcRowi3[k+1]);
                __m256d dCk1dYj03v = load256d(&dcdY[k*8+8]);     
                ddNi0dtByVdYj03v = _mm256_fmadd_pd(ddNi0dtByVdck1v,dCk1dYj03v,ddNi0dtByVdYj03v);
                ddNi1dtByVdYj03v = _mm256_fmadd_pd(ddNi1dtByVdck1v,dCk1dYj03v,ddNi1dtByVdYj03v);
                ddNi2dtByVdYj03v = _mm256_fmadd_pd(ddNi2dtByVdck1v,dCk1dYj03v,ddNi2dtByVdYj03v);
                ddNi3dtByVdYj03v = _mm256_fmadd_pd(ddNi3dtByVdck1v,dCk1dYj03v,ddNi3dtByVdYj03v);
                __m256d dCk1dYj47v = load256d(&dcdY[k*8+12]);     
                ddNi0dtByVdYj47v = _mm256_fmadd_pd(ddNi0dtByVdck1v,dCk1dYj47v,ddNi0dtByVdYj47v);
                ddNi1dtByVdYj47v = _mm256_fmadd_pd(ddNi1dtByVdck1v,dCk1dYj47v,ddNi1dtByVdYj47v);
                ddNi2dtByVdYj47v = _mm256_fmadd_pd(ddNi2dtByVdck1v,dCk1dYj47v,ddNi2dtByVdYj47v);
                ddNi3dtByVdYj47v = _mm256_fmadd_pd(ddNi3dtByVdck1v,dCk1dYj47v,ddNi3dtByVdYj47v);
            }

             __m256d WiByrhoM_0 = _mm256_set1_pd(Wi0ByrhoM_);
            __m256d dYidtv0 = _mm256_set1_pd(dYi0dt);
            __m256d ddYi0dtdYj03v = _mm256_mul_pd(rhoMbyRhoJ03,dYidtv0);
            ddYi0dtdYj03v = _mm256_fmadd_pd(WiByrhoM_0,ddNi0dtByVdYj03v,ddYi0dtdYj03v);
            __m256d ddYi0dtdYj47v = _mm256_mul_pd(rhoMbyRhoJ47,dYidtv0);
            ddYi0dtdYj47v = _mm256_fmadd_pd(WiByrhoM_0,ddNi0dtByVdYj47v,ddYi0dtdYj47v);            


            __m256d WiByrhoM_1 = _mm256_set1_pd(Wi1ByrhoM_);
            __m256d dYidtv1 = _mm256_set1_pd(dYi1dt);
            __m256d ddYi1dtdYj03v = _mm256_mul_pd(rhoMbyRhoJ03,dYidtv1);
            ddYi1dtdYj03v = _mm256_fmadd_pd(WiByrhoM_1,ddNi1dtByVdYj03v,ddYi1dtdYj03v);
            __m256d ddYi1dtdYj47v = _mm256_mul_pd(rhoMbyRhoJ47,dYidtv1);
            ddYi1dtdYj47v = _mm256_fmadd_pd(WiByrhoM_1,ddNi1dtByVdYj47v,ddYi1dtdYj47v);


            __m256d WiByrhoM_2 = _mm256_set1_pd(Wi2ByrhoM_);
            __m256d dYidtv2 = _mm256_set1_pd(dYi2dt);
            __m256d ddYi2dtdYj03v = _mm256_mul_pd(rhoMbyRhoJ03,dYidtv2);
            ddYi2dtdYj03v = _mm256_fmadd_pd(WiByrhoM_2,ddNi2dtByVdYj03v,ddYi2dtdYj03v);
            __m256d ddYi2dtdYj47v = _mm256_mul_pd(rhoMbyRhoJ47,dYidtv2);
            ddYi2dtdYj47v = _mm256_fmadd_pd(WiByrhoM_2,ddNi2dtByVdYj47v,ddYi2dtdYj47v);

            __m256d WiByrhoM_3 = _mm256_set1_pd(Wi3ByrhoM_);
            __m256d dYidtv3 = _mm256_set1_pd(dYi3dt);
            __m256d ddYi3dtdYj03v = _mm256_mul_pd(rhoMbyRhoJ03,dYidtv3);
            ddYi3dtdYj03v = _mm256_fmadd_pd(WiByrhoM_3,ddNi3dtByVdYj03v,ddYi3dtdYj03v);
            __m256d ddYi3dtdYj47v = _mm256_mul_pd(rhoMbyRhoJ47,dYidtv3);
            ddYi3dtdYj47v = _mm256_fmadd_pd(WiByrhoM_3,ddNi3dtByVdYj47v,ddYi3dtdYj47v);

            store256d(&Jac[(i+0)*(this->alignN)+ j+0],ddYi0dtdYj03v);
            store256d(&Jac[(i+1)*(this->alignN)+ j+0],ddYi1dtdYj03v);
            store256d(&Jac[(i+2)*(this->alignN)+ j+0],ddYi2dtdYj03v);
            store256d(&Jac[(i+3)*(this->alignN)+ j+0],ddYi3dtdYj03v);
            store256d(&Jac[(i+0)*(this->alignN)+ j+4],ddYi0dtdYj47v);
            store256d(&Jac[(i+1)*(this->alignN)+ j+4],ddYi1dtdYj47v);
            store256d(&Jac[(i+2)*(this->alignN)+ j+4],ddYi2dtdYj47v);
            store256d(&Jac[(i+3)*(this->alignN)+ j+4],ddYi3dtdYj47v);
        }
        // do the remain 2 species
        {
            int i = this->nSpecie()-2;
            const double Wi0ByrhoM_ = WiByrhoM[i+0];
            const double Wi1ByrhoM_ = WiByrhoM[i+1];
            const double dYi0dt = dPhidt[i+0]*Wi0ByrhoM_;
            const double dYi1dt = dPhidt[i+1]*Wi1ByrhoM_;
            __m256d ddNi0dtByVdYj03v = _mm256_setzero_pd();
            __m256d ddNi1dtByVdYj03v = _mm256_setzero_pd();
            __m256d ddNi0dtByVdYj47v = _mm256_setzero_pd();
            __m256d ddNi1dtByVdYj47v = _mm256_setzero_pd();
            const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*this->alignN];
            const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*this->alignN];

            for (int k=0; k<this->nSpecie()-2; k=k+4)
            {
                __m256d ddNi0dtByVdck0v = _mm256_set1_pd(JcRowi0[k+0]);
                __m256d ddNi1dtByVdck0v = _mm256_set1_pd(JcRowi1[k+0]);

                __m256d dCk0dYj03v = load256d(&dcdY[k*8]);
                ddNi0dtByVdYj03v = _mm256_fmadd_pd(ddNi0dtByVdck0v,dCk0dYj03v,ddNi0dtByVdYj03v);
                ddNi1dtByVdYj03v = _mm256_fmadd_pd(ddNi1dtByVdck0v,dCk0dYj03v,ddNi1dtByVdYj03v);

                __m256d dCk0dYj47v = load256d(&dcdY[k*8+4]);
                ddNi0dtByVdYj47v = _mm256_fmadd_pd(ddNi0dtByVdck0v,dCk0dYj47v,ddNi0dtByVdYj47v);
                ddNi1dtByVdYj47v = _mm256_fmadd_pd(ddNi1dtByVdck0v,dCk0dYj47v,ddNi1dtByVdYj47v);
                
                __m256d ddNi0dtByVdck1v = _mm256_set1_pd(JcRowi0[k+1]);
                __m256d ddNi1dtByVdck1v = _mm256_set1_pd(JcRowi1[k+1]);

                __m256d dCk1dYj03v = load256d(&dcdY[k*8+8]);     
                ddNi0dtByVdYj03v = _mm256_fmadd_pd(ddNi0dtByVdck1v,dCk1dYj03v,ddNi0dtByVdYj03v);
                ddNi1dtByVdYj03v = _mm256_fmadd_pd(ddNi1dtByVdck1v,dCk1dYj03v,ddNi1dtByVdYj03v);

                __m256d dCk1dYj47v = load256d(&dcdY[k*8+12]);     
                ddNi0dtByVdYj47v = _mm256_fmadd_pd(ddNi0dtByVdck1v,dCk1dYj47v,ddNi0dtByVdYj47v);
                ddNi1dtByVdYj47v = _mm256_fmadd_pd(ddNi1dtByVdck1v,dCk1dYj47v,ddNi1dtByVdYj47v);
                
                __m256d ddNi0dtByVdck2v = _mm256_set1_pd(JcRowi0[k+2]);
                __m256d ddNi1dtByVdck2v = _mm256_set1_pd(JcRowi1[k+2]);

                __m256d dCk2dYj03v = load256d(&dcdY[k*8+16]);
                ddNi0dtByVdYj03v = _mm256_fmadd_pd(ddNi0dtByVdck2v,dCk2dYj03v,ddNi0dtByVdYj03v);
                ddNi1dtByVdYj03v = _mm256_fmadd_pd(ddNi1dtByVdck2v,dCk2dYj03v,ddNi1dtByVdYj03v);

                __m256d dCk2dYj47v = load256d(&dcdY[k*8+20]);
                ddNi0dtByVdYj47v = _mm256_fmadd_pd(ddNi0dtByVdck2v,dCk2dYj47v,ddNi0dtByVdYj47v);
                ddNi1dtByVdYj47v = _mm256_fmadd_pd(ddNi1dtByVdck2v,dCk2dYj47v,ddNi1dtByVdYj47v);

                
                __m256d ddNi0dtByVdck3v = _mm256_set1_pd(JcRowi0[k+3]);
                __m256d ddNi1dtByVdck3v = _mm256_set1_pd(JcRowi1[k+3]);

                __m256d dCk3dYj03v = load256d(&dcdY[k*8+24]);              
                ddNi0dtByVdYj03v = _mm256_fmadd_pd(ddNi0dtByVdck3v,dCk3dYj03v,ddNi0dtByVdYj03v);
                ddNi1dtByVdYj03v = _mm256_fmadd_pd(ddNi1dtByVdck3v,dCk3dYj03v,ddNi1dtByVdYj03v);

                __m256d dCk3dYj47v = load256d(&dcdY[k*8+28]);              
                ddNi0dtByVdYj47v = _mm256_fmadd_pd(ddNi0dtByVdck3v,dCk3dYj47v,ddNi0dtByVdYj47v);
                ddNi1dtByVdYj47v = _mm256_fmadd_pd(ddNi1dtByVdck3v,dCk3dYj47v,ddNi1dtByVdYj47v);         
            }
            // do the remain 2 species
            {
                int k = nSpecie()-2; 
                __m256d ddNi0dtByVdck0v = _mm256_set1_pd(JcRowi0[k+0]);
                __m256d ddNi1dtByVdck0v = _mm256_set1_pd(JcRowi1[k+0]);

                __m256d dCk0dYj03v = load256d(&dcdY[k*8]);
                ddNi0dtByVdYj03v = _mm256_fmadd_pd(ddNi0dtByVdck0v,dCk0dYj03v,ddNi0dtByVdYj03v);
                ddNi1dtByVdYj03v = _mm256_fmadd_pd(ddNi1dtByVdck0v,dCk0dYj03v,ddNi1dtByVdYj03v);

                __m256d dCk0dYj47v = load256d(&dcdY[k*8+4]);
                ddNi0dtByVdYj47v = _mm256_fmadd_pd(ddNi0dtByVdck0v,dCk0dYj47v,ddNi0dtByVdYj47v);
                ddNi1dtByVdYj47v = _mm256_fmadd_pd(ddNi1dtByVdck0v,dCk0dYj47v,ddNi1dtByVdYj47v);

                __m256d ddNi0dtByVdck1v = _mm256_set1_pd(JcRowi0[k+1]);
                __m256d ddNi1dtByVdck1v = _mm256_set1_pd(JcRowi1[k+1]);

                __m256d dCk1dYj03v = load256d(&dcdY[k*8+8]);     
                ddNi0dtByVdYj03v = _mm256_fmadd_pd(ddNi0dtByVdck1v,dCk1dYj03v,ddNi0dtByVdYj03v);
                ddNi1dtByVdYj03v = _mm256_fmadd_pd(ddNi1dtByVdck1v,dCk1dYj03v,ddNi1dtByVdYj03v);

                __m256d dCk1dYj47v = load256d(&dcdY[k*8+12]);     
                ddNi0dtByVdYj47v = _mm256_fmadd_pd(ddNi0dtByVdck1v,dCk1dYj47v,ddNi0dtByVdYj47v);
                ddNi1dtByVdYj47v = _mm256_fmadd_pd(ddNi1dtByVdck1v,dCk1dYj47v,ddNi1dtByVdYj47v);
            }

             __m256d WiByrhoM_0 = _mm256_set1_pd(Wi0ByrhoM_);
            __m256d dYidtv0 = _mm256_set1_pd(dYi0dt);
            __m256d ddYi0dtdYj03v = _mm256_mul_pd(rhoMbyRhoJ03,dYidtv0);
            ddYi0dtdYj03v = _mm256_fmadd_pd(WiByrhoM_0,ddNi0dtByVdYj03v,ddYi0dtdYj03v);
            __m256d ddYi0dtdYj47v = _mm256_mul_pd(rhoMbyRhoJ47,dYidtv0);
            ddYi0dtdYj47v = _mm256_fmadd_pd(WiByrhoM_0,ddNi0dtByVdYj47v,ddYi0dtdYj47v);            


            __m256d WiByrhoM_1 = _mm256_set1_pd(Wi1ByrhoM_);
            __m256d dYidtv1 = _mm256_set1_pd(dYi1dt);
            __m256d ddYi1dtdYj03v = _mm256_mul_pd(rhoMbyRhoJ03,dYidtv1);
            ddYi1dtdYj03v = _mm256_fmadd_pd(WiByrhoM_1,ddNi1dtByVdYj03v,ddYi1dtdYj03v);
            __m256d ddYi1dtdYj47v = _mm256_mul_pd(rhoMbyRhoJ47,dYidtv1);
            ddYi1dtdYj47v = _mm256_fmadd_pd(WiByrhoM_1,ddNi1dtByVdYj47v,ddYi1dtdYj47v);

            store256d(&Jac[(i+0)*(this->alignN)+ j+0],ddYi0dtdYj03v);
            store256d(&Jac[(i+1)*(this->alignN)+ j+0],ddYi1dtdYj03v);

            store256d(&Jac[(i+0)*(this->alignN)+ j+4],ddYi0dtdYj47v);
            store256d(&Jac[(i+1)*(this->alignN)+ j+4],ddYi1dtdYj47v);
        }
    }
    // Do the last 2 species
    {
        int j0 = this->nSpecie()-2;
        int j1 = this->nSpecie()-1;

        __m128d rhoMvjv = load128d(&rhoMByRhoi[j0]);
        for(int i=0; i<this->nSpecie()-2; i=i+4)
        {
            const double rhoMByWi0 = rhoM*invW[i+0];
            const double rhoMByWi1 = rhoM*invW[i+1];
            const double rhoMByWi2 = rhoM*invW[i+2];
            const double rhoMByWi3 = rhoM*invW[i+3];

            __m128d Arr0 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-rhoMByWi0),rhoMvjv),_mm_set1_pd(Phi[i+0]));
            store128d(&dcdY[i*4+0],Arr0);
            __m128d Arr1 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-rhoMByWi1),rhoMvjv),_mm_set1_pd(Phi[i+1]));
            store128d(&dcdY[i*4+4],Arr1);
            __m128d Arr2 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-rhoMByWi2),rhoMvjv),_mm_set1_pd(Phi[i+2]));
            store128d(&dcdY[i*4+8],Arr2); 
            __m128d Arr3 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-rhoMByWi3),rhoMvjv),_mm_set1_pd(Phi[i+3]));
            store128d(&dcdY[i*4+12],Arr3);                
        }
        {
            int i = this->nSpecie()-2;
            const double rhoMByWi0 = rhoM*invW[i+0];
            const double rhoMByWi1 = rhoM*invW[i+1];
            __m128d Arr0 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-rhoMByWi0),rhoMvjv),_mm_set1_pd(Phi[i+0]));
            store128d(&dcdY[i*4+0],Arr0);  
            __m128d Arr1 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-rhoMByWi1),rhoMvjv),_mm_set1_pd(Phi[i+1]));
            store128d(&dcdY[i*4+4],Arr1);             
        }
        dcdY[j0*4+0] = dcdY[j0*4+0] + rhoM*invW[j0];
        dcdY[j1*4+1] = dcdY[j1*4+1] + rhoM*invW[j1];

        for(int i=0; i<this->nSpecie()-2; i=i+4)
        {
            const double WiByrhoM_0 = WiByrhoM[i+0];
            const double WiByrhoM_1 = WiByrhoM[i+1];
            const double WiByrhoM_2 = WiByrhoM[i+2];
            const double WiByrhoM_3 = WiByrhoM[i+3];

            double dYidt0 = dPhidt[i+0]*WiByrhoM_0;
            double dYidt1 = dPhidt[i+1]*WiByrhoM_1;
            double dYidt2 = dPhidt[i+2]*WiByrhoM_2;
            double dYidt3 = dPhidt[i+3]*WiByrhoM_3;

            const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*(this->alignN)];
            const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*(this->alignN)];
            const double* __restrict__ JcRowi2 = &ddNdtByVdcT[(i+2)*(this->alignN)];
            const double* __restrict__ JcRowi3 = &ddNdtByVdcT[(i+3)*(this->alignN)];

            __m128d ddNi0dtByVdYjv = _mm_setzero_pd();
            __m128d ddNi1dtByVdYjv = _mm_setzero_pd();
            __m128d ddNi2dtByVdYjv = _mm_setzero_pd();
            __m128d ddNi3dtByVdYjv = _mm_setzero_pd();
            for (int k=0; k<this->nSpecie()-2; k=k+4)
            {
                const double ddNi0dtByVdck0 = JcRowi0[k+0];
                const double ddNi0dtByVdck1 = JcRowi0[k+1];
                const double ddNi0dtByVdck2 = JcRowi0[k+2];
                const double ddNi0dtByVdck3 = JcRowi0[k+3];

                const double ddNi1dtByVdck0 = JcRowi1[k+0];
                const double ddNi1dtByVdck1 = JcRowi1[k+1];
                const double ddNi1dtByVdck2 = JcRowi1[k+2];
                const double ddNi1dtByVdck3 = JcRowi1[k+3];

                const double ddNi2dtByVdck0 = JcRowi2[k+0];
                const double ddNi2dtByVdck1 = JcRowi2[k+1];
                const double ddNi2dtByVdck2 = JcRowi2[k+2];
                const double ddNi2dtByVdck3 = JcRowi2[k+3];

                const double ddNi3dtByVdck0 = JcRowi3[k+0];
                const double ddNi3dtByVdck1 = JcRowi3[k+1];
                const double ddNi3dtByVdck2 = JcRowi3[k+2];
                const double ddNi3dtByVdck3 = JcRowi3[k+3];


                __m128d dCk0dYj = load128d(&dcdY[k*4+0]);
                __m128d dCk1dYj = load128d(&dcdY[k*4+4]);
                __m128d dCk2dYj = load128d(&dcdY[k*4+8]);
                __m128d dCk3dYj = load128d(&dcdY[k*4+12]);

                ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck0),dCk0dYj,ddNi0dtByVdYjv);
                ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck1),dCk1dYj,ddNi0dtByVdYjv);
                ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck2),dCk2dYj,ddNi0dtByVdYjv);
                ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck3),dCk3dYj,ddNi0dtByVdYjv);

                ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck0),dCk0dYj,ddNi1dtByVdYjv);
                ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck1),dCk1dYj,ddNi1dtByVdYjv);
                ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck2),dCk2dYj,ddNi1dtByVdYjv);
                ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck3),dCk3dYj,ddNi1dtByVdYjv);   
                
                ddNi2dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi2dtByVdck0),dCk0dYj,ddNi2dtByVdYjv);
                ddNi2dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi2dtByVdck1),dCk1dYj,ddNi2dtByVdYjv);
                ddNi2dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi2dtByVdck2),dCk2dYj,ddNi2dtByVdYjv);
                ddNi2dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi2dtByVdck3),dCk3dYj,ddNi2dtByVdYjv);

                ddNi3dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi3dtByVdck0),dCk0dYj,ddNi3dtByVdYjv);
                ddNi3dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi3dtByVdck1),dCk1dYj,ddNi3dtByVdYjv);
                ddNi3dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi3dtByVdck2),dCk2dYj,ddNi3dtByVdYjv);
                ddNi3dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi3dtByVdck3),dCk3dYj,ddNi3dtByVdYjv);            
            }
            {
                int k = this->nSpecie()-2;
                const double ddNi0dtByVdck0 = JcRowi0[k+0];
                const double ddNi0dtByVdck1 = JcRowi0[k+1];

                const double ddNi1dtByVdck0 = JcRowi1[k+0];
                const double ddNi1dtByVdck1 = JcRowi1[k+1];  

                const double ddNi2dtByVdck0 = JcRowi2[k+0];
                const double ddNi2dtByVdck1 = JcRowi2[k+1];  

                const double ddNi3dtByVdck0 = JcRowi3[k+0];
                const double ddNi3dtByVdck1 = JcRowi3[k+1];
                __m128d dCk0dYj0 = load128d(&dcdY[k*4+0]);
                __m128d dCk1dYj0 = load128d(&dcdY[k*4+4]);
                ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck0),dCk0dYj0,ddNi0dtByVdYjv);
                ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck1),dCk1dYj0,ddNi0dtByVdYjv);

                ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck0),dCk0dYj0,ddNi1dtByVdYjv);
                ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck1),dCk1dYj0,ddNi1dtByVdYjv);
                
                ddNi2dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi2dtByVdck0),dCk0dYj0,ddNi2dtByVdYjv);
                ddNi2dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi2dtByVdck1),dCk1dYj0,ddNi2dtByVdYjv);
                
                ddNi3dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi3dtByVdck0),dCk0dYj0,ddNi3dtByVdYjv);
                ddNi3dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi3dtByVdck1),dCk1dYj0,ddNi3dtByVdYjv);
            }
            __m128d Wi0ByrhoMv = _mm_set1_pd(WiByrhoM_0);
            __m128d Wi1ByrhoMv = _mm_set1_pd(WiByrhoM_1);
            __m128d Wi2ByrhoMv = _mm_set1_pd(WiByrhoM_2);
            __m128d Wi3ByrhoMv = _mm_set1_pd(WiByrhoM_3);

            __m128d rhoMByrhoi = load128d(&rhoMByRhoi[j0]);

            __m128d dYi0dtv = _mm_set1_pd(dYidt0);
            __m128d dYi1dtv = _mm_set1_pd(dYidt1);
            __m128d dYi2dtv = _mm_set1_pd(dYidt2);
            __m128d dYi3dtv = _mm_set1_pd(dYidt3);

            __m128d result0 = _mm_mul_pd(Wi0ByrhoMv,ddNi0dtByVdYjv);
            result0 = _mm_fmadd_pd(rhoMByrhoi,dYi0dtv,result0);
            store128d(&Jac[(i+0)*(this->alignN) + j0],result0);

            __m128d result1 = _mm_mul_pd(Wi1ByrhoMv,ddNi1dtByVdYjv);
            result1 = _mm_fmadd_pd(rhoMByrhoi,dYi1dtv,result1);
            store128d(&Jac[(i+1)*(this->alignN) + j0],result1);

            __m128d result2 = _mm_mul_pd(Wi2ByrhoMv,ddNi2dtByVdYjv);
            result2 = _mm_fmadd_pd(rhoMByrhoi,dYi2dtv,result2);
            store128d(&Jac[(i+2)*(this->alignN) + j0],result2);

            __m128d result3 = _mm_mul_pd(Wi3ByrhoMv,ddNi3dtByVdYjv);
            result3 = _mm_fmadd_pd(rhoMByrhoi,dYi3dtv,result3);
            store128d(&Jac[(i+3)*(this->alignN) + j0],result3);               
        }

        {
            int i = nSpecie() - 2;

            const double WiByrhoM_0 = WiByrhoM[i+0];
            const double WiByrhoM_1 = WiByrhoM[i+1];

            double dYidt0 = dPhidt[i+0]*WiByrhoM_0;
            double dYidt1 = dPhidt[i+1]*WiByrhoM_1;

            const double* __restrict__ JcRowi0 = &ddNdtByVdcT[(i+0)*(this->alignN)];
            const double* __restrict__ JcRowi1 = &ddNdtByVdcT[(i+1)*(this->alignN)];
            __m128d ddNi0dtByVdYjv = _mm_setzero_pd();
            __m128d ddNi1dtByVdYjv = _mm_setzero_pd();
            for (int k=0; k<this->nSpecie()-2; k=k+4)
            {
                const double ddNi0dtByVdck0 = JcRowi0[k+0];
                const double ddNi0dtByVdck1 = JcRowi0[k+1];
                const double ddNi0dtByVdck2 = JcRowi0[k+2];
                const double ddNi0dtByVdck3 = JcRowi0[k+3];
                const double ddNi1dtByVdck0 = JcRowi1[k+0];
                const double ddNi1dtByVdck1 = JcRowi1[k+1];
                const double ddNi1dtByVdck2 = JcRowi1[k+2];
                const double ddNi1dtByVdck3 = JcRowi1[k+3]; 

                __m128d dCk0dYj = load128d(&dcdY[k*4+0]);
                __m128d dCk1dYj = load128d(&dcdY[k*4+4]);
                __m128d dCk2dYj = load128d(&dcdY[k*4+8]);
                __m128d dCk3dYj = load128d(&dcdY[k*4+12]);

                ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck0),dCk0dYj,ddNi0dtByVdYjv);
                ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck1),dCk1dYj,ddNi0dtByVdYjv);
                ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck2),dCk2dYj,ddNi0dtByVdYjv);
                ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck3),dCk3dYj,ddNi0dtByVdYjv);

                ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck0),dCk0dYj,ddNi1dtByVdYjv);
                ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck1),dCk1dYj,ddNi1dtByVdYjv);
                ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck2),dCk2dYj,ddNi1dtByVdYjv);
                ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck3),dCk3dYj,ddNi1dtByVdYjv);            
            }
            {
                int k = this->nSpecie()-2;
                const double ddNi0dtByVdck0 = JcRowi0[k+0];
                const double ddNi0dtByVdck1 = JcRowi0[k+1];

                const double ddNi1dtByVdck0 = JcRowi1[k+0];
                const double ddNi1dtByVdck1 = JcRowi1[k+1];
                __m128d dCk0dYj0 = load128d(&dcdY[k*4+0]);
                __m128d dCk1dYj0 = load128d(&dcdY[k*4+4]);
                ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck0),dCk0dYj0,ddNi0dtByVdYjv);
                ddNi0dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi0dtByVdck1),dCk1dYj0,ddNi0dtByVdYjv);

                ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck0),dCk0dYj0,ddNi1dtByVdYjv);
                ddNi1dtByVdYjv = _mm_fmadd_pd(_mm_set1_pd(ddNi1dtByVdck1),dCk1dYj0,ddNi1dtByVdYjv);            
            }
            __m128d Wi0ByrhoMv = _mm_set1_pd(WiByrhoM_0);
            __m128d Wi1ByrhoMv = _mm_set1_pd(WiByrhoM_1);
            __m128d rhiMvjv = load128d(&rhoMByRhoi[j0]);
            __m128d dYi0dtv = _mm_set1_pd(dYidt0);
            __m128d dYi1dtv = _mm_set1_pd(dYidt1);
            __m128d result0 = _mm_mul_pd(Wi0ByrhoMv,ddNi0dtByVdYjv);
            __m128d result1 = _mm_mul_pd(Wi1ByrhoMv,ddNi1dtByVdYjv);
            result0 = _mm_fmadd_pd(rhiMvjv,dYi0dtv,result0);
            result1 = _mm_fmadd_pd(rhiMvjv,dYi1dtv,result1);
            store128d(&Jac[(i+0)*(this->alignN) + j0],result0);
            store128d(&Jac[(i+1)*(this->alignN) + j0],result1);
        }  
    }

}