/*---------------------------------------------------------------------------*\
  Description
      Implementation of the idealGas class.

      All thermo kernels (JacobianThermoYT / JacobianThermoNT /
      DerivativeThermoYT / DerivativeThermoNT / negGstdByRT) are written with
      AVX-2 (256-bit) SIMD intrinsics and operate on 4 species per group.
      The coefficients are read from the SoA (column-major) storage filled by
      createIdealGasFromFoamDict. For every call the temperature selects one
      of three regimes (computed once at the top of each kernel):
        - T <= TcommonMin  -> all species use the low-temperature fit (pure L)
        - T >= TcommonMax  -> all species use the high-temperature fit (pure H)
        - otherwise        -> per-species selection (T < Tcommon[i]) via
                               vblendvpd (mixed)
      The trailing nSpecies%%4 species are handled with masked loads / masked
      stores (load256d0/01/012, store256d0/01/012) so the padded lanes are
      neither read (zero-filled) nor written.

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. Standard C++ library headers
//---------------------------------
#include <cstring>

//---------------------------------
// 2. FastChemistry headers
//---------------------------------
#include "idealGas.H"
#include "vec_math_avx2.H"
//=============================================================================//

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

FastChemistry::idealGas::idealGas(const int n)
:
    nSpecies(n),
    HCoeffs(n),
    LCoeffs(n),
    Tlow(n),
    Thigh(n),
    Tcommon(n),
    TcommonMin(1e100),
    TcommonMax(0),
    TlowMin(1e100),
    TlowMax(1),
    ThighMin(1e100),
    ThighMax(1),
    bufferSoA(nullptr),
    buffer(nullptr),
    W(nullptr),
    invW(nullptr),
    Hf(nullptr),
    logT(0),
    invT(0),
    sqrT(0),
    logP(0),
    rhoM(0),
    vM(0),
    Mw(0)
{
    int alignSpecies = ((nSpecies+3)/4)*4;
    if (posix_memalign(reinterpret_cast<void**>(&this->buffer), 32, alignSpecies*4*sizeof(double)))
    {
        throw std::bad_alloc();
    }
    memset(this->buffer, 0, alignSpecies*4*sizeof(double));
    this->W = &this->buffer[0];
    this->invW = &this->buffer[alignSpecies*1];
    this->Hf = &this->buffer[alignSpecies*2];
    this->invWRu = &this->buffer[alignSpecies*3];

    // SoA coefficient storage: 15 segments of alignSpecies64 doubles.
    // Each segment is a multiple of 8 doubles (64 bytes), so every segment
    // boundary stays 64-byte aligned with the 64-byte aligned base pointer.
    // Layout: [L0..L6 | H0..H6 | Tcommon]
    const int alignSpecies64 = ((nSpecies+7)/8)*8;
    if (posix_memalign(reinterpret_cast<void**>(&this->bufferSoA), 64, alignSpecies64*15*sizeof(double)))
    {
        throw std::bad_alloc();
    }
    memset(this->bufferSoA, 0, alignSpecies64*15*sizeof(double));
    for(int j = 0; j < 7; j++)
    {
        this->LCoeffsSoA[j] = &this->bufferSoA[j*alignSpecies64];
        this->HCoeffsSoA[j] = &this->bufferSoA[(7+j)*alignSpecies64];
    }
    this->TcommonSoA = &this->bufferSoA[14*alignSpecies64];

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FastChemistry::idealGas::~idealGas()
{
    this->W=nullptr;
    this->invW=nullptr;
    this->Hf=nullptr;
    this->invWRu=nullptr;
    if(this->buffer!=nullptr)
    {
        free(this->buffer);
    }
    this->LCoeffsSoA[0]=nullptr;
    this->HCoeffsSoA[0]=nullptr;
    this->TcommonSoA=nullptr;
    if(this->bufferSoA!=nullptr)
    {
        free(this->bufferSoA);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void FastChemistry::idealGas::JacobianThermoYT
(
    double p,
    double T_,
    double* __restrict__ Phi,
    double* __restrict__ concentration,
    double* __restrict__ negGstdByRT,
    double* __restrict__ negGstdByRT2,
    double* __restrict__ dBdT,
    double* __restrict__ dCpdT,
    double* __restrict__ Cp,
    double* __restrict__ Ha,
    double* __restrict__ rhoMvj,
    double* __restrict__ WiByrhoM
)const 
{
    this->rhoM = 0;
    double MW = 0;    

    unsigned int remain = this->nSpecies%4;
    const int regime = (T_ <= this->TcommonMin) ? 0 : (T_ >= this->TcommonMax) ? 1 : 2;
    __m256d MWv = _mm256_setzero_pd();
    __m256d rhoMv = _mm256_setzero_pd();
    __m256d vT = _mm256_set1_pd(T_);
    __m256d vInvT = _mm256_set1_pd(this->invT);
    __m256d vlogT1 = _mm256_set1_pd(this->logT-1);
    for(unsigned int i = 0; i < this->nSpecies-remain;i=i+4)
    {
        unsigned int i0 = i+0;

        __m256d invWv = load256d(&invW[i0]);
        __m256d YTpv = load256d(&Phi[i0]);

        rhoMv = _mm256_add_pd(_mm256_mul_pd(YTpv,invWv),rhoMv);
        MWv = _mm256_fmadd_pd(YTpv,invWv,MWv);

        __m256d A0,A1,A2,A3,A4,A5,A6;

        if (regime == 2)
        {
            const __m256d mask = _mm256_cmp_pd(vT,load256d(&this->TcommonSoA[i]),_CMP_LT_OQ);
            A0 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[0][i]), load256d(&this->LCoeffsSoA[0][i]), mask);
            A1 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[1][i]), load256d(&this->LCoeffsSoA[1][i]), mask);
            A2 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[2][i]), load256d(&this->LCoeffsSoA[2][i]), mask);
            A3 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[3][i]), load256d(&this->LCoeffsSoA[3][i]), mask);
            A4 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[4][i]), load256d(&this->LCoeffsSoA[4][i]), mask);
            A5 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[5][i]), load256d(&this->LCoeffsSoA[5][i]), mask);
            A6 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[6][i]), load256d(&this->LCoeffsSoA[6][i]), mask);
        }
        else if(regime == 1)
        {
            A0 = load256d(&this->HCoeffsSoA[0][i]);
            A1 = load256d(&this->HCoeffsSoA[1][i]);
            A2 = load256d(&this->HCoeffsSoA[2][i]);
            A3 = load256d(&this->HCoeffsSoA[3][i]);
            A4 = load256d(&this->HCoeffsSoA[4][i]);
            A5 = load256d(&this->HCoeffsSoA[5][i]);
            A6 = load256d(&this->HCoeffsSoA[6][i]);
        }
        else
        {
            A0 = load256d(&this->LCoeffsSoA[0][i]);
            A1 = load256d(&this->LCoeffsSoA[1][i]);
            A2 = load256d(&this->LCoeffsSoA[2][i]);
            A3 = load256d(&this->LCoeffsSoA[3][i]);
            A4 = load256d(&this->LCoeffsSoA[4][i]);
            A5 = load256d(&this->LCoeffsSoA[5][i]);
            A6 = load256d(&this->LCoeffsSoA[6][i]);
        }

        __m256d vdBdT = _mm256_fmadd_pd(A4*0.2 ,vT,A3*0.25);
        vdBdT = _mm256_fmadd_pd(vdBdT,vT,A2*(1.0/3.0));      
        vdBdT = _mm256_fmadd_pd(vdBdT,vT,A1*0.5);
        vdBdT = _mm256_fmadd_pd(_mm256_fmadd_pd(A5,vInvT,A0),vInvT,vdBdT);
        store256d(&dBdT[i0],vdBdT);

        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,vT,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,vT,_mm256_mul_pd(A5,vInvT));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,vlogT1,vExpNegGstdByRT);
        store256d(&negGstdByRT[i0],vExpNegGstdByRT);
        store256d(&negGstdByRT2[i0],vExpNegGstdByRT);

        //__m256d RuInvW = _mm256_mul_pd(invWv,Ruv);
        __m256d RuInvWv = load256d(&this->invWRu[i0]);
        __m256d vdCpdT = _mm256_fmadd_pd(vT,4*A4,3*A3);
        vdCpdT = _mm256_fmadd_pd(vT,vdCpdT,2*A2);  
        vdCpdT = _mm256_fmadd_pd(vT,vdCpdT,A1);  
        vdCpdT = _mm256_mul_pd(RuInvWv,vdCpdT);
        store256d(&dCpdT[i0],vdCpdT);   

        __m256d vCp = _mm256_fmadd_pd(A4 ,vT,A3);
        vCp = _mm256_fmadd_pd(vCp,vT,A2);
        vCp = _mm256_fmadd_pd(vCp,vT,A1);
        vCp = _mm256_fmadd_pd(vCp,vT,A0);
        vCp = _mm256_mul_pd(RuInvWv,vCp);
        store256d(&Cp[i0],vCp);   

        __m256d vHa = _mm256_fmadd_pd(A4 ,vT*0.2,A3*0.25);
        vHa = _mm256_fmadd_pd(vHa,vT,A2*(1.0/3.0));
        vHa = _mm256_fmadd_pd(vHa,vT,A1*0.5);
        vHa = _mm256_fmadd_pd(vHa,vT,A0);
        vHa = _mm256_fmadd_pd(vHa,vT,A5);
        vHa = _mm256_mul_pd(RuInvWv,vHa);
        store256d(&Ha[i0],vHa);   

    }
    if(remain)
    {
        unsigned int i = this->nSpecies-remain;
        __m256d invWv = load256d(&invW[i]);
        __m256d YTpv = load256d(&Phi[i]);

        rhoMv = _mm256_add_pd(_mm256_mul_pd(YTpv,invWv),rhoMv);
        MWv = _mm256_fmadd_pd(YTpv,invWv,MWv);

        __m256d A0,A1,A2,A3,A4,A5,A6;

        if (regime == 2)
        {
            const __m256d mask = _mm256_cmp_pd(vT,load256d(&this->TcommonSoA[i]),_CMP_LT_OQ);
            A0 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[0][i]), load256d(&this->LCoeffsSoA[0][i]), mask);
            A1 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[1][i]), load256d(&this->LCoeffsSoA[1][i]), mask);
            A2 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[2][i]), load256d(&this->LCoeffsSoA[2][i]), mask);
            A3 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[3][i]), load256d(&this->LCoeffsSoA[3][i]), mask);
            A4 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[4][i]), load256d(&this->LCoeffsSoA[4][i]), mask);
            A5 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[5][i]), load256d(&this->LCoeffsSoA[5][i]), mask);
            A6 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[6][i]), load256d(&this->LCoeffsSoA[6][i]), mask);
        }
        else if(regime == 1)
        {
            A0 = load256d(&this->HCoeffsSoA[0][i]);
            A1 = load256d(&this->HCoeffsSoA[1][i]);
            A2 = load256d(&this->HCoeffsSoA[2][i]);
            A3 = load256d(&this->HCoeffsSoA[3][i]);
            A4 = load256d(&this->HCoeffsSoA[4][i]);
            A5 = load256d(&this->HCoeffsSoA[5][i]);
            A6 = load256d(&this->HCoeffsSoA[6][i]);
        }
        else
        {
            A0 = load256d(&this->LCoeffsSoA[0][i]);
            A1 = load256d(&this->LCoeffsSoA[1][i]);
            A2 = load256d(&this->LCoeffsSoA[2][i]);
            A3 = load256d(&this->LCoeffsSoA[3][i]);
            A4 = load256d(&this->LCoeffsSoA[4][i]);
            A5 = load256d(&this->LCoeffsSoA[5][i]);
            A6 = load256d(&this->LCoeffsSoA[6][i]);
        }
        __m256d vdBdT = _mm256_fmadd_pd(A4*0.2 ,vT,A3*0.25);
        vdBdT = _mm256_fmadd_pd(vdBdT,vT,A2*(1.0/3.0));      
        vdBdT = _mm256_fmadd_pd(vdBdT,vT,A1*0.5);
        vdBdT = _mm256_fmadd_pd(_mm256_fmadd_pd(A5,vInvT,A0),vInvT,vdBdT);

        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,vT,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,vT,_mm256_mul_pd(A5,vInvT));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,vlogT1,vExpNegGstdByRT);

        //__m256d RuInvW = _mm256_mul_pd(invWv,Ruv);
        __m256d RuInvWv = load256d(&this->invWRu[i]);
        __m256d vdCpdT = _mm256_fmadd_pd(vT,4*A4,3*A3);
        vdCpdT = _mm256_fmadd_pd(vT,vdCpdT,2*A2);  
        vdCpdT = _mm256_fmadd_pd(vT,vdCpdT,A1);  
        vdCpdT = _mm256_mul_pd(RuInvWv,vdCpdT);

        __m256d vCp = _mm256_fmadd_pd(A4 ,vT,A3);
        vCp = _mm256_fmadd_pd(vCp,vT,A2);
        vCp = _mm256_fmadd_pd(vCp,vT,A1);
        vCp = _mm256_fmadd_pd(vCp,vT,A0);
        vCp = _mm256_mul_pd(RuInvWv,vCp);

        __m256d vHa = _mm256_fmadd_pd(A4 ,vT*0.2,A3*0.25);
        vHa = _mm256_fmadd_pd(vHa,vT,A2*(1.0/3.0));
        vHa = _mm256_fmadd_pd(vHa,vT,A1*0.5);
        vHa = _mm256_fmadd_pd(vHa,vT,A0);
        vHa = _mm256_fmadd_pd(vHa,vT,A5);
        vHa = _mm256_mul_pd(RuInvWv,vHa);

        if(remain==1)
        {
            store256d0(&dBdT[i],vdBdT);
            store256d0(&negGstdByRT[i],vExpNegGstdByRT);
            store256d0(&negGstdByRT2[i],vExpNegGstdByRT);
            store256d0(&dCpdT[i],vdCpdT);
            store256d0(&Cp[i],vCp);
            store256d0(&Ha[i],vHa);
        }
        else if(remain==2)
        {
            store256d01(&dBdT[i],vdBdT);
            store256d01(&negGstdByRT[i],vExpNegGstdByRT);
            store256d01(&negGstdByRT2[i],vExpNegGstdByRT);
            store256d01(&dCpdT[i],vdCpdT);
            store256d01(&Cp[i],vCp);
            store256d01(&Ha[i],vHa);
        }
        else if(remain==3)
        {
            store256d012(&dBdT[i],vdBdT);
            store256d012(&negGstdByRT[i],vExpNegGstdByRT);
            store256d012(&negGstdByRT2[i],vExpNegGstdByRT);
            store256d012(&dCpdT[i],vdCpdT);
            store256d012(&Cp[i],vCp);
            store256d012(&Ha[i],vHa);
        }
    }

    MW = MW + hsum4(MWv);
    MW = 1/MW;
    dBdT[this->nSpecies] = MW;
    this->rhoM = this->rhoM + hsum4(rhoMv);
    this->rhoM = this->rhoM*(this->Ru*T_)/p;

    const double invRhoM = this->rhoM;
    this->rhoM = 1/this->rhoM;

    __m256d ArrCpM_ = _mm256_setzero_pd();
    __m256d ArrdCpMdT_ = _mm256_setzero_pd() ;
    __m256d MWvv = _mm256_set1_pd(MW);
    __m256d invrhoMv = _mm256_set1_pd(invRhoM);
    __m256d rhoMvv = _mm256_set1_pd(this->rhoM);
    __m256d zerov = _mm256_setzero_pd();
    for(unsigned int i = 0; i < this->nSpecies-remain;i=i+4)
    {
        __m256d invW_ = load256d(&this->invW[i]);
        __m256d rhoMvj_ = _mm256_mul_pd(MWvv,invW_);
        store256d(&rhoMvj[i],rhoMvj_);
        __m256d Wv= load256d(&this->W[i]);
        __m256d WiByrhoM_ = _mm256_mul_pd(invrhoMv,Wv);
        store256d(&WiByrhoM[i],WiByrhoM_);
        __m256d Phi03v = load256d(&Phi[i]);

        __m256d Cv = _mm256_mul_pd(_mm256_mul_pd(rhoMvv,Phi03v),invW_);
        Cv = _mm256_max_pd(Cv,zerov);
        store256d(&concentration[i],Cv);

        __m256d Cp_ = load256d(&Cp[i]);
        ArrCpM_ = _mm256_fmadd_pd(Phi03v,Cp_,ArrCpM_);
        __m256d dCpdT_ = load256d(&dCpdT[i]);
        ArrdCpMdT_ = _mm256_fmadd_pd(Phi03v,dCpdT_,ArrdCpMdT_);
    }

    if(remain==1)
    {
        unsigned int i = this->nSpecies-1;

        __m256d invW_ = load256d0(&this->invW[i]);
        __m256d rhoMvj_ = mul256d(MWvv,invW_);
        store256d0(&rhoMvj[i],rhoMvj_);
        __m256d Wv= load256d0(&this->W[i]);
        __m256d WiByrhoM_ = mul256d(invrhoMv,Wv);
        store256d0(&WiByrhoM[i],WiByrhoM_);
        __m256d Phi03v = load256d0(&Phi[i]);

        __m256d Cv = mul256d(mul256d(rhoMvv,Phi03v),invW_);
        Cv = _mm256_max_pd(Cv,zerov);
        store256d0(&concentration[i],Cv);

        __m256d Cp_ = load256d0(&Cp[i]);
        ArrCpM_ = fmadd256d(Phi03v,Cp_,ArrCpM_);
        __m256d dCpdT_ = load256d0(&dCpdT[i]);
        ArrdCpMdT_ = fmadd256d(Phi03v,dCpdT_,ArrdCpMdT_);
    }
    else if(remain==2)
    {
        unsigned int i = this->nSpecies-2;
        __m256d invW_ = load256d01(&this->invW[i]);
        __m256d rhoMvj_ = mul256d(MWvv,invW_);
        store256d01(&rhoMvj[i],rhoMvj_);
        __m256d Wv= load256d01(&this->W[i]);
        __m256d WiByrhoM_ = mul256d(invrhoMv,Wv);
        store256d01(&WiByrhoM[i],WiByrhoM_);
        __m256d Phi03v = load256d01(&Phi[i]);

        __m256d Cv = mul256d(mul256d(rhoMvv,Phi03v),invW_);
        Cv = _mm256_max_pd(Cv,zerov);
        store256d01(&concentration[i],Cv);

        __m256d Cp_ = load256d01(&Cp[i]);
        ArrCpM_ = fmadd256d(Phi03v,Cp_,ArrCpM_);
        __m256d dCpdT_ = load256d01(&dCpdT[i]);
        ArrdCpMdT_ = fmadd256d(Phi03v,dCpdT_,ArrdCpMdT_);
    }
    else if(remain==3)
    {
        unsigned int i = this->nSpecies-3;

        __m256d invW_ = load256d(&this->invW[i+0]);
        invW_ = _mm256_blend_pd(invW_,zerov,0b1000);
        __m256d rhoMvj_ = _mm256_mul_pd(MWvv,invW_);
        store256d(&rhoMvj[i],rhoMvj_);
        __m256d Wv= load256d(&this->W[i+0]);
        __m256d WiByrhoM_ = _mm256_mul_pd(invrhoMv,Wv);
        store256d(&WiByrhoM[i],WiByrhoM_);

        __m256d Phi03v = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1000);

        __m256d Cv = _mm256_mul_pd(_mm256_mul_pd(rhoMvv,Phi03v),invW_);
        Cv = _mm256_max_pd(Cv,(zerov));
        store256d(&concentration[i],Cv);

        __m256d Cp_ = _mm256_blend_pd(load256d(&Cp[i]),zerov,0b1000);
        ArrCpM_ = _mm256_fmadd_pd(Phi03v,Cp_,ArrCpM_);
        __m256d dCpdT_ = _mm256_blend_pd(load256d(&dCpdT[i]),zerov,0b1000);
        ArrdCpMdT_ = _mm256_fmadd_pd(Phi03v,dCpdT_,ArrdCpMdT_);
    }

    Cp[this->nSpecies] = hsum4(ArrCpM_);
    dCpdT[this->nSpecies] = hsum4(ArrdCpMdT_);  
}
void FastChemistry::idealGas::JacobianThermoNT
(
    double p,
    double T_,
    double* __restrict__ Phi,
    double* __restrict__ concentration,
    double* __restrict__ negGstdByRT,
    double* __restrict__ negGstdByRT2,
    double* __restrict__ dBdT,
    double* __restrict__ dCpdT,
    double* __restrict__ Cp,
    double* __restrict__ Ha
)const 
{
    this->V = this->Ntot*this->Ru*T_*this->invp;
    this->invV = p/(this->Ntot*this->Ru*T_);
    unsigned int remain = this->nSpecies%4;
    const int regime = (T_ <= this->TcommonMin) ? 0 : (T_ >= this->TcommonMax) ? 1 : 2;
    __m256d vT = _mm256_set1_pd(T_);
    __m256d vInvT = _mm256_set1_pd(this->invT);
    __m256d vlogT1 = _mm256_set1_pd(this->logT-1);
    __m256d Ruv = _mm256_set1_pd(this->Ru);
    for(unsigned int i = 0; i < this->nSpecies-remain;i=i+4)
    {
        __m256d A0,A1,A2,A3,A4,A5,A6;

        if (regime == 2)
        {
            const __m256d mask = _mm256_cmp_pd(vT,load256d(&this->TcommonSoA[i]),_CMP_LT_OQ);
            A0 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[0][i]), load256d(&this->LCoeffsSoA[0][i]), mask);
            A1 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[1][i]), load256d(&this->LCoeffsSoA[1][i]), mask);
            A2 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[2][i]), load256d(&this->LCoeffsSoA[2][i]), mask);
            A3 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[3][i]), load256d(&this->LCoeffsSoA[3][i]), mask);
            A4 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[4][i]), load256d(&this->LCoeffsSoA[4][i]), mask);
            A5 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[5][i]), load256d(&this->LCoeffsSoA[5][i]), mask);
            A6 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[6][i]), load256d(&this->LCoeffsSoA[6][i]), mask);
        }
        else if(regime == 1)//T>max(Tcommon)
        {
            A0 = load256d(&this->HCoeffsSoA[0][i]);
            A1 = load256d(&this->HCoeffsSoA[1][i]);
            A2 = load256d(&this->HCoeffsSoA[2][i]);
            A3 = load256d(&this->HCoeffsSoA[3][i]);
            A4 = load256d(&this->HCoeffsSoA[4][i]);
            A5 = load256d(&this->HCoeffsSoA[5][i]);
            A6 = load256d(&this->HCoeffsSoA[6][i]);
        }
        else//T<min(Tcommon)
        {
            A0 = load256d(&this->LCoeffsSoA[0][i]);
            A1 = load256d(&this->LCoeffsSoA[1][i]);
            A2 = load256d(&this->LCoeffsSoA[2][i]);
            A3 = load256d(&this->LCoeffsSoA[3][i]);
            A4 = load256d(&this->LCoeffsSoA[4][i]);
            A5 = load256d(&this->LCoeffsSoA[5][i]);
            A6 = load256d(&this->LCoeffsSoA[6][i]);
        }

        __m256d dBdTv = _mm256_fmadd_pd(A4*0.2 ,vT,A3*0.25);
        dBdTv = _mm256_fmadd_pd(dBdTv,vT,A2*(1.0/3.0));      
        dBdTv = _mm256_fmadd_pd(dBdTv,vT,A1*0.5);
        dBdTv = _mm256_fmadd_pd(_mm256_fmadd_pd(A5,vInvT,A0),vInvT,dBdTv);
        store256d(&dBdT[i],dBdTv);

        __m256d ExpNegGstdByRTv = _mm256_fmadd_pd(A4*0.05,vT,A3*(1.0/12.0));
        ExpNegGstdByRTv = _mm256_fmadd_pd(ExpNegGstdByRTv,vT,A2*(1.0/6.0));
        ExpNegGstdByRTv = _mm256_fmadd_pd(ExpNegGstdByRTv,vT,A1*0.5);  
        ExpNegGstdByRTv = _mm256_fmsub_pd(ExpNegGstdByRTv,vT,_mm256_mul_pd(A5,vInvT));    
        ExpNegGstdByRTv = _mm256_add_pd(ExpNegGstdByRTv,A6);
        ExpNegGstdByRTv = _mm256_fmadd_pd(A0,vlogT1,ExpNegGstdByRTv);
        store256d(&negGstdByRT[i],ExpNegGstdByRTv);
        store256d(&negGstdByRT2[i],ExpNegGstdByRTv);

        __m256d dCpdTv = _mm256_fmadd_pd(vT,4*A4,3*A3);
        dCpdTv = _mm256_fmadd_pd(vT,dCpdTv,2*A2);  
        dCpdTv = _mm256_fmadd_pd(vT,dCpdTv,A1);  
        dCpdTv = _mm256_mul_pd(Ruv,dCpdTv);
        store256d(&dCpdT[i],dCpdTv);   

        __m256d Cpv = _mm256_fmadd_pd(A4 ,vT,A3);
        Cpv = _mm256_fmadd_pd(Cpv,vT,A2);
        Cpv = _mm256_fmadd_pd(Cpv,vT,A1);
        Cpv = _mm256_fmadd_pd(Cpv,vT,A0);
        Cpv = _mm256_mul_pd(Ruv,Cpv);
        store256d(&Cp[i],Cpv);   

        __m256d Hav = _mm256_fmadd_pd(A4 ,vT*0.2,A3*0.25);
        Hav = _mm256_fmadd_pd(Hav,vT,A2*(1.0/3.0));
        Hav = _mm256_fmadd_pd(Hav,vT,A1*0.5);
        Hav = _mm256_fmadd_pd(Hav,vT,A0);
        Hav = _mm256_fmadd_pd(Hav,vT,A5);
        Hav = _mm256_mul_pd(Ruv,Hav);
        store256d(&Ha[i],Hav);   
    }
    if(remain)
    {
        unsigned int i = this->nSpecies-remain;
        __m256d A0,A1,A2,A3,A4,A5,A6;
        if (regime == 2)
        {
            const __m256d mask = _mm256_cmp_pd(vT,load256d(&this->TcommonSoA[i]),_CMP_LT_OQ);
            A0 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[0][i]), load256d(&this->LCoeffsSoA[0][i]), mask);
            A1 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[1][i]), load256d(&this->LCoeffsSoA[1][i]), mask);
            A2 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[2][i]), load256d(&this->LCoeffsSoA[2][i]), mask);
            A3 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[3][i]), load256d(&this->LCoeffsSoA[3][i]), mask);
            A4 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[4][i]), load256d(&this->LCoeffsSoA[4][i]), mask);
            A5 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[5][i]), load256d(&this->LCoeffsSoA[5][i]), mask);
            A6 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[6][i]), load256d(&this->LCoeffsSoA[6][i]), mask);
        }
        else if(regime == 1)
        {
            A0 = load256d(&this->HCoeffsSoA[0][i]);
            A1 = load256d(&this->HCoeffsSoA[1][i]);
            A2 = load256d(&this->HCoeffsSoA[2][i]);
            A3 = load256d(&this->HCoeffsSoA[3][i]);
            A4 = load256d(&this->HCoeffsSoA[4][i]);
            A5 = load256d(&this->HCoeffsSoA[5][i]);
            A6 = load256d(&this->HCoeffsSoA[6][i]);
        }
        else
        {
            A0 = load256d(&this->LCoeffsSoA[0][i]);
            A1 = load256d(&this->LCoeffsSoA[1][i]);
            A2 = load256d(&this->LCoeffsSoA[2][i]);
            A3 = load256d(&this->LCoeffsSoA[3][i]);
            A4 = load256d(&this->LCoeffsSoA[4][i]);
            A5 = load256d(&this->LCoeffsSoA[5][i]);
            A6 = load256d(&this->LCoeffsSoA[6][i]);
        }
        __m256d vdBdT = _mm256_fmadd_pd(A4*0.2 ,vT,A3*0.25);
        vdBdT = _mm256_fmadd_pd(vdBdT,vT,A2*(1.0/3.0));      
        vdBdT = _mm256_fmadd_pd(vdBdT,vT,A1*0.5);
        vdBdT = _mm256_fmadd_pd(_mm256_fmadd_pd(A5,vInvT,A0),vInvT,vdBdT);
        

        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,vT,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,vT,_mm256_mul_pd(A5,vInvT));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,vlogT1,vExpNegGstdByRT);

        __m256d vdCpdT = _mm256_fmadd_pd(vT,4*A4,3*A3);
        vdCpdT = _mm256_fmadd_pd(vT,vdCpdT,2*A2);  
        vdCpdT = _mm256_fmadd_pd(vT,vdCpdT,A1);  
        vdCpdT = _mm256_mul_pd(Ruv,vdCpdT);
 

        __m256d vCp = _mm256_fmadd_pd(A4 ,vT,A3);
        vCp = _mm256_fmadd_pd(vCp,vT,A2);
        vCp = _mm256_fmadd_pd(vCp,vT,A1);
        vCp = _mm256_fmadd_pd(vCp,vT,A0);
        vCp = _mm256_mul_pd(Ruv,vCp);

        __m256d vHa = _mm256_fmadd_pd(A4 ,vT*0.2,A3*0.25);
        vHa = _mm256_fmadd_pd(vHa,vT,A2*(1.0/3.0));
        vHa = _mm256_fmadd_pd(vHa,vT,A1*0.5);
        vHa = _mm256_fmadd_pd(vHa,vT,A0);
        vHa = _mm256_fmadd_pd(vHa,vT,A5);
        vHa = _mm256_mul_pd(Ruv,vHa);
  
        if(remain==1)
        {
            store256d0(&dBdT[i],vdBdT);
            store256d0(&negGstdByRT[i],vExpNegGstdByRT);
            store256d0(&negGstdByRT2[i],vExpNegGstdByRT);
            store256d0(&dCpdT[i],vdCpdT);  
            store256d0(&Cp[i],vCp);   
            store256d0(&Ha[i],vHa); 
        }
        else if(remain==2)
        {
            store256d01(&dBdT[i],vdBdT);
            store256d01(&negGstdByRT[i],vExpNegGstdByRT);
            store256d01(&negGstdByRT2[i],vExpNegGstdByRT);
            store256d01(&dCpdT[i],vdCpdT);  
            store256d01(&Cp[i],vCp);   
            store256d01(&Ha[i],vHa); 
        }
        else
        {
            store256d012(&dBdT[i],vdBdT);
            store256d012(&negGstdByRT[i],vExpNegGstdByRT);
            store256d012(&negGstdByRT2[i],vExpNegGstdByRT);
            store256d012(&dCpdT[i],vdCpdT);  
            store256d012(&Cp[i],vCp);   
            store256d012(&Ha[i],vHa); 
        }
    }

    __m256d cptotv = _mm256_setzero_pd();
    __m256d dcptotdTv = _mm256_setzero_pd() ;

    __m256d invVv = _mm256_set1_pd(this->invV);
    __m256d zerov = _mm256_setzero_pd();
    for(unsigned int i = 0; i < this->nSpecies-remain;i=i+4)
    {
        __m256d Phi03v = load256d(&Phi[i]);
        __m256d Cv = _mm256_mul_pd(invVv,Phi03v);
        store256d(&concentration[i],Cv);
        __m256d Cp_ = load256d(&Cp[i]);
        cptotv = _mm256_fmadd_pd(Phi03v,Cp_,cptotv);
        __m256d dCpdT_ = load256d(&dCpdT[i]);
        dcptotdTv = _mm256_fmadd_pd(Phi03v,dCpdT_,dcptotdTv);
    }

    if(remain==1)
    {
        unsigned int i = this->nSpecies-1;
        __m256d Phi03v = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1110);
        __m256d Cv = _mm256_mul_pd(invVv,Phi03v);
        store256d0(&concentration[i],Cv);
        __m256d Cp_ = _mm256_blend_pd(load256d(&Cp[i]),zerov,0b1110);
        cptotv = _mm256_fmadd_pd(Phi03v,Cp_,cptotv);
        __m256d dCpdT_ = _mm256_blend_pd(load256d(&dCpdT[i]),zerov,0b1110);
        dcptotdTv = _mm256_fmadd_pd(Phi03v,dCpdT_,dcptotdTv);
    }
    else if(remain==2)
    {
        unsigned int i = this->nSpecies-2;
        __m256d Phi03v = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1100);
        __m256d Cv = _mm256_mul_pd(invVv,Phi03v);
        store256d01(&concentration[i],Cv);
        __m256d Cp_ = _mm256_blend_pd(load256d(&Cp[i]),zerov,0b1100);
        cptotv = _mm256_fmadd_pd(Phi03v,Cp_,cptotv);
        __m256d dCpdT_ = _mm256_blend_pd(load256d(&dCpdT[i]),zerov,0b1100);
        dcptotdTv = _mm256_fmadd_pd(Phi03v,dCpdT_,dcptotdTv);
    }
    else if(remain==3)
    {
        unsigned int i = this->nSpecies-3;
        __m256d Phi03v = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1000);
        __m256d Cv = _mm256_mul_pd(invVv,Phi03v);
        store256d012(&concentration[i],Cv);
        __m256d Cp_ = _mm256_blend_pd(load256d(&Cp[i]),zerov,0b1000);
        cptotv = _mm256_fmadd_pd(Phi03v,Cp_,cptotv);
        __m256d dCpdT_ = _mm256_blend_pd(load256d(&dCpdT[i]),zerov,0b1000);
        dcptotdTv = _mm256_fmadd_pd(Phi03v,dCpdT_,dcptotdTv);
    }

    Cp[this->nSpecies] = hsum4(cptotv);
    dCpdT[this->nSpecies] =  hsum4(dcptotdTv);  
}

void FastChemistry::idealGas::DerivativeThermoYT
(
    double T_,
    const double p,
    double* __restrict__ Phi,
    double* __restrict__ concentration,
    double* __restrict__ Cp,
    double* __restrict__ Ha,
    double* __restrict__ negGstdByRT,
    double* __restrict__ negGstdByRT2
)const 
{
    int remain = this->nSpecies%4;
    const int regime = (T_ <= this->TcommonMin) ? 0 : (T_ >= this->TcommonMax) ? 1 : 2;
    for (int i=0; i<this->nSpecies; i++)
    {
        Phi[i] = std::max(Phi[i], 0.0);
    }
    
    this->rhoM = 0;

    double RuTByP = this->Ru*T_/p;
    __m256d RuTByPv = _mm256_set1_pd(RuTByP);
    __m256d rhoMv = _mm256_setzero_pd();

    for (int i=0; i<this->nSpecies-remain; i=i+4)
    {
        __m256d YTpv = load256d(&Phi[i]);
        __m256d invWv = load256d(&this->invW[i]);
        rhoMv = _mm256_fmadd_pd(_mm256_mul_pd(YTpv,invWv),RuTByPv,rhoMv);
    }
    if(remain==1)
    {
        int i = this->nSpecies-1;
        __m256d zerov = _mm256_setzero_pd();
        __m256d YTpv = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1110);
        __m256d invWv = _mm256_blend_pd(load256d(&this->invW[i]),zerov,0b1110);
        rhoMv = _mm256_fmadd_pd(_mm256_mul_pd(YTpv,invWv),RuTByPv,rhoMv);
    }
    else if(remain==2)
    {
        int i = this->nSpecies-2;
        __m256d zerov = _mm256_setzero_pd();
        __m256d YTpv = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1100);
        __m256d invWv = _mm256_blend_pd(load256d(&this->invW[i]),zerov,0b1100);
        rhoMv = _mm256_fmadd_pd(_mm256_mul_pd(YTpv,invWv),RuTByPv,rhoMv);
    }
    else if(remain==3)
    {
        int i = this->nSpecies-3;
        __m256d zerov = _mm256_setzero_pd();
        __m256d YTpv = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1000);
        __m256d invWv = _mm256_blend_pd(load256d(&this->invW[i]),zerov,0b1000);
        rhoMv = _mm256_fmadd_pd(_mm256_mul_pd(YTpv,invWv),RuTByPv,rhoMv);
    }

    this->rhoM = hsum4(rhoMv);
    this->vM = this->rhoM;
    this->rhoM = 1/this->rhoM;

    __m256d rhoMvv = _mm256_set1_pd(this->rhoM);

    // Compute molar concentration for each species [kmol/m^3]
    __m256d zerov = _mm256_setzero_pd();
    for (int i=0; i<this->nSpecies-remain; i=i+4)
    {
        __m256d invWv = load256d(&this->invW[i]);
        __m256d Phiv = load256d(&Phi[i]);
        __m256d cv = _mm256_mul_pd(rhoMvv,_mm256_mul_pd(Phiv,invWv));
        cv = _mm256_max_pd(cv,zerov);
        store256d(&concentration[i],cv);
    }
    if(remain==1)
    {
        int i = this->nSpecies-1;
        __m256d invWv = _mm256_blend_pd(load256d(&this->invW[i]),zerov,0b1110);
        __m256d Phiv = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1110);
        __m256d cv = _mm256_mul_pd(rhoMvv,_mm256_mul_pd(Phiv,invWv));
        cv = _mm256_max_pd(cv,zerov);
        store256d0(&concentration[i],cv);
    }
    else if(remain==2)
    {
        int i = this->nSpecies-2;
        __m256d invWv = _mm256_blend_pd(load256d(&this->invW[i]),zerov,0b1100);
        __m256d Phiv = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1100);
        __m256d cv = _mm256_mul_pd(rhoMvv,_mm256_mul_pd(Phiv,invWv));
        cv = _mm256_max_pd(cv,zerov);
        store256d01(&concentration[i],cv);
    }
    else if(remain==3)
    {
        int i = this->nSpecies-3;
        __m256d invWv = _mm256_blend_pd(load256d(&this->invW[i]),zerov,0b1000);
        __m256d Phiv = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1000);
        __m256d cv = _mm256_mul_pd(rhoMvv,_mm256_mul_pd(Phiv,invWv));
        cv = _mm256_max_pd(cv,zerov);
        store256d012(&concentration[i],cv);
    }
    
    __m256d vT = _mm256_set1_pd(T_);
    __m256d vInvT = _mm256_set1_pd(this->invT);
    __m256d Ruv = _mm256_set1_pd(this->Ru);
    __m256d logT1v = _mm256_set1_pd(this->logT-1);
    __m256d Cpmv = _mm256_setzero_pd();

    for(int i = 0; i < this->nSpecies-remain;i=i+4)
    {
        __m256d A0,A1,A2,A3,A4,A5,A6;
        if (regime == 2)
        {
            const __m256d mask = _mm256_cmp_pd(vT,load256d(&this->TcommonSoA[i]),_CMP_LT_OQ);
            A0 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[0][i]), load256d(&this->LCoeffsSoA[0][i]), mask);
            A1 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[1][i]), load256d(&this->LCoeffsSoA[1][i]), mask);
            A2 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[2][i]), load256d(&this->LCoeffsSoA[2][i]), mask);
            A3 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[3][i]), load256d(&this->LCoeffsSoA[3][i]), mask);
            A4 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[4][i]), load256d(&this->LCoeffsSoA[4][i]), mask);
            A5 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[5][i]), load256d(&this->LCoeffsSoA[5][i]), mask);
            A6 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[6][i]), load256d(&this->LCoeffsSoA[6][i]), mask);
        }
        else if(regime == 1)
        {
            A0 = load256d(&this->HCoeffsSoA[0][i]);
            A1 = load256d(&this->HCoeffsSoA[1][i]);
            A2 = load256d(&this->HCoeffsSoA[2][i]);
            A3 = load256d(&this->HCoeffsSoA[3][i]);
            A4 = load256d(&this->HCoeffsSoA[4][i]);
            A5 = load256d(&this->HCoeffsSoA[5][i]);
            A6 = load256d(&this->HCoeffsSoA[6][i]);
        }
        else
        {
            A0 = load256d(&this->LCoeffsSoA[0][i]);
            A1 = load256d(&this->LCoeffsSoA[1][i]);
            A2 = load256d(&this->LCoeffsSoA[2][i]);
            A3 = load256d(&this->LCoeffsSoA[3][i]);
            A4 = load256d(&this->LCoeffsSoA[4][i]);
            A5 = load256d(&this->LCoeffsSoA[5][i]);
            A6 = load256d(&this->LCoeffsSoA[6][i]);
        }

        __m256d RuInvW = load256d(&this->invW[i+0]);
        RuInvW = _mm256_mul_pd(RuInvW,Ruv);

        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,vT,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,vT,_mm256_mul_pd(A5,vInvT));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,logT1v,vExpNegGstdByRT);
        store256d(&negGstdByRT[i+0],vExpNegGstdByRT);
        store256d(&negGstdByRT2[i+0],vExpNegGstdByRT);

        __m256d vCp = _mm256_fmadd_pd(A4 ,vT,A3);
        vCp = _mm256_fmadd_pd(vCp,vT,A2);
        vCp = _mm256_fmadd_pd(vCp,vT,A1);
        vCp = _mm256_fmadd_pd(vCp,vT,A0);
        vCp = _mm256_mul_pd(RuInvW,vCp);
        store256d(&Cp[i+0],vCp);
        __m256d Yv = load256d(&Phi[i+0]);
        Cpmv = _mm256_fmadd_pd(Yv,vCp,Cpmv);

        __m256d vHa = _mm256_fmadd_pd(A4 ,vT*0.2,A3*0.25);
        vHa = _mm256_fmadd_pd(vHa,vT,A2*(1.0/3.0));
        vHa = _mm256_fmadd_pd(vHa,vT,A1*0.5);
        vHa = _mm256_fmadd_pd(vHa,vT,A0);
        vHa = _mm256_fmadd_pd(vHa,vT,A5);
        vHa = _mm256_mul_pd(RuInvW,vHa);
        store256d(&Ha[i+0],vHa);
    }
    if(remain)
    {
        int i = this->nSpecies-remain;
        __m256d A0,A1,A2,A3,A4,A5,A6;
        if (regime == 2)
        {
            const __m256d mask = _mm256_cmp_pd(vT,load256d(&this->TcommonSoA[i]),_CMP_LT_OQ);
            A0 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[0][i]), load256d(&this->LCoeffsSoA[0][i]), mask);
            A1 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[1][i]), load256d(&this->LCoeffsSoA[1][i]), mask);
            A2 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[2][i]), load256d(&this->LCoeffsSoA[2][i]), mask);
            A3 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[3][i]), load256d(&this->LCoeffsSoA[3][i]), mask);
            A4 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[4][i]), load256d(&this->LCoeffsSoA[4][i]), mask);
            A5 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[5][i]), load256d(&this->LCoeffsSoA[5][i]), mask);
            A6 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[6][i]), load256d(&this->LCoeffsSoA[6][i]), mask);
        }
        else if(regime == 1)
        {
            A0 = load256d(&this->HCoeffsSoA[0][i]);
            A1 = load256d(&this->HCoeffsSoA[1][i]);
            A2 = load256d(&this->HCoeffsSoA[2][i]);
            A3 = load256d(&this->HCoeffsSoA[3][i]);
            A4 = load256d(&this->HCoeffsSoA[4][i]);
            A5 = load256d(&this->HCoeffsSoA[5][i]);
            A6 = load256d(&this->HCoeffsSoA[6][i]);
        }
        else
        {
            A0 = load256d(&this->LCoeffsSoA[0][i]);
            A1 = load256d(&this->LCoeffsSoA[1][i]);
            A2 = load256d(&this->LCoeffsSoA[2][i]);
            A3 = load256d(&this->LCoeffsSoA[3][i]);
            A4 = load256d(&this->LCoeffsSoA[4][i]);
            A5 = load256d(&this->LCoeffsSoA[5][i]);
            A6 = load256d(&this->LCoeffsSoA[6][i]);
        }

        __m256d RuInvW = load256d(&this->invW[i+0]);
        RuInvW = _mm256_mul_pd(RuInvW,Ruv);
        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,vT,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,vT,_mm256_mul_pd(A5,vInvT));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,logT1v,vExpNegGstdByRT);

        __m256d vCp = _mm256_fmadd_pd(A4 ,vT,A3);
        vCp = _mm256_fmadd_pd(vCp,vT,A2);
        vCp = _mm256_fmadd_pd(vCp,vT,A1);
        vCp = _mm256_fmadd_pd(vCp,vT,A0);
        vCp = _mm256_mul_pd(RuInvW,vCp);

        __m256d vHa = _mm256_fmadd_pd(A4 ,vT*0.2,A3*0.25);
        vHa = _mm256_fmadd_pd(vHa,vT,A2*(1.0/3.0));
        vHa = _mm256_fmadd_pd(vHa,vT,A1*0.5);
        vHa = _mm256_fmadd_pd(vHa,vT,A0);
        vHa = _mm256_fmadd_pd(vHa,vT,A5);
        vHa = _mm256_mul_pd(RuInvW,vHa);

        if(remain==1)
        {
            __m256d Yv = _mm256_blend_pd(load256d(&Phi[i+0]),_mm256_setzero_pd(),0b1110);
            Cpmv = _mm256_fmadd_pd(Yv,vCp,Cpmv);
            store256d0(&negGstdByRT[i+0],vExpNegGstdByRT);
            store256d0(&negGstdByRT2[i+0],vExpNegGstdByRT);
            store256d0(&Cp[i],vCp);
            store256d0(&Ha[i],vHa);
        }
        else if(remain==2)
        {
            __m256d Yv = _mm256_blend_pd(load256d(&Phi[i+0]),_mm256_setzero_pd(),0b1100);
            Cpmv = _mm256_fmadd_pd(Yv,vCp,Cpmv);
            store256d01(&negGstdByRT[i+0],vExpNegGstdByRT);
            store256d01(&negGstdByRT2[i+0],vExpNegGstdByRT);
            store256d01(&Cp[i],vCp);
            store256d01(&Ha[i],vHa);
        }
        else if(remain==3)
        {
            __m256d Yv = _mm256_blend_pd(load256d(&Phi[i+0]),_mm256_setzero_pd(),0b1000);
            Cpmv = _mm256_fmadd_pd(Yv,vCp,Cpmv);
            store256d012(&negGstdByRT[i+0],vExpNegGstdByRT);
            store256d012(&negGstdByRT2[i+0],vExpNegGstdByRT);
            store256d012(&Cp[i],vCp);
            store256d012(&Ha[i],vHa);
        }
    }

    Cp[this->nSpecies] = hsum4(Cpmv);
}

void FastChemistry::idealGas::DerivativeThermoNT
(
    double T_,
    const double p,
    double* __restrict__ Phi,
    double* __restrict__ concentration,
    double* __restrict__ cp,
    double* __restrict__ ha,
    double* __restrict__ negGstdByRT,
    double* __restrict__ negGstdByRT2
)const 
{

    int remain = this->nSpecies%4;
    const int regime = (T_ <= this->TcommonMin) ? 0 : (T_ >= this->TcommonMax) ? 1 : 2;

    
    this->V = this->Ntot*this->Ru*T_*this->invp;
    this->invV = p/(this->Ntot*this->Ru*T_);
    // Compute molar concentration for each species [kmol/m^3]
    
    __m256d zerov = _mm256_setzero_pd();
    __m256d invVv = _mm256_set1_pd(invV);
    for (int i=0; i<this->nSpecies-remain; i=i+4)
    {
        __m256d Nv = load256d(&Phi[i]);
        __m256d cv = _mm256_mul_pd(Nv,invVv);
        store256d(&concentration[i],cv);
    }
    if(remain==1)
    {
        int i = this->nSpecies-1;
        __m256d Nv = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1110);
        __m256d cv = _mm256_mul_pd(invVv,Nv);
        store256d0(&concentration[i],cv);
    }
    else if(remain==2)
    {
        int i = this->nSpecies-2;
        __m256d Nv = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1100);
        __m256d cv = _mm256_mul_pd(invVv,Nv);
        store256d01(&concentration[i],cv);
    }
    else if(remain==3)
    {
        int i = this->nSpecies-3;
        __m256d Nv = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1000);
        __m256d cv = _mm256_mul_pd(invVv,Nv);
        store256d012(&concentration[i],cv);
    }
    

    __m256d Tv = _mm256_set1_pd(T_);
    __m256d invTv = _mm256_set1_pd(this->invT);
    __m256d Ruv = _mm256_set1_pd(this->Ru);
    __m256d logT1v = _mm256_set1_pd(this->logT-1);
    __m256d cptotv = _mm256_setzero_pd();

    for(int i = 0; i < this->nSpecies-remain;i=i+4)
    {
        __m256d A0,A1,A2,A3,A4,A5,A6;
        if (regime == 2)
        {
            const __m256d mask = _mm256_cmp_pd(Tv,load256d(&this->TcommonSoA[i]),_CMP_LT_OQ);
            A0 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[0][i]), load256d(&this->LCoeffsSoA[0][i]), mask);
            A1 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[1][i]), load256d(&this->LCoeffsSoA[1][i]), mask);
            A2 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[2][i]), load256d(&this->LCoeffsSoA[2][i]), mask);
            A3 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[3][i]), load256d(&this->LCoeffsSoA[3][i]), mask);
            A4 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[4][i]), load256d(&this->LCoeffsSoA[4][i]), mask);
            A5 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[5][i]), load256d(&this->LCoeffsSoA[5][i]), mask);
            A6 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[6][i]), load256d(&this->LCoeffsSoA[6][i]), mask);
        }
        else if(regime == 1)
        {
            A0 = load256d(&this->HCoeffsSoA[0][i]);
            A1 = load256d(&this->HCoeffsSoA[1][i]);
            A2 = load256d(&this->HCoeffsSoA[2][i]);
            A3 = load256d(&this->HCoeffsSoA[3][i]);
            A4 = load256d(&this->HCoeffsSoA[4][i]);
            A5 = load256d(&this->HCoeffsSoA[5][i]);
            A6 = load256d(&this->HCoeffsSoA[6][i]);
        }
        else
        {
            A0 = load256d(&this->LCoeffsSoA[0][i]);
            A1 = load256d(&this->LCoeffsSoA[1][i]);
            A2 = load256d(&this->LCoeffsSoA[2][i]);
            A3 = load256d(&this->LCoeffsSoA[3][i]);
            A4 = load256d(&this->LCoeffsSoA[4][i]);
            A5 = load256d(&this->LCoeffsSoA[5][i]);
            A6 = load256d(&this->LCoeffsSoA[6][i]);
        }

        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,Tv,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,Tv,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,Tv,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,Tv,_mm256_mul_pd(A5,invTv));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,logT1v,vExpNegGstdByRT);
        store256d(&negGstdByRT[i+0],vExpNegGstdByRT);
        store256d(&negGstdByRT2[i+0],vExpNegGstdByRT);

        __m256d cpv = _mm256_fmadd_pd(A4 ,Tv,A3);
        cpv = _mm256_fmadd_pd(cpv,Tv,A2);
        cpv = _mm256_fmadd_pd(cpv,Tv,A1);
        cpv = _mm256_fmadd_pd(cpv,Tv,A0);
        cpv = _mm256_mul_pd(Ruv,cpv);
        store256d(&cp[i+0],cpv);
        __m256d Nv = load256d(&Phi[i+0]);
        cptotv = _mm256_fmadd_pd(Nv,cpv,cptotv);

        __m256d hav = _mm256_fmadd_pd(A4 ,Tv*0.2,A3*0.25);
        hav = _mm256_fmadd_pd(hav,Tv,A2*(1.0/3.0));
        hav = _mm256_fmadd_pd(hav,Tv,A1*0.5);
        hav = _mm256_fmadd_pd(hav,Tv,A0);
        hav = _mm256_fmadd_pd(hav,Tv,A5);
        hav = _mm256_mul_pd(Ruv,hav);
        store256d(&ha[i+0],hav);
    }
    if(remain)
    {
        int i = this->nSpecies-remain;
        __m256d A0,A1,A2,A3,A4,A5,A6;
        if (regime == 2)
        {
            const __m256d mask = _mm256_cmp_pd(Tv,load256d(&this->TcommonSoA[i]),_CMP_LT_OQ);
            A0 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[0][i]), load256d(&this->LCoeffsSoA[0][i]), mask);
            A1 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[1][i]), load256d(&this->LCoeffsSoA[1][i]), mask);
            A2 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[2][i]), load256d(&this->LCoeffsSoA[2][i]), mask);
            A3 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[3][i]), load256d(&this->LCoeffsSoA[3][i]), mask);
            A4 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[4][i]), load256d(&this->LCoeffsSoA[4][i]), mask);
            A5 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[5][i]), load256d(&this->LCoeffsSoA[5][i]), mask);
            A6 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[6][i]), load256d(&this->LCoeffsSoA[6][i]), mask);
        }
        else if(regime == 1)
        {
            A0 = load256d(&this->HCoeffsSoA[0][i]);
            A1 = load256d(&this->HCoeffsSoA[1][i]);
            A2 = load256d(&this->HCoeffsSoA[2][i]);
            A3 = load256d(&this->HCoeffsSoA[3][i]);
            A4 = load256d(&this->HCoeffsSoA[4][i]);
            A5 = load256d(&this->HCoeffsSoA[5][i]);
            A6 = load256d(&this->HCoeffsSoA[6][i]);
        }
        else
        {
            A0 = load256d(&this->LCoeffsSoA[0][i]);
            A1 = load256d(&this->LCoeffsSoA[1][i]);
            A2 = load256d(&this->LCoeffsSoA[2][i]);
            A3 = load256d(&this->LCoeffsSoA[3][i]);
            A4 = load256d(&this->LCoeffsSoA[4][i]);
            A5 = load256d(&this->LCoeffsSoA[5][i]);
            A6 = load256d(&this->LCoeffsSoA[6][i]);
        }

        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,Tv,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,Tv,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,Tv,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,Tv,_mm256_mul_pd(A5,invTv));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,logT1v,vExpNegGstdByRT);

        __m256d cpv = _mm256_fmadd_pd(A4 ,Tv,A3);
        cpv = _mm256_fmadd_pd(cpv,Tv,A2);
        cpv = _mm256_fmadd_pd(cpv,Tv,A1);
        cpv = _mm256_fmadd_pd(cpv,Tv,A0);
        cpv = _mm256_mul_pd(Ruv,cpv);

        __m256d hav = _mm256_fmadd_pd(A4 ,Tv*0.2,A3*0.25);
        hav = _mm256_fmadd_pd(hav,Tv,A2*(1.0/3.0));
        hav = _mm256_fmadd_pd(hav,Tv,A1*0.5);
        hav = _mm256_fmadd_pd(hav,Tv,A0);
        hav = _mm256_fmadd_pd(hav,Tv,A5);
        hav = _mm256_mul_pd(Ruv,hav);
        if(remain==1)
        {
            __m256d Nv = load256d0(&Phi[i+0]);
            cptotv = _mm256_fmadd_pd(Nv,cpv,cptotv);
            store256d0(&negGstdByRT[i+0],vExpNegGstdByRT);
            store256d0(&negGstdByRT2[i+0],vExpNegGstdByRT);
            store256d0(&cp[i+0],cpv);
            store256d0(&ha[i+0],hav);
        }
        else if(remain==2)
        {
            __m256d Nv = load256d01(&Phi[i+0]);
            cptotv = _mm256_fmadd_pd(Nv,cpv,cptotv);
            store256d01(&negGstdByRT[i+0],vExpNegGstdByRT);
            store256d01(&negGstdByRT2[i+0],vExpNegGstdByRT);
            store256d01(&cp[i+0],cpv);
            store256d01(&ha[i+0],hav);
        }
        else if(remain==3)
        {
            __m256d Nv = load256d012(&Phi[i+0]);
            cptotv = _mm256_fmadd_pd(Nv,cpv,cptotv);
            store256d012(&negGstdByRT[i+0],vExpNegGstdByRT);
            store256d012(&negGstdByRT2[i+0],vExpNegGstdByRT);
            store256d012(&cp[i+0],cpv);
            store256d012(&ha[i+0],hav);
        }
    }
    cp[this->nSpecies] = hsum4(cptotv);
}

void FastChemistry::idealGas::negGstdByRT
(
    double T_,
    double* __restrict__ negGstdByRT,
    double* __restrict__ negGstdByRT2
)const 
{
    int remain = this->nSpecies%4;
    const int regime = (T_ <= this->TcommonMin) ? 0 : (T_ >= this->TcommonMax) ? 1 : 2;
    __m256d Tv = _mm256_set1_pd(T_);
    __m256d vInvT = _mm256_set1_pd(invT);
    __m256d Ruv = _mm256_set1_pd(this->Ru);
    __m256d logT1v = _mm256_set1_pd(this->logT-1);
    for(int i = 0; i < this->nSpecies-remain;i=i+4)
    {
        __m256d A0,A1,A2,A3,A4,A5,A6;
        if (regime == 2)
        {
            const __m256d mask = _mm256_cmp_pd(Tv,load256d(&this->TcommonSoA[i]),_CMP_LT_OQ);
            A0 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[0][i]), load256d(&this->LCoeffsSoA[0][i]), mask);
            A1 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[1][i]), load256d(&this->LCoeffsSoA[1][i]), mask);
            A2 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[2][i]), load256d(&this->LCoeffsSoA[2][i]), mask);
            A3 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[3][i]), load256d(&this->LCoeffsSoA[3][i]), mask);
            A4 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[4][i]), load256d(&this->LCoeffsSoA[4][i]), mask);
            A5 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[5][i]), load256d(&this->LCoeffsSoA[5][i]), mask);
            A6 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[6][i]), load256d(&this->LCoeffsSoA[6][i]), mask);
        }
        else if(regime == 1)
        {
            A0 = load256d(&this->HCoeffsSoA[0][i]);
            A1 = load256d(&this->HCoeffsSoA[1][i]);
            A2 = load256d(&this->HCoeffsSoA[2][i]);
            A3 = load256d(&this->HCoeffsSoA[3][i]);
            A4 = load256d(&this->HCoeffsSoA[4][i]);
            A5 = load256d(&this->HCoeffsSoA[5][i]);
            A6 = load256d(&this->HCoeffsSoA[6][i]);
        }
        else
        {
            A0 = load256d(&this->LCoeffsSoA[0][i]);
            A1 = load256d(&this->LCoeffsSoA[1][i]);
            A2 = load256d(&this->LCoeffsSoA[2][i]);
            A3 = load256d(&this->LCoeffsSoA[3][i]);
            A4 = load256d(&this->LCoeffsSoA[4][i]);
            A5 = load256d(&this->LCoeffsSoA[5][i]);
            A6 = load256d(&this->LCoeffsSoA[6][i]);
        }

        __m256d RuInvW = load256d(&this->invW[i+0]);
        RuInvW = _mm256_mul_pd(RuInvW,Ruv);

        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,Tv,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,Tv,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,Tv,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,Tv,_mm256_mul_pd(A5,vInvT));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,logT1v,vExpNegGstdByRT);
        store256d(&negGstdByRT[i+0],vExpNegGstdByRT);
        store256d(&negGstdByRT2[i+0],vExpNegGstdByRT);
    }
    if(remain)
    {
        int i = this->nSpecies-remain;
        __m256d A0,A1,A2,A3,A4,A5,A6;
        if (regime == 2)
        {
            const __m256d mask = _mm256_cmp_pd(Tv,load256d(&this->TcommonSoA[i]),_CMP_LT_OQ);
            A0 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[0][i]), load256d(&this->LCoeffsSoA[0][i]), mask);
            A1 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[1][i]), load256d(&this->LCoeffsSoA[1][i]), mask);
            A2 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[2][i]), load256d(&this->LCoeffsSoA[2][i]), mask);
            A3 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[3][i]), load256d(&this->LCoeffsSoA[3][i]), mask);
            A4 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[4][i]), load256d(&this->LCoeffsSoA[4][i]), mask);
            A5 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[5][i]), load256d(&this->LCoeffsSoA[5][i]), mask);
            A6 = _mm256_blendv_pd(load256d(&this->HCoeffsSoA[6][i]), load256d(&this->LCoeffsSoA[6][i]), mask);
        }
        else if(regime == 1)
        {
            A0 = load256d(&this->HCoeffsSoA[0][i]);
            A1 = load256d(&this->HCoeffsSoA[1][i]);
            A2 = load256d(&this->HCoeffsSoA[2][i]);
            A3 = load256d(&this->HCoeffsSoA[3][i]);
            A4 = load256d(&this->HCoeffsSoA[4][i]);
            A5 = load256d(&this->HCoeffsSoA[5][i]);
            A6 = load256d(&this->HCoeffsSoA[6][i]);
        }
        else
        {
            A0 = load256d(&this->LCoeffsSoA[0][i]);
            A1 = load256d(&this->LCoeffsSoA[1][i]);
            A2 = load256d(&this->LCoeffsSoA[2][i]);
            A3 = load256d(&this->LCoeffsSoA[3][i]);
            A4 = load256d(&this->LCoeffsSoA[4][i]);
            A5 = load256d(&this->LCoeffsSoA[5][i]);
            A6 = load256d(&this->LCoeffsSoA[6][i]);
        }

        __m256d RuInvW = load256d(&this->invW[i+0]);
        RuInvW = _mm256_mul_pd(RuInvW,Ruv);

        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,Tv,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,Tv,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,Tv,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,Tv,_mm256_mul_pd(A5,vInvT));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,logT1v,vExpNegGstdByRT);
        if(remain==1)
        {
            store256d0(&negGstdByRT[i+0],vExpNegGstdByRT);
            store256d0(&negGstdByRT2[i+0],vExpNegGstdByRT);
        }
        else if(remain==2)
        {
            store256d01(&negGstdByRT[i+0],vExpNegGstdByRT);
            store256d01(&negGstdByRT2[i+0],vExpNegGstdByRT);
        }
        else if(remain==3)
        {
            store256d012(&negGstdByRT[i+0],vExpNegGstdByRT);
            store256d012(&negGstdByRT2[i+0],vExpNegGstdByRT);
        }
    }
}

// * * * * * * * * * * * * * * * * Member Functions (SoA)  * * * * * * * * * //

