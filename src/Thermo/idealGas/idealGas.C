#include "idealGas.H"
#include <cstring>


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
    ThighMax(1),
    PtrCoeffs(n),
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
    if (posix_memalign(reinterpret_cast<void**>(&this->buffer), 32, alignSpecies*3*sizeof(double)))
    {
        throw std::bad_alloc();
    }
    memset(this->buffer, 0, alignSpecies*3*sizeof(double));
    this->W = &this->buffer[0];
    this->invW = &this->buffer[alignSpecies*1];
    this->Hf = &this->buffer[alignSpecies*2];
}

FastChemistry::idealGas::idealGas(const idealGas&& gas)
:
    nSpecies(gas.nSpecies),
    HCoeffs(std::move(gas.HCoeffs)),
    LCoeffs(std::move(gas.LCoeffs)),
    Tlow(std::move(gas.Tlow)),
    Thigh(std::move(gas.Thigh)),
    Tcommon(std::move(gas.Tcommon)),
    TcommonMin(gas.TcommonMin),
    TcommonMax(gas.TcommonMax),
    TlowMin(gas.TlowMin),
    ThighMax(gas.ThighMax),
    PtrCoeffs(std::move(gas.PtrCoeffs)),
    buffer(gas.buffer),
    W(gas.W),
    invW(gas.invW),
    Hf(gas.Hf),
    logT(gas.logT),
    invT(gas.invT),
    sqrT(gas.sqrT),
    logP(gas.logP),
    rhoM(gas.rhoM),
    vM(gas.vM),
    Mw(gas.Mw)
{
    gas.buffer = nullptr;
    gas.W = nullptr;
    gas.invW = nullptr;
    gas.invW = nullptr;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FastChemistry::idealGas::~idealGas()
{
    this->W=nullptr;
    this->invW=nullptr;
    this->Hf=nullptr;

    if(this->buffer!=nullptr)
    {
        free(this->buffer);
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

    this->setPtrCoeffs(T_);
    this->rhoM = 0;
    double MW = 0;    

    const double RuTByP = (this->Ru*T_)/p;
    unsigned int remain = this->nSpecies%4;
    __m256d MWv = _mm256_setzero_pd();
    __m256d rhoMv = _mm256_setzero_pd();
    __m256d vT = _mm256_set1_pd(T_);
    __m256d vInvT = _mm256_set1_pd(this->invT);
    __m256d vlogT1 = _mm256_set1_pd(this->logT-1);
    __m256d Ruv = _mm256_set1_pd(this->Ru);
    for(unsigned int i = 0; i < this->nSpecies-remain;i=i+4)
    {
        unsigned int i0 = i+0;
        unsigned int i1 = i+1;
        unsigned int i2 = i+2;
        unsigned int i3 = i+3;
        std::array<double,7>& a = *this->PtrCoeffs[i0];
        std::array<double,7>& b = *this->PtrCoeffs[i1];
        std::array<double,7>& c = *this->PtrCoeffs[i2];
        std::array<double,7>& d = *this->PtrCoeffs[i3];

        __m256d invWv = load256d(&invW[i0]);
        __m256d YTpv = load256d(&Phi[i0]);
        __m256d RuTByPv = _mm256_set1_pd(RuTByP);

        rhoMv = _mm256_fmadd_pd(RuTByPv,_mm256_mul_pd(YTpv,invWv),rhoMv);
        MWv = _mm256_fmadd_pd(YTpv,invWv,MWv);

        __m256d A0 = _mm256_setr_pd(a[0],b[0],c[0],d[0]);
        __m256d A1 = _mm256_setr_pd(a[1],b[1],c[1],d[1]);
        __m256d A2 = _mm256_setr_pd(a[2],b[2],c[2],d[2]);
        __m256d A3 = _mm256_setr_pd(a[3],b[3],c[3],d[3]);
        __m256d A4 = _mm256_setr_pd(a[4],b[4],c[4],d[4]);
        __m256d A5 = _mm256_setr_pd(a[5],b[5],c[5],d[5]);
        __m256d A6 = _mm256_setr_pd(a[6],b[6],c[6],d[6]);


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
        //store256d(&ExpNegGstdByRT[i0],vExpNegGstdByRT);
        store256d(&negGstdByRT[i0],vExpNegGstdByRT);
        store256d(&negGstdByRT2[i0],vExpNegGstdByRT);

        __m256d RuInvW = _mm256_mul_pd(invWv,Ruv);
        __m256d vdCpdT = _mm256_fmadd_pd(vT,4*A4,3*A3);
        vdCpdT = _mm256_fmadd_pd(vT,vdCpdT,2*A2);  
        vdCpdT = _mm256_fmadd_pd(vT,vdCpdT,A1);  
        vdCpdT = _mm256_mul_pd(RuInvW,vdCpdT);
        store256d(&dCpdT[i0],vdCpdT);   

        __m256d vCp = _mm256_fmadd_pd(A4 ,vT,A3);
        vCp = _mm256_fmadd_pd(vCp,vT,A2);
        vCp = _mm256_fmadd_pd(vCp,vT,A1);
        vCp = _mm256_fmadd_pd(vCp,vT,A0);
        vCp = _mm256_mul_pd(RuInvW,vCp);
        store256d(&Cp[i0],vCp);   

        __m256d vHa = _mm256_fmadd_pd(A4 ,vT*0.2,A3*0.25);
        vHa = _mm256_fmadd_pd(vHa,vT,A2*(1.0/3.0));
        vHa = _mm256_fmadd_pd(vHa,vT,A1*0.5);
        vHa = _mm256_fmadd_pd(vHa,vT,A0);
        vHa = _mm256_fmadd_pd(vHa,vT,A5);
        vHa = _mm256_mul_pd(RuInvW,vHa);
        store256d(&Ha[i0],vHa);   

    }
    if(remain==1)
    {
        unsigned int i = this->nSpecies -1;

        const std::array<double,7>& a = *this->PtrCoeffs[i];
        this->rhoM += Phi[i]*RuTByP*invW[i];
        MW += Phi[i]*this->invW[i];
        dBdT[i] = (a[0]+a[5]*this->invT)*this->invT + a[1]*0.5 + T_*(a[2]*(1.0/3.0)+T_*(a[3]*0.25+a[4]*0.2*T_));
        double negGstdByRT_ = ((((a[4]*0.05*T_+a[3]*(1.0/12.0))*T_+a[2]*(1.0/6.0))*T_+a[1]*0.5)*T_-a[5]*this->invT+a[6]+(this->logT-1)*a[0]);   
        //ExpNegGstdByRT[i] = negGstdByRT_;
        negGstdByRT[i] = negGstdByRT_;
        negGstdByRT2[i] = negGstdByRT_;
        dCpdT[i] = ( a[1] + T_*(2*a[2]+T_*(3*a[3]+T_*4*a[4])) )*this->Ru*this->invW[i]; 
        Cp[i] = ((((a[4]*T_+a[3])*T_+a[2])*T_+a[1])*T_+a[0])*this->Ru*this->invW[i];
        Ha[i] = (((((a[4]*T_*0.2+a[3]*0.25)*T_+a[2]*(1.0/3.0))*T_+a[1]*0.5)*T_+a[0])*T_ +a[5])*this->Ru*this->invW[i]; 
    }
    else if(remain==2)
    {

        unsigned int i = this->nSpecies - 2;

        std::array<double,7>& a = *this->PtrCoeffs[i+0];
        std::array<double,7>& b = *this->PtrCoeffs[i+1];
        this->rhoM += Phi[i+0]*RuTByP*invW[i+0];
        this->rhoM += Phi[i+1]*RuTByP*invW[i+1];
        MW += Phi[i+0]*this->invW[i+0];
        MW += Phi[i+1]*this->invW[i+1];
        
        __m128d invWv = _mm_load_pd(&invW[i]);
        __m128d A0 = _mm_setr_pd(a[0],b[0]);
        __m128d A1 = _mm_setr_pd(a[1],b[1]);
        __m128d A2 = _mm_setr_pd(a[2],b[2]);
        __m128d A3 = _mm_setr_pd(a[3],b[3]);
        __m128d A4 = _mm_setr_pd(a[4],b[4]);
        __m128d A5 = _mm_setr_pd(a[5],b[5]);
        __m128d A6 = _mm_setr_pd(a[6],b[6]);


        __m128d vdBdT = _mm_fmadd_pd(A4*0.2 ,_mm256_castpd256_pd128(vT),A3*0.25);
        vdBdT = _mm_fmadd_pd(vdBdT,_mm256_castpd256_pd128(vT),A2*(1.0/3.0));      
        vdBdT = _mm_fmadd_pd(vdBdT,_mm256_castpd256_pd128(vT),A1*0.5);
        vdBdT = _mm_fmadd_pd(_mm_fmadd_pd(A5,_mm256_castpd256_pd128(vInvT),A0),_mm256_castpd256_pd128(vInvT),vdBdT);

        store128d(&dBdT[i],vdBdT);

        __m128d vExpNegGstdByRT = _mm_fmadd_pd(A4*0.05,_mm256_castpd256_pd128(vT),A3*(1.0/12.0));
        vExpNegGstdByRT = _mm_fmadd_pd(vExpNegGstdByRT,_mm256_castpd256_pd128(vT),A2*(1.0/6.0));
        vExpNegGstdByRT = _mm_fmadd_pd(vExpNegGstdByRT,_mm256_castpd256_pd128(vT),A1*0.5);  
        vExpNegGstdByRT = _mm_fmsub_pd(vExpNegGstdByRT,_mm256_castpd256_pd128(vT),_mm_mul_pd(A5,_mm256_castpd256_pd128(vInvT)));    
        vExpNegGstdByRT = _mm_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm_fmadd_pd(A0,_mm256_castpd256_pd128(vlogT1),vExpNegGstdByRT);
        //store128d(&ExpNegGstdByRT[i],vExpNegGstdByRT);
        store128d(&negGstdByRT[i],vExpNegGstdByRT);
        store128d(&negGstdByRT2[i],vExpNegGstdByRT);

        __m128d RuInvW = _mm_mul_pd(invWv,_mm256_castpd256_pd128(Ruv));
        __m128d vdCpdT = _mm_fmadd_pd(_mm256_castpd256_pd128(vT),4*A4,3*A3);
        vdCpdT = _mm_fmadd_pd(_mm256_castpd256_pd128(vT),vdCpdT,2*A2);  
        vdCpdT = _mm_fmadd_pd(_mm256_castpd256_pd128(vT),vdCpdT,A1);  
        vdCpdT = _mm_mul_pd(RuInvW,vdCpdT);
        store128d(&dCpdT[i],vdCpdT);

        __m128d vCp = _mm_fmadd_pd(A4 ,_mm256_castpd256_pd128(vT),A3);
        vCp = _mm_fmadd_pd(vCp,_mm256_castpd256_pd128(vT),A2);
        vCp = _mm_fmadd_pd(vCp,_mm256_castpd256_pd128(vT),A1);
        vCp = _mm_fmadd_pd(vCp,_mm256_castpd256_pd128(vT),A0);
        vCp = _mm_mul_pd(RuInvW,vCp);
        store128d(&Cp[i],vCp);

        __m128d vHa = _mm_fmadd_pd(A4 ,_mm256_castpd256_pd128(vT)*0.2,A3*0.25);
        vHa = _mm_fmadd_pd(vHa,_mm256_castpd256_pd128(vT),A2*(1.0/3.0));
        vHa = _mm_fmadd_pd(vHa,_mm256_castpd256_pd128(vT),A1*0.5);
        vHa = _mm_fmadd_pd(vHa,_mm256_castpd256_pd128(vT),A0);
        vHa = _mm_fmadd_pd(vHa,_mm256_castpd256_pd128(vT),A5);
        vHa = _mm_mul_pd(RuInvW,vHa);
        store128d(&Ha[i],vHa);
    }
    else if(remain==3)
    {
        unsigned int i = this->nSpecies - 3;

        std::array<double,7>& a = *this->PtrCoeffs[i+0];
        std::array<double,7>& b = *this->PtrCoeffs[i+1];
        std::array<double,7>& c = *this->PtrCoeffs[i+2];
        this->rhoM += Phi[i+0]*RuTByP*invW[i+0];
        this->rhoM += Phi[i+1]*RuTByP*invW[i+1];
        this->rhoM += Phi[i+2]*RuTByP*invW[i+2];    

        MW += Phi[i+0]*this->invW[i+0];
        MW += Phi[i+1]*this->invW[i+1];
        MW += Phi[i+2]*this->invW[i+2];

        __m256d invWv = _mm256_setr_pd(invW[i+0],invW[i+1],invW[i+2],1);
        __m256d A0 = _mm256_setr_pd(a[0],b[0],c[0],0);
        __m256d A1 = _mm256_setr_pd(a[1],b[1],c[1],0);
        __m256d A2 = _mm256_setr_pd(a[2],b[2],c[2],0);
        __m256d A3 = _mm256_setr_pd(a[3],b[3],c[3],0);
        __m256d A4 = _mm256_setr_pd(a[4],b[4],c[4],0);
        __m256d A5 = _mm256_setr_pd(a[5],b[5],c[5],0);
        __m256d A6 = _mm256_setr_pd(a[6],b[6],c[6],0);


        __m256d vdBdT = _mm256_fmadd_pd(A4*0.2 ,vT,A3*0.25);
        vdBdT = _mm256_fmadd_pd(vdBdT,vT,A2*(1.0/3.0));      
        vdBdT = _mm256_fmadd_pd(vdBdT,vT,A1*0.5);
        vdBdT = _mm256_fmadd_pd(_mm256_fmadd_pd(A5,vInvT,A0),vInvT,vdBdT);
        store256d(&dBdT[i],vdBdT);

        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,vT,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,vT,_mm256_mul_pd(A5,vInvT));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,vlogT1,vExpNegGstdByRT);
        //store256d(&ExpNegGstdByRT[i],vExpNegGstdByRT);
        store256d(&negGstdByRT[i],vExpNegGstdByRT);
        store256d(&negGstdByRT2[i],vExpNegGstdByRT);

        __m256d RuInvW = _mm256_mul_pd(invWv,Ruv);
        __m256d vdCpdT = _mm256_fmadd_pd(vT,4*A4,3*A3);
        vdCpdT = _mm256_fmadd_pd(vT,vdCpdT,2*A2);  
        vdCpdT = _mm256_fmadd_pd(vT,vdCpdT,A1);  
        vdCpdT = _mm256_mul_pd(RuInvW,vdCpdT);
        store256d(&dCpdT[i],vdCpdT);


        __m256d vCp = _mm256_fmadd_pd(A4 ,vT,A3);
        vCp = _mm256_fmadd_pd(vCp,vT,A2);
        vCp = _mm256_fmadd_pd(vCp,vT,A1);
        vCp = _mm256_fmadd_pd(vCp,vT,A0);
        vCp = _mm256_mul_pd(RuInvW,vCp);
        store256d(&Cp[i],vCp);
         

        __m256d vHa = _mm256_fmadd_pd(A4 ,vT*0.2,A3*0.25);
        vHa = _mm256_fmadd_pd(vHa,vT,A2*(1.0/3.0));
        vHa = _mm256_fmadd_pd(vHa,vT,A1*0.5);
        vHa = _mm256_fmadd_pd(vHa,vT,A0);
        vHa = _mm256_fmadd_pd(vHa,vT,A5);
        vHa = _mm256_mul_pd(RuInvW,vHa);
        store256d(&Ha[i],vHa);
    }
    MW = MW + hsum4(MWv);
    MW = 1/MW;
    dBdT[this->nSpecies] = MW;
    this->rhoM = this->rhoM + hsum4(rhoMv);
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


    double CpM = 0;
    double dCpMdT = 0;
    if(remain==1)
    {
        unsigned int i = this->nSpecies-1;
        rhoMvj[i] = MW*this->invW[i];
        WiByrhoM[i] = this->W[i]*invRhoM;
        concentration[i] = std::max(this->rhoM*this->invW[i]*Phi[i],0.0);
        CpM += Phi[i]*Cp[i];
        dCpMdT += Phi[i]*dCpdT[i];
    }
    else if(remain==2)
    {
        unsigned int i = this->nSpecies-2;

        __m128d invW_ = load128d(&this->invW[i]);
        __m128d rhoMvj_ = _mm_mul_pd(_mm256_castpd256_pd128(MWvv),invW_);
        store128d(&rhoMvj[i],rhoMvj_);

        __m128d Wv= load128d(&this->W[i]);
        __m128d WiByrhoM_ = _mm_mul_pd(_mm256_castpd256_pd128(invrhoMv),Wv);
        store128d(&WiByrhoM[i],WiByrhoM_);

        __m128d Phi01v = load128d(&Phi[i]);
        __m128d Cv = _mm_mul_pd(_mm_mul_pd(_mm256_castpd256_pd128(rhoMvv),Phi01v),invW_);
        Cv = _mm_max_pd(Cv,_mm256_castpd256_pd128(zerov));
        store128d(&concentration[i],Cv);

        __m128d Cp_ = load128d(&Cp[i]);
        ArrCpM_ = _mm256_fmadd_pd(_mm256_zextpd128_pd256(Phi01v),_mm256_zextpd128_pd256(Cp_),ArrCpM_);

        __m128d dCpdT_ = load128d(&dCpdT[i]);
        ArrdCpMdT_ = _mm256_fmadd_pd(_mm256_zextpd128_pd256(Phi01v),_mm256_zextpd128_pd256(dCpdT_),ArrdCpMdT_);
    }
    else if(remain==3)
    {
        unsigned int i = this->nSpecies-3;
        __m256d zerov = _mm256_setzero_pd();

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

    Cp[this->nSpecies] = CpM + hsum4(ArrCpM_);
    dCpdT[this->nSpecies] = dCpMdT + hsum4(ArrdCpMdT_);  
}

/*void FastChemistry::idealGas::CpHaNegGstdByRT
(
    double T,
    double* __restrict__ Cp,
    double* __restrict__ Ha,
    double* __restrict__ negGstdByRT
)const 
{
    T = T<TlowMin?TlowMin:T;
    T = T>ThighMax?ThighMax:T;

    this->setPtrCoeffs(T);
    int remain = this->nSpecies%4;
    __m256d vT = _mm256_set1_pd(T);
    __m256d vInvT = _mm256_set1_pd(invT);
    __m256d Ruv = _mm256_set1_pd(this->Ru);
    __m256d logT1v = _mm256_set1_pd(this->logT-1);
    for(int i = 0; i < this->nSpecies-remain;i=i+4)
    {
        const std::array<double,7>& a = *this->PtrCoeffs[i+0];
        const std::array<double,7>& b = *this->PtrCoeffs[i+1];
        const std::array<double,7>& c = *this->PtrCoeffs[i+2];
        const std::array<double,7>& d = *this->PtrCoeffs[i+3]; 
        __m256d A0 = _mm256_setr_pd(a[0],b[0],c[0],d[0]);
        __m256d A1 = _mm256_setr_pd(a[1],b[1],c[1],d[1]);
        __m256d A2 = _mm256_setr_pd(a[2],b[2],c[2],d[2]);
        __m256d A3 = _mm256_setr_pd(a[3],b[3],c[3],d[3]);
        __m256d A4 = _mm256_setr_pd(a[4],b[4],c[4],d[4]);
        __m256d A5 = _mm256_setr_pd(a[5],b[5],c[5],d[5]);
        __m256d A6 = _mm256_setr_pd(a[6],b[6],c[6],d[6]);

        __m256d RuInvW = load256d(&this->invW[i+0]);
        RuInvW = _mm256_mul_pd(RuInvW,Ruv);

        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,vT,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,vT,_mm256_mul_pd(A5,vInvT));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,logT1v,vExpNegGstdByRT);
        store256d(&negGstdByRT[i+0],vExpNegGstdByRT);

        __m256d vCp = _mm256_fmadd_pd(A4 ,vT,A3);
        vCp = _mm256_fmadd_pd(vCp,vT,A2);
        vCp = _mm256_fmadd_pd(vCp,vT,A1);
        vCp = _mm256_fmadd_pd(vCp,vT,A0);
        vCp = _mm256_mul_pd(RuInvW,vCp);
        store256d(&Cp[i+0],vCp);

        __m256d vHa = _mm256_fmadd_pd(A4 ,vT*0.2,A3*0.25);
        vHa = _mm256_fmadd_pd(vHa,vT,A2*(1.0/3.0));
        vHa = _mm256_fmadd_pd(vHa,vT,A1*0.5);
        vHa = _mm256_fmadd_pd(vHa,vT,A0);
        vHa = _mm256_fmadd_pd(vHa,vT,A5);
        vHa = _mm256_mul_pd(RuInvW,vHa);
        store256d(&Ha[i+0],vHa);
    }
    if(remain ==1)
    {
        int i = this->nSpecies-1;
        const std::array<double,7>& a = *this->PtrCoeffs[i];        
        Cp[i] = ((((a[4]*T+a[3])*T+a[2])*T+a[1])*T+a[0])*this->Ru*this->invW[i];
        Ha[i] = (((((a[4]*T*0.2+a[3]*0.25)*T+a[2]*(1.0/3.0))*T+a[1]*0.5)*T+a[0])*T +a[5])*this->Ru*this->invW[i]; 
        double negGstdByRT_ = (((a[4]*0.05*T+a[3]*(1.0/12.0))*T+a[2]*(1.0/6.0))*T+a[1]*0.5)*T-a[5]*invT+a[6]+(this->logT-1)*a[0];
        negGstdByRT[i] = negGstdByRT_;
    }
    else if(remain==2)
    {
        int i = this->nSpecies-2;

        const std::array<double,7>& a = *this->PtrCoeffs[i+0];
        const std::array<double,7>& b = *this->PtrCoeffs[i+1];

        __m128d A0 = _mm_setr_pd(a[0],b[0]);
        __m128d A1 = _mm_setr_pd(a[1],b[1]);
        __m128d A2 = _mm_setr_pd(a[2],b[2]);
        __m128d A3 = _mm_setr_pd(a[3],b[3]);
        __m128d A4 = _mm_setr_pd(a[4],b[4]);
        __m128d A5 = _mm_setr_pd(a[5],b[5]);
        __m128d A6 = _mm_setr_pd(a[6],b[6]);

        __m128d RuInvW = load128d(&this->invW[i+0]);
        RuInvW = _mm_mul_pd(RuInvW,_mm256_castpd256_pd128(Ruv));

        __m128d vExpNegGstdByRT = _mm_fmadd_pd(A4*0.05,_mm256_castpd256_pd128(vT),A3*(1.0/12.0));
        vExpNegGstdByRT = _mm_fmadd_pd(vExpNegGstdByRT,_mm256_castpd256_pd128(vT),A2*(1.0/6.0));
        vExpNegGstdByRT = _mm_fmadd_pd(vExpNegGstdByRT,_mm256_castpd256_pd128(vT),A1*0.5);  
        vExpNegGstdByRT = _mm_fmsub_pd(vExpNegGstdByRT,_mm256_castpd256_pd128(vT),_mm_mul_pd(A5,_mm256_castpd256_pd128(vInvT)));    
        vExpNegGstdByRT = _mm_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm_fmadd_pd(A0,_mm256_castpd256_pd128(logT1v),vExpNegGstdByRT);
        store128d(&negGstdByRT[i+0],vExpNegGstdByRT);

        __m128d vCp = _mm_fmadd_pd(A4,_mm256_castpd256_pd128(vT),A3);
        vCp = _mm_fmadd_pd(vCp,_mm256_castpd256_pd128(vT),A2);
        vCp = _mm_fmadd_pd(vCp,_mm256_castpd256_pd128(vT),A1);
        vCp = _mm_fmadd_pd(vCp,_mm256_castpd256_pd128(vT),A0);
        vCp = _mm_mul_pd(RuInvW,vCp);
        store128d(&Cp[i+0],vCp);

        __m128d vHa = _mm_fmadd_pd(A4 ,_mm256_castpd256_pd128(vT)*0.2,A3*0.25);
        vHa = _mm_fmadd_pd(vHa,_mm256_castpd256_pd128(vT),A2*(1.0/3.0));
        vHa = _mm_fmadd_pd(vHa,_mm256_castpd256_pd128(vT),A1*0.5);
        vHa = _mm_fmadd_pd(vHa,_mm256_castpd256_pd128(vT),A0);
        vHa = _mm_fmadd_pd(vHa,_mm256_castpd256_pd128(vT),A5);
        vHa = _mm_mul_pd(RuInvW,vHa);
        store128d(&Ha[i+0],vHa);
    }
    else if(remain==3)
    {
        int i = this->nSpecies-3;
        const std::array<double,7>& a = *this->PtrCoeffs[i+0];
        const std::array<double,7>& b = *this->PtrCoeffs[i+1];
        const std::array<double,7>& c = *this->PtrCoeffs[i+2];

        __m256d A0 = _mm256_setr_pd(a[0],b[0],c[0],0);
        __m256d A1 = _mm256_setr_pd(a[1],b[1],c[1],0);
        __m256d A2 = _mm256_setr_pd(a[2],b[2],c[2],0);
        __m256d A3 = _mm256_setr_pd(a[3],b[3],c[3],0);
        __m256d A4 = _mm256_setr_pd(a[4],b[4],c[4],0);
        __m256d A5 = _mm256_setr_pd(a[5],b[5],c[5],0);
        __m256d A6 = _mm256_setr_pd(a[6],b[6],c[6],0);

        __m256d zerov = _mm256_setzero_pd();
        __m256d RuInvW = _mm256_blend_pd(load256d(&this->invW[i+0]),zerov,0b1000);
        RuInvW = _mm256_mul_pd(RuInvW,Ruv);

        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,vT,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,vT,_mm256_mul_pd(A5,vInvT));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,logT1v,vExpNegGstdByRT);

        store256d(&negGstdByRT[i+0],vExpNegGstdByRT);

        __m256d vCp = _mm256_fmadd_pd(A4 ,vT,A3);
        vCp = _mm256_fmadd_pd(vCp,vT,A2);
        vCp = _mm256_fmadd_pd(vCp,vT,A1);
        vCp = _mm256_fmadd_pd(vCp,vT,A0);
        vCp = _mm256_mul_pd(RuInvW,vCp);
        store256d(&Cp[i+0],vCp);

        __m256d vHa = _mm256_fmadd_pd(A4 ,vT*0.2,A3*0.25);
        vHa = _mm256_fmadd_pd(vHa,vT,A2*(1.0/3.0));
        vHa = _mm256_fmadd_pd(vHa,vT,A1*0.5);
        vHa = _mm256_fmadd_pd(vHa,vT,A0);
        vHa = _mm256_fmadd_pd(vHa,vT,A5);
        vHa = _mm256_mul_pd(RuInvW,vHa);
        store256d(&Ha[i+0],vHa);
    }
}*/

/*void FastChemistry::idealGas::rhoMc
(
    double T,
    const double p,
    double* __restrict__ Phi,
    double* __restrict__ c
)const 
{
    T = T<TlowMin?TlowMin:T;
    T = T>ThighMax?ThighMax:T;

    int remain = this->nSpecies%4;
    for (int i=0; i<this->nSpecies; i++)
    {
        Phi[i] = std::max(Phi[i], 0.0);
    }
    
    this->rhoM = 0;

    double RuTByP = this->Ru*T/p;
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
        this->rhoM += Phi[i]*this->invW[i]*RuTByP;    
    }
    else if(remain==2)
    {
        int i = this->nSpecies-2;
        __m256d YTpv = _mm256_zextpd128_pd256(load128d(&Phi[i]));
        __m256d invWv = _mm256_zextpd128_pd256(load128d(&this->invW[i+0]));
        rhoMv = _mm256_fmadd_pd(_mm256_mul_pd(YTpv,invWv),RuTByPv,rhoMv);
    }
    else
    {
        int i = this->nSpecies-3;
        __m256d zerov = _mm256_setzero_pd();
        __m256d YTpv = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1000);
        __m256d invWv = _mm256_blend_pd(load256d(&this->invW[i]),zerov,0b1000);
        rhoMv = _mm256_fmadd_pd(_mm256_mul_pd(YTpv,invWv),RuTByPv,rhoMv);
    }

    this->rhoM += this->hsum4(rhoMv);
    this->vM = this->rhoM;
    this->rhoM = 1/this->rhoM;

    __m256d rhoMvv = _mm256_set1_pd(this->rhoM);

    for (int i=0; i<this->nSpecies-remain; i=i+4)
    {
        //c[i] = rhoM*this->invW[i]*Phi[i];
        __m256d invWv = load256d(&this->invW[i]);
        __m256d Phiv = load256d(&Phi[i]);
        __m256d cv = _mm256_mul_pd(rhoMvv,_mm256_mul_pd(Phiv,invWv));
        store256d(&c[i],cv);
    }
    if(remain==1)
    {
        int i = this->nSpecies-1;
        c[i] = rhoM*this->invW[i]*Phi[i];
    }
    else if(remain==2)
    {
        int i = this->nSpecies-2;
        __m128d invWv = load128d(&this->invW[i]);
        __m128d Phiv = load128d(&Phi[i]);
        __m128d cv = _mm_mul_pd(_mm256_castpd256_pd128(rhoMvv),_mm_mul_pd(Phiv,invWv));
        store128d(&c[i],cv);
    }
    else
    {
        int i = this->nSpecies-3;
        __m256d zerov = _mm256_setzero_pd();
        __m256d invWv = _mm256_blend_pd(load256d(&this->invW[i]),zerov,0b1000);
        __m256d Phiv = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1000);
        __m256d cv = _mm256_mul_pd(rhoMvv,_mm256_mul_pd(Phiv,invWv));
        store256d(&c[i],cv);
    }

}*/

void FastChemistry::idealGas::DerivativeThermoYT
(
    double T,
    const double p,
    double* __restrict__ Phi,
    double* __restrict__ concentration,
    double* __restrict__ Cp,
    double* __restrict__ Ha,
    double* __restrict__ negGstdByRT,
    double* __restrict__ negGstdByRT2
)const 
{
    //T = T<TlowMin?TlowMin:T;
    //T = T>ThighMax?ThighMax:T;

    this->setPtrCoeffs(T);


    int remain = this->nSpecies%4;
    for (int i=0; i<this->nSpecies; i++)
    {
        Phi[i] = std::max(Phi[i], 0.0);
    }
    
    this->rhoM = 0;

    double RuTByP = this->Ru*T/p;
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
        this->rhoM += Phi[i]*this->invW[i]*RuTByP;    
    }
    else if(remain==2)
    {
        int i = this->nSpecies-2;
        __m256d YTpv = _mm256_zextpd128_pd256(load128d(&Phi[i]));
        __m256d invWv = _mm256_zextpd128_pd256(load128d(&this->invW[i+0]));
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

    this->rhoM += this->hsum4(rhoMv);
    this->vM = this->rhoM;
    this->rhoM = 1/this->rhoM;

    __m256d rhoMvv = _mm256_set1_pd(this->rhoM);

    __m256d zerov = _mm256_setzero_pd();
    for (int i=0; i<this->nSpecies-remain; i=i+4)
    {
        //concentration[i] = rhoM*this->invW[i]*Phi[i];
        __m256d invWv = load256d(&this->invW[i]);
        __m256d Phiv = load256d(&Phi[i]);
        __m256d cv = _mm256_mul_pd(rhoMvv,_mm256_mul_pd(Phiv,invWv));
        cv = _mm256_max_pd(cv,zerov);
        store256d(&concentration[i],cv);
    }
    if(remain==1)
    {
        int i = this->nSpecies-1;
        concentration[i] = std::max(rhoM*this->invW[i]*Phi[i],0.0);
    }
    else if(remain==2)
    {
        int i = this->nSpecies-2;
        __m128d invWv = load128d(&this->invW[i]);
        __m128d Phiv = load128d(&Phi[i]);
        __m128d cv = _mm_mul_pd(_mm256_castpd256_pd128(rhoMvv),_mm_mul_pd(Phiv,invWv));
        cv = _mm_max_pd(cv,_mm256_castpd256_pd128(zerov));
        store128d(&concentration[i],cv);
    }
    else if(remain==3)
    {
        int i = this->nSpecies-3;
        __m256d zerov = _mm256_setzero_pd();
        __m256d invWv = _mm256_blend_pd(load256d(&this->invW[i]),zerov,0b1000);
        __m256d Phiv = _mm256_blend_pd(load256d(&Phi[i]),zerov,0b1000);
        __m256d cv = _mm256_mul_pd(rhoMvv,_mm256_mul_pd(Phiv,invWv));
        cv = _mm256_max_pd(cv,zerov);
        store256d(&concentration[i],cv);
    }

    __m256d vT = _mm256_set1_pd(T);
    __m256d vInvT = _mm256_set1_pd(this->invT);
    __m256d Ruv = _mm256_set1_pd(this->Ru);
    __m256d logT1v = _mm256_set1_pd(this->logT-1);
    __m256d Cpmv = _mm256_setzero_pd();
    double Cpm = 0;
    for(int i = 0; i < this->nSpecies-remain;i=i+4)
    {
        const std::array<double,7>& a = *this->PtrCoeffs[i+0];
        const std::array<double,7>& b = *this->PtrCoeffs[i+1];
        const std::array<double,7>& c = *this->PtrCoeffs[i+2];
        const std::array<double,7>& d = *this->PtrCoeffs[i+3]; 
        __m256d A0 = _mm256_setr_pd(a[0],b[0],c[0],d[0]);
        __m256d A1 = _mm256_setr_pd(a[1],b[1],c[1],d[1]);
        __m256d A2 = _mm256_setr_pd(a[2],b[2],c[2],d[2]);
        __m256d A3 = _mm256_setr_pd(a[3],b[3],c[3],d[3]);
        __m256d A4 = _mm256_setr_pd(a[4],b[4],c[4],d[4]);
        __m256d A5 = _mm256_setr_pd(a[5],b[5],c[5],d[5]);
        __m256d A6 = _mm256_setr_pd(a[6],b[6],c[6],d[6]);

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
    if(remain ==1)
    {
        int i = this->nSpecies-1;
        const std::array<double,7>& a = *this->PtrCoeffs[i];        
        Cp[i] = ((((a[4]*T+a[3])*T+a[2])*T+a[1])*T+a[0])*this->Ru*this->invW[i];
        Ha[i] = (((((a[4]*T*0.2+a[3]*0.25)*T+a[2]*(1.0/3.0))*T+a[1]*0.5)*T+a[0])*T +a[5])*this->Ru*this->invW[i]; 
        double negGstdByRT_ = (((a[4]*0.05*T+a[3]*(1.0/12.0))*T+a[2]*(1.0/6.0))*T+a[1]*0.5)*T-a[5]*invT+a[6]+(this->logT-1)*a[0];
        negGstdByRT[i] = negGstdByRT_;
        negGstdByRT2[i] = negGstdByRT_;
        Cpm = Cpm + Cp[i]*Phi[i];
    }
    else if(remain==2)
    {
        int i = this->nSpecies-2;

        const std::array<double,7>& a = *this->PtrCoeffs[i+0];
        const std::array<double,7>& b = *this->PtrCoeffs[i+1];

        __m128d A0 = _mm_setr_pd(a[0],b[0]);
        __m128d A1 = _mm_setr_pd(a[1],b[1]);
        __m128d A2 = _mm_setr_pd(a[2],b[2]);
        __m128d A3 = _mm_setr_pd(a[3],b[3]);
        __m128d A4 = _mm_setr_pd(a[4],b[4]);
        __m128d A5 = _mm_setr_pd(a[5],b[5]);
        __m128d A6 = _mm_setr_pd(a[6],b[6]);

        __m128d RuInvW = load128d(&this->invW[i+0]);
        RuInvW = _mm_mul_pd(RuInvW,_mm256_castpd256_pd128(Ruv));

        __m128d vExpNegGstdByRT = _mm_fmadd_pd(A4*0.05,_mm256_castpd256_pd128(vT),A3*(1.0/12.0));
        vExpNegGstdByRT = _mm_fmadd_pd(vExpNegGstdByRT,_mm256_castpd256_pd128(vT),A2*(1.0/6.0));
        vExpNegGstdByRT = _mm_fmadd_pd(vExpNegGstdByRT,_mm256_castpd256_pd128(vT),A1*0.5);  
        vExpNegGstdByRT = _mm_fmsub_pd(vExpNegGstdByRT,_mm256_castpd256_pd128(vT),_mm_mul_pd(A5,_mm256_castpd256_pd128(vInvT)));    
        vExpNegGstdByRT = _mm_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm_fmadd_pd(A0,_mm256_castpd256_pd128(logT1v),vExpNegGstdByRT);
        store128d(&negGstdByRT[i+0],vExpNegGstdByRT);
        store128d(&negGstdByRT2[i+0],vExpNegGstdByRT);

        __m128d vCp = _mm_fmadd_pd(A4,_mm256_castpd256_pd128(vT),A3);
        vCp = _mm_fmadd_pd(vCp,_mm256_castpd256_pd128(vT),A2);
        vCp = _mm_fmadd_pd(vCp,_mm256_castpd256_pd128(vT),A1);
        vCp = _mm_fmadd_pd(vCp,_mm256_castpd256_pd128(vT),A0);
        vCp = _mm_mul_pd(RuInvW,vCp);
        store128d(&Cp[i+0],vCp);
        __m256d Yv = _mm256_zextpd128_pd256(load128d(&Phi[i+0]));
        Cpmv = _mm256_fmadd_pd(Yv,_mm256_zextpd128_pd256(vCp),Cpmv);

        __m128d vHa = _mm_fmadd_pd(A4 ,_mm256_castpd256_pd128(vT)*0.2,A3*0.25);
        vHa = _mm_fmadd_pd(vHa,_mm256_castpd256_pd128(vT),A2*(1.0/3.0));
        vHa = _mm_fmadd_pd(vHa,_mm256_castpd256_pd128(vT),A1*0.5);
        vHa = _mm_fmadd_pd(vHa,_mm256_castpd256_pd128(vT),A0);
        vHa = _mm_fmadd_pd(vHa,_mm256_castpd256_pd128(vT),A5);
        vHa = _mm_mul_pd(RuInvW,vHa);
        store128d(&Ha[i+0],vHa);
    }
    else if(remain==3)
    {
        int i = this->nSpecies-3;
        const std::array<double,7>& a = *this->PtrCoeffs[i+0];
        const std::array<double,7>& b = *this->PtrCoeffs[i+1];
        const std::array<double,7>& c = *this->PtrCoeffs[i+2];

        __m256d A0 = _mm256_setr_pd(a[0],b[0],c[0],0);
        __m256d A1 = _mm256_setr_pd(a[1],b[1],c[1],0);
        __m256d A2 = _mm256_setr_pd(a[2],b[2],c[2],0);
        __m256d A3 = _mm256_setr_pd(a[3],b[3],c[3],0);
        __m256d A4 = _mm256_setr_pd(a[4],b[4],c[4],0);
        __m256d A5 = _mm256_setr_pd(a[5],b[5],c[5],0);
        __m256d A6 = _mm256_setr_pd(a[6],b[6],c[6],0);

        __m256d zerov = _mm256_setzero_pd();
        __m256d RuInvW = _mm256_blend_pd(load256d(&this->invW[i+0]),zerov,0b1000);
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
        __m256d Yv = _mm256_blend_pd(load256d(&Phi[i+0]),_mm256_setzero_pd(),0b1000);
        Cpmv = _mm256_fmadd_pd(Yv,vCp,Cpmv);

        __m256d vHa = _mm256_fmadd_pd(A4 ,vT*0.2,A3*0.25);
        vHa = _mm256_fmadd_pd(vHa,vT,A2*(1.0/3.0));
        vHa = _mm256_fmadd_pd(vHa,vT,A1*0.5);
        vHa = _mm256_fmadd_pd(vHa,vT,A0);
        vHa = _mm256_fmadd_pd(vHa,vT,A5);
        vHa = _mm256_mul_pd(RuInvW,vHa);
        store256d(&Ha[i+0],vHa);
    }
    Cpm = Cpm + this->hsum4(Cpmv);
    Cp[this->nSpecies] = Cpm;
}

void FastChemistry::idealGas::negGstdByRT
(
    double T,
    double* __restrict__ negGstdByRT,
    double* __restrict__ negGstdByRT2
)const 
{


    this->setPtrCoeffs(T);
    int remain = this->nSpecies%4;
    __m256d vT = _mm256_set1_pd(T);
    __m256d vInvT = _mm256_set1_pd(invT);
    __m256d Ruv = _mm256_set1_pd(this->Ru);
    __m256d logT1v = _mm256_set1_pd(this->logT-1);
    for(int i = 0; i < this->nSpecies-remain;i=i+4)
    {
        const std::array<double,7>& a = *this->PtrCoeffs[i+0];
        const std::array<double,7>& b = *this->PtrCoeffs[i+1];
        const std::array<double,7>& c = *this->PtrCoeffs[i+2];
        const std::array<double,7>& d = *this->PtrCoeffs[i+3]; 
        __m256d A0 = _mm256_setr_pd(a[0],b[0],c[0],d[0]);
        __m256d A1 = _mm256_setr_pd(a[1],b[1],c[1],d[1]);
        __m256d A2 = _mm256_setr_pd(a[2],b[2],c[2],d[2]);
        __m256d A3 = _mm256_setr_pd(a[3],b[3],c[3],d[3]);
        __m256d A4 = _mm256_setr_pd(a[4],b[4],c[4],d[4]);
        __m256d A5 = _mm256_setr_pd(a[5],b[5],c[5],d[5]);
        __m256d A6 = _mm256_setr_pd(a[6],b[6],c[6],d[6]);

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
    }
    if(remain ==1)
    {
        int i = this->nSpecies-1;
        const std::array<double,7>& a = *this->PtrCoeffs[i];        
        double negGstdByRT_ = (((a[4]*0.05*T+a[3]*(1.0/12.0))*T+a[2]*(1.0/6.0))*T+a[1]*0.5)*T-a[5]*invT+a[6]+(this->logT-1)*a[0];
        negGstdByRT[i] = negGstdByRT_;
        negGstdByRT2[i] = negGstdByRT_;
    }
    else if(remain==2)
    {
        int i = this->nSpecies-2;

        const std::array<double,7>& a = *this->PtrCoeffs[i+0];
        const std::array<double,7>& b = *this->PtrCoeffs[i+1];

        __m128d A0 = _mm_setr_pd(a[0],b[0]);
        __m128d A1 = _mm_setr_pd(a[1],b[1]);
        __m128d A2 = _mm_setr_pd(a[2],b[2]);
        __m128d A3 = _mm_setr_pd(a[3],b[3]);
        __m128d A4 = _mm_setr_pd(a[4],b[4]);
        __m128d A5 = _mm_setr_pd(a[5],b[5]);
        __m128d A6 = _mm_setr_pd(a[6],b[6]);

        __m128d RuInvW = load128d(&this->invW[i+0]);
        RuInvW = _mm_mul_pd(RuInvW,_mm256_castpd256_pd128(Ruv));

        __m128d vExpNegGstdByRT = _mm_fmadd_pd(A4*0.05,_mm256_castpd256_pd128(vT),A3*(1.0/12.0));
        vExpNegGstdByRT = _mm_fmadd_pd(vExpNegGstdByRT,_mm256_castpd256_pd128(vT),A2*(1.0/6.0));
        vExpNegGstdByRT = _mm_fmadd_pd(vExpNegGstdByRT,_mm256_castpd256_pd128(vT),A1*0.5);  
        vExpNegGstdByRT = _mm_fmsub_pd(vExpNegGstdByRT,_mm256_castpd256_pd128(vT),_mm_mul_pd(A5,_mm256_castpd256_pd128(vInvT)));    
        vExpNegGstdByRT = _mm_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm_fmadd_pd(A0,_mm256_castpd256_pd128(logT1v),vExpNegGstdByRT);
        store128d(&negGstdByRT[i+0],vExpNegGstdByRT);
        store128d(&negGstdByRT2[i+0],vExpNegGstdByRT);
    }
    else if(remain==3)
    {
        int i = this->nSpecies-3;
        const std::array<double,7>& a = *this->PtrCoeffs[i+0];
        const std::array<double,7>& b = *this->PtrCoeffs[i+1];
        const std::array<double,7>& c = *this->PtrCoeffs[i+2];

        __m256d A0 = _mm256_setr_pd(a[0],b[0],c[0],0);
        __m256d A1 = _mm256_setr_pd(a[1],b[1],c[1],0);
        __m256d A2 = _mm256_setr_pd(a[2],b[2],c[2],0);
        __m256d A3 = _mm256_setr_pd(a[3],b[3],c[3],0);
        __m256d A4 = _mm256_setr_pd(a[4],b[4],c[4],0);
        __m256d A5 = _mm256_setr_pd(a[5],b[5],c[5],0);
        __m256d A6 = _mm256_setr_pd(a[6],b[6],c[6],0);

        __m256d zerov = _mm256_setzero_pd();
        __m256d RuInvW = _mm256_blend_pd(load256d(&this->invW[i+0]),zerov,0b1000);
        RuInvW = _mm256_mul_pd(RuInvW,Ruv);

        __m256d vExpNegGstdByRT = _mm256_fmadd_pd(A4*0.05,vT,A3*(1.0/12.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A2*(1.0/6.0));
        vExpNegGstdByRT = _mm256_fmadd_pd(vExpNegGstdByRT,vT,A1*0.5);  
        vExpNegGstdByRT = _mm256_fmsub_pd(vExpNegGstdByRT,vT,_mm256_mul_pd(A5,vInvT));    
        vExpNegGstdByRT = _mm256_add_pd(vExpNegGstdByRT,A6);
        vExpNegGstdByRT = _mm256_fmadd_pd(A0,logT1v,vExpNegGstdByRT);
        store256d(&negGstdByRT[i+0],vExpNegGstdByRT);
        store256d(&negGstdByRT2[i+0],vExpNegGstdByRT);
    }
}

inline void FastChemistry::idealGas::setPtrCoeffs(double T) const
{
    if(T<this->TcommonMin)
    {
        for(int i = 0; i < this->nSpecies;i++)
        {
            this->PtrCoeffs[i] = &this->LCoeffs[i];
        }
    }
    else if(T>this->TcommonMax)
    {
        for(int i = 0; i < this->nSpecies;i++)
        {
            this->PtrCoeffs[i] = &this->HCoeffs[i];
        }               
    }
    else 
    {
        for(int i = 0; i < this->nSpecies;i++)
        {
            if(T < this->Tcommon[i])
            {
                this->PtrCoeffs[i] = &this->LCoeffs[i];              
            }
            else
            {
                this->PtrCoeffs[i] = &this->HCoeffs[i];                
            }
        }
    }    
}

inline double FastChemistry::idealGas::hsum4(__m256d v) const
{
    __m128d lo = _mm256_castpd256_pd128(v);
    __m128d hi = _mm256_extractf128_pd(v,1);
    lo = _mm_add_pd(lo,hi);               // [v0+v2, v1+v3]
    __m128d sh = _mm_unpackhi_pd(lo,lo);  // [v1+v3, v1+v3]
    __m128d s = _mm_add_sd(lo, sh);       // [sum, …]
    return _mm_cvtsd_f64(s);
}

inline __m256d FastChemistry::idealGas::hsum4x4(__m256d sum0, __m256d sum1,__m256d sum2, __m256d sum3) const
{
    __m256d t0 = _mm256_hadd_pd(sum0, sum1);  // [a0+a1, b0+b1, a2+a3, b2+b3]
    __m256d t1 = _mm256_hadd_pd(sum2, sum3);  // [c0+c1, d0+d1, c2+c3, d2+d3]
    __m256d t2 = _mm256_permute4x64_pd (t0, _MM_SHUFFLE(3,1,2,0)); // [a0+a1, a2+a3, b0+b1,  b2+b3]
    __m256d t3 = _mm256_permute4x64_pd(t1, _MM_SHUFFLE(3,1,2,0)); // [c0+c1, c2+c3, d0+d1, d2+d3]
    __m256d t4 = _mm256_hadd_pd(t2, t3); // [a0+a1+a2+a3, c0+c1+c2+c3, b0+b1+b2+b3, d0+d1+d2+d3]
    __m256d t5 = _mm256_permute4x64_pd(t4, _MM_SHUFFLE(3,1,2,0));
    return t5;
}

inline void FastChemistry::idealGas::transpose4x4_pd(__m256d& v0,__m256d& v1,__m256d& v2,__m256d& v3)const
{
    __m256d t0 = _mm256_unpacklo_pd(v0, v1); // a0 b0 a1 b1
    __m256d t1 = _mm256_unpackhi_pd(v0, v1); // a2 b2 a3 b3
    __m256d t2 = _mm256_unpacklo_pd(v2, v3); // c0 d0 c1 d1
    __m256d t3 = _mm256_unpackhi_pd(v2, v3); // c2 d2 c3 d3

    v0 = _mm256_permute2f128_pd(t0, t2, 0x20); // [a0 b0 c0 d0]
    v1 = _mm256_permute2f128_pd(t1, t3, 0x20); // [a1 b1 c1 d1]
    v2 = _mm256_permute2f128_pd(t0, t2, 0x31); // [a2 b2 c2 d2]
    v3 = _mm256_permute2f128_pd(t1, t3, 0x31); // [a3 b3 c3 d3]
}
        
inline double FastChemistry::idealGas::get_elem0(__m256d vec) const
{
    return _mm256_cvtsd_f64(vec); 
}

inline double FastChemistry::idealGas::get_elem1(__m256d vec) const
{
    __m128d low = _mm256_castpd256_pd128(vec); 
    return _mm_cvtsd_f64(_mm_unpackhi_pd(low, low));
}

inline double FastChemistry::idealGas::get_elem2(__m256d vec) const
{
    __m128d high = _mm256_extractf128_pd(vec, 1); 
    return _mm_cvtsd_f64(high);
}

inline double FastChemistry::idealGas::get_elem3(__m256d vec) const
{
    __m128d high = _mm256_extractf128_pd(vec, 1); 
    return _mm_cvtsd_f64(_mm_unpackhi_pd(high, high));
}