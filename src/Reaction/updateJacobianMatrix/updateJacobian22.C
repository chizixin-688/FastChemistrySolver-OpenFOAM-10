#include "OptReaction.H"
#include <immintrin.h>  

void  OptReaction::updateJacobian22
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = reaction22index.size();

    std::size_t lhsIndex = 0;

    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->reaction22index[k];

 
        const unsigned int sl0 = lhsSpeciesIndex1D22[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D22[lhsIndex+1];
        
        const unsigned int sr0 = rhsSpeciesIndex1D22[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D22[rhsIndex+1];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];
        double Kr = 0;
        double dKrdT = 0;
        double invKc = 0;


        if(this->isIrreversible[i]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
            double Kc = Kp;  
            Kc = Kc > KcLimiter?Kc:KcLimiter;
            const double sumVdBdT = (dBdT[sr0] + dBdT[sr1]) - (dBdT[sl0] + dBdT[sl1]);
            const double dKcdTByKc = sumVdBdT;    
            invKc = 1.0/Kc;    
            Kr = Kf*invKc;
            dKrdT = (dKfdT*invKc - (Kc > KcLimiter ? Kr*dKcdTByKc : 0));     

            const double dCrdC0 = C[sr1];
            const double dCrdC1 = C[sr0];
            ddNdtByVdcTp[sl0*(this->alignN)+ sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sl1*(this->alignN)+ sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sl0*(this->alignN)+ sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sl1*(this->alignN)+ sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sr0*(this->alignN)+ sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr1*(this->alignN)+ sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr0*(this->alignN)+ sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sr1*(this->alignN)+ sr1] += (-Kr*dCrdC1);        
        }
        else if(this->isIrreversible[i]==2)
        {
            Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
            dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
            invKc = Kr/Kf;
            const double dCrdC0 = C[sr1];
            const double dCrdC1 = C[sr0]; 
            ddNdtByVdcTp[sl0*(this->alignN)+ sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sl1*(this->alignN)+ sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sl0*(this->alignN)+ sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sl1*(this->alignN)+ sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sr0*(this->alignN)+ sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr1*(this->alignN)+ sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr0*(this->alignN)+ sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sr1*(this->alignN)+ sr1] += (-Kr*dCrdC1);
        }

        const double CF = C[sl0]*C[sl1];
        const double CR = C[sr0]*C[sr1];

        const double q = (Kf*CF) - (Kr*CR);
        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;

        const double dqdT = (dKfdT*CF)-(dKrdT*CR);
        ddNdtByVdcTp[sl0*(this->alignN)+ (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sl1*(this->alignN)+ (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN)+ (this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr1*(this->alignN)+ (this->nSpecies)] += dqdT;

        const double dCfdC0 = C[sl1];
        const double dCfdC1 = C[sl0];
        ddNdtByVdcTp[sl0*(this->alignN)+ sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl1*(this->alignN)+ sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl0*(this->alignN)+ sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sl1*(this->alignN)+ sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sr0*(this->alignN)+ sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr1*(this->alignN)+ sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr0*(this->alignN)+ sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sr1*(this->alignN)+ sl1] += Kf*dCfdC1;     

        if(i>=this->Ikf[2] && i <this->Ikf[6])
        { 
            const unsigned int k = i - this->Ikf[2];

            __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[k]);
            __m256d CF_ = _mm256_set1_pd(CF);            
            __m256d CR_ = _mm256_set1_pd(CR);  
            double dKrdC = 0;
            dKrdC = (this->dKfdC_[k]*invKc);
            __m256d dKrdC_ = _mm256_set1_pd(dKrdC);

            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[k*this->AlignSpecies+j]);

                __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0*(this->alignN)+j],sl0v);

                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1*(this->alignN)+j],sl1v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1*(this->alignN)+j],sr1v);
            } 
        }

        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 2;

    }
}   


void  OptReaction::updateJacobian22ReversibleReaction
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = reversibleReaction22index.size();

    std::size_t remain = end%4;

    std::size_t lhsIndex = 0;

    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->reversibleReaction22index[k+0];
        const unsigned int i1 = this->reversibleReaction22index[k+1];
        const unsigned int i2 = this->reversibleReaction22index[k+2];
        const unsigned int i3 = this->reversibleReaction22index[k+3];

 
        const unsigned int sl0a = lhsSpeciesIndex1D22RR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D22RR[lhsIndex+1];
        
        const unsigned int sr0a = rhsSpeciesIndex1D22RR[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D22RR[rhsIndex+1];

        double Kf0 = this->Kf_[i0];
        double dKfdT0 = this->dKfdT_[i0];
        double Kr0 = 0;
        double dKrdT0 = 0;
        double invKc0 = 0;

        
            const double Kp0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
            double Kc0 = Kp0;  
            Kc0 = Kc0 > KcLimiter?Kc0:KcLimiter;
            const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a]) - (dBdT[sl0a] + dBdT[sl1a]);
            const double dKcdTByKc0 = sumVdBdT0;    
            invKc0 = 1.0/Kc0;    
            Kr0 = Kf0*invKc0;
            dKrdT0 = (dKfdT0*invKc0 - (Kc0 > KcLimiter ? Kr0*dKcdTByKc0 : 0));     

            const double dCrdC0a = C[sr1a];
            const double dCrdC1a = C[sr0a];
            ddNdtByVdcTp[sl0a*(this->alignN)+ sr0a] -= (-Kr0*dCrdC0a);
            ddNdtByVdcTp[sl1a*(this->alignN)+ sr0a] -= (-Kr0*dCrdC0a);
            ddNdtByVdcTp[sl0a*(this->alignN)+ sr1a] -= (-Kr0*dCrdC1a);
            ddNdtByVdcTp[sl1a*(this->alignN)+ sr1a] -= (-Kr0*dCrdC1a);
            ddNdtByVdcTp[sr0a*(this->alignN)+ sr0a] += (-Kr0*dCrdC0a);
            ddNdtByVdcTp[sr1a*(this->alignN)+ sr0a] += (-Kr0*dCrdC0a);
            ddNdtByVdcTp[sr0a*(this->alignN)+ sr1a] += (-Kr0*dCrdC1a);
            ddNdtByVdcTp[sr1a*(this->alignN)+ sr1a] += (-Kr0*dCrdC1a);        
        


        const double CFa = C[sl0a]*C[sl1a];
        const double CRa = C[sr0a]*C[sr1a];

        const double qa = (Kf0*CFa) - (Kr0*CRa);
        dNdtByV[sl0a] = dNdtByV[sl0a] - qa;
        dNdtByV[sl1a] = dNdtByV[sl1a] - qa;
        dNdtByV[sr0a] = dNdtByV[sr0a] + qa;
        dNdtByV[sr1a] = dNdtByV[sr1a] + qa;

        const double dqdTa = (dKfdT0*CFa)-(dKrdT0*CRa);
        ddNdtByVdcTp[sl0a*(this->alignN)+ (this->nSpecies)] -= dqdTa;
        ddNdtByVdcTp[sl1a*(this->alignN)+ (this->nSpecies)] -= dqdTa;
        ddNdtByVdcTp[sr0a*(this->alignN)+ (this->nSpecies)] += dqdTa;
        ddNdtByVdcTp[sr1a*(this->alignN)+ (this->nSpecies)] += dqdTa;

        const double dCfdC0a = C[sl1a];
        const double dCfdC1a = C[sl0a];
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl1a] += Kf0*dCfdC1a;     

        if(i0>=this->Ikf[2] && i0 <this->Ikf[6])
        { 
            const unsigned int k = i0 - this->Ikf[2];

            __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[k]);
            __m256d CF_ = _mm256_set1_pd(CFa);            
            __m256d CR_ = _mm256_set1_pd(CRa);  
            double dKrdC = 0;
            dKrdC = (this->dKfdC_[k]*invKc0);
            __m256d dKrdC_ = _mm256_set1_pd(dKrdC);

            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[k*this->AlignSpecies+j]);

                __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);

                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j],sr0v);

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j],sr1v);
            } 
        }

        const unsigned int sl0b = lhsSpeciesIndex1D22RR[lhsIndex+2];
        const unsigned int sl1b = lhsSpeciesIndex1D22RR[lhsIndex+3];
        
        const unsigned int sr0b = rhsSpeciesIndex1D22RR[rhsIndex+2];
        const unsigned int sr1b = rhsSpeciesIndex1D22RR[rhsIndex+3];

        double Kf1 = this->Kf_[i1];
        double dKfdT1 = this->dKfdT_[i1];
        double Kr1 = 0;
        double dKrdT1 = 0;
        double invKc1 = 0;

        
            const double Kp1 = (ExpNegGbyRT[sr0b]*ExpNegGbyRT[sr1b])/(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]);
            double Kc1 = Kp1;  
            Kc1 = Kc1 > KcLimiter?Kc1:KcLimiter;
            const double sumVdBdT1 = (dBdT[sr0b] + dBdT[sr1b]) - (dBdT[sl0b] + dBdT[sl1b]);
            const double dKcdTByKc1 = sumVdBdT1;    
            invKc1 = 1.0/Kc1;    
            Kr1 = Kf1*invKc1;
            dKrdT1 = (dKfdT1*invKc1 - (Kc1 > KcLimiter ? Kr1*dKcdTByKc1 : 0));     

            const double dCrdC0b = C[sr1b];
            const double dCrdC1b = C[sr0b];
            ddNdtByVdcTp[sl0b*(this->alignN)+ sr0b] -= (-Kr1*dCrdC0b);
            ddNdtByVdcTp[sl1b*(this->alignN)+ sr0b] -= (-Kr1*dCrdC0b);
            ddNdtByVdcTp[sl0b*(this->alignN)+ sr1b] -= (-Kr1*dCrdC1b);
            ddNdtByVdcTp[sl1b*(this->alignN)+ sr1b] -= (-Kr1*dCrdC1b);
            ddNdtByVdcTp[sr0b*(this->alignN)+ sr0b] += (-Kr1*dCrdC0b);
            ddNdtByVdcTp[sr1b*(this->alignN)+ sr0b] += (-Kr1*dCrdC0b);
            ddNdtByVdcTp[sr0b*(this->alignN)+ sr1b] += (-Kr1*dCrdC1b);
            ddNdtByVdcTp[sr1b*(this->alignN)+ sr1b] += (-Kr1*dCrdC1b);        
        


        const double CFb = C[sl0b]*C[sl1b];
        const double CRb = C[sr0b]*C[sr1b];

        const double qb = (Kf1*CFb) - (Kr1*CRb);
        dNdtByV[sl0b] = dNdtByV[sl0b] - qb;
        dNdtByV[sl1b] = dNdtByV[sl1b] - qb;
        dNdtByV[sr0b] = dNdtByV[sr0b] + qb;
        dNdtByV[sr1b] = dNdtByV[sr1b] + qb;

        const double dqdTb = (dKfdT1*CFb)-(dKrdT1*CRb);
        ddNdtByVdcTp[sl0b*(this->alignN)+ (this->nSpecies)] -= dqdTb;
        ddNdtByVdcTp[sl1b*(this->alignN)+ (this->nSpecies)] -= dqdTb;
        ddNdtByVdcTp[sr0b*(this->alignN)+ (this->nSpecies)] += dqdTb;
        ddNdtByVdcTp[sr1b*(this->alignN)+ (this->nSpecies)] += dqdTb;

        const double dCfdC0b = C[sl1b];
        const double dCfdC1b = C[sl0b];
        ddNdtByVdcTp[sl0b*(this->alignN)+ sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl1b*(this->alignN)+ sl0b] -= Kf1*dCfdC0b;
        ddNdtByVdcTp[sl0b*(this->alignN)+ sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sl1b*(this->alignN)+ sl1b] -= Kf1*dCfdC1b;
        ddNdtByVdcTp[sr0b*(this->alignN)+ sl0b] += Kf1*dCfdC0b;
        ddNdtByVdcTp[sr1b*(this->alignN)+ sl0b] += Kf1*dCfdC0b;
        ddNdtByVdcTp[sr0b*(this->alignN)+ sl1b] += Kf1*dCfdC1b;
        ddNdtByVdcTp[sr1b*(this->alignN)+ sl1b] += Kf1*dCfdC1b;     

        if(i1>=this->Ikf[2] && i1 <this->Ikf[6])
        { 
            const unsigned int k = i1 - this->Ikf[2];

            __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[k]);
            __m256d CF_ = _mm256_set1_pd(CFb);            
            __m256d CR_ = _mm256_set1_pd(CRb);  
            double dKrdC = 0;
            dKrdC = (this->dKfdC_[k]*invKc1);
            __m256d dKrdC_ = _mm256_set1_pd(dKrdC);

            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[k*this->AlignSpecies+j]);

                __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0b*(this->alignN)+j],sl0v);

                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1b*(this->alignN)+j],sl1v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0b*(this->alignN)+j],sr0v);

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1b*(this->alignN)+j],sr1v);
            } 
        }


        const unsigned int sl0c = lhsSpeciesIndex1D22RR[lhsIndex+4];
        const unsigned int sl1c = lhsSpeciesIndex1D22RR[lhsIndex+5];
        
        const unsigned int sr0c = rhsSpeciesIndex1D22RR[rhsIndex+4];
        const unsigned int sr1c = rhsSpeciesIndex1D22RR[rhsIndex+5];

        double Kf2 = this->Kf_[i2];
        double dKfdT2 = this->dKfdT_[i2];
        double Kr2 = 0;
        double dKrdT2 = 0;
        double invKc2 = 0;

        
            const double Kp2 = (ExpNegGbyRT[sr0c]*ExpNegGbyRT[sr1c])/(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]);
            double Kc2 = Kp2;  
            Kc2 = Kc2 > KcLimiter?Kc2:KcLimiter;
            const double sumVdBdT2 = (dBdT[sr0c] + dBdT[sr1c]) - (dBdT[sl0c] + dBdT[sl1c]);
            const double dKcdTByKc2 = sumVdBdT2;    
            invKc2 = 1.0/Kc2;    
            Kr2 = Kf2*invKc2;
            dKrdT2 = (dKfdT2*invKc2 - (Kc2 > KcLimiter ? Kr2*dKcdTByKc2 : 0));     

            const double dCrdC0c = C[sr1c];
            const double dCrdC1c = C[sr0c];
            ddNdtByVdcTp[sl0c*(this->alignN)+ sr0c] -= (-Kr2*dCrdC0c);
            ddNdtByVdcTp[sl1c*(this->alignN)+ sr0c] -= (-Kr2*dCrdC0c);
            ddNdtByVdcTp[sl0c*(this->alignN)+ sr1c] -= (-Kr2*dCrdC1c);
            ddNdtByVdcTp[sl1c*(this->alignN)+ sr1c] -= (-Kr2*dCrdC1c);
            ddNdtByVdcTp[sr0c*(this->alignN)+ sr0c] += (-Kr2*dCrdC0c);
            ddNdtByVdcTp[sr1c*(this->alignN)+ sr0c] += (-Kr2*dCrdC0c);
            ddNdtByVdcTp[sr0c*(this->alignN)+ sr1c] += (-Kr2*dCrdC1c);
            ddNdtByVdcTp[sr1c*(this->alignN)+ sr1c] += (-Kr2*dCrdC1c);        
        


        const double CFc = C[sl0c]*C[sl1c];
        const double CRc = C[sr0c]*C[sr1c];

        const double qc = (Kf2*CFc) - (Kr2*CRc);
        dNdtByV[sl0c] = dNdtByV[sl0c] - qc;
        dNdtByV[sl1c] = dNdtByV[sl1c] - qc;
        dNdtByV[sr0c] = dNdtByV[sr0c] + qc;
        dNdtByV[sr1c] = dNdtByV[sr1c] + qc;

        const double dqdTc = (dKfdT2*CFc)-(dKrdT2*CRc);
        ddNdtByVdcTp[sl0c*(this->alignN)+ (this->nSpecies)] -= dqdTc;
        ddNdtByVdcTp[sl1c*(this->alignN)+ (this->nSpecies)] -= dqdTc;
        ddNdtByVdcTp[sr0c*(this->alignN)+ (this->nSpecies)] += dqdTc;
        ddNdtByVdcTp[sr1c*(this->alignN)+ (this->nSpecies)] += dqdTc;

        const double dCfdC0c = C[sl1c];
        const double dCfdC1c = C[sl0c];
        ddNdtByVdcTp[sl0c*(this->alignN)+ sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl1c*(this->alignN)+ sl0c] -= Kf2*dCfdC0c;
        ddNdtByVdcTp[sl0c*(this->alignN)+ sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sl1c*(this->alignN)+ sl1c] -= Kf2*dCfdC1c;
        ddNdtByVdcTp[sr0c*(this->alignN)+ sl0c] += Kf2*dCfdC0c;
        ddNdtByVdcTp[sr1c*(this->alignN)+ sl0c] += Kf2*dCfdC0c;
        ddNdtByVdcTp[sr0c*(this->alignN)+ sl1c] += Kf2*dCfdC1c;
        ddNdtByVdcTp[sr1c*(this->alignN)+ sl1c] += Kf2*dCfdC1c;     

        if(i2>=this->Ikf[2] && i2 <this->Ikf[6])
        { 
            const unsigned int k = i2 - this->Ikf[2];

            __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[k]);
            __m256d CF_ = _mm256_set1_pd(CFc);            
            __m256d CR_ = _mm256_set1_pd(CRc);  
            double dKrdC = 0;
            dKrdC = (this->dKfdC_[k]*invKc2);
            __m256d dKrdC_ = _mm256_set1_pd(dKrdC);

            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[k*this->AlignSpecies+j]);

                __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0c*(this->alignN)+j],sl0v);

                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1c*(this->alignN)+j],sl1v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0c*(this->alignN)+j],sr0v);

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1c*(this->alignN)+j],sr1v);
            } 
        }


        const unsigned int sl0d = lhsSpeciesIndex1D22RR[lhsIndex+6];
        const unsigned int sl1d = lhsSpeciesIndex1D22RR[lhsIndex+7];
        
        const unsigned int sr0d = rhsSpeciesIndex1D22RR[rhsIndex+6];
        const unsigned int sr1d = rhsSpeciesIndex1D22RR[rhsIndex+7];

        double Kf3 = this->Kf_[i3];
        double dKfdT3 = this->dKfdT_[i3];
        double Kr3 = 0;
        double dKrdT3 = 0;
        double invKc3 = 0;

        
            const double Kp3 = (ExpNegGbyRT[sr0d]*ExpNegGbyRT[sr1d])/(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]);
            double Kc3 = Kp3;  
            Kc3 = Kc3 > KcLimiter?Kc3:KcLimiter;
            const double sumVdBdT3 = (dBdT[sr0d] + dBdT[sr1d]) - (dBdT[sl0d] + dBdT[sl1d]);
            const double dKcdTByKc3 = sumVdBdT3;    
            invKc3 = 1.0/Kc3;    
            Kr3 = Kf3*invKc3;
            dKrdT3 = (dKfdT3*invKc3 - (Kc3 > KcLimiter ? Kr3*dKcdTByKc3 : 0));     

            const double dCrdC0d = C[sr1d];
            const double dCrdC1d = C[sr0d];
            ddNdtByVdcTp[sl0d*(this->alignN)+ sr0d] -= (-Kr3*dCrdC0d);
            ddNdtByVdcTp[sl1d*(this->alignN)+ sr0d] -= (-Kr3*dCrdC0d);
            ddNdtByVdcTp[sl0d*(this->alignN)+ sr1d] -= (-Kr3*dCrdC1d);
            ddNdtByVdcTp[sl1d*(this->alignN)+ sr1d] -= (-Kr3*dCrdC1d);
            ddNdtByVdcTp[sr0d*(this->alignN)+ sr0d] += (-Kr3*dCrdC0d);
            ddNdtByVdcTp[sr1d*(this->alignN)+ sr0d] += (-Kr3*dCrdC0d);
            ddNdtByVdcTp[sr0d*(this->alignN)+ sr1d] += (-Kr3*dCrdC1d);
            ddNdtByVdcTp[sr1d*(this->alignN)+ sr1d] += (-Kr3*dCrdC1d);        
        


        const double CFd = C[sl0d]*C[sl1d];
        const double CRd = C[sr0d]*C[sr1d];

        const double qd = (Kf3*CFd) - (Kr3*CRd);
        dNdtByV[sl0d] = dNdtByV[sl0d] - qd;
        dNdtByV[sl1d] = dNdtByV[sl1d] - qd;
        dNdtByV[sr0d] = dNdtByV[sr0d] + qd;
        dNdtByV[sr1d] = dNdtByV[sr1d] + qd;

        const double dqdTd = (dKfdT3*CFd)-(dKrdT3*CRd);
        ddNdtByVdcTp[sl0d*(this->alignN)+ (this->nSpecies)] -= dqdTd;
        ddNdtByVdcTp[sl1d*(this->alignN)+ (this->nSpecies)] -= dqdTd;
        ddNdtByVdcTp[sr0d*(this->alignN)+ (this->nSpecies)] += dqdTd;
        ddNdtByVdcTp[sr1d*(this->alignN)+ (this->nSpecies)] += dqdTd;

        const double dCfdC0d = C[sl1d];
        const double dCfdC1d = C[sl0d];
        ddNdtByVdcTp[sl0d*(this->alignN)+ sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl1d*(this->alignN)+ sl0d] -= Kf3*dCfdC0d;
        ddNdtByVdcTp[sl0d*(this->alignN)+ sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sl1d*(this->alignN)+ sl1d] -= Kf3*dCfdC1d;
        ddNdtByVdcTp[sr0d*(this->alignN)+ sl0d] += Kf3*dCfdC0d;
        ddNdtByVdcTp[sr1d*(this->alignN)+ sl0d] += Kf3*dCfdC0d;
        ddNdtByVdcTp[sr0d*(this->alignN)+ sl1d] += Kf3*dCfdC1d;
        ddNdtByVdcTp[sr1d*(this->alignN)+ sl1d] += Kf3*dCfdC1d;     

        if(i3>=this->Ikf[2] && i3 <this->Ikf[6])
        { 
            const unsigned int k = i3 - this->Ikf[2];

            __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[k]);
            __m256d CF_ = _mm256_set1_pd(CFd);            
            __m256d CR_ = _mm256_set1_pd(CRd);  
            double dKrdC = 0;
            dKrdC = (this->dKfdC_[k]*invKc3);
            __m256d dKrdC_ = _mm256_set1_pd(dKrdC);

            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[k*this->AlignSpecies+j]);

                __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0d*(this->alignN)+j],sl0v);

                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1d*(this->alignN)+j],sl1v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0d*(this->alignN)+j],sr0v);

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1d*(this->alignN)+j],sr1v);
            } 
        }
        lhsIndex = lhsIndex + 8;
        rhsIndex = rhsIndex + 8;
    }
    for(std::size_t k=end-remain; k<end; k=k+1)
    {
        const unsigned int i0 = this->reversibleReaction22index[k+0];


 
        const unsigned int sl0a = lhsSpeciesIndex1D22RR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D22RR[lhsIndex+1];
        
        const unsigned int sr0a = rhsSpeciesIndex1D22RR[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D22RR[rhsIndex+1];

        double Kf0 = this->Kf_[i0];
        double dKfdT0 = this->dKfdT_[i0];
        double Kr0 = 0;
        double dKrdT0 = 0;
        double invKc0 = 0;

        
            const double Kp0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
            double Kc0 = Kp0;  
            Kc0 = Kc0 > KcLimiter?Kc0:KcLimiter;
            const double sumVdBdT0 = (dBdT[sr0a] + dBdT[sr1a]) - (dBdT[sl0a] + dBdT[sl1a]);
            const double dKcdTByKc0 = sumVdBdT0;    
            invKc0 = 1.0/Kc0;    
            Kr0 = Kf0*invKc0;
            dKrdT0 = (dKfdT0*invKc0 - (Kc0 > KcLimiter ? Kr0*dKcdTByKc0 : 0));     

            const double dCrdC0a = C[sr1a];
            const double dCrdC1a = C[sr0a];
            ddNdtByVdcTp[sl0a*(this->alignN)+ sr0a] -= (-Kr0*dCrdC0a);
            ddNdtByVdcTp[sl1a*(this->alignN)+ sr0a] -= (-Kr0*dCrdC0a);
            ddNdtByVdcTp[sl0a*(this->alignN)+ sr1a] -= (-Kr0*dCrdC1a);
            ddNdtByVdcTp[sl1a*(this->alignN)+ sr1a] -= (-Kr0*dCrdC1a);
            ddNdtByVdcTp[sr0a*(this->alignN)+ sr0a] += (-Kr0*dCrdC0a);
            ddNdtByVdcTp[sr1a*(this->alignN)+ sr0a] += (-Kr0*dCrdC0a);
            ddNdtByVdcTp[sr0a*(this->alignN)+ sr1a] += (-Kr0*dCrdC1a);
            ddNdtByVdcTp[sr1a*(this->alignN)+ sr1a] += (-Kr0*dCrdC1a);        
        


        const double CFa = C[sl0a]*C[sl1a];
        const double CRa = C[sr0a]*C[sr1a];

        const double qa = (Kf0*CFa) - (Kr0*CRa);
        dNdtByV[sl0a] = dNdtByV[sl0a] - qa;
        dNdtByV[sl1a] = dNdtByV[sl1a] - qa;
        dNdtByV[sr0a] = dNdtByV[sr0a] + qa;
        dNdtByV[sr1a] = dNdtByV[sr1a] + qa;

        const double dqdTa = (dKfdT0*CFa)-(dKrdT0*CRa);
        ddNdtByVdcTp[sl0a*(this->alignN)+ (this->nSpecies)] -= dqdTa;
        ddNdtByVdcTp[sl1a*(this->alignN)+ (this->nSpecies)] -= dqdTa;
        ddNdtByVdcTp[sr0a*(this->alignN)+ (this->nSpecies)] += dqdTa;
        ddNdtByVdcTp[sr1a*(this->alignN)+ (this->nSpecies)] += dqdTa;

        const double dCfdC0a = C[sl1a];
        const double dCfdC1a = C[sl0a];
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl0a] -= Kf0*dCfdC0a;
        ddNdtByVdcTp[sl0a*(this->alignN)+ sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sl1a*(this->alignN)+ sl1a] -= Kf0*dCfdC1a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl0a] += Kf0*dCfdC0a;
        ddNdtByVdcTp[sr0a*(this->alignN)+ sl1a] += Kf0*dCfdC1a;
        ddNdtByVdcTp[sr1a*(this->alignN)+ sl1a] += Kf0*dCfdC1a;     

        if(i0>=this->Ikf[2] && i0 <this->Ikf[6])
        { 
            const unsigned int k = i0 - this->Ikf[2];

            __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[k]);
            __m256d CF_ = _mm256_set1_pd(CFa);            
            __m256d CR_ = _mm256_set1_pd(CRa);  
            double dKrdC = 0;
            dKrdC = (this->dKfdC_[k]*invKc0);
            __m256d dKrdC_ = _mm256_set1_pd(dKrdC);

            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[k*this->AlignSpecies+j]);

                __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0a*(this->alignN)+j],sl0v);

                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1a*(this->alignN)+j],sl1v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0a*(this->alignN)+j],sr0v);

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1a*(this->alignN)+j],sr1v);
            } 
        }
        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 2;
    }
}   

void  OptReaction::updateJacobian22IrreversibleReaction
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = irreversibleReaction22index.size();

    std::size_t lhsIndex = 0;

    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->irreversibleReaction22index[k];

 
        const unsigned int sl0 = lhsSpeciesIndex1D22IR[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D22IR[lhsIndex+1];
        
        const unsigned int sr0 = rhsSpeciesIndex1D22IR[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D22IR[rhsIndex+1];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];
        //double Kr = 0;
        //double dKrdT = 0;
        //double invKc = 0;


        const double CF = C[sl0]*C[sl1];
        const double CR = C[sr0]*C[sr1];

        const double q = (Kf*CF);
        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;

        const double dqdT = (dKfdT*CF);
        ddNdtByVdcTp[sl0*(this->alignN)+ (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sl1*(this->alignN)+ (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN)+ (this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr1*(this->alignN)+ (this->nSpecies)] += dqdT;

        const double dCfdC0 = C[sl1];
        const double dCfdC1 = C[sl0];
        ddNdtByVdcTp[sl0*(this->alignN)+ sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl1*(this->alignN)+ sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl0*(this->alignN)+ sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sl1*(this->alignN)+ sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sr0*(this->alignN)+ sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr1*(this->alignN)+ sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr0*(this->alignN)+ sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sr1*(this->alignN)+ sl1] += Kf*dCfdC1;     

        if(i>=this->Ikf[2] && i <this->Ikf[6])
        { 
            const unsigned int k = i - this->Ikf[2];

            __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[k]);
            __m256d CF_ = _mm256_set1_pd(CF);            
            __m256d CR_ = _mm256_set1_pd(CR);  
            //double dKrdC = 0;
            //dKrdC = (0);
            __m256d dKrdC_ = _mm256_set1_pd(0);

            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[k*this->AlignSpecies+j]);

                __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0*(this->alignN)+j],sl0v);

                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1*(this->alignN)+j],sl1v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1*(this->alignN)+j],sr1v);
            } 
        }

        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 2;
    }
}   

void  OptReaction::updateJacobian22NonEquilibriumReaction
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = nonEquilibriumReaction22index.size();

    std::size_t lhsIndex = 0;

    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->nonEquilibriumReaction22index[k];

 
        const unsigned int sl0 = lhsSpeciesIndex1D22NER[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D22NER[lhsIndex+1];
        
        const unsigned int sr0 = rhsSpeciesIndex1D22NER[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D22NER[rhsIndex+1];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];
        double Kr = 0;
        double dKrdT = 0;
        double invKc = 0;

        {
            Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
            dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
            invKc = Kr/Kf;
            const double dCrdC0 = C[sr1];
            const double dCrdC1 = C[sr0]; 
            ddNdtByVdcTp[sl0*(this->alignN)+ sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sl1*(this->alignN)+ sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sl0*(this->alignN)+ sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sl1*(this->alignN)+ sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sr0*(this->alignN)+ sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr1*(this->alignN)+ sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr0*(this->alignN)+ sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sr1*(this->alignN)+ sr1] += (-Kr*dCrdC1);
        }

        const double CF = C[sl0]*C[sl1];
        const double CR = C[sr0]*C[sr1];

        const double q = (Kf*CF) - (Kr*CR);
        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;

        const double dqdT = (dKfdT*CF)-(dKrdT*CR);
        ddNdtByVdcTp[sl0*(this->alignN)+ (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sl1*(this->alignN)+ (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN)+ (this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr1*(this->alignN)+ (this->nSpecies)] += dqdT;

        const double dCfdC0 = C[sl1];
        const double dCfdC1 = C[sl0];
        ddNdtByVdcTp[sl0*(this->alignN)+ sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl1*(this->alignN)+ sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl0*(this->alignN)+ sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sl1*(this->alignN)+ sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sr0*(this->alignN)+ sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr1*(this->alignN)+ sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr0*(this->alignN)+ sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sr1*(this->alignN)+ sl1] += Kf*dCfdC1;     

        if(i>=this->Ikf[2] && i <this->Ikf[6])
        { 
            const unsigned int k = i - this->Ikf[2];

            __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[k]);
            __m256d CF_ = _mm256_set1_pd(CF);            
            __m256d CR_ = _mm256_set1_pd(CR);  
            double dKrdC = 0;
            dKrdC = (this->dKfdC_[k]*invKc);
            __m256d dKrdC_ = _mm256_set1_pd(dKrdC);

            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[k*this->AlignSpecies+j]);

                __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0*(this->alignN)+j],sl0v);

                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1*(this->alignN)+j],sl1v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1*(this->alignN)+j],sr1v);
            } 
        }

        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 2;
    }
}   
