#include "OptReaction.H"
#include <immintrin.h>  

/*void  FastChemistry::OptReaction::updateJacobian11
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = reaction11index.size();
    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->reaction11index[k];

        const unsigned int sl0 = lhsSpeciesIndex1D11[lhsIndex+0];
        
        const unsigned int sr0 = rhsSpeciesIndex1D11[rhsIndex+0];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];
        double Kr = 0;
        double dKrdT = 0;
        double invKc = 0;

        if(this->isIrreversible[i]==0)
        {
            const double Kp = ExpNegGbyRT[sr0]/ExpNegGbyRT[sl0];
            const double Kc = Kp > KcLimiter?Kp:KcLimiter;
            const double sumVdBdT = (dBdT[sr0] ) - (dBdT[sl0] );
            const double dKcdTByKc = sumVdBdT;
            invKc = 1.0/Kc;        
            Kr = Kf*invKc;
            dKrdT = (dKfdT*invKc - (Kc > KcLimiter ? Kr*dKcdTByKc : 0));
            const double dCrdC0 = 1;
            ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
        }
        else if(this->isIrreversible[i]==2)
        {
            Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
            dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
            invKc = Kr/Kf;
            const double dCrdC0 = 1;
            ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);        
        }
        const double CF = C[sl0];
        const double CR = C[sr0];

        const double q = (Kf*CF) - (Kr*CR);
        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;

        const double dqdT = (dKfdT*CF)-(dKrdT*CR);
        ddNdtByVdcTp[sl0*(this->alignN)+(this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN)+(this->nSpecies)] += dqdT;

        const double dCfdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;

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

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);            
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);                
            }   
        }
        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 1;
    }
}   */


void  FastChemistry::OptReaction::JF11RR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = reversibleReaction11index.size();
    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->reversibleReaction11index[k];

        const unsigned int sl0 = lhsSpeciesIndex1D11RR[lhsIndex+0];
        
        const unsigned int sr0 = rhsSpeciesIndex1D11RR[rhsIndex+0];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];
        double Kr = 0;
        double dKrdT = 0;
        double invKc = 0;


        const double Kp = ExpNegGbyRT[sr0]/ExpNegGbyRT[sl0];
        const double Kc = Kp > KcLimiter?Kp:KcLimiter;
        const double sumVdBdT = (dBdT[sr0] ) - (dBdT[sl0] );
        const double dKcdTByKc = sumVdBdT;
        invKc = 1.0/Kc;        
        Kr = Kf*invKc;
        dKrdT = (dKfdT*invKc - (Kc > KcLimiter ? Kr*dKcdTByKc : 0));
        const double dCrdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
        

        const double CF = C[sl0];
        const double CR = C[sr0];

        const double q = (Kf*CF) - (Kr*CR);
        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;

        const double dqdT = (dKfdT*CF)-(dKrdT*CR);
        ddNdtByVdcTp[sl0*(this->alignN)+(this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN)+(this->nSpecies)] += dqdT;

        const double dCfdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;

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

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);            
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);                
            }   
        }
        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 1;
    }
}   



void  FastChemistry::OptReaction::JF11IR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = irreversibleReaction11index.size();
    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->irreversibleReaction11index[k];

        const unsigned int sl0 = lhsSpeciesIndex1D11IR[lhsIndex+0];
        
        const unsigned int sr0 = rhsSpeciesIndex1D11IR[rhsIndex+0];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];


        const double CF = C[sl0];
        const double CR = C[sr0];

        const double q = (Kf*CF);
        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;

        const double dqdT = (dKfdT*CF);
        ddNdtByVdcTp[sl0*(this->alignN)+(this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN)+(this->nSpecies)] += dqdT;

        const double dCfdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;

        if(i>=this->Ikf[2] && i <this->Ikf[6])
        { 
            const unsigned int k = i - this->Ikf[2];

            __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[k]);
            __m256d CF_ = _mm256_set1_pd(CF);            
            __m256d CR_ = _mm256_set1_pd(CR);

            __m256d dKrdC_ = _mm256_set1_pd(0);
            for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
            {
                __m256d dMdC = load256d(&this->ThirdBodyFactor1D[k*this->AlignSpecies+j]);
                __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0*(this->alignN)+j],sl0v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);            
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);                
            }   
        }
        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 1;
    }
}   


void  FastChemistry::OptReaction::JF11NER
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = nonEquilibriumReaction11index.size();
    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->nonEquilibriumReaction11index[k];

        const unsigned int sl0 = lhsSpeciesIndex1D11NER[lhsIndex+0];
        
        const unsigned int sr0 = rhsSpeciesIndex1D11NER[rhsIndex+0];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];
        double Kr = 0;
        double dKrdT = 0;
        double invKc = 0;



        
        Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
        dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
        invKc = Kr/Kf;
        const double dCrdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);        
        
        const double CF = C[sl0];
        const double CR = C[sr0];

        const double q = (Kf*CF) - (Kr*CR);
        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;

        const double dqdT = (dKfdT*CF)-(dKrdT*CR);
        ddNdtByVdcTp[sl0*(this->alignN)+(this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN)+(this->nSpecies)] += dqdT;

        const double dCfdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;

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

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);            
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);                
            }   
        }
        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 1;
    }
}   