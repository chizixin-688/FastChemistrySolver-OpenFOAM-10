#include "OptReaction.H"
#include <immintrin.h>  

/*void  FastChemistry::OptReaction::updateJacobian33
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = reaction33index.size();

    std::size_t lhsIndex = 0;

    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->reaction33index[k];

 
        const unsigned int sl0 = lhsSpeciesIndex1D33[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D33[lhsIndex+1];
        const unsigned int sl2 = lhsSpeciesIndex1D33[lhsIndex+2];

        const unsigned int sr0 = rhsSpeciesIndex1D33[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D33[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D33[rhsIndex+2];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];
        double Kr = 0;
        double dKrdT = 0;
        double invKc = 0;

        if(this->isIrreversible[i]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
            double Kc = Kp;
            Kc = Kc > KcLimiter?Kc:KcLimiter;  
            invKc = 1.0/Kc;
            Kr = Kf*invKc;
            const double sumVdBdT = (dBdT[sr0] + dBdT[sr1] + dBdT[sr2]) - (dBdT[sl0] + dBdT[sl1] + dBdT[sl2]);
            const double dKcdTByKc = sumVdBdT;    
            dKrdT = (dKfdT*invKc - (Kc > KcLimiter ? Kr*dKcdTByKc : 0));  
            
            const double dCrdC0 = C[sr1]*C[sr2];
            const double dCrdC1 = C[sr0]*C[sr2];
            const double dCrdC2 = C[sr0]*C[sr1];  
            ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sl1*(this->alignN)+sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sl2*(this->alignN)+sr0] -= (-Kr*dCrdC0);    
            ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr2*(this->alignN)+sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sl1*(this->alignN)+sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sl2*(this->alignN)+sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sr2*(this->alignN)+sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sl0*(this->alignN)+sr2] -= (-Kr*dCrdC2);
            ddNdtByVdcTp[sl1*(this->alignN)+sr2] -= (-Kr*dCrdC2);
            ddNdtByVdcTp[sl2*(this->alignN)+sr2] -= (-Kr*dCrdC2);
            ddNdtByVdcTp[sr0*(this->alignN)+sr2] += (-Kr*dCrdC2);
            ddNdtByVdcTp[sr1*(this->alignN)+sr2] += (-Kr*dCrdC2);
            ddNdtByVdcTp[sr2*(this->alignN)+sr2] += (-Kr*dCrdC2); 
        }
        else if(this->isIrreversible[i]==2)
        {
            Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
            dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
            invKc = Kr/Kf;
            const double dCrdC0 = C[sr1]*C[sr2];
            const double dCrdC1 = C[sr0]*C[sr2];
            const double dCrdC2 = C[sr0]*C[sr1];  

            ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sl1*(this->alignN)+sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sl2*(this->alignN)+sr0] -= (-Kr*dCrdC0);    
            ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr2*(this->alignN)+sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sl1*(this->alignN)+sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sl2*(this->alignN)+sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sr2*(this->alignN)+sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sl0*(this->alignN)+sr2] -= (-Kr*dCrdC2);
            ddNdtByVdcTp[sl1*(this->alignN)+sr2] -= (-Kr*dCrdC2);
            ddNdtByVdcTp[sl2*(this->alignN)+sr2] -= (-Kr*dCrdC2);
            ddNdtByVdcTp[sr0*(this->alignN)+sr2] += (-Kr*dCrdC2);
            ddNdtByVdcTp[sr1*(this->alignN)+sr2] += (-Kr*dCrdC2);
            ddNdtByVdcTp[sr2*(this->alignN)+sr2] += (-Kr*dCrdC2); 
        }

        const double CF = C[sl0]*C[sl1]*C[sl2];
        const double CR = C[sr0]*C[sr1]*C[sr2];
        
        const double q = (Kf*CF) - (Kr*CR);
        dNdtByV[sl0] = dNdtByV[sl0] - q; 
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        const double dqdT = (dKfdT*CF)-(dKrdT*CR);
        ddNdtByVdcTp[sl0*(this->alignN) + (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sl1*(this->alignN) + (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sl2*(this->alignN) + (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN) + (this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr1*(this->alignN) + (this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr2*(this->alignN) + (this->nSpecies)] += dqdT;
        
        const double dCfdC0 = C[sl1]*C[sl2];
        const double dCfdC1 = C[sl0]*C[sl2];
        const double dCfdC2 = C[sl0]*C[sl1];
        ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl1*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl2*(this->alignN)+sl0] -= Kf*dCfdC0;    
        ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr1*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr2*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sl0*(this->alignN)+sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sl1*(this->alignN)+sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sl2*(this->alignN)+sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sr0*(this->alignN)+sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sr1*(this->alignN)+sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sr2*(this->alignN)+sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sl0*(this->alignN)+sl2] -= Kf*dCfdC2;
        ddNdtByVdcTp[sl1*(this->alignN)+sl2] -= Kf*dCfdC2;
        ddNdtByVdcTp[sl2*(this->alignN)+sl2] -= Kf*dCfdC2;
        ddNdtByVdcTp[sr0*(this->alignN)+sl2] += Kf*dCfdC2;
        ddNdtByVdcTp[sr1*(this->alignN)+sl2] += Kf*dCfdC2;
        ddNdtByVdcTp[sr2*(this->alignN)+sl2] += Kf*dCfdC2;

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
                __m256d WdMdC = _mm256_mul_pd(dKrdC_,CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0*(this->alignN)+j],sl0v);

                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1*(this->alignN)+j],sl1v);

                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2*(this->alignN)+j],sl2v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1*(this->alignN)+j],sr1v);

                __m256d sr2v = load256d(&ddNdtByVdcTp[sr2*(this->alignN)+j]);
                sr2v = _mm256_add_pd(sr2v,WdMdC);
                store256d(&ddNdtByVdcTp[sr2*(this->alignN)+j],sr2v);

            } 
        }  

        lhsIndex = lhsIndex + 3;
        rhsIndex = rhsIndex + 3;
    }
}   */

void  FastChemistry::OptReaction::JF33RR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = this->reversibleReaction33index.size();

    std::size_t lhsIndex = 0;

    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->reversibleReaction33index[k];

 
        const unsigned int sl0 = lhsSpeciesIndex1D33RR[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D33RR[lhsIndex+1];
        const unsigned int sl2 = lhsSpeciesIndex1D33RR[lhsIndex+2];

        const unsigned int sr0 = rhsSpeciesIndex1D33RR[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D33RR[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D33RR[rhsIndex+2];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];
        double Kr = 0;
        double dKrdT = 0;


        /*const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
        double Kc = Kp;
        Kc = Kc > KcLimiter?Kc:KcLimiter;  
        invKc = 1.0/Kc;
        Kr = Kf*invKc;
        const double sumVdBdT = (dBdT[sr0] + dBdT[sr1] + dBdT[sr2]) - (dBdT[sl0] + dBdT[sl1] + dBdT[sl2]);
        const double dKcdTByKc = sumVdBdT;    
        dKrdT = (dKfdT*invKc - (Kc > KcLimiter ? Kr*dKcdTByKc : 0));*/

        const double invKp = (invNegGstdByRT[sr0]*invNegGstdByRT[sr1])*
            (invNegGstdByRT[sr2]*ExpNegGbyRT[sl0])*
            (ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
        double invKc = invKp;

        invKc = std::min(invKc,invKcLimiter);

        Kr = Kf*invKc;
        const double sumVdBdT = (dBdT[sr0] + dBdT[sr1] + dBdT[sr2]) - (dBdT[sl0] + dBdT[sl1] + dBdT[sl2]);
        const double dKcdTByKc = sumVdBdT;    
        dKrdT = (dKfdT*invKc - (invKc < invKcLimiter ? Kr*dKcdTByKc : 0));
            
        const double dCrdC0 = C[sr1]*C[sr2];
        const double dCrdC1 = C[sr0]*C[sr2];
        const double dCrdC2 = C[sr0]*C[sr1];  
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl1*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl2*(this->alignN)+sr0] -= (-Kr*dCrdC0);    
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr2*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sl1*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sl2*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr2*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sl0*(this->alignN)+sr2] -= (-Kr*dCrdC2);
        ddNdtByVdcTp[sl1*(this->alignN)+sr2] -= (-Kr*dCrdC2);
        ddNdtByVdcTp[sl2*(this->alignN)+sr2] -= (-Kr*dCrdC2);
        ddNdtByVdcTp[sr0*(this->alignN)+sr2] += (-Kr*dCrdC2);
        ddNdtByVdcTp[sr1*(this->alignN)+sr2] += (-Kr*dCrdC2);
        ddNdtByVdcTp[sr2*(this->alignN)+sr2] += (-Kr*dCrdC2); 
        
        const double CF = C[sl0]*C[sl1]*C[sl2];
        const double CR = C[sr0]*C[sr1]*C[sr2];
        
        const double q = (Kf*CF) - (Kr*CR);
        dNdtByV[sl0] = dNdtByV[sl0] - q; 
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        const double dqdT = (dKfdT*CF)-(dKrdT*CR);
        ddNdtByVdcTp[sl0*(this->alignN) + (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sl1*(this->alignN) + (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sl2*(this->alignN) + (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN) + (this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr1*(this->alignN) + (this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr2*(this->alignN) + (this->nSpecies)] += dqdT;
        
        const double dCfdC0 = C[sl1]*C[sl2];
        const double dCfdC1 = C[sl0]*C[sl2];
        const double dCfdC2 = C[sl0]*C[sl1];
        ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl1*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl2*(this->alignN)+sl0] -= Kf*dCfdC0;    
        ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr1*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr2*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sl0*(this->alignN)+sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sl1*(this->alignN)+sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sl2*(this->alignN)+sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sr0*(this->alignN)+sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sr1*(this->alignN)+sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sr2*(this->alignN)+sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sl0*(this->alignN)+sl2] -= Kf*dCfdC2;
        ddNdtByVdcTp[sl1*(this->alignN)+sl2] -= Kf*dCfdC2;
        ddNdtByVdcTp[sl2*(this->alignN)+sl2] -= Kf*dCfdC2;
        ddNdtByVdcTp[sr0*(this->alignN)+sl2] += Kf*dCfdC2;
        ddNdtByVdcTp[sr1*(this->alignN)+sl2] += Kf*dCfdC2;
        ddNdtByVdcTp[sr2*(this->alignN)+sl2] += Kf*dCfdC2;

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
                __m256d WdMdC = _mm256_mul_pd(dKrdC_,CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0*(this->alignN)+j],sl0v);

                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1*(this->alignN)+j],sl1v);

                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2*(this->alignN)+j],sl2v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1*(this->alignN)+j],sr1v);

                __m256d sr2v = load256d(&ddNdtByVdcTp[sr2*(this->alignN)+j]);
                sr2v = _mm256_add_pd(sr2v,WdMdC);
                store256d(&ddNdtByVdcTp[sr2*(this->alignN)+j],sr2v);

            } 
        }  

        lhsIndex = lhsIndex + 3;
        rhsIndex = rhsIndex + 3;
    }
}   

void  FastChemistry::OptReaction::JF33IR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = this->irreversibleReaction33index.size();

    std::size_t lhsIndex = 0;

    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->irreversibleReaction33index[k];

 
        const unsigned int sl0 = lhsSpeciesIndex1D33IR[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D33IR[lhsIndex+1];
        const unsigned int sl2 = lhsSpeciesIndex1D33IR[lhsIndex+2];

        const unsigned int sr0 = rhsSpeciesIndex1D33IR[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D33IR[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D33IR[rhsIndex+2];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];

        const double CF = C[sl0]*C[sl1]*C[sl2];
        const double CR = C[sr0]*C[sr1]*C[sr2];
        
        const double q = (Kf*CF);
        dNdtByV[sl0] = dNdtByV[sl0] - q; 
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        const double dqdT = (dKfdT*CF);
        ddNdtByVdcTp[sl0*(this->alignN) + (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sl1*(this->alignN) + (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sl2*(this->alignN) + (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN) + (this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr1*(this->alignN) + (this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr2*(this->alignN) + (this->nSpecies)] += dqdT;
        
        const double dCfdC0 = C[sl1]*C[sl2];
        const double dCfdC1 = C[sl0]*C[sl2];
        const double dCfdC2 = C[sl0]*C[sl1];
        ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl1*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl2*(this->alignN)+sl0] -= Kf*dCfdC0;    
        ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr1*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr2*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sl0*(this->alignN)+sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sl1*(this->alignN)+sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sl2*(this->alignN)+sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sr0*(this->alignN)+sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sr1*(this->alignN)+sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sr2*(this->alignN)+sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sl0*(this->alignN)+sl2] -= Kf*dCfdC2;
        ddNdtByVdcTp[sl1*(this->alignN)+sl2] -= Kf*dCfdC2;
        ddNdtByVdcTp[sl2*(this->alignN)+sl2] -= Kf*dCfdC2;
        ddNdtByVdcTp[sr0*(this->alignN)+sl2] += Kf*dCfdC2;
        ddNdtByVdcTp[sr1*(this->alignN)+sl2] += Kf*dCfdC2;
        ddNdtByVdcTp[sr2*(this->alignN)+sl2] += Kf*dCfdC2;

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
                __m256d WdMdC = _mm256_mul_pd(dKrdC_,CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0*(this->alignN)+j],sl0v);

                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1*(this->alignN)+j],sl1v);

                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2*(this->alignN)+j],sl2v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1*(this->alignN)+j],sr1v);

                __m256d sr2v = load256d(&ddNdtByVdcTp[sr2*(this->alignN)+j]);
                sr2v = _mm256_add_pd(sr2v,WdMdC);
                store256d(&ddNdtByVdcTp[sr2*(this->alignN)+j],sr2v);

            } 
        }  

        lhsIndex = lhsIndex + 3;
        rhsIndex = rhsIndex + 3;
    }
}   


void  FastChemistry::OptReaction::JF33NER
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = this->nonEquilibriumReaction33index.size();

    std::size_t lhsIndex = 0;

    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->nonEquilibriumReaction33index[k];

 
        const unsigned int sl0 = lhsSpeciesIndex1D33NER[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D33NER[lhsIndex+1];
        const unsigned int sl2 = lhsSpeciesIndex1D33NER[lhsIndex+2];

        const unsigned int sr0 = rhsSpeciesIndex1D33NER[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D33NER[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D33NER[rhsIndex+2];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];
        double Kr = 0;
        double dKrdT = 0;
        double invKc = 0;


        
        Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
        dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
        invKc = Kr/Kf;
        const double dCrdC0 = C[sr1]*C[sr2];
        const double dCrdC1 = C[sr0]*C[sr2];
        const double dCrdC2 = C[sr0]*C[sr1];  

        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl1*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl2*(this->alignN)+sr0] -= (-Kr*dCrdC0);    
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr2*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sl1*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sl2*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr2*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sl0*(this->alignN)+sr2] -= (-Kr*dCrdC2);
        ddNdtByVdcTp[sl1*(this->alignN)+sr2] -= (-Kr*dCrdC2);
        ddNdtByVdcTp[sl2*(this->alignN)+sr2] -= (-Kr*dCrdC2);
        ddNdtByVdcTp[sr0*(this->alignN)+sr2] += (-Kr*dCrdC2);
        ddNdtByVdcTp[sr1*(this->alignN)+sr2] += (-Kr*dCrdC2);
        ddNdtByVdcTp[sr2*(this->alignN)+sr2] += (-Kr*dCrdC2); 
        

        const double CF = C[sl0]*C[sl1]*C[sl2];
        const double CR = C[sr0]*C[sr1]*C[sr2];
        
        const double q = (Kf*CF) - (Kr*CR);
        dNdtByV[sl0] = dNdtByV[sl0] - q; 
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        const double dqdT = (dKfdT*CF)-(dKrdT*CR);
        ddNdtByVdcTp[sl0*(this->alignN) + (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sl1*(this->alignN) + (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sl2*(this->alignN) + (this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN) + (this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr1*(this->alignN) + (this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr2*(this->alignN) + (this->nSpecies)] += dqdT;
        
        const double dCfdC0 = C[sl1]*C[sl2];
        const double dCfdC1 = C[sl0]*C[sl2];
        const double dCfdC2 = C[sl0]*C[sl1];
        ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl1*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sl2*(this->alignN)+sl0] -= Kf*dCfdC0;    
        ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr1*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr2*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sl0*(this->alignN)+sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sl1*(this->alignN)+sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sl2*(this->alignN)+sl1] -= Kf*dCfdC1;
        ddNdtByVdcTp[sr0*(this->alignN)+sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sr1*(this->alignN)+sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sr2*(this->alignN)+sl1] += Kf*dCfdC1;
        ddNdtByVdcTp[sl0*(this->alignN)+sl2] -= Kf*dCfdC2;
        ddNdtByVdcTp[sl1*(this->alignN)+sl2] -= Kf*dCfdC2;
        ddNdtByVdcTp[sl2*(this->alignN)+sl2] -= Kf*dCfdC2;
        ddNdtByVdcTp[sr0*(this->alignN)+sl2] += Kf*dCfdC2;
        ddNdtByVdcTp[sr1*(this->alignN)+sl2] += Kf*dCfdC2;
        ddNdtByVdcTp[sr2*(this->alignN)+sl2] += Kf*dCfdC2;

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
                __m256d WdMdC = _mm256_mul_pd(dKrdC_,CR_);
                WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

                __m256d sl0v = load256d(&ddNdtByVdcTp[sl0*(this->alignN)+j]);
                sl0v = _mm256_sub_pd(sl0v,WdMdC);
                store256d(&ddNdtByVdcTp[sl0*(this->alignN)+j],sl0v);

                __m256d sl1v = load256d(&ddNdtByVdcTp[sl1*(this->alignN)+j]);
                sl1v = _mm256_sub_pd(sl1v,WdMdC);
                store256d(&ddNdtByVdcTp[sl1*(this->alignN)+j],sl1v);

                __m256d sl2v = load256d(&ddNdtByVdcTp[sl2*(this->alignN)+j]);
                sl2v = _mm256_sub_pd(sl2v,WdMdC);
                store256d(&ddNdtByVdcTp[sl2*(this->alignN)+j],sl2v);

                __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);
                sr0v = _mm256_add_pd(sr0v,WdMdC);
                store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1*(this->alignN)+j],sr1v);

                __m256d sr2v = load256d(&ddNdtByVdcTp[sr2*(this->alignN)+j]);
                sr2v = _mm256_add_pd(sr2v,WdMdC);
                store256d(&ddNdtByVdcTp[sr2*(this->alignN)+j],sr2v);

            } 
        }  

        lhsIndex = lhsIndex + 3;
        rhsIndex = rhsIndex + 3;
    }
}   