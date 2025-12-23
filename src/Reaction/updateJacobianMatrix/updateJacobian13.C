#include "OptReaction.H"
#include <immintrin.h>  

/*void  FastChemistry::OptReaction::updateJacobian13
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = reaction13index.size();

    std::size_t lhsIndex = 0;

    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->reaction13index[k];

        const unsigned int sl0 = lhsSpeciesIndex1D13[lhsIndex+0];
        
        const unsigned int sr0 = rhsSpeciesIndex1D13[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D13[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D13[rhsIndex+2];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];
        double Kr = 0;
        double dKrdT = 0;
        double invKc = 0;

        if(this->isIrreversible[i]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]);
            double Kc = Kp*this->Pow_pByRT_SumVki[4];  
            Kc = Kc > KcLimiter?Kc:KcLimiter;
            invKc = 1.0/Kc;
            const double sumVdBdT = (dBdT[sr0] + dBdT[sr1] + dBdT[sr2]) - (dBdT[sl0] );
            const double dKcdTByKc = sumVdBdT - 1*invT;
            Kr = Kf*invKc;
            dKrdT = (dKfdT*invKc - (Kc > KcLimiter ? Kr*dKcdTByKc : 0));

            const double dCrdC0 = C[sr1]*C[sr2];
            const double dCrdC1 = C[sr0]*C[sr2];
            const double dCrdC2 = C[sr0]*C[sr1];
            ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
            ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr2*(this->alignN)+sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sr2*(this->alignN)+sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sl0*(this->alignN)+sr2] -= (-Kr*dCrdC2);
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
            ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sr2*(this->alignN)+sr0] += (-Kr*dCrdC0);
            ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
            ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sr2*(this->alignN)+sr1] += (-Kr*dCrdC1);
            ddNdtByVdcTp[sl0*(this->alignN)+sr2] -= (-Kr*dCrdC2);
            ddNdtByVdcTp[sr0*(this->alignN)+sr2] += (-Kr*dCrdC2);
            ddNdtByVdcTp[sr1*(this->alignN)+sr2] += (-Kr*dCrdC2);
            ddNdtByVdcTp[sr2*(this->alignN)+sr2] += (-Kr*dCrdC2);
        }

        const double CF = C[sl0];
        const double CR = C[sr0]*C[sr1]*C[sr2];

        const double q = (Kf*CF) - (Kr*CR);
        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        const double dqdT = (dKfdT*CF)-(dKrdT*CR);
        ddNdtByVdcTp[sl0*(this->alignN)+(this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN)+(this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr1*(this->alignN)+(this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr2*(this->alignN)+(this->nSpecies)] += dqdT;

        const double dCfdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr1*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr2*(this->alignN)+sl0] += Kf*dCfdC0;

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

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1*(this->alignN)+j],sr1v);

                __m256d sr2v = load256d(&ddNdtByVdcTp[sr2*(this->alignN)+j]);
                sr2v = _mm256_add_pd(sr2v,WdMdC);
                store256d(&ddNdtByVdcTp[sr2*(this->alignN)+j],sr2v);            
            } 
        }
        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 3;
    }
}   */


void  FastChemistry::OptReaction::JF13RR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = reversibleReaction13index.size();

    std::size_t lhsIndex = 0;

    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->reversibleReaction13index[k];

        const unsigned int sl0 = lhsSpeciesIndex1D13RR[lhsIndex+0];
        
        const unsigned int sr0 = rhsSpeciesIndex1D13RR[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D13RR[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D13RR[rhsIndex+2];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];
        double Kr = 0;
        double dKrdT = 0;


        /*const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]);
        double Kc = Kp*this->Pow_pByRT_SumVki[4];  
        Kc = Kc > KcLimiter?Kc:KcLimiter;
        invKc = 1.0/Kc;
        const double sumVdBdT = (dBdT[sr0] + dBdT[sr1] + dBdT[sr2]) - (dBdT[sl0] );
        const double dKcdTByKc = sumVdBdT - 1*invT;
        Kr = Kf*invKc;
        dKrdT = (dKfdT*invKc - (Kc > KcLimiter ? Kr*dKcdTByKc : 0));*/

        const double invKp = (invNegGstdByRT[sr0]*invNegGstdByRT[sr1])*(invNegGstdByRT[sr2]*ExpNegGbyRT[sl0]);
        double invKc = invKp*this->Pow_pByRT_SumVki[0];  
        invKc = std::min(invKc,invKcLimiter);
        const double sumVdBdT = (dBdT[sr0] + dBdT[sr1] + dBdT[sr2]) - (dBdT[sl0] );
        const double dKcdTByKc = sumVdBdT - 1*invT;
        Kr = Kf*invKc;
        dKrdT = (dKfdT*invKc - (invKc < invKcLimiter ? Kr*dKcdTByKc : 0));

        const double dCrdC0 = C[sr1]*C[sr2];
        const double dCrdC1 = C[sr0]*C[sr2];
        const double dCrdC2 = C[sr0]*C[sr1];
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr2*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr2*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sl0*(this->alignN)+sr2] -= (-Kr*dCrdC2);
        ddNdtByVdcTp[sr0*(this->alignN)+sr2] += (-Kr*dCrdC2);
        ddNdtByVdcTp[sr1*(this->alignN)+sr2] += (-Kr*dCrdC2);
        ddNdtByVdcTp[sr2*(this->alignN)+sr2] += (-Kr*dCrdC2);
        


        const double CF = C[sl0];
        const double CR = C[sr0]*C[sr1]*C[sr2];

        const double q = (Kf*CF) - (Kr*CR);
        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        const double dqdT = (dKfdT*CF)-(dKrdT*CR);
        ddNdtByVdcTp[sl0*(this->alignN)+(this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN)+(this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr1*(this->alignN)+(this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr2*(this->alignN)+(this->nSpecies)] += dqdT;

        const double dCfdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr1*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr2*(this->alignN)+sl0] += Kf*dCfdC0;

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

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1*(this->alignN)+j],sr1v);

                __m256d sr2v = load256d(&ddNdtByVdcTp[sr2*(this->alignN)+j]);
                sr2v = _mm256_add_pd(sr2v,WdMdC);
                store256d(&ddNdtByVdcTp[sr2*(this->alignN)+j],sr2v);            
            } 
        }
        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 3;
    }
}   


void  FastChemistry::OptReaction::JF13IR
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = irreversibleReaction13index.size();

    std::size_t lhsIndex = 0;

    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->irreversibleReaction13index[k];

        const unsigned int sl0 = lhsSpeciesIndex1D13IR[lhsIndex+0];
        
        const unsigned int sr0 = rhsSpeciesIndex1D13IR[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D13IR[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D13IR[rhsIndex+2];

        double Kf = this->Kf_[i];
        double dKfdT = this->dKfdT_[i];

        const double CF = C[sl0];
        const double CR = C[sr0]*C[sr1]*C[sr2];

        const double q = (Kf*CF);
        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        const double dqdT = (dKfdT*CF);
        ddNdtByVdcTp[sl0*(this->alignN)+(this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN)+(this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr1*(this->alignN)+(this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr2*(this->alignN)+(this->nSpecies)] += dqdT;

        const double dCfdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr1*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr2*(this->alignN)+sl0] += Kf*dCfdC0;

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

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1*(this->alignN)+j],sr1v);

                __m256d sr2v = load256d(&ddNdtByVdcTp[sr2*(this->alignN)+j]);
                sr2v = _mm256_add_pd(sr2v,WdMdC);
                store256d(&ddNdtByVdcTp[sr2*(this->alignN)+j],sr2v);            
            } 
        }
        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 3;
    }
}   



void  FastChemistry::OptReaction::JF13NER
(
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    std::size_t end = nonEquilibriumReaction13index.size();

    std::size_t lhsIndex = 0;

    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end; k++)
    {
        const unsigned int i = this->nonEquilibriumReaction13index[k];

        const unsigned int sl0 = lhsSpeciesIndex1D13NER[lhsIndex+0];
        
        const unsigned int sr0 = rhsSpeciesIndex1D13NER[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D13NER[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D13NER[rhsIndex+2];

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
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr2*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr2*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sl0*(this->alignN)+sr2] -= (-Kr*dCrdC2);
        ddNdtByVdcTp[sr0*(this->alignN)+sr2] += (-Kr*dCrdC2);
        ddNdtByVdcTp[sr1*(this->alignN)+sr2] += (-Kr*dCrdC2);
        ddNdtByVdcTp[sr2*(this->alignN)+sr2] += (-Kr*dCrdC2);
        

        const double CF = C[sl0];
        const double CR = C[sr0]*C[sr1]*C[sr2];

        const double q = (Kf*CF) - (Kr*CR);
        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        const double dqdT = (dKfdT*CF)-(dKrdT*CR);
        ddNdtByVdcTp[sl0*(this->alignN)+(this->nSpecies)] -= dqdT;
        ddNdtByVdcTp[sr0*(this->alignN)+(this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr1*(this->alignN)+(this->nSpecies)] += dqdT;
        ddNdtByVdcTp[sr2*(this->alignN)+(this->nSpecies)] += dqdT;

        const double dCfdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
        ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr1*(this->alignN)+sl0] += Kf*dCfdC0;
        ddNdtByVdcTp[sr2*(this->alignN)+sl0] += Kf*dCfdC0;

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

                __m256d sr1v = load256d(&ddNdtByVdcTp[sr1*(this->alignN)+j]);
                sr1v = _mm256_add_pd(sr1v,WdMdC);
                store256d(&ddNdtByVdcTp[sr1*(this->alignN)+j],sr1v);

                __m256d sr2v = load256d(&ddNdtByVdcTp[sr2*(this->alignN)+j]);
                sr2v = _mm256_add_pd(sr2v,WdMdC);
                store256d(&ddNdtByVdcTp[sr2*(this->alignN)+j],sr2v);            
            } 
        }
        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 3;
    }
}   