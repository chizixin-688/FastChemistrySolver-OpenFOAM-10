#include "OptReaction.H"

/*void 
FastChemistry::OptReaction::update13Reaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reaction13index.size();
    //std::size_t remain = end%4;

    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    for(std::size_t i=0; i<end; i++)
    {
        const unsigned int KfIndex = this->reaction13index[i];

        const unsigned int sl0 = lhsSpeciesIndex1D13[lhsIndex+0];

        const unsigned int sr0 = rhsSpeciesIndex1D13[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D13[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D13[rhsIndex+2];


        double Kr = 0;
        double Kf = Kf_[KfIndex];
        if(this->isIrreversible[KfIndex]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]);
            double Kc = Kp*this->Pow_pByRT_SumVki[4];
            Kc = std::max(Kc,KcLimiter);
            Kr = Kf/Kc;         
        }
        else if(this->isIrreversible[KfIndex]==2)
        {
            unsigned int l = KfIndex - this->Ikf[1] + this->Ikf[9];
            Kr = this->Kf_[l];
        }

        const double CF = c[sl0];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 3;
    }
}*/

void 
FastChemistry::OptReaction::RF13RR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{

    std::size_t end = reversibleReaction13index.size();


    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;

    std::size_t remain = end%4;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->reversibleReaction13index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D13RR[lhsIndex+0];

        const unsigned int sr0a = rhsSpeciesIndex1D13RR[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D13RR[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D13RR[rhsIndex+2];

        double Kf0 = Kf_[i0];
        const double Kp0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a]*ExpNegGbyRT[sr2a])*(invNegGstdByRT[sl0a]);
        double Kc0 = Kp0*this->Pow_pByRT_SumVki[4];
        Kc0 = std::max(Kc0,KcLimiter);
        const double Kr0 = Kf0/Kc0;         
        
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;


        const unsigned int i1 = this->reversibleReaction13index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D13RR[lhsIndex+1];

        const unsigned int sr0b = rhsSpeciesIndex1D13RR[rhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D13RR[rhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D13RR[rhsIndex+5];

        double Kf1 = Kf_[i1];
        const double Kp1 = (ExpNegGbyRT[sr0b]*ExpNegGbyRT[sr1b]*ExpNegGbyRT[sr2b])*(invNegGstdByRT[sl0b]);
        double Kc1 = Kp1*this->Pow_pByRT_SumVki[4];
        Kc1 = std::max(Kc1,KcLimiter);
        const double Kr1 = Kf1/Kc1;         
        
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;


        const unsigned int i2 = this->reversibleReaction13index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D13RR[lhsIndex+2];

        const unsigned int sr0c = rhsSpeciesIndex1D13RR[rhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D13RR[rhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D13RR[rhsIndex+8];

        double Kf2 = Kf_[i2];
        const double Kp2 = (ExpNegGbyRT[sr0c]*ExpNegGbyRT[sr1c]*ExpNegGbyRT[sr2c])*(invNegGstdByRT[sl0c]);
        double Kc2 = Kp2*this->Pow_pByRT_SumVki[4];
        Kc2 = std::max(Kc2,KcLimiter);
        const double Kr2 = Kf2/Kc2;         
        
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;


        const unsigned int i3 = this->reversibleReaction13index[k+3];

        const unsigned int sl0d = lhsSpeciesIndex1D13RR[lhsIndex+3];

        const unsigned int sr0d = rhsSpeciesIndex1D13RR[rhsIndex+9];
        const unsigned int sr1d = rhsSpeciesIndex1D13RR[rhsIndex+10];
        const unsigned int sr2d = rhsSpeciesIndex1D13RR[rhsIndex+11];

        double Kf3 = Kf_[i3];
        const double Kp3 = (ExpNegGbyRT[sr0d]*ExpNegGbyRT[sr1d]*ExpNegGbyRT[sr2d])*(invNegGstdByRT[sl0d]);
        double Kc3 = Kp3*this->Pow_pByRT_SumVki[4];
        Kc3 = std::max(Kc3,KcLimiter);
        const double Kr3 = Kf3/Kc3;         
        
        const double CF3 = c[sl0d];
        const double CR3 = c[sr0d]*c[sr1d]*c[sr2d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;

        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
        dNdtByV[sr2d] = dNdtByV[sr2d] + q3;



        lhsIndex = lhsIndex + 4;
        rhsIndex = rhsIndex + 12;
    }
    if(remain==1)
    {
        std::size_t k = end-1;
        const unsigned int i0 = this->reversibleReaction13index[k];

        const unsigned int sl0a = lhsSpeciesIndex1D13RR[lhsIndex+0];

        const unsigned int sr0a = rhsSpeciesIndex1D13RR[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D13RR[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D13RR[rhsIndex+2];



        double Kf0 = Kf_[i0];
        const double Kp0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a]*ExpNegGbyRT[sr2a])*(invNegGstdByRT[sl0a]);
        double Kc0 = Kp0*this->Pow_pByRT_SumVki[4];
        Kc0 = std::max(Kc0,KcLimiter);
        const double Kr0 = Kf0/Kc0;         
        
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 3;
    }
    else if(remain==2)
    {
        std::size_t k = end-2;
        const unsigned int i0 = this->reversibleReaction13index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D13RR[lhsIndex+0];

        const unsigned int sr0a = rhsSpeciesIndex1D13RR[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D13RR[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D13RR[rhsIndex+2];

        double Kf0 = Kf_[i0];
        const double Kp0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a]*ExpNegGbyRT[sr2a])*(invNegGstdByRT[sl0a]);
        double Kc0 = Kp0*this->Pow_pByRT_SumVki[4];
        Kc0 = std::max(Kc0,KcLimiter);
        const double Kr0 = Kf0/Kc0;         
        
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;


        const unsigned int i1 = this->reversibleReaction13index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D13RR[lhsIndex+1];

        const unsigned int sr0b = rhsSpeciesIndex1D13RR[rhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D13RR[rhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D13RR[rhsIndex+5];

        double Kf1 = Kf_[i1];
        const double Kp1 = (ExpNegGbyRT[sr0b]*ExpNegGbyRT[sr1b]*ExpNegGbyRT[sr2b])*(invNegGstdByRT[sl0b]);
        double Kc1 = Kp1*this->Pow_pByRT_SumVki[4];
        Kc1 = std::max(Kc1,KcLimiter);
        const double Kr1 = Kf1/Kc1;         
        
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;


        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 6;
    }
    else if(remain==3)
    {
        std::size_t k = end-3;
        const unsigned int i0 = this->reversibleReaction13index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D13RR[lhsIndex+0];

        const unsigned int sr0a = rhsSpeciesIndex1D13RR[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D13RR[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D13RR[rhsIndex+2];

        double Kf0 = Kf_[i0];
        const double Kp0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a]*ExpNegGbyRT[sr2a])*(invNegGstdByRT[sl0a]);
        double Kc0 = Kp0*this->Pow_pByRT_SumVki[4];
        Kc0 = std::max(Kc0,KcLimiter);
        const double Kr0 = Kf0/Kc0;         
        
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;


        const unsigned int i1 = this->reversibleReaction13index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D13RR[lhsIndex+1];

        const unsigned int sr0b = rhsSpeciesIndex1D13RR[rhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D13RR[rhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D13RR[rhsIndex+5];

        double Kf1 = Kf_[i1];
        const double Kp1 = (ExpNegGbyRT[sr0b]*ExpNegGbyRT[sr1b]*ExpNegGbyRT[sr2b])*(invNegGstdByRT[sl0b]);
        double Kc1 = Kp1*this->Pow_pByRT_SumVki[4];
        Kc1 = std::max(Kc1,KcLimiter);
        const double Kr1 = Kf1/Kc1;         
        
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;


        const unsigned int i2 = this->reversibleReaction13index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D13RR[lhsIndex+2];

        const unsigned int sr0c = rhsSpeciesIndex1D13RR[rhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D13RR[rhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D13RR[rhsIndex+8];

        double Kf2 = Kf_[i2];
        const double Kp2 = (ExpNegGbyRT[sr0c]*ExpNegGbyRT[sr1c]*ExpNegGbyRT[sr2c])*(invNegGstdByRT[sl0c]);
        double Kc2 = Kp2*this->Pow_pByRT_SumVki[4];
        Kc2 = std::max(Kc2,KcLimiter);
        const double Kr2 = Kf2/Kc2;         
        
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;

        lhsIndex = lhsIndex + 3;
        rhsIndex = rhsIndex + 9;
    }


}


void 
FastChemistry::OptReaction::RF13IR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{

    std::size_t end = this->irreversibleReaction13index.size();


    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;

    std::size_t remain = end%4;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->irreversibleReaction13index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D13IR[lhsIndex+0];

        const unsigned int sr0a = rhsSpeciesIndex1D13IR[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D13IR[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D13IR[rhsIndex+2];

        double Kf0 = Kf_[i0];

        
        const double CF0 = c[sl0a];
        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;


        const unsigned int i1 = this->irreversibleReaction13index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D13IR[lhsIndex+1];

        const unsigned int sr0b = rhsSpeciesIndex1D13IR[rhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D13IR[rhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D13IR[rhsIndex+5];

        double Kf1 = Kf_[i1];
 
        
        const double CF1 = c[sl0b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;


        const unsigned int i2 = this->irreversibleReaction13index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D13IR[lhsIndex+2];

        const unsigned int sr0c = rhsSpeciesIndex1D13IR[rhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D13IR[rhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D13IR[rhsIndex+8];

        double Kf2 = Kf_[i2];

        const double CF2 = c[sl0c];

        const double q2 = (Kf2*CF2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;


        const unsigned int i3 = this->irreversibleReaction13index[k+3];

        const unsigned int sl0d = lhsSpeciesIndex1D13IR[lhsIndex+3];

        const unsigned int sr0d = rhsSpeciesIndex1D13IR[rhsIndex+9];
        const unsigned int sr1d = rhsSpeciesIndex1D13IR[rhsIndex+10];
        const unsigned int sr2d = rhsSpeciesIndex1D13IR[rhsIndex+11];

        double Kf3 = Kf_[i3];

        
        const double CF3 = c[sl0d];

        const double q3 = (Kf3*CF3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;

        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
        dNdtByV[sr2d] = dNdtByV[sr2d] + q3;



        lhsIndex = lhsIndex + 4;
        rhsIndex = rhsIndex + 12;
    }
    if(remain==1)
    {
        std::size_t k = end-1;
        const unsigned int i0 = this->irreversibleReaction13index[k];

        const unsigned int sl0a = lhsSpeciesIndex1D13IR[lhsIndex+0];

        const unsigned int sr0a = rhsSpeciesIndex1D13IR[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D13IR[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D13IR[rhsIndex+2];



        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 3;
    }
    else if(remain==2)
    {
        std::size_t k = end-2;
        const unsigned int i0 = this->irreversibleReaction13index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D13IR[lhsIndex+0];

        const unsigned int sr0a = rhsSpeciesIndex1D13IR[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D13IR[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D13IR[rhsIndex+2];

        double Kf0 = Kf_[i0];

        
        const double CF0 = c[sl0a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;


        const unsigned int i1 = this->irreversibleReaction13index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D13IR[lhsIndex+1];

        const unsigned int sr0b = rhsSpeciesIndex1D13IR[rhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D13IR[rhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D13IR[rhsIndex+5];

        double Kf1 = Kf_[i1];
 
        
        const double CF1 = c[sl0b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;


        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 6;
    }
    else if(remain==3)
    {
        std::size_t k = end-3;
        const unsigned int i0 = this->irreversibleReaction13index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D13IR[lhsIndex+0];

        const unsigned int sr0a = rhsSpeciesIndex1D13IR[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D13IR[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D13IR[rhsIndex+2];

        double Kf0 = Kf_[i0];

        
        const double CF0 = c[sl0a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;


        const unsigned int i1 = this->irreversibleReaction13index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D13IR[lhsIndex+1];

        const unsigned int sr0b = rhsSpeciesIndex1D13IR[rhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D13IR[rhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D13IR[rhsIndex+5];

        double Kf1 = Kf_[i1];

        
        const double CF1 = c[sl0b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;


        const unsigned int i2 = this->irreversibleReaction13index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D13IR[lhsIndex+2];

        const unsigned int sr0c = rhsSpeciesIndex1D13IR[rhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D13IR[rhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D13IR[rhsIndex+8];

        double Kf2 = Kf_[i2];

        
        const double CF2 = c[sl0c];

        const double q2 = (Kf2*CF2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;

        lhsIndex = lhsIndex + 3;
        rhsIndex = rhsIndex + 9;

    }
}
void 
FastChemistry::OptReaction::RF13NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = this->nonEquilibriumReaction13index.size();


    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;

    std::size_t remain = end%4;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->nonEquilibriumReaction13index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D13NER[lhsIndex+0];

        const unsigned int sr0a = rhsSpeciesIndex1D13NER[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D13NER[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D13NER[rhsIndex+2];

        double Kf0 = Kf_[i0];
    
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];

        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;


        const unsigned int i1 = this->nonEquilibriumReaction13index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D13NER[lhsIndex+1];

        const unsigned int sr0b = rhsSpeciesIndex1D13NER[rhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D13NER[rhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D13NER[rhsIndex+5];

        double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;


        const unsigned int i2 = this->nonEquilibriumReaction13index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D13NER[lhsIndex+2];

        const unsigned int sr0c = rhsSpeciesIndex1D13NER[rhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D13NER[rhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D13NER[rhsIndex+8];

        double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];  
        
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;


        const unsigned int i3 = this->nonEquilibriumReaction13index[k+3];

        const unsigned int sl0d = lhsSpeciesIndex1D13NER[lhsIndex+3];

        const unsigned int sr0d = rhsSpeciesIndex1D13NER[rhsIndex+9];
        const unsigned int sr1d = rhsSpeciesIndex1D13NER[rhsIndex+10];
        const unsigned int sr2d = rhsSpeciesIndex1D13NER[rhsIndex+11];

        double Kf3 = Kf_[i3];
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        const double Kr3 = this->Kf_[l3];  
        
        const double CF3 = c[sl0d];
        const double CR3 = c[sr0d]*c[sr1d]*c[sr2d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;

        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
        dNdtByV[sr2d] = dNdtByV[sr2d] + q3;



        lhsIndex = lhsIndex + 4;
        rhsIndex = rhsIndex + 12;
    }
    if(remain==1)
    {
        std::size_t k = end-1;
        const unsigned int i0 = this->nonEquilibriumReaction13index[k];

        const unsigned int sl0a = lhsSpeciesIndex1D13NER[lhsIndex+0];

        const unsigned int sr0a = rhsSpeciesIndex1D13NER[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D13NER[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D13NER[rhsIndex+2];



        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];  
        
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        lhsIndex = lhsIndex + 1;
        rhsIndex = rhsIndex + 3;
    }
    else if(remain==2)
    {
        std::size_t k = end-2;
        const unsigned int i0 = this->nonEquilibriumReaction13index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D13NER[lhsIndex+0];

        const unsigned int sr0a = rhsSpeciesIndex1D13NER[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D13NER[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D13NER[rhsIndex+2];

        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];
        
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;


        const unsigned int i1 = this->nonEquilibriumReaction13index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D13NER[lhsIndex+1];

        const unsigned int sr0b = rhsSpeciesIndex1D13NER[rhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D13NER[rhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D13NER[rhsIndex+5];

        double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];
        
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;


        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 6;
    }
    else if(remain==3)
    {
        std::size_t k = end-3;
        const unsigned int i0 = this->nonEquilibriumReaction13index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D13NER[lhsIndex+0];

        const unsigned int sr0a = rhsSpeciesIndex1D13NER[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D13NER[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D13NER[rhsIndex+2];

        double Kf0 = Kf_[i0];
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        const double Kr0 = this->Kf_[l0];     
        
        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;


        const unsigned int i1 = this->nonEquilibriumReaction13index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D13NER[lhsIndex+1];

        const unsigned int sr0b = rhsSpeciesIndex1D13NER[rhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D13NER[rhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D13NER[rhsIndex+5];

        double Kf1 = Kf_[i1];
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        const double Kr1 = this->Kf_[l1];   
        
        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;


        const unsigned int i2 = this->nonEquilibriumReaction13index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D13NER[lhsIndex+2];

        const unsigned int sr0c = rhsSpeciesIndex1D13NER[rhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D13NER[rhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D13NER[rhsIndex+8];

        double Kf2 = Kf_[i2];
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        const double Kr2 = this->Kf_[l2];
        
        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;

        lhsIndex = lhsIndex + 3;
        rhsIndex = rhsIndex + 9;
    }
}