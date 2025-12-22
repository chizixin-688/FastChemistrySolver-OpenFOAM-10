#include "OptReaction.H"

void 
OptReaction::update33Reaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reaction33index.size();
    //std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    std::size_t lhsIndex = 0;
    for(std::size_t i=0; i<end; i++)
    {
        const unsigned int KfIndex = this->reaction33index[i];

        const unsigned int sl0 = lhsSpeciesIndex1D33[rhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D33[rhsIndex+1];
        const unsigned int sl2 = lhsSpeciesIndex1D33[rhsIndex+2];

        const unsigned int sr0 = rhsSpeciesIndex1D33[lhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D33[lhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D33[lhsIndex+2];

        double Kr = 0;
        double Kf = Kf_[KfIndex];
        if(this->isIrreversible[KfIndex]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
            double Kc = Kp;
            Kc = std::max(Kc,KcLimiter);
            Kr = Kf/Kc;         
        }
        else if(this->isIrreversible[KfIndex]==2)
        {
            unsigned int l = KfIndex - this->Ikf[1] + this->Ikf[9];
            Kr = this->Kf_[l];
        }

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        rhsIndex = rhsIndex + 3;
        lhsIndex = lhsIndex + 3;
    }
}


void 
OptReaction::RF33RR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = this->reversibleReaction33index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    std::size_t lhsIndex = 0;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->reversibleReaction33index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D33RR[rhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D33RR[rhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D33RR[rhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D33RR[lhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D33RR[lhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D33RR[lhsIndex+2];


        double Kf0 = Kf_[i0];
        const double Kp0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a]*ExpNegGbyRT[sr2a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a]);
        double Kc0 = Kp0;
        Kc0 = std::max(Kc0,KcLimiter);
        double Kr0 = Kf0/Kc0;         
        
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;
        
        const unsigned int i1 = this->reversibleReaction33index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D33RR[rhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D33RR[rhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D33RR[rhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D33RR[lhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D33RR[lhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D33RR[lhsIndex+5];


        double Kf1 = Kf_[i1];
        const double Kp1 = (ExpNegGbyRT[sr0b]*ExpNegGbyRT[sr1b]*ExpNegGbyRT[sr2b])/(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b]);
        double Kc1 = Kp1;
        Kc1 = std::max(Kc1,KcLimiter);
        double Kr1 = Kf1/Kc1;         
        
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;
        
        
        const unsigned int i2 = this->reversibleReaction33index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D33RR[rhsIndex+6];
        const unsigned int sl1c = lhsSpeciesIndex1D33RR[rhsIndex+7];
        const unsigned int sl2c = lhsSpeciesIndex1D33RR[rhsIndex+8];

        const unsigned int sr0c = rhsSpeciesIndex1D33RR[lhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D33RR[lhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D33RR[lhsIndex+8];


        double Kf2 = Kf_[i2];
        const double Kp2 = (ExpNegGbyRT[sr0c]*ExpNegGbyRT[sr1c]*ExpNegGbyRT[sr2c])/(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c]);
        double Kc2 = Kp2;
        Kc2 = std::max(Kc2,KcLimiter);
        double Kr2 = Kf2/Kc2;         
        
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;
    
        const unsigned int i3 = this->reversibleReaction33index[k+3];

        const unsigned int sl0d = lhsSpeciesIndex1D33RR[rhsIndex+9];
        const unsigned int sl1d = lhsSpeciesIndex1D33RR[rhsIndex+10];
        const unsigned int sl2d = lhsSpeciesIndex1D33RR[rhsIndex+11];

        const unsigned int sr0d = rhsSpeciesIndex1D33RR[lhsIndex+9];
        const unsigned int sr1d = rhsSpeciesIndex1D33RR[lhsIndex+10];
        const unsigned int sr2d = rhsSpeciesIndex1D33RR[lhsIndex+11];


        double Kf3 = Kf_[i3];
        const double Kp3 = (ExpNegGbyRT[sr0d]*ExpNegGbyRT[sr1d]*ExpNegGbyRT[sr2d])/(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]*ExpNegGbyRT[sl2d]);
        double Kc3 = Kp3;
        Kc3 = std::max(Kc3,KcLimiter);
        double Kr3 = Kf3/Kc3;
        
        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d]*c[sr1d]*c[sr2d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;

        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
        dNdtByV[sr2d] = dNdtByV[sr2d] + q3;
        
        rhsIndex = rhsIndex + 12;
        lhsIndex = lhsIndex + 12;
    }
    if(remain==1)
    {
        std::size_t k = end-1;
        const unsigned int i0 = this->reversibleReaction33index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D33RR[rhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D33RR[rhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D33RR[rhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D33RR[lhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D33RR[lhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D33RR[lhsIndex+2];


        double Kf0 = Kf_[i0];
        const double Kp0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a]*ExpNegGbyRT[sr2a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a]);
        double Kc0 = Kp0;
        Kc0 = std::max(Kc0,KcLimiter);
        double Kr0 = Kf0/Kc0;         
        
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;
        
        rhsIndex = rhsIndex + 3;
        lhsIndex = lhsIndex + 3;
    }
    else if(remain==2)
    {
        std::size_t k = end-2;
        const unsigned int i0 = this->reversibleReaction33index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D33RR[rhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D33RR[rhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D33RR[rhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D33RR[lhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D33RR[lhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D33RR[lhsIndex+2];


        double Kf0 = Kf_[i0];
        const double Kp0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a]*ExpNegGbyRT[sr2a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a]);
        double Kc0 = Kp0;
        Kc0 = std::max(Kc0,KcLimiter);
        double Kr0 = Kf0/Kc0;         
        
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;
        
        const unsigned int i1 = this->reversibleReaction33index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D33RR[rhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D33RR[rhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D33RR[rhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D33RR[lhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D33RR[lhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D33RR[lhsIndex+5];


        double Kf1 = Kf_[i1];
        const double Kp1 = (ExpNegGbyRT[sr0b]*ExpNegGbyRT[sr1b]*ExpNegGbyRT[sr2b])/(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b]);
        double Kc1 = Kp1;
        Kc1 = std::max(Kc1,KcLimiter);
        double Kr1 = Kf1/Kc1;         
        
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;

        rhsIndex = rhsIndex + 6;
        lhsIndex = lhsIndex + 6;
    }
    else if(remain==3)
    {
        std::size_t k = end-3;
        const unsigned int i0 = this->reversibleReaction33index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D33RR[rhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D33RR[rhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D33RR[rhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D33RR[lhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D33RR[lhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D33RR[lhsIndex+2];


        double Kf0 = Kf_[i0];
        const double Kp0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a]*ExpNegGbyRT[sr2a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a]);
        double Kc0 = Kp0;
        Kc0 = std::max(Kc0,KcLimiter);
        double Kr0 = Kf0/Kc0;         
        
        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;
        
        const unsigned int i1 = this->reversibleReaction33index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D33RR[rhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D33RR[rhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D33RR[rhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D33RR[lhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D33RR[lhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D33RR[lhsIndex+5];


        double Kf1 = Kf_[i1];
        const double Kp1 = (ExpNegGbyRT[sr0b]*ExpNegGbyRT[sr1b]*ExpNegGbyRT[sr2b])/(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b]);
        double Kc1 = Kp1;
        Kc1 = std::max(Kc1,KcLimiter);
        double Kr1 = Kf1/Kc1;         
        
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;
        
        
        const unsigned int i2 = this->reversibleReaction33index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D33RR[rhsIndex+6];
        const unsigned int sl1c = lhsSpeciesIndex1D33RR[rhsIndex+7];
        const unsigned int sl2c = lhsSpeciesIndex1D33RR[rhsIndex+8];

        const unsigned int sr0c = rhsSpeciesIndex1D33RR[lhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D33RR[lhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D33RR[lhsIndex+8];


        double Kf2 = Kf_[i2];
        const double Kp2 = (ExpNegGbyRT[sr0c]*ExpNegGbyRT[sr1c]*ExpNegGbyRT[sr2c])/(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c]);
        double Kc2 = Kp2;
        Kc2 = std::max(Kc2,KcLimiter);
        double Kr2 = Kf2/Kc2;         
        
        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;

        rhsIndex = rhsIndex + 9;
        lhsIndex = lhsIndex + 9;
    }
}


void 
OptReaction::RF33IR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = this->irreversibleReaction33index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    std::size_t lhsIndex = 0;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->irreversibleReaction33index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D33IR[rhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D33IR[rhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D33IR[rhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D33IR[lhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D33IR[lhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D33IR[lhsIndex+2];


        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;
        
        const unsigned int i1 = this->irreversibleReaction33index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D33IR[rhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D33IR[rhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D33IR[rhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D33IR[lhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D33IR[lhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D33IR[lhsIndex+5];


        double Kf1 = Kf_[i1];

        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;
        
        
        const unsigned int i2 = this->irreversibleReaction33index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D33IR[rhsIndex+6];
        const unsigned int sl1c = lhsSpeciesIndex1D33IR[rhsIndex+7];
        const unsigned int sl2c = lhsSpeciesIndex1D33IR[rhsIndex+8];

        const unsigned int sr0c = rhsSpeciesIndex1D33IR[lhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D33IR[lhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D33IR[lhsIndex+8];


        double Kf2 = Kf_[i2];

        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];

        const double q2 = (Kf2*CF2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;
        
        
        const unsigned int i3 = this->irreversibleReaction33index[k+3];

        const unsigned int sl0d = lhsSpeciesIndex1D33IR[rhsIndex+9];
        const unsigned int sl1d = lhsSpeciesIndex1D33IR[rhsIndex+10];
        const unsigned int sl2d = lhsSpeciesIndex1D33IR[rhsIndex+11];

        const unsigned int sr0d = rhsSpeciesIndex1D33IR[lhsIndex+9];
        const unsigned int sr1d = rhsSpeciesIndex1D33IR[lhsIndex+10];
        const unsigned int sr2d = rhsSpeciesIndex1D33IR[lhsIndex+11];


        double Kf3 = Kf_[i3];

        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];

        const double q3 = (Kf3*CF3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;

        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
        dNdtByV[sr2d] = dNdtByV[sr2d] + q3;
        

        rhsIndex = rhsIndex + 12;
        lhsIndex = lhsIndex + 12;
    }
    if(remain==1)
    {
        std::size_t k = end-1;
        const unsigned int i0 = this->irreversibleReaction33index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D33IR[rhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D33IR[rhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D33IR[rhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D33IR[lhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D33IR[lhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D33IR[lhsIndex+2];


        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;
    
        
        rhsIndex = rhsIndex + 3;
        lhsIndex = lhsIndex + 3;
    }
    else if(remain==2)
    {
        std::size_t k = end-2;

        const unsigned int i0 = this->irreversibleReaction33index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D33IR[rhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D33IR[rhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D33IR[rhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D33IR[lhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D33IR[lhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D33IR[lhsIndex+2];


        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;
        
        const unsigned int i1 = this->irreversibleReaction33index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D33IR[rhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D33IR[rhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D33IR[rhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D33IR[lhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D33IR[lhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D33IR[lhsIndex+5];


        double Kf1 = Kf_[i1];

        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;
        
        rhsIndex = rhsIndex + 6;
        lhsIndex = lhsIndex + 6;
    }
    else if(remain==3)
    {
        std::size_t k = end-3;
        const unsigned int i0 = this->irreversibleReaction33index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D33IR[rhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D33IR[rhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D33IR[rhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D33IR[lhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D33IR[lhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D33IR[lhsIndex+2];


        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;
        
        const unsigned int i1 = this->irreversibleReaction33index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D33IR[rhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D33IR[rhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D33IR[rhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D33IR[lhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D33IR[lhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D33IR[lhsIndex+5];


        double Kf1 = Kf_[i1];

        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;
        
        
        const unsigned int i2 = this->irreversibleReaction33index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D33IR[rhsIndex+6];
        const unsigned int sl1c = lhsSpeciesIndex1D33IR[rhsIndex+7];
        const unsigned int sl2c = lhsSpeciesIndex1D33IR[rhsIndex+8];

        const unsigned int sr0c = rhsSpeciesIndex1D33IR[lhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D33IR[lhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D33IR[lhsIndex+8];


        double Kf2 = Kf_[i2];

        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];

        const double q2 = (Kf2*CF2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;
        

        rhsIndex = rhsIndex + 9;
        lhsIndex = lhsIndex + 9;
    }
}



void 
OptReaction::RF33NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reaction33index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    std::size_t lhsIndex = 0;
    for(std::size_t k=0; k<end; k++)
    {

        {
        const unsigned int i0 = this->reaction33index[k+0];

        const unsigned int sl0 = lhsSpeciesIndex1D33[rhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D33[rhsIndex+1];
        const unsigned int sl2 = lhsSpeciesIndex1D33[rhsIndex+2];

        const unsigned int sr0 = rhsSpeciesIndex1D33[lhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D33[lhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D33[lhsIndex+2];

        double Kr = 0;
        double Kf = Kf_[i0];


        
        unsigned int l = i0 - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[l];
        

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        }
        {
        const unsigned int i0 = this->reaction33index[k+1];

        const unsigned int sl0 = lhsSpeciesIndex1D33[rhsIndex+3];
        const unsigned int sl1 = lhsSpeciesIndex1D33[rhsIndex+4];
        const unsigned int sl2 = lhsSpeciesIndex1D33[rhsIndex+5];

        const unsigned int sr0 = rhsSpeciesIndex1D33[lhsIndex+3];
        const unsigned int sr1 = rhsSpeciesIndex1D33[lhsIndex+4];
        const unsigned int sr2 = rhsSpeciesIndex1D33[lhsIndex+5];

        double Kr = 0;
        double Kf = Kf_[i0];


        
        unsigned int l = i0 - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[l];
        

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        }
        {
        const unsigned int i0 = this->reaction33index[k+2];

        const unsigned int sl0 = lhsSpeciesIndex1D33[rhsIndex+6];
        const unsigned int sl1 = lhsSpeciesIndex1D33[rhsIndex+7];
        const unsigned int sl2 = lhsSpeciesIndex1D33[rhsIndex+8];

        const unsigned int sr0 = rhsSpeciesIndex1D33[lhsIndex+6];
        const unsigned int sr1 = rhsSpeciesIndex1D33[lhsIndex+7];
        const unsigned int sr2 = rhsSpeciesIndex1D33[lhsIndex+8];

        double Kr = 0;
        double Kf = Kf_[i0];


        
        unsigned int l = i0 - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[l];
        

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        }
        {
        const unsigned int i0 = this->reaction33index[k+3];

        const unsigned int sl0 = lhsSpeciesIndex1D33[rhsIndex+9];
        const unsigned int sl1 = lhsSpeciesIndex1D33[rhsIndex+10];
        const unsigned int sl2 = lhsSpeciesIndex1D33[rhsIndex+11];

        const unsigned int sr0 = rhsSpeciesIndex1D33[lhsIndex+9];
        const unsigned int sr1 = rhsSpeciesIndex1D33[lhsIndex+10];
        const unsigned int sr2 = rhsSpeciesIndex1D33[lhsIndex+11];

        double Kr = 0;
        double Kf = Kf_[i0];


        
        unsigned int l = i0 - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[l];
        

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        }

        rhsIndex = rhsIndex + 12;
        lhsIndex = lhsIndex + 12;
    }
    if(remain==1)
    {
        std::size_t k = end-1;
        {
        const unsigned int i0 = this->reaction33index[k+0];

        const unsigned int sl0 = lhsSpeciesIndex1D33[rhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D33[rhsIndex+1];
        const unsigned int sl2 = lhsSpeciesIndex1D33[rhsIndex+2];

        const unsigned int sr0 = rhsSpeciesIndex1D33[lhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D33[lhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D33[lhsIndex+2];

        double Kr = 0;
        double Kf = Kf_[i0];


        
        unsigned int l = i0 - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[l];
        

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        }
        rhsIndex = rhsIndex + 3;
        lhsIndex = lhsIndex + 3;
    }
    else if(remain==2)
    {
        std::size_t k = end-2;
        {
        const unsigned int i0 = this->reaction33index[k+0];

        const unsigned int sl0 = lhsSpeciesIndex1D33[rhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D33[rhsIndex+1];
        const unsigned int sl2 = lhsSpeciesIndex1D33[rhsIndex+2];

        const unsigned int sr0 = rhsSpeciesIndex1D33[lhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D33[lhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D33[lhsIndex+2];

        double Kr = 0;
        double Kf = Kf_[i0];


        
        unsigned int l = i0 - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[l];
        

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        }
        {
        const unsigned int i0 = this->reaction33index[k+1];

        const unsigned int sl0 = lhsSpeciesIndex1D33[rhsIndex+3];
        const unsigned int sl1 = lhsSpeciesIndex1D33[rhsIndex+4];
        const unsigned int sl2 = lhsSpeciesIndex1D33[rhsIndex+5];

        const unsigned int sr0 = rhsSpeciesIndex1D33[lhsIndex+3];
        const unsigned int sr1 = rhsSpeciesIndex1D33[lhsIndex+4];
        const unsigned int sr2 = rhsSpeciesIndex1D33[lhsIndex+5];

        double Kr = 0;
        double Kf = Kf_[i0];


        
        unsigned int l = i0 - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[l];
        

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        }
        rhsIndex = rhsIndex + 6;
        lhsIndex = lhsIndex + 6;
    }
    else if(remain==3)
    {
        std::size_t k = end-3;
        {
        const unsigned int i0 = this->reaction33index[k+0];

        const unsigned int sl0 = lhsSpeciesIndex1D33[rhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D33[rhsIndex+1];
        const unsigned int sl2 = lhsSpeciesIndex1D33[rhsIndex+2];

        const unsigned int sr0 = rhsSpeciesIndex1D33[lhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D33[lhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D33[lhsIndex+2];

        double Kr = 0;
        double Kf = Kf_[i0];


        
        unsigned int l = i0 - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[l];
        

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        }
        {
        const unsigned int i0 = this->reaction33index[k+1];

        const unsigned int sl0 = lhsSpeciesIndex1D33[rhsIndex+3];
        const unsigned int sl1 = lhsSpeciesIndex1D33[rhsIndex+4];
        const unsigned int sl2 = lhsSpeciesIndex1D33[rhsIndex+5];

        const unsigned int sr0 = rhsSpeciesIndex1D33[lhsIndex+3];
        const unsigned int sr1 = rhsSpeciesIndex1D33[lhsIndex+4];
        const unsigned int sr2 = rhsSpeciesIndex1D33[lhsIndex+5];

        double Kr = 0;
        double Kf = Kf_[i0];


        
        unsigned int l = i0 - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[l];
        

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        }
        {
        const unsigned int i0 = this->reaction33index[k+2];

        const unsigned int sl0 = lhsSpeciesIndex1D33[rhsIndex+6];
        const unsigned int sl1 = lhsSpeciesIndex1D33[rhsIndex+7];
        const unsigned int sl2 = lhsSpeciesIndex1D33[rhsIndex+8];

        const unsigned int sr0 = rhsSpeciesIndex1D33[lhsIndex+6];
        const unsigned int sr1 = rhsSpeciesIndex1D33[lhsIndex+7];
        const unsigned int sr2 = rhsSpeciesIndex1D33[lhsIndex+8];

        double Kr = 0;
        double Kf = Kf_[i0];


        
        unsigned int l = i0 - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[l];
        

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        }
        rhsIndex = rhsIndex + 9;
        lhsIndex = lhsIndex + 9;
    }
}