#include "OptReaction.H"

void 
OptReaction::update31Reaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reaction31index.size();
    //std::size_t remain = end%4;

    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    
    for(std::size_t i=0; i<end; i++)
    {
        const unsigned int KfIndex = this->reaction31index[i];

        const unsigned int sl0 = lhsSpeciesIndex1D31[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D31[lhsIndex+1];
        const unsigned int sl2 = lhsSpeciesIndex1D31[lhsIndex+2];

        const unsigned int sr0 = rhsSpeciesIndex1D31[rhsIndex+0];


        double Kr = 0;
        double Kf = Kf_[KfIndex];
        if(this->isIrreversible[KfIndex]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
            double Kc = Kp*this->Pow_pByRT_SumVki[0];
            Kc = std::max(Kc,KcLimiter);
            Kr = Kf/Kc;         
        }
        else if(this->isIrreversible[KfIndex]==2)
        {
            unsigned int l = KfIndex - this->Ikf[1] + this->Ikf[9];
            Kr = this->Kf_[l];
        }

        const double CF = c[sl0]*c[sl1]*c[sl2];
        const double CR = c[sr0];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sl2] = dNdtByV[sl2] - q;

        dNdtByV[sr0] = dNdtByV[sr0] + q;


        lhsIndex = lhsIndex + 3;
        rhsIndex = rhsIndex + 1;

    }
}


void 
OptReaction::RF31RR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reversibleReaction31index.size();
    std::size_t remain = end%4;

    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->reversibleReaction31index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D31RR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D31RR[lhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D31RR[lhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D31RR[rhsIndex+0];


        double Kr0 = 0;
        double Kf0 = Kf_[i0];

        
        const double Kp0 = (ExpNegGbyRT[sr0a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a]);
        double Kc0 = Kp0*this->Pow_pByRT_SumVki[0];
        Kc0 = std::max(Kc0,KcLimiter);
        Kr0 = Kf0/Kc0;
        

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        
        const unsigned int i1 = this->reversibleReaction31index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D31RR[lhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D31RR[lhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D31RR[lhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D31RR[rhsIndex+1];


        double Kr1 = 0;
        double Kf1 = Kf_[i1];

        
        const double Kp1 = (ExpNegGbyRT[sr0b])/(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b]);
        double Kc1 = Kp1*this->Pow_pByRT_SumVki[0];
        Kc1 = std::max(Kc1,KcLimiter);
        Kr1 = Kf1/Kc1;
        

        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        
        
        const unsigned int i2 = this->reversibleReaction31index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D31RR[lhsIndex+6];
        const unsigned int sl1c = lhsSpeciesIndex1D31RR[lhsIndex+7];
        const unsigned int sl2c = lhsSpeciesIndex1D31RR[lhsIndex+8];

        const unsigned int sr0c = rhsSpeciesIndex1D31RR[rhsIndex+2];


        double Kr2 = 0;
        double Kf2 = Kf_[i2];

        
        const double Kp2 = (ExpNegGbyRT[sr0c])/(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c]);
        double Kc2 = Kp2*this->Pow_pByRT_SumVki[0];
        Kc2 = std::max(Kc2,KcLimiter);
        Kr2 = Kf2/Kc2;
        

        

        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        
        
        const unsigned int i3 = this->reversibleReaction31index[k+3];

        const unsigned int sl0d = lhsSpeciesIndex1D31RR[lhsIndex+9];
        const unsigned int sl1d = lhsSpeciesIndex1D31RR[lhsIndex+10];
        const unsigned int sl2d = lhsSpeciesIndex1D31RR[lhsIndex+11];

        const unsigned int sr0d = rhsSpeciesIndex1D31RR[rhsIndex+3];


        double Kr3 = 0;
        double Kf3 = Kf_[i3];

        
        const double Kp3 = (ExpNegGbyRT[sr0d])/(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]*ExpNegGbyRT[sl2d]);
        double Kc3 = Kp3*this->Pow_pByRT_SumVki[0];
        Kc3 = std::max(Kc3,KcLimiter);
        Kr3 = Kf3/Kc3;
        

        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;

        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        

        lhsIndex = lhsIndex + 12;
        rhsIndex = rhsIndex + 4;
    }
    if(remain==1)
    {
        std::size_t k = end-1;
        const unsigned int i0 = this->reversibleReaction31index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D31RR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D31RR[lhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D31RR[lhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D31RR[rhsIndex+0];


        double Kr0 = 0;
        double Kf0 = Kf_[i0];

        
        const double Kp0 = (ExpNegGbyRT[sr0a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a]);
        double Kc0 = Kp0*this->Pow_pByRT_SumVki[0];
        Kc0 = std::max(Kc0,KcLimiter);
        Kr0 = Kf0/Kc0;
        

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;


        lhsIndex = lhsIndex + 3;
        rhsIndex = rhsIndex + 1;
    }
    else if(remain==2)
    {
        std::size_t k = end-2;
        const unsigned int i0 = this->reversibleReaction31index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D31RR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D31RR[lhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D31RR[lhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D31RR[rhsIndex+0];


        double Kr0 = 0;
        double Kf0 = Kf_[i0];

        
        const double Kp0 = (ExpNegGbyRT[sr0a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a]);
        double Kc0 = Kp0*this->Pow_pByRT_SumVki[0];
        Kc0 = std::max(Kc0,KcLimiter);
        Kr0 = Kf0/Kc0;
        

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        
        const unsigned int i1 = this->reversibleReaction31index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D31RR[lhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D31RR[lhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D31RR[lhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D31RR[rhsIndex+1];


        double Kr1 = 0;
        double Kf1 = Kf_[i1];

        
        const double Kp1 = (ExpNegGbyRT[sr0b])/(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b]);
        double Kc1 = Kp1*this->Pow_pByRT_SumVki[0];
        Kc1 = std::max(Kc1,KcLimiter);
        Kr1 = Kf1/Kc1;
        

        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
    

        lhsIndex = lhsIndex + 6;
        rhsIndex = rhsIndex + 2;
    }
    else if(remain==3)
    {
        std::size_t k = end-3;
        const unsigned int i0 = this->reversibleReaction31index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D31RR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D31RR[lhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D31RR[lhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D31RR[rhsIndex+0];


        double Kr0 = 0;
        double Kf0 = Kf_[i0];

        
        const double Kp0 = (ExpNegGbyRT[sr0a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a]);
        double Kc0 = Kp0*this->Pow_pByRT_SumVki[0];
        Kc0 = std::max(Kc0,KcLimiter);
        Kr0 = Kf0/Kc0;
        

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        
        const unsigned int i1 = this->reversibleReaction31index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D31RR[lhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D31RR[lhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D31RR[lhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D31RR[rhsIndex+1];


        double Kr1 = 0;
        double Kf1 = Kf_[i1];

        
        const double Kp1 = (ExpNegGbyRT[sr0b])/(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]*ExpNegGbyRT[sl2b]);
        double Kc1 = Kp1*this->Pow_pByRT_SumVki[0];
        Kc1 = std::max(Kc1,KcLimiter);
        Kr1 = Kf1/Kc1;
        

        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        
        
        const unsigned int i2 = this->reversibleReaction31index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D31RR[lhsIndex+6];
        const unsigned int sl1c = lhsSpeciesIndex1D31RR[lhsIndex+7];
        const unsigned int sl2c = lhsSpeciesIndex1D31RR[lhsIndex+8];

        const unsigned int sr0c = rhsSpeciesIndex1D31RR[rhsIndex+2];


        double Kr2 = 0;
        double Kf2 = Kf_[i2];

        
        const double Kp2 = (ExpNegGbyRT[sr0c])/(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]*ExpNegGbyRT[sl2c]);
        double Kc2 = Kp2*this->Pow_pByRT_SumVki[0];
        Kc2 = std::max(Kc2,KcLimiter);
        Kr2 = Kf2/Kc2;
        

        

        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
    
        lhsIndex = lhsIndex + 9;
        rhsIndex = rhsIndex + 3;
    }
}


void 
OptReaction::RF31IR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = irreversibleReaction31index.size();
    std::size_t remain = end%4;

    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->irreversibleReaction31index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D31IR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D31IR[lhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D31IR[lhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D31IR[rhsIndex+0];


        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double q0 = (Kf0*CF0) ;

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        
        const unsigned int i1 = this->irreversibleReaction31index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D31IR[lhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D31IR[lhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D31IR[lhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D31IR[rhsIndex+1];



        double Kf1 = Kf_[i1];
        
        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        
        
        const unsigned int i2 = this->irreversibleReaction31index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D31IR[lhsIndex+6];
        const unsigned int sl1c = lhsSpeciesIndex1D31IR[lhsIndex+7];
        const unsigned int sl2c = lhsSpeciesIndex1D31IR[lhsIndex+8];

        const unsigned int sr0c = rhsSpeciesIndex1D31IR[rhsIndex+2];



        double Kf2 = Kf_[i2];

        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];

        const double q2 = (Kf2*CF2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        
        
        const unsigned int i3 = this->irreversibleReaction31index[k+3];

        const unsigned int sl0d = lhsSpeciesIndex1D31IR[lhsIndex+9];
        const unsigned int sl1d = lhsSpeciesIndex1D31IR[lhsIndex+10];
        const unsigned int sl2d = lhsSpeciesIndex1D31IR[lhsIndex+11];

        const unsigned int sr0d = rhsSpeciesIndex1D31IR[rhsIndex+3];



        double Kf3 = Kf_[i3];

        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];

        const double q3 = (Kf3*CF3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;

        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        

        lhsIndex = lhsIndex + 12;
        rhsIndex = rhsIndex + 4;
    }
    if(remain==1)
    {
        std::size_t k = end-1;
        const unsigned int i0 = this->irreversibleReaction31index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D31IR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D31IR[lhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D31IR[lhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D31IR[rhsIndex+0];

        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;


        lhsIndex = lhsIndex + 3;
        rhsIndex = rhsIndex + 1;
    }
    else if(remain==2)
    {
        std::size_t k = end-2;
        const unsigned int i0 = this->irreversibleReaction31index[k+1];

        const unsigned int sl0a = lhsSpeciesIndex1D31IR[lhsIndex+3];
        const unsigned int sl1a = lhsSpeciesIndex1D31IR[lhsIndex+4];
        const unsigned int sl2a = lhsSpeciesIndex1D31IR[lhsIndex+5];

        const unsigned int sr0a = rhsSpeciesIndex1D31IR[rhsIndex+1];


        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];

        const double q0 = (Kf0*CF0) ;

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        
        const unsigned int i1 = this->irreversibleReaction31index[k+2];

        const unsigned int sl0b = lhsSpeciesIndex1D31IR[lhsIndex+6];
        const unsigned int sl1b = lhsSpeciesIndex1D31IR[lhsIndex+7];
        const unsigned int sl2b = lhsSpeciesIndex1D31IR[lhsIndex+8];

        const unsigned int sr0b = rhsSpeciesIndex1D31IR[rhsIndex+2];

        double Kf1 = Kf_[i1];


        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
    

        lhsIndex = lhsIndex + 6;
        rhsIndex = rhsIndex + 2;
    }
    else if(remain==3)
    {
        std::size_t k = end-3;
        const unsigned int i0 = this->irreversibleReaction31index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D31IR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D31IR[lhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D31IR[lhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D31IR[rhsIndex+0];


        double Kr0 = 0;
        double Kf0 = Kf_[i0];

        
        const double Kp0 = (ExpNegGbyRT[sr0a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]*ExpNegGbyRT[sl2a]);
        double Kc0 = Kp0*this->Pow_pByRT_SumVki[0];
        Kc0 = std::max(Kc0,KcLimiter);
        Kr0 = Kf0/Kc0;
        

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        
        const unsigned int i1 = this->irreversibleReaction31index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D31IR[lhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D31IR[lhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D31IR[lhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D31IR[rhsIndex+1];


        double Kf1 = Kf_[i1];

        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        
        
        const unsigned int i2 = this->irreversibleReaction31index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D31IR[lhsIndex+6];
        const unsigned int sl1c = lhsSpeciesIndex1D31IR[lhsIndex+7];
        const unsigned int sl2c = lhsSpeciesIndex1D31IR[lhsIndex+8];

        const unsigned int sr0c = rhsSpeciesIndex1D31IR[rhsIndex+2];


        double Kf2 = Kf_[i2];

        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];

        const double q2 = (Kf2*CF2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
    
        lhsIndex = lhsIndex + 9;
        rhsIndex = rhsIndex + 3;
    }
}



void 
OptReaction::RF31NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = nonEquilibriumReaction31index.size();
    std::size_t remain = end%4;

    std::size_t lhsIndex = 0;
    std::size_t rhsIndex = 0;
    
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->nonEquilibriumReaction31index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D31NER[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D31NER[lhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D31NER[lhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D31NER[rhsIndex+0];

        double Kf0 = Kf_[i0];

        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        double Kr0 = this->Kf_[l0];

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        
        const unsigned int i1 = this->nonEquilibriumReaction31index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D31NER[lhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D31NER[lhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D31NER[lhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D31NER[rhsIndex+1];

        double Kf1 = Kf_[i1];

        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        double Kr1 = this->Kf_[l1];

        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        
        
        const unsigned int i2 = this->nonEquilibriumReaction31index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D31NER[lhsIndex+6];
        const unsigned int sl1c = lhsSpeciesIndex1D31NER[lhsIndex+7];
        const unsigned int sl2c = lhsSpeciesIndex1D31NER[lhsIndex+8];

        const unsigned int sr0c = rhsSpeciesIndex1D31NER[rhsIndex+2];

        double Kf2 = Kf_[i2];

        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        double Kr2 = this->Kf_[l2];

        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        
        
        const unsigned int i3 = this->nonEquilibriumReaction31index[k+3];

        const unsigned int sl0d = lhsSpeciesIndex1D31NER[lhsIndex+9];
        const unsigned int sl1d = lhsSpeciesIndex1D31NER[lhsIndex+10];
        const unsigned int sl2d = lhsSpeciesIndex1D31NER[lhsIndex+11];

        const unsigned int sr0d = rhsSpeciesIndex1D31NER[rhsIndex+3];

        double Kf3 = Kf_[i3];

        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        double Kr3 = this->Kf_[l3];

        const double CF3 = c[sl0d]*c[sl1d]*c[sl2d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sl2d] = dNdtByV[sl2d] - q3;

        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        


        lhsIndex = lhsIndex + 12;
        rhsIndex = rhsIndex + 4;

    }
    if(remain==1)
    {
        std::size_t k = end-1;
        const unsigned int i0 = this->nonEquilibriumReaction31index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D31NER[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D31NER[lhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D31NER[lhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D31NER[rhsIndex+0];

        double Kf0 = Kf_[i0];

        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        double Kr0 = this->Kf_[l0];

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        lhsIndex = lhsIndex + 3;
        rhsIndex = rhsIndex + 1;
    }
    else if(remain==2)
    {
        std::size_t k = end-2;
        const unsigned int i0 = this->nonEquilibriumReaction31index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D31NER[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D31NER[lhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D31NER[lhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D31NER[rhsIndex+0];

        double Kf0 = Kf_[i0];

        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        double Kr0 = this->Kf_[l0];

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        
        const unsigned int i1 = this->nonEquilibriumReaction31index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D31NER[lhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D31NER[lhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D31NER[lhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D31NER[rhsIndex+1];

        double Kf1 = Kf_[i1];

        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        double Kr1 = this->Kf_[l1];

        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;

        lhsIndex = lhsIndex + 6;
        rhsIndex = rhsIndex + 2;
    }
    else if(remain==3)
    {
        std::size_t k = end-3;
        const unsigned int i0 = this->nonEquilibriumReaction31index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D31NER[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D31NER[lhsIndex+1];
        const unsigned int sl2a = lhsSpeciesIndex1D31NER[lhsIndex+2];

        const unsigned int sr0a = rhsSpeciesIndex1D31NER[rhsIndex+0];



        double Kf0 = Kf_[i0];

        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        double Kr0 = this->Kf_[l0];

        const double CF0 = c[sl0a]*c[sl1a]*c[sl2a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sl2a] = dNdtByV[sl2a] - q0;

        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        
        const unsigned int i1 = this->nonEquilibriumReaction31index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D31NER[lhsIndex+3];
        const unsigned int sl1b = lhsSpeciesIndex1D31NER[lhsIndex+4];
        const unsigned int sl2b = lhsSpeciesIndex1D31NER[lhsIndex+5];

        const unsigned int sr0b = rhsSpeciesIndex1D31NER[rhsIndex+1];

        double Kf1 = Kf_[i1];

        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        double  Kr1 = this->Kf_[l1];

        const double CF1 = c[sl0b]*c[sl1b]*c[sl2b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sl2b] = dNdtByV[sl2b] - q1;

        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        
        
        const unsigned int i2 = this->nonEquilibriumReaction31index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D31NER[lhsIndex+6];
        const unsigned int sl1c = lhsSpeciesIndex1D31NER[lhsIndex+7];
        const unsigned int sl2c = lhsSpeciesIndex1D31NER[lhsIndex+8];

        const unsigned int sr0c = rhsSpeciesIndex1D31NER[rhsIndex+2];



        double Kf2 = Kf_[i2];

        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        double Kr2 = this->Kf_[l2];

        const double CF2 = c[sl0c]*c[sl1c]*c[sl2c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sl2c] = dNdtByV[sl2c] - q2;

        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;

        lhsIndex = lhsIndex + 9;
        rhsIndex = rhsIndex + 3;
    }
}