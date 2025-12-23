#include "OptReaction.H"

/*void 
FastChemistry::OptReaction::update11Reaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reaction11index.size();
    //std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    for(std::size_t i=0; i<end; i++)
    {
        const unsigned int KfIndex = this->reaction11index[i];

        const unsigned int sl0 = lhsSpeciesIndex1D11[rhsIndex+0];
        
        const unsigned int sr0 = rhsSpeciesIndex1D11[rhsIndex+0];


        double Kr = 0;
        double Kf = Kf_[KfIndex];
        if(this->isIrreversible[KfIndex]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0])/(ExpNegGbyRT[sl0]);
            double Kc = Kp;
            Kc = std::max(Kc,KcLimiter);
            Kr = Kf/Kc;         
        }
        else if(this->isIrreversible[KfIndex]==2)
        {
            unsigned int l = KfIndex - this->Ikf[1] + this->Ikf[9];
            Kr = this->Kf_[l];
        }

        const double CF = c[sl0];
        const double CR = c[sr0];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q; 
        dNdtByV[sr0] = dNdtByV[sr0] + q;


        rhsIndex = rhsIndex+1;
    }
}*/

void 
FastChemistry::OptReaction::RF11RR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reversibleReaction11index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        
        const unsigned int i0 = this->reversibleReaction11index[k+0];

        {
        const unsigned int sl0a = lhsSpeciesIndex1D11RR[rhsIndex+0];
        
        const unsigned int sr0a = rhsSpeciesIndex1D11RR[rhsIndex+0];

        double Kf0 = Kf_[i0];

        //const double Kp0 = (ExpNegGbyRT[sr0a])*(invNegGstdByRT[sl0a]);
        //double Kc0 = Kp0;
        //Kc0 = std::max(Kc0,KcLimiter);
        //double Kr0 = Kf0/Kc0;

        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0;
        invKc0 = std::min(invKc0,invKcLimiter);
        double Kr0 = Kf0*invKc0;

        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        }
        {
        const unsigned int i1 = this->reversibleReaction11index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D11RR[rhsIndex+1];
        
        const unsigned int sr0b = rhsSpeciesIndex1D11RR[rhsIndex+1];

        double Kf1 = Kf_[i1];

        //const double Kp1 = (ExpNegGbyRT[sr0b])*(invNegGstdByRT[sl0b]);
        //double Kc1 = Kp1;
        //Kc1 = std::max(Kc1,KcLimiter);
        //double Kr1 = Kf1/Kc1;

        const double invKp1 = (invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1;
        invKc1 = std::min(invKc1,invKcLimiter);
        double Kr1 = Kf1*invKc1;

        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        }

        
        
        const unsigned int i2 = this->reversibleReaction11index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D11RR[rhsIndex+2];
        
        const unsigned int sr0c = rhsSpeciesIndex1D11RR[rhsIndex+2];

        double Kf2 = Kf_[i2];

        //const double Kp2 = (ExpNegGbyRT[sr0c])*(invNegGstdByRT[sl0c]);
        //double Kc2 = Kp2;
        //Kc2 = std::max(Kc2,KcLimiter);
        //double Kr2 = Kf2/Kc2;
        const double invKp2 = (invNegGstdByRT[sr0c])*(ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2;
        invKc2 = std::min(invKc2,invKcLimiter);
        double Kr2 = Kf2*invKc2;

        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        
        
        const unsigned int i3 = this->reversibleReaction11index[k+3];

        const unsigned int sl0d = lhsSpeciesIndex1D11RR[rhsIndex+3];
        
        const unsigned int sr0d = rhsSpeciesIndex1D11RR[rhsIndex+3];

        double Kf3 = Kf_[i3];
        const double Kp3 = (ExpNegGbyRT[sr0d])*(invNegGstdByRT[sl0d]);
        double Kc3 = Kp3;
        Kc3 = std::max(Kc3,KcLimiter);
        double Kr3 = Kf3/Kc3;

        const double CF3 = c[sl0d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        

        rhsIndex = rhsIndex+4;
    }
    if(remain==1)
    {
        std::size_t k = end-1;
        const unsigned int i0 = this->reversibleReaction11index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D11RR[rhsIndex+0];
        
        const unsigned int sr0a = rhsSpeciesIndex1D11RR[rhsIndex+0];

        double Kf0 = Kf_[i0];

        //const double Kp0 = (ExpNegGbyRT[sr0a])*(invNegGstdByRT[sl0a]);
        //double Kc0 = Kp0;
        //Kc0 = std::max(Kc0,KcLimiter);
        //double Kr0 = Kf0/Kc0;
        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0;
        invKc0 = std::min(invKc0,invKcLimiter);
        double Kr0 = Kf0*invKc0;


        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        
        rhsIndex = rhsIndex+1;
    }
    else if(remain==2)
    {
        std::size_t k = end-2;
        const unsigned int i0 = this->reversibleReaction11index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D11RR[rhsIndex+0];
        
        const unsigned int sr0a = rhsSpeciesIndex1D11RR[rhsIndex+0];

        double Kf0 = Kf_[i0];

        //const double Kp0 = (ExpNegGbyRT[sr0a])*(invNegGstdByRT[sl0a]);
        //double Kc0 = Kp0;
        //Kc0 = std::max(Kc0,KcLimiter);
        //double Kr0 = Kf0/Kc0;
        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0;
        invKc0 = std::min(invKc0,invKcLimiter);
        double Kr0 = Kf0*invKc0;

        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        
        const unsigned int i1 = this->reversibleReaction11index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D11RR[rhsIndex+1];
        
        const unsigned int sr0b = rhsSpeciesIndex1D11RR[rhsIndex+1];

        double Kf1 = Kf_[i1];

        //const double Kp1 = (ExpNegGbyRT[sr0b])*(invNegGstdByRT[sl0b]);
        //double Kc1 = Kp1;
        //Kc1 = std::max(Kc1,KcLimiter);
        //double Kr1 = Kf1/Kc1;
        const double invKp1 = (invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1;
        invKc1 = std::min(invKc1,invKcLimiter);
        double Kr1 = Kf1*invKc1;

        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        

        rhsIndex = rhsIndex+2;
    }
    else if(remain==3)
    {
        std::size_t k = end-3;
        const unsigned int i0 = this->reversibleReaction11index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D11RR[rhsIndex+0];
        
        const unsigned int sr0a = rhsSpeciesIndex1D11RR[rhsIndex+0];

        double Kf0 = Kf_[i0];

        //const double Kp0 = (ExpNegGbyRT[sr0a])*(invNegGstdByRT[sl0a]);
        //double Kc0 = Kp0;
        //Kc0 = std::max(Kc0,KcLimiter);
        //double Kr0 = Kf0/Kc0;
        const double invKp0 = (invNegGstdByRT[sr0a])*(ExpNegGbyRT[sl0a]);
        double invKc0 = invKp0;
        invKc0 = std::min(invKc0,invKcLimiter);
        double Kr0 = Kf0*invKc0;

        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        
        const unsigned int i1 = this->reversibleReaction11index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D11RR[rhsIndex+1];
        
        const unsigned int sr0b = rhsSpeciesIndex1D11RR[rhsIndex+1];

        double Kf1 = Kf_[i1];

        //const double Kp1 = (ExpNegGbyRT[sr0b])*(invNegGstdByRT[sl0b]);
        //double Kc1 = Kp1;
        //Kc1 = std::max(Kc1,KcLimiter);
        //double Kr1 = Kf1/Kc1;
        const double invKp1 = (invNegGstdByRT[sr0b])*(ExpNegGbyRT[sl0b]);
        double invKc1 = invKp1;
        invKc1 = std::min(invKc1,invKcLimiter);
        double Kr1 = Kf1*invKc1;

        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        
        
        const unsigned int i2 = this->reversibleReaction11index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D11RR[rhsIndex+2];
        
        const unsigned int sr0c = rhsSpeciesIndex1D11RR[rhsIndex+2];

        double Kf2 = Kf_[i2];

        //const double Kp2 = (ExpNegGbyRT[sr0c])*(invNegGstdByRT[sl0c]);
        //double Kc2 = Kp2;
        //Kc2 = std::max(Kc2,KcLimiter);
        //double Kr2 = Kf2/Kc2;
        const double invKp2 = (invNegGstdByRT[sr0c])*(ExpNegGbyRT[sl0c]);
        double invKc2 = invKp2;
        invKc2 = std::min(invKc2,invKcLimiter);
        double Kr2 = Kf2*invKc2;

        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;

        rhsIndex = rhsIndex+3;
    }
}

void 
FastChemistry::OptReaction::RF11IR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = irreversibleReaction11index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->irreversibleReaction11index[k];

        const unsigned int sl0a = lhsSpeciesIndex1D11IR[rhsIndex+0];
        
        const unsigned int sr0a = rhsSpeciesIndex1D11IR[rhsIndex+0];

        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        
        const unsigned int i1 = this->irreversibleReaction11index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D11IR[rhsIndex+1];
        
        const unsigned int sr0b = rhsSpeciesIndex1D11IR[rhsIndex+1];

        double Kf1 = Kf_[i1];

        const double CF1 = c[sl0b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        
        
        const unsigned int i2 = this->irreversibleReaction11index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D11IR[rhsIndex+2];
        
        const unsigned int sr0c = rhsSpeciesIndex1D11IR[rhsIndex+2];

        double Kf2 = Kf_[i2];

        const double CF2 = c[sl0c];

        const double q2 = (Kf2*CF2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        
        
        const unsigned int i3 = this->irreversibleReaction11index[k+3];

        const unsigned int sl0d = lhsSpeciesIndex1D11IR[rhsIndex+3];
        
        const unsigned int sr0d = rhsSpeciesIndex1D11IR[rhsIndex+3];

        double Kf3 = Kf_[i3];

        const double CF3 = c[sl0d];

        const double q3 = (Kf3*CF3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        
        rhsIndex = rhsIndex+4;
    }
    if(remain==1)
    {
        std::size_t k = end-1;
        const unsigned int i0 = this->irreversibleReaction11index[k];

        const unsigned int sl0a = lhsSpeciesIndex1D11IR[rhsIndex+0];
        
        const unsigned int sr0a = rhsSpeciesIndex1D11IR[rhsIndex+0];

        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        rhsIndex = rhsIndex+1;
    }
    else if(remain==2)
    {
        std::size_t k = end-2;
        const unsigned int i0 = this->irreversibleReaction11index[k];

        const unsigned int sl0a = lhsSpeciesIndex1D11IR[rhsIndex+0];
        
        const unsigned int sr0a = rhsSpeciesIndex1D11IR[rhsIndex+0];

        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        
        const unsigned int i1 = this->irreversibleReaction11index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D11IR[rhsIndex+1];
        
        const unsigned int sr0b = rhsSpeciesIndex1D11IR[rhsIndex+1];

        double Kf1 = Kf_[i1];

        const double CF1 = c[sl0b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;

        rhsIndex = rhsIndex+4;
    }
    else if(remain==3)
    {
        std::size_t k = end-3;
        const unsigned int i0 = this->irreversibleReaction11index[k];

        const unsigned int sl0a = lhsSpeciesIndex1D11IR[rhsIndex+0];
        
        const unsigned int sr0a = rhsSpeciesIndex1D11IR[rhsIndex+0];

        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a];

        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;

        
        const unsigned int i1 = this->irreversibleReaction11index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D11IR[rhsIndex+1];
        
        const unsigned int sr0b = rhsSpeciesIndex1D11IR[rhsIndex+1];

        double Kf1 = Kf_[i1];

        const double CF1 = c[sl0b];

        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        
        
        const unsigned int i2 = this->irreversibleReaction11index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D11IR[rhsIndex+2];
        
        const unsigned int sr0c = rhsSpeciesIndex1D11IR[rhsIndex+2];

        double Kf2 = Kf_[i2];

        const double CF2 = c[sl0c];

        const double q2 = (Kf2*CF2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        
        rhsIndex = rhsIndex+4;
    }
}

void 
FastChemistry::OptReaction::RF11NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = nonEquilibriumReaction11index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->nonEquilibriumReaction11index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D11NER[rhsIndex+0];
        
        const unsigned int sr0a = rhsSpeciesIndex1D11NER[rhsIndex+0];

        double Kf0 = Kf_[i0];

        
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        double Kr0 = this->Kf_[l0];
        

        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        
        const unsigned int i1 = this->nonEquilibriumReaction11index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D11NER[rhsIndex+1];
        
        const unsigned int sr0b = rhsSpeciesIndex1D11NER[rhsIndex+1];

        double Kf1 = Kf_[i1];

        
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        double Kr1 = this->Kf_[l1];
        

        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        
        
        const unsigned int i2 = this->nonEquilibriumReaction11index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D11NER[rhsIndex+2];
        
        const unsigned int sr0c = rhsSpeciesIndex1D11NER[rhsIndex+2];

        double Kf2 = Kf_[i2];

        
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        double Kr2 = this->Kf_[l2];
        

        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        
        
        const unsigned int i3 = this->nonEquilibriumReaction11index[k+3];

        const unsigned int sl0d = lhsSpeciesIndex1D11NER[rhsIndex+3];
        
        const unsigned int sr0d = rhsSpeciesIndex1D11NER[rhsIndex+3];

        double Kf3 = Kf_[i3];

        
        unsigned int l3 = i3 - this->Ikf[1] + this->Ikf[9];
        double Kr3 = this->Kf_[l3];
        

        const double CF3 = c[sl0d];
        const double CR3 = c[sr0d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        

        rhsIndex = rhsIndex+4;
    }
    if(remain==1)
    {
        std::size_t k = end-1;
        const unsigned int i0 = this->nonEquilibriumReaction11index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D11NER[rhsIndex+0];
        
        const unsigned int sr0a = rhsSpeciesIndex1D11NER[rhsIndex+0];

        double Kf0 = Kf_[i0];

        
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        double Kr0 = this->Kf_[l0];
        

        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        
        rhsIndex = rhsIndex+1;
    }
    else if(remain==2)
    {
        std::size_t k = end-2;
        const unsigned int i0 = this->nonEquilibriumReaction11index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D11NER[rhsIndex+0];
        
        const unsigned int sr0a = rhsSpeciesIndex1D11NER[rhsIndex+0];

        double Kf0 = Kf_[i0];

        
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        double Kr0 = this->Kf_[l0];
        

        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        
        const unsigned int i1 = this->nonEquilibriumReaction11index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D11NER[rhsIndex+1];
        
        const unsigned int sr0b = rhsSpeciesIndex1D11NER[rhsIndex+1];

        double Kf1 = Kf_[i1];

        
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        double Kr1 = this->Kf_[l1];
        

        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;

        rhsIndex = rhsIndex+2;
    }
    else if(remain==3)
    {
        std::size_t k = end-3;
        const unsigned int i0 = this->nonEquilibriumReaction11index[k+0];

        const unsigned int sl0a = lhsSpeciesIndex1D11NER[rhsIndex+0];
        
        const unsigned int sr0a = rhsSpeciesIndex1D11NER[rhsIndex+0];

        double Kf0 = Kf_[i0];

        
        unsigned int l0 = i0 - this->Ikf[1] + this->Ikf[9];
        double Kr0 = this->Kf_[l0];
        

        const double CF0 = c[sl0a];
        const double CR0 = c[sr0a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        
        const unsigned int i1 = this->nonEquilibriumReaction11index[k+1];

        const unsigned int sl0b = lhsSpeciesIndex1D11NER[rhsIndex+1];
        
        const unsigned int sr0b = rhsSpeciesIndex1D11NER[rhsIndex+1];

        double Kf1 = Kf_[i1];

        
        unsigned int l1 = i1 - this->Ikf[1] + this->Ikf[9];
        double Kr1 = this->Kf_[l1];
        

        const double CF1 = c[sl0b];
        const double CR1 = c[sr0b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        
        
        const unsigned int i2 = this->nonEquilibriumReaction11index[k+2];

        const unsigned int sl0c = lhsSpeciesIndex1D11NER[rhsIndex+2];
        
        const unsigned int sr0c = rhsSpeciesIndex1D11NER[rhsIndex+2];

        double Kf2 = Kf_[i2];

        
        unsigned int l2 = i2 - this->Ikf[1] + this->Ikf[9];
        double Kr2 = this->Kf_[l2];
        

        const double CF2 = c[sl0c];
        const double CR2 = c[sr0c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        
    
        rhsIndex = rhsIndex+3;
    }
}