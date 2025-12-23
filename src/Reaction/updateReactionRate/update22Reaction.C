#include "OptReaction.H"

/*void 
FastChemistry::OptReaction::update22Reaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reaction22index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    for(std::size_t i=0; i<end-remain; i=i+4)
    {

        const unsigned int KfIndex0 = this->reaction22index[i+0];
        const unsigned int KfIndex1 = this->reaction22index[i+1];
        const unsigned int KfIndex2 = this->reaction22index[i+2];
        const unsigned int KfIndex3 = this->reaction22index[i+3];
        {
            const unsigned int sl0a = lhsSpeciesIndex1D22[rhsIndex+0];
            const unsigned int sl1a = lhsSpeciesIndex1D22[rhsIndex+1];
            
            const unsigned int sr0a = rhsSpeciesIndex1D22[rhsIndex+0];
            const unsigned int sr1a = rhsSpeciesIndex1D22[rhsIndex+1];

            double Kr0 = 0;
            double Kf0 = Kf_[KfIndex0];
            if(this->isIrreversible[KfIndex0]==0)
            {
                const double Kp = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
                double Kc = Kp;
                Kc = std::max(Kc,KcLimiter);
                Kr0 = Kf0/Kc;         
            }
            else if(this->isIrreversible[KfIndex0]==2)
            {
                unsigned int l = KfIndex0 - this->Ikf[1] + this->Ikf[9];
                Kr0 = this->Kf_[l];
            }

            const double CF0 = c[sl0a]*c[sl1a];
            const double CR0 = c[sr0a]*c[sr1a];
            const double q0 = (Kf0*CF0) - (Kr0*CR0);

            dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
            dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
            dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
            dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

            //rhsIndex = rhsIndex+2;

            const unsigned int sl0b = lhsSpeciesIndex1D22[rhsIndex+2];
            const unsigned int sl1b = lhsSpeciesIndex1D22[rhsIndex+3];
            
            const unsigned int sr0b = rhsSpeciesIndex1D22[rhsIndex+2];
            const unsigned int sr1b = rhsSpeciesIndex1D22[rhsIndex+3];

            double Kr1 = 0;
            double Kf1 = Kf_[KfIndex1];
            if(this->isIrreversible[KfIndex1]==0)
            {
                const double Kp = (ExpNegGbyRT[sr0b]*ExpNegGbyRT[sr1b])/(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]);
                double Kc = Kp;
                Kc = std::max(Kc,KcLimiter);
                Kr1 = Kf1/Kc;         
            }
            else if(this->isIrreversible[KfIndex1]==2)
            {
                unsigned int l = KfIndex1 - this->Ikf[1] + this->Ikf[9];
                Kr1 = this->Kf_[l];
            }

            const double CF1 = c[sl0b]*c[sl1b];
            const double CR1 = c[sr0b]*c[sr1b];
            const double q1 = (Kf1*CF1) - (Kr1*CR1);

            dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
            dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
            dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
            dNdtByV[sr1b] = dNdtByV[sr1b] + q1;


            const unsigned int sl0c = lhsSpeciesIndex1D22[rhsIndex+4];
            const unsigned int sl1c = lhsSpeciesIndex1D22[rhsIndex+5];
            
            const unsigned int sr0c = rhsSpeciesIndex1D22[rhsIndex+4];
            const unsigned int sr1c = rhsSpeciesIndex1D22[rhsIndex+5];

            double Kr2 = 0;
            double Kf2 = Kf_[KfIndex2];
            if(this->isIrreversible[KfIndex2]==0)
            {
                const double Kp = (ExpNegGbyRT[sr0c]*ExpNegGbyRT[sr1c])/(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]);
                double Kc = Kp;
                Kc = std::max(Kc,KcLimiter);
                Kr2 = Kf2/Kc;         
            }
            else if(this->isIrreversible[KfIndex2]==2)
            {
                unsigned int l = KfIndex2 - this->Ikf[1] + this->Ikf[9];
                Kr2 = this->Kf_[l];
            }

            const double CF2 = c[sl0c]*c[sl1c];
            const double CR2 = c[sr0c]*c[sr1c];
            const double q2 = (Kf2*CF2) - (Kr2*CR2);

            dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
            dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
            dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
            dNdtByV[sr1c] = dNdtByV[sr1c] + q2;



            const unsigned int sl0d = lhsSpeciesIndex1D22[rhsIndex+6];
            const unsigned int sl1d = lhsSpeciesIndex1D22[rhsIndex+7];
            
            const unsigned int sr0d = rhsSpeciesIndex1D22[rhsIndex+6];
            const unsigned int sr1d = rhsSpeciesIndex1D22[rhsIndex+7];

            double Kr3 = 0;
            double Kf3 = Kf_[KfIndex3];
            if(this->isIrreversible[KfIndex3]==0)
            {
                const double Kp = (ExpNegGbyRT[sr0d]*ExpNegGbyRT[sr1d])/(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]);
                double Kc = Kp;
                Kc = std::max(Kc,KcLimiter);
                Kr3 = Kf3/Kc;         
            }
            else if(this->isIrreversible[KfIndex3]==2)
            {
                unsigned int l = KfIndex3 - this->Ikf[1] + this->Ikf[9];
                Kr3 = this->Kf_[l];
            }

            const double CF3 = c[sl0d]*c[sl1d];
            const double CR3 = c[sr0d]*c[sr1d];
            const double q3 = (Kf3*CF3) - (Kr3*CR3);

            dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
            dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
            dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
            dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

            rhsIndex = rhsIndex+8;
        }
    }
    for(std::size_t i=end-remain; i<end; i=i+1)
    {
        const unsigned int KfIndex3 = this->reaction22index[i+0];
        {
            const unsigned int sl0d = lhsSpeciesIndex1D22[rhsIndex+0];
            const unsigned int sl1d = lhsSpeciesIndex1D22[rhsIndex+1];
            
            const unsigned int sr0d = rhsSpeciesIndex1D22[rhsIndex+0];
            const unsigned int sr1d = rhsSpeciesIndex1D22[rhsIndex+1];

            double Kr3 = 0;
            double Kf3 = Kf_[KfIndex3];
            if(this->isIrreversible[KfIndex3]==0)
            {
                const double Kp = (ExpNegGbyRT[sr0d]*ExpNegGbyRT[sr1d])/(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]);
                double Kc = Kp;
                Kc = std::max(Kc,KcLimiter);
                Kr3 = Kf3/Kc;         
            }
            else if(this->isIrreversible[KfIndex3]==2)
            {
                unsigned int l = KfIndex3 - this->Ikf[1] + this->Ikf[9];
                Kr3 = this->Kf_[l];
            }

            const double CF3 = c[sl0d]*c[sl1d];
            const double CR3 = c[sr0d]*c[sr1d];
            const double q3 = (Kf3*CF3) - (Kr3*CR3);

            dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
            dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
            dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
            dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

            rhsIndex = rhsIndex+2;
        }
    }
}*/


void 
FastChemistry::OptReaction::RF22RR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reversibleReaction22index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;

    for(std::size_t k=0; k<end-remain; k=k+4)
    {

        const unsigned int i0 = this->reversibleReaction22index[k+0];
        const unsigned int i1 = this->reversibleReaction22index[k+1];
        const unsigned int i2 = this->reversibleReaction22index[k+2];
        const unsigned int i3 = this->reversibleReaction22index[k+3];

        {
            const unsigned int sl0a = lhsSpeciesIndex1D22RR[rhsIndex+0];
            const unsigned int sl1a = lhsSpeciesIndex1D22RR[rhsIndex+1];
            
            const unsigned int sr0a = rhsSpeciesIndex1D22RR[rhsIndex+0];
            const unsigned int sr1a = rhsSpeciesIndex1D22RR[rhsIndex+1];

            double Kr0 = 0;
            double Kf0 = Kf_[i0];

            //double Kc0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
            //Kc0 = std::max(Kc0,KcLimiter);
            //Kr0 = Kf0/Kc0; 

            //double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])/(ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a]);
            double invKc0 = (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a])*(invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]);
            invKc0 = std::min(invKc0,invKcLimiter);
            Kr0 = Kf0*invKc0; 

            const double CF0 = c[sl0a]*c[sl1a];
            const double CR0 = c[sr0a]*c[sr1a];
            const double q0 = (Kf0*CF0) - (Kr0*CR0);

            dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
            dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
            dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
            dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        }



        {
            const unsigned int sl0b = lhsSpeciesIndex1D22RR[rhsIndex+2];
            const unsigned int sl1b = lhsSpeciesIndex1D22RR[rhsIndex+3];
            
            const unsigned int sr0b = rhsSpeciesIndex1D22RR[rhsIndex+2];
            const unsigned int sr1b = rhsSpeciesIndex1D22RR[rhsIndex+3];

            double Kr1 = 0;
            double Kf1 = Kf_[i1];

            //double Kc1 = (ExpNegGbyRT[sr0b]*ExpNegGbyRT[sr1b])/(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]);
            //Kc1 = std::max(Kc1,KcLimiter);
            //Kr1 = Kf1/Kc1;
            //double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b])/(ExpNegGbyRT[sr0b]*ExpNegGbyRT[sr1b]);
            double invKc1 = (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b])*(invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]);
            invKc1 = std::min(invKc1,invKcLimiter);
            Kr1 = Kf1*invKc1; 
            

            const double CF1 = c[sl0b]*c[sl1b];
            const double CR1 = c[sr0b]*c[sr1b];
            const double q1 = (Kf1*CF1) - (Kr1*CR1);

            dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
            dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
            dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
            dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        }
        {
            const unsigned int sl0c = lhsSpeciesIndex1D22RR[rhsIndex+4];
            const unsigned int sl1c = lhsSpeciesIndex1D22RR[rhsIndex+5];
            
            const unsigned int sr0c = rhsSpeciesIndex1D22RR[rhsIndex+4];
            const unsigned int sr1c = rhsSpeciesIndex1D22RR[rhsIndex+5];

            double Kr2 = 0;
            double Kf2 = Kf_[i2];

            //double Kc2 = (ExpNegGbyRT[sr0c]*ExpNegGbyRT[sr1c])/(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]);
            //Kc2 = std::max(Kc2,KcLimiter);
            //Kr2 = Kf2/Kc2;

            //double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c])/(ExpNegGbyRT[sr0c]*ExpNegGbyRT[sr1c]);
            double invKc2 = (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c])*(invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c]);
            invKc2 = std::min(invKc2,invKcLimiter);
            Kr2 = Kf2*invKc2; 
            
            const double CF2 = c[sl0c]*c[sl1c];
            const double CR2 = c[sr0c]*c[sr1c];
            const double q2 = (Kf2*CF2) - (Kr2*CR2);

            dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
            dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
            dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
            dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

        }
        {
            const unsigned int sl0d = lhsSpeciesIndex1D22RR[rhsIndex+6];
            const unsigned int sl1d = lhsSpeciesIndex1D22RR[rhsIndex+7];
            
            const unsigned int sr0d = rhsSpeciesIndex1D22RR[rhsIndex+6];
            const unsigned int sr1d = rhsSpeciesIndex1D22RR[rhsIndex+7];

            double Kr3 = 0;
            double Kf3 = Kf_[i3];
            
            //double Kc3 = (ExpNegGbyRT[sr0d]*ExpNegGbyRT[sr1d])/(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]);
            //Kc3 = std::max(Kc3,KcLimiter);
            //Kr3 = Kf3/Kc3;
            //double invKc3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d])/(ExpNegGbyRT[sr0d]*ExpNegGbyRT[sr1d]);
            double invKc3 = (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d])*(invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d]);
            invKc3 = std::min(invKc3,invKcLimiter);
            Kr3 = Kf3*invKc3; 

            const double CF3 = c[sl0d]*c[sl1d];
            const double CR3 = c[sr0d]*c[sr1d];
            const double q3 = (Kf3*CF3) - (Kr3*CR3);

            dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
            dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
            dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
            dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
        }


            rhsIndex = rhsIndex+8;
    }
    for(std::size_t k=end-remain; k<end; k=k+1)
    {
        const unsigned int i = this->reversibleReaction22index[k+0];
        {
            const unsigned int sl0 = lhsSpeciesIndex1D22RR[rhsIndex+0];
            const unsigned int sl1 = lhsSpeciesIndex1D22RR[rhsIndex+1];
            
            const unsigned int sr0 = rhsSpeciesIndex1D22RR[rhsIndex+0];
            const unsigned int sr1 = rhsSpeciesIndex1D22RR[rhsIndex+1];

            double Kr = 0;
            double Kf = Kf_[i];

            //double Kc = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
            //Kc = std::max(Kc,KcLimiter);
            //Kr = Kf/Kc;  

            //double invKc = (ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1])/(ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]);
            double invKc = (ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1])*(invNegGstdByRT[sr0]*invNegGstdByRT[sr1]);
            invKc = std::min(invKc,invKcLimiter);
            Kr = Kf*invKc;      

            const double CF = c[sl0]*c[sl1];
            const double CR = c[sr0]*c[sr1];
            const double q = (Kf*CF) - (Kr*CR);

            dNdtByV[sl0] = dNdtByV[sl0] - q; 
            dNdtByV[sl1] = dNdtByV[sl1] - q;
            dNdtByV[sr0] = dNdtByV[sr0] + q;
            dNdtByV[sr1] = dNdtByV[sr1] + q;

            rhsIndex = rhsIndex+2;
        }
    }
}

void 
FastChemistry::OptReaction::RF22IR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{

    std::size_t end = irreversibleReaction22index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    for(std::size_t i=0; i<end-remain; i=i+4)
    {

        const unsigned int KfIndex0 = this->irreversibleReaction22index[i+0];
        const unsigned int KfIndex1 = this->irreversibleReaction22index[i+1];
        const unsigned int KfIndex2 = this->irreversibleReaction22index[i+2];
        const unsigned int KfIndex3 = this->irreversibleReaction22index[i+3];
        {
            const unsigned int sl0a = lhsSpeciesIndex1D22IR[rhsIndex+0];
            const unsigned int sl1a = lhsSpeciesIndex1D22IR[rhsIndex+1];
            
            const unsigned int sr0a = rhsSpeciesIndex1D22IR[rhsIndex+0];
            const unsigned int sr1a = rhsSpeciesIndex1D22IR[rhsIndex+1];

            double Kr0 = 0;
            double Kf0 = Kf_[KfIndex0];

            const double CF0 = c[sl0a]*c[sl1a];
            const double CR0 = c[sr0a]*c[sr1a];
            const double q0 = (Kf0*CF0) - (Kr0*CR0);

            dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
            dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
            dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
            dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

            //rhsIndex = rhsIndex+2;

            const unsigned int sl0b = lhsSpeciesIndex1D22IR[rhsIndex+2];
            const unsigned int sl1b = lhsSpeciesIndex1D22IR[rhsIndex+3];
            
            const unsigned int sr0b = rhsSpeciesIndex1D22IR[rhsIndex+2];
            const unsigned int sr1b = rhsSpeciesIndex1D22IR[rhsIndex+3];

            double Kr1 = 0;
            double Kf1 = Kf_[KfIndex1];

            const double CF1 = c[sl0b]*c[sl1b];
            const double CR1 = c[sr0b]*c[sr1b];
            const double q1 = (Kf1*CF1) - (Kr1*CR1);

            dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
            dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
            dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
            dNdtByV[sr1b] = dNdtByV[sr1b] + q1;


            const unsigned int sl0c = lhsSpeciesIndex1D22IR[rhsIndex+4];
            const unsigned int sl1c = lhsSpeciesIndex1D22IR[rhsIndex+5];
            
            const unsigned int sr0c = rhsSpeciesIndex1D22IR[rhsIndex+4];
            const unsigned int sr1c = rhsSpeciesIndex1D22IR[rhsIndex+5];

            double Kr2 = 0;
            double Kf2 = Kf_[KfIndex2];


            const double CF2 = c[sl0c]*c[sl1c];
            const double CR2 = c[sr0c]*c[sr1c];
            const double q2 = (Kf2*CF2) - (Kr2*CR2);

            dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
            dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
            dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
            dNdtByV[sr1c] = dNdtByV[sr1c] + q2;



            const unsigned int sl0d = lhsSpeciesIndex1D22IR[rhsIndex+6];
            const unsigned int sl1d = lhsSpeciesIndex1D22IR[rhsIndex+7];
            
            const unsigned int sr0d = rhsSpeciesIndex1D22IR[rhsIndex+6];
            const unsigned int sr1d = rhsSpeciesIndex1D22IR[rhsIndex+7];

            double Kr3 = 0;
            double Kf3 = Kf_[KfIndex3];


            const double CF3 = c[sl0d]*c[sl1d];
            const double CR3 = c[sr0d]*c[sr1d];
            const double q3 = (Kf3*CF3) - (Kr3*CR3);

            dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
            dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
            dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
            dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

            rhsIndex = rhsIndex+8;
        }
    }
    for(std::size_t i=end-remain; i<end; i=i+1)
    {
        const unsigned int KfIndex3 = this->irreversibleReaction22index[i+0];
        {
            const unsigned int sl0d = lhsSpeciesIndex1D22IR[rhsIndex+0];
            const unsigned int sl1d = lhsSpeciesIndex1D22IR[rhsIndex+1];
            
            const unsigned int sr0d = rhsSpeciesIndex1D22IR[rhsIndex+0];
            const unsigned int sr1d = rhsSpeciesIndex1D22IR[rhsIndex+1];

            double Kr3 = 0;
            double Kf3 = Kf_[KfIndex3];

            const double CF3 = c[sl0d]*c[sl1d];
            const double CR3 = c[sr0d]*c[sr1d];
            const double q3 = (Kf3*CF3) - (Kr3*CR3);

            dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
            dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
            dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
            dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

            rhsIndex = rhsIndex+2;
        }
    }
}


void 
FastChemistry::OptReaction::RF22NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{

    std::size_t end = nonEquilibriumReaction22index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    for(std::size_t i=0; i<end-remain; i=i+4)
    {

        const unsigned int KfIndex0 = this->nonEquilibriumReaction22index[i+0];
        const unsigned int KfIndex1 = this->nonEquilibriumReaction22index[i+1];
        const unsigned int KfIndex2 = this->nonEquilibriumReaction22index[i+2];
        const unsigned int KfIndex3 = this->nonEquilibriumReaction22index[i+3];
        {
            const unsigned int sl0a = lhsSpeciesIndex1D22NER[rhsIndex+0];
            const unsigned int sl1a = lhsSpeciesIndex1D22NER[rhsIndex+1];
            
            const unsigned int sr0a = rhsSpeciesIndex1D22NER[rhsIndex+0];
            const unsigned int sr1a = rhsSpeciesIndex1D22NER[rhsIndex+1];

            double Kr0 = 0;
            double Kf0 = Kf_[KfIndex0];
            {
            unsigned int l = KfIndex0 - this->Ikf[1] + this->Ikf[9];
            Kr0 = this->Kf_[l];
            }
            

            const double CF0 = c[sl0a]*c[sl1a];
            const double CR0 = c[sr0a]*c[sr1a];
            const double q0 = (Kf0*CF0) - (Kr0*CR0);

            dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
            dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
            dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
            dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

            //rhsIndex = rhsIndex+2;

            const unsigned int sl0b = lhsSpeciesIndex1D22NER[rhsIndex+2];
            const unsigned int sl1b = lhsSpeciesIndex1D22NER[rhsIndex+3];
            
            const unsigned int sr0b = rhsSpeciesIndex1D22NER[rhsIndex+2];
            const unsigned int sr1b = rhsSpeciesIndex1D22NER[rhsIndex+3];

            double Kr1 = 0;
            double Kf1 = Kf_[KfIndex1];


            {
            unsigned int l = KfIndex1 - this->Ikf[1] + this->Ikf[9];
            Kr1 = this->Kf_[l];
            }

            

            const double CF1 = c[sl0b]*c[sl1b];
            const double CR1 = c[sr0b]*c[sr1b];
            const double q1 = (Kf1*CF1) - (Kr1*CR1);

            dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
            dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
            dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
            dNdtByV[sr1b] = dNdtByV[sr1b] + q1;


            const unsigned int sl0c = lhsSpeciesIndex1D22NER[rhsIndex+4];
            const unsigned int sl1c = lhsSpeciesIndex1D22NER[rhsIndex+5];
            
            const unsigned int sr0c = rhsSpeciesIndex1D22NER[rhsIndex+4];
            const unsigned int sr1c = rhsSpeciesIndex1D22NER[rhsIndex+5];

            double Kr2 = 0;
            double Kf2 = Kf_[KfIndex2];


            {
            unsigned int l = KfIndex2 - this->Ikf[1] + this->Ikf[9];
            Kr2 = this->Kf_[l];
            }
            

            const double CF2 = c[sl0c]*c[sl1c];
            const double CR2 = c[sr0c]*c[sr1c];
            const double q2 = (Kf2*CF2) - (Kr2*CR2);

            dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
            dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
            dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
            dNdtByV[sr1c] = dNdtByV[sr1c] + q2;



            const unsigned int sl0d = lhsSpeciesIndex1D22NER[rhsIndex+6];
            const unsigned int sl1d = lhsSpeciesIndex1D22NER[rhsIndex+7];
            
            const unsigned int sr0d = rhsSpeciesIndex1D22NER[rhsIndex+6];
            const unsigned int sr1d = rhsSpeciesIndex1D22NER[rhsIndex+7];

            double Kr3 = 0;
            double Kf3 = Kf_[KfIndex3];
            {
            unsigned int l = KfIndex3 - this->Ikf[1] + this->Ikf[9];
            Kr3 = this->Kf_[l];
            }


            const double CF3 = c[sl0d]*c[sl1d];
            const double CR3 = c[sr0d]*c[sr1d];
            const double q3 = (Kf3*CF3) - (Kr3*CR3);

            dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
            dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
            dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
            dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

            rhsIndex = rhsIndex+8;
        }
    }
    for(std::size_t i=end-remain; i<end; i=i+1)
    {
        const unsigned int KfIndex3 = this->nonEquilibriumReaction22index[i+0];
        {
            const unsigned int sl0d = lhsSpeciesIndex1D22NER[rhsIndex+0];
            const unsigned int sl1d = lhsSpeciesIndex1D22NER[rhsIndex+1];
            
            const unsigned int sr0d = rhsSpeciesIndex1D22NER[rhsIndex+0];
            const unsigned int sr1d = rhsSpeciesIndex1D22NER[rhsIndex+1];

            double Kr3 = 0;
            double Kf3 = Kf_[KfIndex3];


            
            unsigned int l = KfIndex3 - this->Ikf[1] + this->Ikf[9];
            Kr3 = this->Kf_[l];
            

            const double CF3 = c[sl0d]*c[sl1d];
            const double CR3 = c[sr0d]*c[sr1d];
            const double q3 = (Kf3*CF3) - (Kr3*CR3);

            dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
            dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
            dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
            dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

            rhsIndex = rhsIndex+2;
        }
    }
}