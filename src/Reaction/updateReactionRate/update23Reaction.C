#include "OptReaction.H"

/*void 
FastChemistry::OptReaction::update23Reaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reaction23index.size();
    //std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    std::size_t lhsIndex = 0;
    //std::vector<double> debugRR;
    //std::vector<double> debugIR;
    //std::vector<double> debugNER;
    for(std::size_t i=0; i<end; i++)
    {
        const unsigned int KfIndex = this->reaction23index[i];

        const unsigned int sl0 = lhsSpeciesIndex1D23[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D23[lhsIndex+1];

        const unsigned int sr0 = rhsSpeciesIndex1D23[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D23[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D23[rhsIndex+2];

        double Kr = 0;
        double Kf = Kf_[KfIndex];
        if(this->isIrreversible[KfIndex]==0)
        {
            const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
            double Kc = Kp*this->Pow_pByRT_SumVki[3];
            Kc = std::max(Kc,KcLimiter);
            Kr = Kf/Kc;         
        }
        else if(this->isIrreversible[KfIndex]==2)
        {
            unsigned int l = KfIndex - this->Ikf[1] + this->Ikf[9];
            Kr = this->Kf_[l];
        }

        const double CF = c[sl0]*c[sl1];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;

        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 3;

    }
}*/

void 
FastChemistry::OptReaction::RF23RR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reversibleReaction23index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    std::size_t lhsIndex = 0;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->reversibleReaction23index[k+0];
        const unsigned int i1 = this->reversibleReaction23index[k+1];
        const unsigned int i2 = this->reversibleReaction23index[k+2];
        const unsigned int i3 = this->reversibleReaction23index[k+3];

        const unsigned int sl0a = lhsSpeciesIndex1D23RR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D23RR[lhsIndex+1];

        const unsigned int sr0a = rhsSpeciesIndex1D23RR[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D23RR[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D23RR[rhsIndex+2];

        double Kr0 = 0;
        double Kf0 = Kf_[i0];

        
            //const double Kp0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a]*ExpNegGbyRT[sr2a])*
            //    (invNegGstdByRT[sl0a]*invNegGstdByRT[sl1a]);
            //double Kc0 = Kp0*this->Pow_pByRT_SumVki[3];
            //Kc0 = std::max(Kc0,KcLimiter);
            //Kr0 = Kf0/Kc0;

            const double invKp0 = (invNegGstdByRT[sr0a]*invNegGstdByRT[sr1a]*invNegGstdByRT[sr2a])*
                (ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
            double invKc0 = invKp0*this->Pow_pByRT_SumVki[1];
            invKc0 = std::min(invKc0, invKcLimiter);
            Kr0 = Kf0*invKc0;
        

        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        //lhsIndex = lhsIndex + 2;
        //rhsIndex = rhsIndex + 3;



        const unsigned int sl0b = lhsSpeciesIndex1D23RR[lhsIndex+2];
        const unsigned int sl1b = lhsSpeciesIndex1D23RR[lhsIndex+3];

        const unsigned int sr0b = rhsSpeciesIndex1D23RR[rhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D23RR[rhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D23RR[rhsIndex+5];

        double Kr1 = 0;
        double Kf1 = Kf_[i1];

        
            //const double Kp1 = (ExpNegGbyRT[sr0b]*ExpNegGbyRT[sr1b]*ExpNegGbyRT[sr2b])*
            //    (invNegGstdByRT[sl0b]*invNegGstdByRT[sl1b]);
            //double Kc1 = Kp1*this->Pow_pByRT_SumVki[3];
            //Kc1 = std::max(Kc1,KcLimiter);
            //Kr1 = Kf1/Kc1;
            const double invKp1 = (invNegGstdByRT[sr0b]*invNegGstdByRT[sr1b]*invNegGstdByRT[sr2b])*
                (ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]);
            double invKc1 = invKp1*this->Pow_pByRT_SumVki[1];
            invKc1 = std::min(invKc1,invKcLimiter);
            Kr1 = Kf1*invKc1;

        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;



        const unsigned int sl0c = lhsSpeciesIndex1D23RR[lhsIndex+4];
        const unsigned int sl1c = lhsSpeciesIndex1D23RR[lhsIndex+5];

        const unsigned int sr0c = rhsSpeciesIndex1D23RR[rhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D23RR[rhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D23RR[rhsIndex+8];

        double Kr2 = 0;
        double Kf2 = Kf_[i2];

        
            //const double Kp2 = (ExpNegGbyRT[sr0c]*ExpNegGbyRT[sr1c]*ExpNegGbyRT[sr2c])*
            //    (invNegGstdByRT[sl0c]*invNegGstdByRT[sl1c]);
            //double Kc2 = Kp2*this->Pow_pByRT_SumVki[3];
            //Kc2 = std::max(Kc2,KcLimiter);
            //Kr2 = Kf2/Kc2;
            const double invKp2 = (invNegGstdByRT[sr0c]*invNegGstdByRT[sr1c]*invNegGstdByRT[sr2c])*
                (ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]);
            double invKc2 = invKp2*this->Pow_pByRT_SumVki[1];
            invKc2 = std::min(invKc2,invKcLimiter);
            Kr2 = Kf2*invKc2;
        

        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;




        const unsigned int sl0d = lhsSpeciesIndex1D23RR[lhsIndex+6];
        const unsigned int sl1d = lhsSpeciesIndex1D23RR[lhsIndex+7];

        const unsigned int sr0d = rhsSpeciesIndex1D23RR[rhsIndex+9];
        const unsigned int sr1d = rhsSpeciesIndex1D23RR[rhsIndex+10];
        const unsigned int sr2d = rhsSpeciesIndex1D23RR[rhsIndex+11];

        double Kr3 = 0;
        double Kf3 = Kf_[i3];

        
            //const double Kp3 = (ExpNegGbyRT[sr0d]*ExpNegGbyRT[sr1d]*ExpNegGbyRT[sr2d])*
            //    (invNegGstdByRT[sl0d]*invNegGstdByRT[sl1d]);
            //double Kc3 = Kp3*this->Pow_pByRT_SumVki[3];
            //Kc3 = std::max(Kc3,KcLimiter);
            //Kr3 = Kf3/Kc3;
            const double invKp3 = (invNegGstdByRT[sr0d]*invNegGstdByRT[sr1d]*invNegGstdByRT[sr2d])*
                (ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]);
            double invKc3 = invKp3*this->Pow_pByRT_SumVki[1];
            invKc3 = std::min(invKc3,invKcLimiter);
            Kr3 = Kf3*invKc3;
        

        const double CF3 = c[sl0d]*c[sl1d];
        const double CR3 = c[sr0d]*c[sr1d]*c[sr2d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
        dNdtByV[sr2d] = dNdtByV[sr2d] + q3;



        lhsIndex = lhsIndex + 8;
        rhsIndex = rhsIndex + 12;
        //std::cout<<"RR "<<q0<<std::endl;
        //std::cout<<"RR "<<q1<<std::endl;
        //std::cout<<"RR "<<q2<<std::endl;
        //std::cout<<"RR "<<q3<<std::endl;
    }
    for(std::size_t k=end-remain; k<end; k=k+1)
    {
        const unsigned int i = this->reversibleReaction23index[k+0];

        const unsigned int sl0 = lhsSpeciesIndex1D23RR[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D23RR[lhsIndex+1];

        const unsigned int sr0 = rhsSpeciesIndex1D23RR[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D23RR[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D23RR[rhsIndex+2];

        double Kr = 0;
        double Kf = Kf_[i];

        
            //const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])*
            //    (invNegGstdByRT[sl0]*invNegGstdByRT[sl1]);
            //double Kc = Kp*this->Pow_pByRT_SumVki[3];
            //Kc = std::max(Kc,KcLimiter);
            //Kr = Kf/Kc;         
            const double invKp = (invNegGstdByRT[sr0]*invNegGstdByRT[sr1]*invNegGstdByRT[sr2])*
                (ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
            double invKc = invKp*this->Pow_pByRT_SumVki[1];
            invKc = std::min(invKc,invKcLimiter);
            Kr = Kf*invKc;    

        const double CF = c[sl0]*c[sl1];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 3;
        //std::cout<<"RR "<<q<<std::endl;

    }
}

void 
FastChemistry::OptReaction::RF23IR
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = irreversibleReaction23index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    std::size_t lhsIndex = 0;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->irreversibleReaction23index[k+0];
        const unsigned int i1 = this->irreversibleReaction23index[k+1];
        const unsigned int i2 = this->irreversibleReaction23index[k+2];
        const unsigned int i3 = this->irreversibleReaction23index[k+3];

        const unsigned int sl0a = lhsSpeciesIndex1D23IR[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D23IR[lhsIndex+1];

        const unsigned int sr0a = rhsSpeciesIndex1D23IR[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D23IR[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D23IR[rhsIndex+2];

        double Kf0 = Kf_[i0];

        const double CF0 = c[sl0a]*c[sl1a];
        const double q0 = (Kf0*CF0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;




        const unsigned int sl0b = lhsSpeciesIndex1D23IR[lhsIndex+2];
        const unsigned int sl1b = lhsSpeciesIndex1D23IR[lhsIndex+3];

        const unsigned int sr0b = rhsSpeciesIndex1D23IR[rhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D23IR[rhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D23IR[rhsIndex+5];

        double Kf1 = Kf_[i1];

        const double CF1 = c[sl0b]*c[sl1b];
        const double q1 = (Kf1*CF1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;



        const unsigned int sl0c = lhsSpeciesIndex1D23IR[lhsIndex+4];
        const unsigned int sl1c = lhsSpeciesIndex1D23IR[lhsIndex+5];

        const unsigned int sr0c = rhsSpeciesIndex1D23IR[rhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D23IR[rhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D23IR[rhsIndex+8];

        double Kf2 = Kf_[i2];

        const double CF2 = c[sl0c]*c[sl1c];
        const double q2 = (Kf2*CF2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;




        const unsigned int sl0d = lhsSpeciesIndex1D23IR[lhsIndex+6];
        const unsigned int sl1d = lhsSpeciesIndex1D23IR[lhsIndex+7];

        const unsigned int sr0d = rhsSpeciesIndex1D23IR[rhsIndex+9];
        const unsigned int sr1d = rhsSpeciesIndex1D23IR[rhsIndex+10];
        const unsigned int sr2d = rhsSpeciesIndex1D23IR[rhsIndex+11];

        double Kf3 = Kf_[i3];

        const double CF3 = c[sl0d]*c[sl1d];

        const double q3 = (Kf3*CF3) ;

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
        dNdtByV[sr2d] = dNdtByV[sr2d] + q3;



        lhsIndex = lhsIndex + 8;
        rhsIndex = rhsIndex + 12;
        //std::cout<<"IR "<<q0<<std::endl;
        //std::cout<<"IR "<<q1<<std::endl;
        //std::cout<<"IR "<<q2<<std::endl;
        //std::cout<<"IR "<<q3<<std::endl;
    }
    for(std::size_t k=end-remain; k<end; k=k+1)
    {
        const unsigned int i = this->irreversibleReaction23index[k+0];

        const unsigned int sl0 = lhsSpeciesIndex1D23IR[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D23IR[lhsIndex+1];

        const unsigned int sr0 = rhsSpeciesIndex1D23IR[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D23IR[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D23IR[rhsIndex+2];

        double Kf = Kf_[i];

        const double CF = c[sl0]*c[sl1];

        const double q = (Kf*CF) ;

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 3;
        //std::cout<<"IR "<<q<<std::endl;
    }
}
void 
FastChemistry::OptReaction::RF23NER
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = nonEquilibriumReaction23index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    std::size_t lhsIndex = 0;
    for(std::size_t k=0; k<end-remain; k=k+4)
    {
        const unsigned int i0 = this->nonEquilibriumReaction23index[k+0];
        const unsigned int i1 = this->nonEquilibriumReaction23index[k+1];
        const unsigned int i2 = this->nonEquilibriumReaction23index[k+2];
        const unsigned int i3 = this->nonEquilibriumReaction23index[k+3];

        const unsigned int sl0a = lhsSpeciesIndex1D23NER[lhsIndex+0];
        const unsigned int sl1a = lhsSpeciesIndex1D23NER[lhsIndex+1];

        const unsigned int sr0a = rhsSpeciesIndex1D23NER[rhsIndex+0];
        const unsigned int sr1a = rhsSpeciesIndex1D23NER[rhsIndex+1];
        const unsigned int sr2a = rhsSpeciesIndex1D23NER[rhsIndex+2];

        double Kr0 = 0;
        double Kf0 = Kf_[i0];

        {
            unsigned int l = i0 - this->Ikf[1] + this->Ikf[9];
            Kr0 = this->Kf_[l];
        }

        const double CF0 = c[sl0a]*c[sl1a];
        const double CR0 = c[sr0a]*c[sr1a]*c[sr2a];
        const double q0 = (Kf0*CF0) - (Kr0*CR0);

        dNdtByV[sl0a] = dNdtByV[sl0a] - q0;
        dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
        dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
        dNdtByV[sr1a] = dNdtByV[sr1a] + q0;
        dNdtByV[sr2a] = dNdtByV[sr2a] + q0;

        //lhsIndex = lhsIndex + 2;
        //rhsIndex = rhsIndex + 3;



        const unsigned int sl0b = lhsSpeciesIndex1D23NER[lhsIndex+2];
        const unsigned int sl1b = lhsSpeciesIndex1D23NER[lhsIndex+3];

        const unsigned int sr0b = rhsSpeciesIndex1D23NER[rhsIndex+3];
        const unsigned int sr1b = rhsSpeciesIndex1D23NER[rhsIndex+4];
        const unsigned int sr2b = rhsSpeciesIndex1D23NER[rhsIndex+5];

        double Kr1 = 0;
        double Kf1 = Kf_[i1];

        {
            unsigned int l = i1 - this->Ikf[1] + this->Ikf[9];
            Kr1 = this->Kf_[l];    
        }

        const double CF1 = c[sl0b]*c[sl1b];
        const double CR1 = c[sr0b]*c[sr1b]*c[sr2b];
        const double q1 = (Kf1*CF1) - (Kr1*CR1);

        dNdtByV[sl0b] = dNdtByV[sl0b] - q1;
        dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
        dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
        dNdtByV[sr1b] = dNdtByV[sr1b] + q1;
        dNdtByV[sr2b] = dNdtByV[sr2b] + q1;



        const unsigned int sl0c = lhsSpeciesIndex1D23NER[lhsIndex+4];
        const unsigned int sl1c = lhsSpeciesIndex1D23NER[lhsIndex+5];

        const unsigned int sr0c = rhsSpeciesIndex1D23NER[rhsIndex+6];
        const unsigned int sr1c = rhsSpeciesIndex1D23NER[rhsIndex+7];
        const unsigned int sr2c = rhsSpeciesIndex1D23NER[rhsIndex+8];

        double Kr2 = 0;
        double Kf2 = Kf_[i2];

        {
            unsigned int l = i2 - this->Ikf[1] + this->Ikf[9];
            Kr2 = this->Kf_[l];
        }

        const double CF2 = c[sl0c]*c[sl1c];
        const double CR2 = c[sr0c]*c[sr1c]*c[sr2c];
        const double q2 = (Kf2*CF2) - (Kr2*CR2);

        dNdtByV[sl0c] = dNdtByV[sl0c] - q2;
        dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
        dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
        dNdtByV[sr1c] = dNdtByV[sr1c] + q2;
        dNdtByV[sr2c] = dNdtByV[sr2c] + q2;




        const unsigned int sl0d = lhsSpeciesIndex1D23NER[lhsIndex+6];
        const unsigned int sl1d = lhsSpeciesIndex1D23NER[lhsIndex+7];

        const unsigned int sr0d = rhsSpeciesIndex1D23NER[rhsIndex+9];
        const unsigned int sr1d = rhsSpeciesIndex1D23NER[rhsIndex+10];
        const unsigned int sr2d = rhsSpeciesIndex1D23NER[rhsIndex+11];

        double Kr3 = 0;
        double Kf3 = Kf_[i3];

        {
            unsigned int l = i3 - this->Ikf[1] + this->Ikf[9];
            Kr3 = this->Kf_[l];
        }

        const double CF3 = c[sl0d]*c[sl1d];
        const double CR3 = c[sr0d]*c[sr1d]*c[sr2d];
        const double q3 = (Kf3*CF3) - (Kr3*CR3);

        dNdtByV[sl0d] = dNdtByV[sl0d] - q3;
        dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
        dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
        dNdtByV[sr1d] = dNdtByV[sr1d] + q3;
        dNdtByV[sr2d] = dNdtByV[sr2d] + q3;


        lhsIndex = lhsIndex + 8;
        rhsIndex = rhsIndex + 12;
        //std::cout<<q0<<std::endl;
        //std::cout<<q1<<std::endl;
        //std::cout<<q2<<std::endl;
        //std::cout<<q3<<std::endl;
    }
    for(std::size_t k=end-remain; k<end; k=k+1)
    {
        const unsigned int i = this->nonEquilibriumReaction23index[k+0];

        const unsigned int sl0 = lhsSpeciesIndex1D23NER[lhsIndex+0];
        const unsigned int sl1 = lhsSpeciesIndex1D23NER[lhsIndex+1];

        const unsigned int sr0 = rhsSpeciesIndex1D23NER[rhsIndex+0];
        const unsigned int sr1 = rhsSpeciesIndex1D23NER[rhsIndex+1];
        const unsigned int sr2 = rhsSpeciesIndex1D23NER[rhsIndex+2];

        double Kr = 0;
        double Kf = Kf_[i];

        {
            unsigned int l = i - this->Ikf[1] + this->Ikf[9];
            Kr = this->Kf_[l];     
        }

        const double CF = c[sl0]*c[sl1];
        const double CR = c[sr0]*c[sr1]*c[sr2];
        const double q = (Kf*CF) - (Kr*CR);

        dNdtByV[sl0] = dNdtByV[sl0] - q;
        dNdtByV[sl1] = dNdtByV[sl1] - q;
        dNdtByV[sr0] = dNdtByV[sr0] + q;
        dNdtByV[sr1] = dNdtByV[sr1] + q;
        dNdtByV[sr2] = dNdtByV[sr2] + q;
        lhsIndex = lhsIndex + 2;
        rhsIndex = rhsIndex + 3;
        //std::cout<<q<<std::endl;
    }
}