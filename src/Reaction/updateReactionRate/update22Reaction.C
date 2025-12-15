#include "OptReaction.H"

void 
OptReaction::update22Reaction
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
}


void 
OptReaction::update22ReversibleReaction
(
    const double* __restrict__ c,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    std::size_t end = reversibleReaction22index.size();
    std::size_t remain = end%4;
    std::size_t rhsIndex = 0;
    __m256d KcLimiterv = _mm256_set1_pd(KcLimiter);
    for(std::size_t i=0; i<end-remain; i=i+4)
    {

        const unsigned int KfIndex0 = this->reversibleReaction22index[i+0];
        const unsigned int KfIndex1 = this->reversibleReaction22index[i+1];
        const unsigned int KfIndex2 = this->reversibleReaction22index[i+2];
        const unsigned int KfIndex3 = this->reversibleReaction22index[i+3];
        const __m256d Kfv = _mm256_setr_pd(Kf_[KfIndex0],Kf_[KfIndex1],Kf_[KfIndex2],Kf_[KfIndex3]);
        {

            /*const unsigned int sl0a = lhsSpeciesIndex1D22RR[rhsIndex+0];
            const unsigned int sl0b = lhsSpeciesIndex1D22RR[rhsIndex+2];
            const unsigned int sl0c = lhsSpeciesIndex1D22RR[rhsIndex+4];
            const unsigned int sl0d = lhsSpeciesIndex1D22RR[rhsIndex+6];
            const __m256d ExpNegGbyRTsl0v = _mm256_setr_pd(ExpNegGbyRT[sl0a],ExpNegGbyRT[sl0b],ExpNegGbyRT[sl0c],ExpNegGbyRT[sl0d]);
            const __m256d Csl0v = _mm256_setr_pd(c[sl0a],c[sl0b],c[sl0c],c[sl0d]);

            const unsigned int sl1a = lhsSpeciesIndex1D22RR[rhsIndex+1];
            const unsigned int sl1b = lhsSpeciesIndex1D22RR[rhsIndex+3];
            const unsigned int sl1c = lhsSpeciesIndex1D22RR[rhsIndex+5];
            const unsigned int sl1d = lhsSpeciesIndex1D22RR[rhsIndex+7];
            const __m256d ExpNegGbyRTsl1v = _mm256_setr_pd(ExpNegGbyRT[sl1a],ExpNegGbyRT[sl1b],ExpNegGbyRT[sl1c],ExpNegGbyRT[sl1d]);
            const __m256d Csl1v = _mm256_setr_pd(c[sl1a],c[sl1b],c[sl1c],c[sl1d]);

            const unsigned int sr0a = rhsSpeciesIndex1D22RR[rhsIndex+0];
            const unsigned int sr0b = rhsSpeciesIndex1D22RR[rhsIndex+2];
            const unsigned int sr0c = rhsSpeciesIndex1D22RR[rhsIndex+4];
            const unsigned int sr0d = rhsSpeciesIndex1D22RR[rhsIndex+6];
            const __m256d ExpNegGbyRTsr0v = _mm256_setr_pd(ExpNegGbyRT[sr0a],ExpNegGbyRT[sr0b],ExpNegGbyRT[sr0c],ExpNegGbyRT[sr0d]);
            const __m256d Csr0v = _mm256_setr_pd(c[sr0a],c[sr0b],c[sr0c],c[sr0d]);

            const unsigned int sr1a = rhsSpeciesIndex1D22RR[rhsIndex+1];
            const unsigned int sr1b = rhsSpeciesIndex1D22RR[rhsIndex+3];
            const unsigned int sr1c = rhsSpeciesIndex1D22RR[rhsIndex+5];
            const unsigned int sr1d = rhsSpeciesIndex1D22RR[rhsIndex+7];
            const __m256d ExpNegGbyRTsr1v = _mm256_setr_pd(ExpNegGbyRT[sr1a],ExpNegGbyRT[sr1b],ExpNegGbyRT[sr1c],ExpNegGbyRT[sr1d]);
            const __m256d Csr1v = _mm256_setr_pd(c[sr1a],c[sr1b],c[sr1c],c[sr1d]);

            const __m256d Kpv = _mm256_div_pd(_mm256_mul_pd(ExpNegGbyRTsr0v,ExpNegGbyRTsr1v),_mm256_mul_pd(ExpNegGbyRTsl0v,ExpNegGbyRTsl1v));
            const __m256d Kcv = _mm256_max_pd(Kpv,KcLimiterv);
            const __m256d Krv = _mm256_div_pd(Kfv,Kcv);
            
            const __m256d CFv = _mm256_mul_pd(Csl0v,Csl1v);
            const __m256d CRv = _mm256_mul_pd(Csr0v,Csr1v);

            const __m256d qv = _mm256_sub_pd(_mm256_mul_pd(Kfv,CFv),_mm256_mul_pd(Krv,CRv));

            const double q0 = this->get_elem0(qv);
            dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
            dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
            dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
            dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

            const double q1 = this->get_elem1(qv);
            dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
            dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
            dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
            dNdtByV[sr1b] = dNdtByV[sr1b] + q1;

            const double q2 = this->get_elem2(qv);
            dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
            dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
            dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
            dNdtByV[sr1c] = dNdtByV[sr1c] + q2;

            const double q3 = this->get_elem3(qv);
            dNdtByV[sl0d] = dNdtByV[sl0d] - q3; 
            dNdtByV[sl1d] = dNdtByV[sl1d] - q3;
            dNdtByV[sr0d] = dNdtByV[sr0d] + q3;
            dNdtByV[sr1d] = dNdtByV[sr1d] + q3;

            rhsIndex = rhsIndex+8;*/

            const unsigned int sl0a = lhsSpeciesIndex1D22RR[rhsIndex+0];
            const unsigned int sl1a = lhsSpeciesIndex1D22RR[rhsIndex+1];
            
            const unsigned int sr0a = rhsSpeciesIndex1D22RR[rhsIndex+0];
            const unsigned int sr1a = rhsSpeciesIndex1D22RR[rhsIndex+1];

            double Kr0 = 0;
            double Kf0 = Kf_[KfIndex0];

            
            const double Kp0 = (ExpNegGbyRT[sr0a]*ExpNegGbyRT[sr1a])/(ExpNegGbyRT[sl0a]*ExpNegGbyRT[sl1a]);
            double Kc0 = Kp0;
            Kc0 = std::max(Kc0,KcLimiter);
            Kr0 = Kf0/Kc0;         
            


            const double CF0 = c[sl0a]*c[sl1a];
            const double CR0 = c[sr0a]*c[sr1a];
            const double q0 = (Kf0*CF0) - (Kr0*CR0);

            dNdtByV[sl0a] = dNdtByV[sl0a] - q0; 
            dNdtByV[sl1a] = dNdtByV[sl1a] - q0;
            dNdtByV[sr0a] = dNdtByV[sr0a] + q0;
            dNdtByV[sr1a] = dNdtByV[sr1a] + q0;

            //rhsIndex = rhsIndex+2;

            const unsigned int sl0b = lhsSpeciesIndex1D22RR[rhsIndex+2];
            const unsigned int sl1b = lhsSpeciesIndex1D22RR[rhsIndex+3];
            
            const unsigned int sr0b = rhsSpeciesIndex1D22RR[rhsIndex+2];
            const unsigned int sr1b = rhsSpeciesIndex1D22RR[rhsIndex+3];

            double Kr1 = 0;
            double Kf1 = Kf_[KfIndex1];

            
            const double Kp1 = (ExpNegGbyRT[sr0b]*ExpNegGbyRT[sr1b])/(ExpNegGbyRT[sl0b]*ExpNegGbyRT[sl1b]);
            double Kc1 = Kp1;
            Kc1 = std::max(Kc1,KcLimiter);
            Kr1 = Kf1/Kc1;         
            


            const double CF1 = c[sl0b]*c[sl1b];
            const double CR1 = c[sr0b]*c[sr1b];
            const double q1 = (Kf1*CF1) - (Kr1*CR1);

            dNdtByV[sl0b] = dNdtByV[sl0b] - q1; 
            dNdtByV[sl1b] = dNdtByV[sl1b] - q1;
            dNdtByV[sr0b] = dNdtByV[sr0b] + q1;
            dNdtByV[sr1b] = dNdtByV[sr1b] + q1;


            const unsigned int sl0c = lhsSpeciesIndex1D22RR[rhsIndex+4];
            const unsigned int sl1c = lhsSpeciesIndex1D22RR[rhsIndex+5];
            
            const unsigned int sr0c = rhsSpeciesIndex1D22RR[rhsIndex+4];
            const unsigned int sr1c = rhsSpeciesIndex1D22RR[rhsIndex+5];

            double Kr2 = 0;
            double Kf2 = Kf_[KfIndex2];

            
            const double Kp2 = (ExpNegGbyRT[sr0c]*ExpNegGbyRT[sr1c])/(ExpNegGbyRT[sl0c]*ExpNegGbyRT[sl1c]);
            double Kc2 = Kp2;
            Kc2 = std::max(Kc2,KcLimiter);
            Kr2 = Kf2/Kc2;         
            


            const double CF2 = c[sl0c]*c[sl1c];
            const double CR2 = c[sr0c]*c[sr1c];
            const double q2 = (Kf2*CF2) - (Kr2*CR2);

            dNdtByV[sl0c] = dNdtByV[sl0c] - q2; 
            dNdtByV[sl1c] = dNdtByV[sl1c] - q2;
            dNdtByV[sr0c] = dNdtByV[sr0c] + q2;
            dNdtByV[sr1c] = dNdtByV[sr1c] + q2;



            const unsigned int sl0d = lhsSpeciesIndex1D22RR[rhsIndex+6];
            const unsigned int sl1d = lhsSpeciesIndex1D22RR[rhsIndex+7];
            
            const unsigned int sr0d = rhsSpeciesIndex1D22RR[rhsIndex+6];
            const unsigned int sr1d = rhsSpeciesIndex1D22RR[rhsIndex+7];

            double Kr3 = 0;
            double Kf3 = Kf_[KfIndex3];

            
            const double Kp3 = (ExpNegGbyRT[sr0d]*ExpNegGbyRT[sr1d])/(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]);
            double Kc3 = Kp3;
            Kc3 = std::max(Kc3,KcLimiter);
            Kr3 = Kf3/Kc3;         
            


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
        const unsigned int KfIndex3 = this->reversibleReaction22index[i+0];
        {
            const unsigned int sl0d = lhsSpeciesIndex1D22RR[rhsIndex+0];
            const unsigned int sl1d = lhsSpeciesIndex1D22RR[rhsIndex+1];
            
            const unsigned int sr0d = rhsSpeciesIndex1D22RR[rhsIndex+0];
            const unsigned int sr1d = rhsSpeciesIndex1D22RR[rhsIndex+1];

            double Kr3 = 0;
            double Kf3 = Kf_[KfIndex3];

            const double Kp = (ExpNegGbyRT[sr0d]*ExpNegGbyRT[sr1d])/(ExpNegGbyRT[sl0d]*ExpNegGbyRT[sl1d]);
            double Kc = Kp;
            Kc = std::max(Kc,KcLimiter);
            Kr3 = Kf3/Kc;         
            
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
OptReaction::update22IrreversibleReaction
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
OptReaction::update22NonEquilibriumReaction
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