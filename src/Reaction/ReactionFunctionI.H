void inline OptReaction::RF11
(
    unsigned int i,
    double Kf,
    const double* __restrict__ C,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    double Kr = 0;
    if(this->isIrreversible[i]==0)
    {
        const double Kp = ExpNegGbyRT[sr0]/ExpNegGbyRT[sl0];
        double Kc = Kp;
        Kc = std::max(Kc,KcLimiter);
        Kr = Kf/Kc;       
    }
    else if(this->isIrreversible[i]==2)
    {
        unsigned int j = i - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[j];
    }

    const double CF = C[sl0];
    const double CR = C[sr0];

    const double q = (Kf*CF) - (Kr*CR);

    dNdtByV[sl0] = dNdtByV[sl0] - q; 
    dNdtByV[sr0] = dNdtByV[sr0] + q;
}

void inline OptReaction::RF12
(
    unsigned int i,
    double Kf,
    const double* __restrict__ C,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];
    double Kr = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1])/(ExpNegGbyRT[sl0]);
        double Kc = Kp*this->Pow_pByRT_SumVki[3];
        Kc = std::max(Kc,KcLimiter);
        Kr = Kf/Kc;       
    }
    else if(this->isIrreversible[i]==2)
    {
        unsigned int j = i - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[j];
    }

    const double CF = C[sl0];
    const double CR = C[sr0]*C[sr1];

    const double q = (Kf*CF) - (Kr*CR);

 
    dNdtByV[sl0] = dNdtByV[sl0] - q; 
    dNdtByV[sr0] = dNdtByV[sr0] + q;
    dNdtByV[sr1] = dNdtByV[sr1] + q;         

}             


void inline OptReaction::RF13
(
    unsigned int i,
    double Kf,
    const double* __restrict__ C,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];
    const unsigned int sr2 = rhsSpeciesIndex1D[rhsOffset[i]+2];

    double Kr = 0;
    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]);
        double Kc = Kp*this->Pow_pByRT_SumVki[4];
        Kc = std::max(Kc,KcLimiter);
        Kr = Kf/Kc;         
    }
    else if(this->isIrreversible[i]==2)
    {
        unsigned int j = i - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[j];
    }

    const double CF = C[sl0];
    const double CR = C[sr0]*C[sr1]*C[sr2];

    const double q = (Kf*CF) - (Kr*CR);

    dNdtByV[sl0] = dNdtByV[sl0] - q;
    dNdtByV[sr0] = dNdtByV[sr0] + q;
    dNdtByV[sr1] = dNdtByV[sr1] + q;
    dNdtByV[sr2] = dNdtByV[sr2] + q; 
}             

void inline OptReaction::RF21
(
    unsigned int i,
    double Kf,
    const double* __restrict__ C,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];


    double Kr = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = ExpNegGbyRT[sr0]/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
        double Kc = Kp*this->Pow_pByRT_SumVki[1];
        Kc = std::max(Kc,KcLimiter);
        Kr = Kf/Kc;        
    }
    else if(this->isIrreversible[i]==2)
    {
        unsigned int j = i - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[j];
    }
    
    const double CF = C[sl0]*C[sl1];
    const double CR = C[sr0];

    const double q = (Kf*CF) - (Kr*CR);



    dNdtByV[sl0] = dNdtByV[sl0] - q; 
    dNdtByV[sl1] = dNdtByV[sl1] - q;
    dNdtByV[sr0] = dNdtByV[sr0] + q;        
}             

void inline OptReaction::RF22
(
    unsigned int i,
    double Kf,
    const double* __restrict__ C,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];

    double Kr = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
        double Kc = Kp;
        Kc = std::max(Kc,KcLimiter);
        Kr = Kf/Kc;         
    }
    else if(this->isIrreversible[i]==2)
    {
        unsigned int j = i - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[j];
    }

    const double CF = C[sl0]*C[sl1];
    const double CR = C[sr0]*C[sr1];


    const double q = (Kf*CF) - (Kr*CR);



    dNdtByV[sl0] = dNdtByV[sl0] - q; 
    dNdtByV[sl1] = dNdtByV[sl1] - q;
    dNdtByV[sr0] = dNdtByV[sr0] + q;
    dNdtByV[sr1] = dNdtByV[sr1] + q; 
}

void inline OptReaction::RF23
(
    unsigned int i,
    double Kf,
    const double* __restrict__ C,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];
    const unsigned int sr2 = rhsSpeciesIndex1D[rhsOffset[i]+2];

    double Kr = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
        double Kc = Kp*this->Pow_pByRT_SumVki[3];
        Kc = std::max(Kc,KcLimiter);
        Kr = Kf/Kc;    
    }
    else if(this->isIrreversible[i]==2)
    {
        unsigned int j = i - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[j];
    }

    const double CF = C[sl0]*C[sl1];
    const double CR = C[sr0]*C[sr1]*C[sr2];
    const double q = (Kf*CF) - (Kr*CR);


    dNdtByV[sl0] = dNdtByV[sl0] - q;
    dNdtByV[sl1] = dNdtByV[sl1] - q;
    dNdtByV[sr0] = dNdtByV[sr0] + q;
    dNdtByV[sr1] = dNdtByV[sr1] + q;
    dNdtByV[sr2] = dNdtByV[sr2] + q;
}

void inline OptReaction::RF31
(
    unsigned int i,
    double Kf,
    const double* __restrict__ C,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{

    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
    const unsigned int sl2 = lhsSpeciesIndex1D[lhsOffset[i]+2];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];

    double Kr = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
        double Kc = Kp*this->Pow_pByRT_SumVki[0];
        Kc = std::max(Kc,KcLimiter);
        Kr = Kf/Kc;         
    }
    else if(this->isIrreversible[i]==2)
    {
        unsigned int j = i - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[j];
    }

    const double CF = C[sl0]*C[sl1]*C[sl2];
    const double CR = C[sr0];

    const double q = (Kf*CF) - (Kr*CR);

    dNdtByV[sl0] = dNdtByV[sl0] - q; 
    dNdtByV[sl1] = dNdtByV[sl1] - q; 
    dNdtByV[sl2] = dNdtByV[sl2] - q; 
    dNdtByV[sr0] = dNdtByV[sr0] + q;
}             

void inline OptReaction::RF32
(
    unsigned int i,
    double Kf,
    const double* __restrict__ C,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
    const unsigned int sl2 = lhsSpeciesIndex1D[lhsOffset[i]+2];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];

    double Kr = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
        double Kc = Kp*this->Pow_pByRT_SumVki[1];
        Kc = std::max(Kc,KcLimiter);
        Kr = Kf/Kc;
    }
    else if(this->isIrreversible[i]==2)
    {
        unsigned int j = i - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[j];
    }

    const double CF = C[sl0]*C[sl1]*C[sl2];
    const double CR = C[sr0]*C[sr1];

    const double q = (Kf*CF) - (Kr*CR);




    
    dNdtByV[sl0] = dNdtByV[sl0] - q; 
    dNdtByV[sl1] = dNdtByV[sl1] - q; 
    dNdtByV[sl2] = dNdtByV[sl2] - q; 
    dNdtByV[sr0] = dNdtByV[sr0] + q;
    dNdtByV[sr1] = dNdtByV[sr1] + q;
}

void inline OptReaction::RF33
(
    unsigned int i,
    double Kf,
    const double* __restrict__ C,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
    const unsigned int sl2 = lhsSpeciesIndex1D[lhsOffset[i]+2];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];
    const unsigned int sr2 = rhsSpeciesIndex1D[rhsOffset[i]+2];
    double Kr = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
        double Kc = Kp;
        Kc = std::max(Kc,KcLimiter);
        Kr = Kf/Kc;
    }
    else if(this->isIrreversible[i]==2)
    {
        unsigned int j = i - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[j];
    }

    const double CF = C[sl0]*C[sl1]*C[sl2];
    const double CR = C[sr0]*C[sr1];

    const double q = (Kf*CF) - (Kr*CR);

    dNdtByV[sl0] = dNdtByV[sl0] - q; 
    dNdtByV[sl1] = dNdtByV[sl1] - q; 
    dNdtByV[sl2] = dNdtByV[sl2] - q; 
    dNdtByV[sr0] = dNdtByV[sr0] + q;
    dNdtByV[sr1] = dNdtByV[sr1] + q;
    dNdtByV[sr2] = dNdtByV[sr2] + q;        
}        


void inline OptReaction::RFGI
(
    unsigned int i,
    double Kf,
    const double* __restrict__ C,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    int sumVki = 0;
    double Kp = 1.0;
    double CR = 1.0;
    double CF = 1.0;
    double Kc_ = 1.0;
    double Kr = 0;
    if(this->isIrreversible[i]==0)
    {
        for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->lhsSpeciesIndex[i][j];
            Kp = Kp / ExpNegGbyRT[si];
            sumVki = sumVki - 1;
            CF = CF * C[si];
        }
        
        for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[i][j];
            Kp = Kp * ExpNegGbyRT[si];
            sumVki = sumVki + 1;
            CR = CR * C[si];
        }

        Kc_ = Kp*this->Pow_pByRT_SumVki_I[sumVki];
        Kc_ = std::max(Kc_,1.49011611938476E-08);
        Kr = Kf/Kc_;
    }
    else if(this->isIrreversible[i]==2)
    {
        unsigned int J = i - this->Ikf[1] + this->Ikf[9];

        for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->lhsSpeciesIndex[i][j];
            CF = CF * C[si];
        }
        
        for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[i][j];
            CR = CR * C[si];
        }
        Kr = this->Kf_[J];
    }
    else
    {
        for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->lhsSpeciesIndex[i][j];
            CF = CF * C[si];
        }
        
        for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[i][j];
            CR = CR * C[si];
        }
    }

    const double q = Kf*CF-Kr*CR;

    for(unsigned int j = 0; j < lhsSpeciesIndex[i].size();j++)
    {
        const unsigned int si = lhsSpeciesIndex[i][j];
        dNdtByV[si] = dNdtByV[si] - q;
    }
    
    for(unsigned int j = 0; j < rhsSpeciesIndex[i].size();j++)
    {
        const unsigned int si = rhsSpeciesIndex[i][j];
        dNdtByV[si] = dNdtByV[si] + q;
    }

}

void inline OptReaction::RFGNI
(
    unsigned int i,
    double Kf,
    const double* __restrict__ C,
    double*  __restrict__ dNdtByV,
    const double*  __restrict__ ExpNegGbyRT
)const noexcept
{
    double CR = 1.0;
    double CF = 1.0;
    double Kc_ = 1.0;
    double Kr = 0;
    if(this->isIrreversible[i]==0) // Reversible reaction
    {
        double Kp = 0;
        double sumVki = 0;
        for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->lhsSpeciesIndex[i][j];
            const double sl = this->lhsStoichCoeff[i][j];
            const double el = this->lhsReactionOrder[i][j];
            Kp += sl*this->negGstdByRT[si];
            sumVki = sumVki - sl;
            CF = CF * (C[si] >= small || el >= 1 ? std::pow(std::max(C[si], 0.0), el) : 0.0);
        }
        
        for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[i][j];
            const double sr = this->rhsStoichCoeff[i][j];
            const double er = this->rhsReactionOrder[i][j];
            Kp -= sr*this->negGstdByRT[si];
            sumVki = sumVki + sr;
            CR = CR * (C[si] >= small || er >= 1 ? std::pow(std::max(C[si], 0.0), er) : 0.0);            
        }
        Kp = std::exp(Kp);
        Kc_ = Kp*std::pow(this->Pstd/(this->Ru*this->T),sumVki);
        Kc_ = std::max(Kc_,1.49011611938476E-08);
        Kr = Kf/Kc_;
    }
    else if(this->isIrreversible[i]==2)
    {
        unsigned int J = i - this->Ikf[1] + this->Ikf[9];
        Kr = this->Kf_[J];

        for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->lhsSpeciesIndex[i][j];
            const double el = this->lhsReactionOrder[i][j];
            CF = CF * (C[si] >= small || el >= 1 ? std::pow(std::max(C[si], 0.0), el) : 0.0);
        }
        
        for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[i][j];
            const double er = this->rhsReactionOrder[i][j];
            CR = CR * (C[si] >= small || er >= 1 ? std::pow(std::max(C[si], 0.0), er) : 0.0);            
        }
    }
    else
    {
        for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->lhsSpeciesIndex[i][j];
            const double el = this->lhsReactionOrder[i][j];
            CF = CF * (C[si] >= small || el >= 1 ? std::pow(std::max(C[si], 0.0), el) : 0.0);
        }
        
        for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[i][j];
            const double er = this->rhsReactionOrder[i][j];
            CR = CR * (C[si] >= small || er >= 1 ? std::pow(std::max(C[si], 0.0), er) : 0.0);            
        }
    }

    const double q = Kf*CF-Kr*CR;

    for(unsigned int j = 0; j < lhsSpeciesIndex[i].size();j++)
    {
        const unsigned int si = lhsSpeciesIndex[i][j];
        const double sl = lhsStoichCoeff[i][j];
        dNdtByV[si] = dNdtByV[si] - sl*q;
    }
    
    for(unsigned int j = 0; j < rhsSpeciesIndex[i].size();j++)
    {
        const unsigned int si = rhsSpeciesIndex[i][j];
        const double sl = rhsStoichCoeff[i][j];        
        dNdtByV[si] = dNdtByV[si] + sl*q;
    }
}             

