void inline OptReaction::JF11
(
    unsigned int i,
    double Kf,
    double dKfdT,
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{
    const unsigned int sl0 = this->lhsSpeciesIndex1D[this->lhsOffset[i]+0];
    const unsigned int sr0 = this->rhsSpeciesIndex1D[this->rhsOffset[i]+0];

    double Kr = 0;
    double dKrdT = 0;
    double invKc = 0;
    if(this->isIrreversible[i]==0)
    {
        const double Kp = ExpNegGbyRT[sr0]/ExpNegGbyRT[sl0];
        const double Kc = Kp > 1.4901171103413047e-8?Kp:1.4901171103413047e-8;
        const double sumVdBdT = (dBdT[sr0] ) - (dBdT[sl0] );
        const double dKcdTByKc = sumVdBdT;
        invKc = 1.0/Kc;        
        Kr = Kf*invKc;
        dKrdT = (dKfdT*invKc - (Kc > 1.4901171103413047e-8 ? Kr*dKcdTByKc : 0));
        const double dCrdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
    }
    else if(this->isIrreversible[i]==2)
    {
        Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
        dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
        invKc = Kr/Kf;
        const double dCrdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);        
    }
    const double CF = C[sl0];
    const double CR = C[sr0];

    const double q = (Kf*CF) - (Kr*CR);
    dNdtByV[sl0] = dNdtByV[sl0] - q;
    dNdtByV[sr0] = dNdtByV[sr0] + q;

    const double dqdT = (dKfdT*CF)-(dKrdT*CR);
    ddNdtByVdcTp[sl0*(this->alignN)+(this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sr0*(this->alignN)+(this->nSpecies)] += dqdT;

    const double dCfdC0 = 1;
    ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
    ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;

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
        }   
    }
}   
void inline OptReaction::JF12
(
    unsigned int i,
    double Kf,
    double dKfdT,
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{
    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];

    double Kr = 0;
    double dKrdT = 0;
    double invKc = 0;
    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1])/(ExpNegGbyRT[sl0]);
        double Kc = Kp*this->Pow_pByRT_SumVki[3];  
        Kc = Kc > 1.4901171103413047e-8?Kc:1.4901171103413047e-8;
        invKc = 1.0/Kc;
        const double sumVdBdT = (dBdT[sr0] + dBdT[sr1]) - (dBdT[sl0] );
        const double dKcdTByKc = sumVdBdT - 1*invT;
        Kr = Kf*invKc;
        dKrdT = (dKfdT*invKc - (Kc > 1.4901171103413047e-8 ? Kr*dKcdTByKc : 0));

        const double dCrdC0 = C[sr1];
        const double dCrdC1 = C[sr0];     
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);
    }
    else if(this->isIrreversible[i]==2)
    {
        Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
        dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
        invKc = Kr/Kf;
        const double dCrdC0 = C[sr1];
        const double dCrdC1 = C[sr0];     
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);
    }


    const double CF = C[sl0];
    const double CR = C[sr0]*C[sr1];

    const double q = (Kf*CF) - (Kr*CR);
    dNdtByV[sl0] = dNdtByV[sl0] - q; 
    dNdtByV[sr0] = dNdtByV[sr0] + q;
    dNdtByV[sr1] = dNdtByV[sr1] + q;   

    const double dqdT = (dKfdT*CF)-(dKrdT*CR);
    ddNdtByVdcTp[sl0*(this->alignN)+(this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sr0*(this->alignN)+(this->nSpecies)] += dqdT;
    ddNdtByVdcTp[sr1*(this->alignN)+(this->nSpecies)] += dqdT;

    const double dCfdC0 = 1; 
    ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
    ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
    ddNdtByVdcTp[sr1*(this->alignN)+sl0] += Kf*dCfdC0;

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
        } 
    }
}             


void inline OptReaction::JF13
(
    unsigned int i,
    double Kf,
    double dKfdT,
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];
    const unsigned int sr2 = rhsSpeciesIndex1D[rhsOffset[i]+2];

    double Kr = 0;
    double dKrdT = 0;
    double invKc = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]);
        double Kc = Kp*this->Pow_pByRT_SumVki[4];  
        Kc = Kc > 1.4901171103413047e-8?Kc:1.4901171103413047e-8;
        invKc = 1.0/Kc;
        Kr = Kf*invKc;
        const double sumVdBdT = (dBdT[sr0] + dBdT[sr1] + dBdT[sr2]) - (dBdT[sl0] );
        const double dKcdTByKc = sumVdBdT - 2*invT;
        dKrdT = (dKfdT*invKc - (Kc > 1.4901171103413047e-8 ? Kr*dKcdTByKc : 0));    

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
}            
void inline OptReaction::JF21
(
    unsigned int i,
    double Kf,
    double dKfdT,
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{
    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];

    double Kr = 0;
    double dKrdT = 0;
    double invKc = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
        double Kc = Kp*this->Pow_pByRT_SumVki[1];  
        Kc = Kc > 1.4901171103413047e-8?Kc:1.4901171103413047e-8;
        const double sumVdBdT = (dBdT[sr0] ) - (dBdT[sl0] + dBdT[sl1]);
        const double dKcdTByKc = sumVdBdT +1*invT;
        invKc = 1.0/Kc;  
        Kr = Kf*invKc;        
        dKrdT = (dKfdT*invKc - (Kc > 1.4901171103413047e-8 ? Kr*dKcdTByKc : 0));  
        
        const double dCrdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl1*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);     
    }
    else if(this->isIrreversible[i]==2)
    {
        Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
        dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
        invKc = Kr/Kf;
        const double dCrdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl1*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
    }

    const double CF = C[sl0]*C[sl1];
    const double CR = C[sr0];

    const double q = (Kf*CF) - (Kr*CR);
    dNdtByV[sl0] = dNdtByV[sl0] - q; 
    dNdtByV[sl1] = dNdtByV[sl1] - q;
    dNdtByV[sr0] = dNdtByV[sr0] + q;

    const double dqdT = (dKfdT*CF)-(dKrdT*CR);
    ddNdtByVdcTp[sl0*(this->alignN)+(this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sl1*(this->alignN)+(this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sr0*(this->alignN)+(this->nSpecies)] += dqdT;
    
    const double dCfdC0 = C[sl1];
    const double dCfdC1 = C[sl0];
    ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
    ddNdtByVdcTp[sl1*(this->alignN)+sl0] -= Kf*dCfdC0;
    ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
    ddNdtByVdcTp[sl0*(this->alignN)+sl1] -= Kf*dCfdC1;
    ddNdtByVdcTp[sl1*(this->alignN)+sl1] -= Kf*dCfdC1;
    ddNdtByVdcTp[sr0*(this->alignN)+sl1] += Kf*dCfdC1;

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

            __m256d sl1v = load256d(&ddNdtByVdcTp[sl1*(this->alignN)+j]);
            sl1v = _mm256_sub_pd(sl1v,WdMdC);
            store256d(&ddNdtByVdcTp[sl1*(this->alignN)+j],sl1v);

            __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);
            sr0v = _mm256_add_pd(sr0v,WdMdC);
            store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);                              
        } 
    }    
}     

void inline OptReaction::JF22
(
    unsigned int i,
    double Kf,
    double dKfdT,
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];


    double Kr = 0;
    double dKrdT = 0;
    double invKc = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
        double Kc = Kp;  
        Kc = Kc > 1.4901171103413047e-8?Kc:1.4901171103413047e-8;
        const double sumVdBdT = (dBdT[sr0] + dBdT[sr1]) - (dBdT[sl0] + dBdT[sl1]);
        const double dKcdTByKc = sumVdBdT;    
        invKc = 1.0/Kc;    
        Kr = Kf*invKc;
        dKrdT = (dKfdT*invKc - (Kc > 1.4901171103413047e-8 ? Kr*dKcdTByKc : 0));     

        const double dCrdC0 = C[sr1];
        const double dCrdC1 = C[sr0];
        ddNdtByVdcTp[sl0*(this->alignN)+ sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl1*(this->alignN)+ sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl0*(this->alignN)+ sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sl1*(this->alignN)+ sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sr0*(this->alignN)+ sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr1*(this->alignN)+ sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+ sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr1*(this->alignN)+ sr1] += (-Kr*dCrdC1);        
    }
    else if(this->isIrreversible[i]==2)
    {
        Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
        dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
        invKc = Kr/Kf;
        const double dCrdC0 = C[sr1];
        const double dCrdC1 = C[sr0]; 
        ddNdtByVdcTp[sl0*(this->alignN)+ sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl1*(this->alignN)+ sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl0*(this->alignN)+ sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sl1*(this->alignN)+ sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sr0*(this->alignN)+ sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr1*(this->alignN)+ sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+ sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr1*(this->alignN)+ sr1] += (-Kr*dCrdC1);
    }

    const double CF = C[sl0]*C[sl1];
    const double CR = C[sr0]*C[sr1];

    const double q = (Kf*CF) - (Kr*CR);
    dNdtByV[sl0] = dNdtByV[sl0] - q;
    dNdtByV[sl1] = dNdtByV[sl1] - q;
    dNdtByV[sr0] = dNdtByV[sr0] + q;
    dNdtByV[sr1] = dNdtByV[sr1] + q;

    const double dqdT = (dKfdT*CF)-(dKrdT*CR);
    ddNdtByVdcTp[sl0*(this->alignN)+ (this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sl1*(this->alignN)+ (this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sr0*(this->alignN)+ (this->nSpecies)] += dqdT;
    ddNdtByVdcTp[sr1*(this->alignN)+ (this->nSpecies)] += dqdT;

    const double dCfdC0 = C[sl1];
    const double dCfdC1 = C[sl0];
    ddNdtByVdcTp[sl0*(this->alignN)+ sl0] -= Kf*dCfdC0;
    ddNdtByVdcTp[sl1*(this->alignN)+ sl0] -= Kf*dCfdC0;
    ddNdtByVdcTp[sl0*(this->alignN)+ sl1] -= Kf*dCfdC1;
    ddNdtByVdcTp[sl1*(this->alignN)+ sl1] -= Kf*dCfdC1;
    ddNdtByVdcTp[sr0*(this->alignN)+ sl0] += Kf*dCfdC0;
    ddNdtByVdcTp[sr1*(this->alignN)+ sl0] += Kf*dCfdC0;
    ddNdtByVdcTp[sr0*(this->alignN)+ sl1] += Kf*dCfdC1;
    ddNdtByVdcTp[sr1*(this->alignN)+ sl1] += Kf*dCfdC1;     

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

            __m256d sl1v = load256d(&ddNdtByVdcTp[sl1*(this->alignN)+j]);
            sl1v = _mm256_sub_pd(sl1v,WdMdC);
            store256d(&ddNdtByVdcTp[sl1*(this->alignN)+j],sl1v);

            __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);
            sr0v = _mm256_add_pd(sr0v,WdMdC);
            store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);

            __m256d sr1v = load256d(&ddNdtByVdcTp[sr1*(this->alignN)+j]);
            sr1v = _mm256_add_pd(sr1v,WdMdC);
            store256d(&ddNdtByVdcTp[sr1*(this->alignN)+j],sr1v);
        } 
    }
}   



void inline OptReaction::JF23
(
    unsigned int i,
    double Kf,
    double dKfdT,
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{
    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];
    const unsigned int sr2 = rhsSpeciesIndex1D[rhsOffset[i]+2];

    double Kr = 0;
    double dKrdT = 0;
    double invKc = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]);
        double Kc = Kp*this->Pow_pByRT_SumVki[3];  
        Kc = Kc > 1.4901171103413047e-8?Kc:1.4901171103413047e-8;
        invKc = 1.0/Kc;
        const double sumVdBdT = (dBdT[sr0] + dBdT[sr1] + dBdT[sr2]) - (dBdT[sl0] + dBdT[sl1] );
        const double dKcdTByKc = sumVdBdT - 1*invT;
        Kr = Kf*invKc;
        dKrdT = (dKfdT*invKc - (Kc > 1.4901171103413047e-8 ? Kr*dKcdTByKc : 0));  

        const double dCrdC0 = C[sr1]*C[sr2];
        const double dCrdC1 = C[sr0]*C[sr2];
        const double dCrdC2 = C[sr0]*C[sr1];
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl1*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr2*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sl1*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr2*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sl0*(this->alignN)+sr2] -= (-Kr*dCrdC2);
        ddNdtByVdcTp[sl1*(this->alignN)+sr2] -= (-Kr*dCrdC2);
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
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr2*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sl1*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr2*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sl0*(this->alignN)+sr2] -= (-Kr*dCrdC2);
        ddNdtByVdcTp[sl1*(this->alignN)+sr2] -= (-Kr*dCrdC2);
        ddNdtByVdcTp[sr0*(this->alignN)+sr2] += (-Kr*dCrdC2);
        ddNdtByVdcTp[sr1*(this->alignN)+sr2] += (-Kr*dCrdC2);
        ddNdtByVdcTp[sr2*(this->alignN)+sr2] += (-Kr*dCrdC2);
    }

    const double CF = C[sl0]*C[sl1];
    const double CR = C[sr0]*C[sr1]*C[sr2];

    const double q = (Kf*CF) - (Kr*CR);
    dNdtByV[sl0] = dNdtByV[sl0] - q;
    dNdtByV[sl1] = dNdtByV[sl1] - q;
    dNdtByV[sr0] = dNdtByV[sr0] + q;
    dNdtByV[sr1] = dNdtByV[sr1] + q;
    dNdtByV[sr2] = dNdtByV[sr2] + q;

    const double dqdT = (dKfdT*CF)-(dKrdT*CR);
    ddNdtByVdcTp[sl0*(this->alignN) + (this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sl1*(this->alignN) + (this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sr0*(this->alignN) + (this->nSpecies)] += dqdT;
    ddNdtByVdcTp[sr1*(this->alignN) + (this->nSpecies)] += dqdT;
    ddNdtByVdcTp[sr2*(this->alignN) + (this->nSpecies)] += dqdT;

    const double dCfdC0 = C[sl1];
    const double dCfdC1 = C[sl0];
    ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
    ddNdtByVdcTp[sl1*(this->alignN)+sl0] -= Kf*dCfdC0;
    ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
    ddNdtByVdcTp[sr1*(this->alignN)+sl0] += Kf*dCfdC0;
    ddNdtByVdcTp[sr2*(this->alignN)+sl0] += Kf*dCfdC0;
    ddNdtByVdcTp[sl0*(this->alignN)+sl1] -= Kf*dCfdC1;
    ddNdtByVdcTp[sl1*(this->alignN)+sl1] -= Kf*dCfdC1;
    ddNdtByVdcTp[sr0*(this->alignN)+sl1] += Kf*dCfdC1;
    ddNdtByVdcTp[sr1*(this->alignN)+sl1] += Kf*dCfdC1;
    ddNdtByVdcTp[sr2*(this->alignN)+sl1] += Kf*dCfdC1;

    if(i>=this->Ikf[2] && i <this->Ikf[6])
    { 
        const unsigned int k = i - this->Ikf[2];

        __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[k]);
        __m256d CF_ = _mm256_set1_pd(CF);            
        __m256d CR_ = _mm256_set1_pd(CR);  
        double dKrdC = 0;
        dKrdC = (this->dKfdC_[k]*invKc);
        __m256d dKrdC_ = _mm256_set1_pd(dKrdC);

        for(unsigned int j = 0; j < AlignSpecies;j=j+4)
        {
            __m256d dMdC = load256d(&this->ThirdBodyFactor1D[k*this->AlignSpecies+j]);
            __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
            WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

            __m256d sl0v = load256d(&ddNdtByVdcTp[sl0*(this->alignN)+j]);
            sl0v = _mm256_sub_pd(sl0v,WdMdC);
            store256d(&ddNdtByVdcTp[sl0*(this->alignN)+j],sl0v);

            __m256d sl1v = load256d(&ddNdtByVdcTp[sl1*(this->alignN)+j]);
            sl1v = _mm256_sub_pd(sl1v,WdMdC);
            store256d(&ddNdtByVdcTp[sl1*(this->alignN)+j],sl1v);

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
}    

void inline OptReaction::JF31
(
    unsigned int i,
    double Kf,
    double dKfdT,
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
    const unsigned int sl2 = lhsSpeciesIndex1D[lhsOffset[i]+2];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];


    double Kr = 0;
    double dKrdT = 0;
    double invKc = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
        double Kc = Kp*this->Pow_pByRT_SumVki[0];
        Kc = Kc > 1.4901171103413047e-8?Kc:1.4901171103413047e-8;  
        invKc = 1.0/Kc;
        Kr = Kf*invKc;
        const double sumVdBdT = (dBdT[sr0]) - (dBdT[sl0] + dBdT[sl1] + dBdT[sl2]);
        const double dKcdTByKc = sumVdBdT + 2*invT;
        dKrdT = (dKfdT*invKc - (Kc > 1.4901171103413047e-8 ? Kr*dKcdTByKc : 0));  
        
        const double dCrdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl1*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl2*(this->alignN)+sr0] -= (-Kr*dCrdC0);    
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);   
    }
    else if(this->isIrreversible[i]==2)
    {
        Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
        dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
        invKc = Kr/Kf;        
        const double dCrdC0 = 1;
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl1*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl2*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
    }

    const double CF = C[sl0]*C[sl1]*C[sl2];
    const double CR = C[sr0];

    const double q = (Kf*CF) - (Kr*CR);
    dNdtByV[sl0] = dNdtByV[sl0] - q; 
    dNdtByV[sl1] = dNdtByV[sl1] - q;
    dNdtByV[sl2] = dNdtByV[sl2] - q;
    dNdtByV[sr0] = dNdtByV[sr0] + q;

    const double dqdT = (dKfdT*CF)-(dKrdT*CR);
    ddNdtByVdcTp[sl0*(this->alignN) + (this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sl1*(this->alignN) + (this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sl2*(this->alignN) + (this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sr0*(this->alignN) + (this->nSpecies)] += dqdT;

    const double dCfdC0 = C[sl1]*C[sl2];
    const double dCfdC1 = C[sl0]*C[sl2];
    const double dCfdC2 = C[sl0]*C[sl1];
    ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
    ddNdtByVdcTp[sl1*(this->alignN)+sl0] -= Kf*dCfdC0;
    ddNdtByVdcTp[sl2*(this->alignN)+sl0] -= Kf*dCfdC0;    
    ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
    ddNdtByVdcTp[sl0*(this->alignN)+sl1] -= Kf*dCfdC1;
    ddNdtByVdcTp[sl1*(this->alignN)+sl1] -= Kf*dCfdC1;
    ddNdtByVdcTp[sl2*(this->alignN)+sl1] -= Kf*dCfdC1;
    ddNdtByVdcTp[sr0*(this->alignN)+sl1] += Kf*dCfdC1;
    ddNdtByVdcTp[sl0*(this->alignN)+sl2] -= Kf*dCfdC2;
    ddNdtByVdcTp[sl1*(this->alignN)+sl2] -= Kf*dCfdC2;
    ddNdtByVdcTp[sl2*(this->alignN)+sl2] -= Kf*dCfdC2;
    ddNdtByVdcTp[sr0*(this->alignN)+sl2] += Kf*dCfdC2;

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


            __m256d sl1v = load256d(&ddNdtByVdcTp[sl1*(this->alignN)+j]);
            sl1v = _mm256_sub_pd(sl1v,WdMdC);
            store256d(&ddNdtByVdcTp[sl1*(this->alignN)+j],sl1v);


            __m256d sl2v = load256d(&ddNdtByVdcTp[sl2*(this->alignN)+j]);
            sl2v = _mm256_sub_pd(sl2v,WdMdC);
            store256d(&ddNdtByVdcTp[sl2*(this->alignN)+j],sl2v);


            __m256d sr0v = load256d(&ddNdtByVdcTp[sr0*(this->alignN)+j]);
            sr0v = _mm256_add_pd(sr0v,WdMdC);
            store256d(&ddNdtByVdcTp[sr0*(this->alignN)+j],sr0v);
        } 
    }  
}

void inline OptReaction::JF32
(
    unsigned int i,
    double Kf,
    double dKfdT,
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{


    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
    const unsigned int sl2 = lhsSpeciesIndex1D[lhsOffset[i]+2];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];


    double Kr = 0;
    double dKrdT = 0;
    double invKc = 0;
    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
        double Kc = Kp*this->Pow_pByRT_SumVki[1];
        Kc = Kc > 1.4901171103413047e-8?Kc:1.4901171103413047e-8;  
        invKc = 1.0/Kc;
        Kr = Kf*invKc;
        const double sumVdBdT = (dBdT[sr0] + dBdT[sr1] ) - (dBdT[sl0] + dBdT[sl1] + dBdT[sl2]);
        const double dKcdTByKc = sumVdBdT + 1*invT;
        dKrdT = (dKfdT/Kc - (Kc > 1.4901171103413047e-8 ? Kr*dKcdTByKc : 0));  
        
        const double dCrdC0 = C[sr1];
        const double dCrdC1 = C[sr0];
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl1*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl2*(this->alignN)+sr0] -= (-Kr*dCrdC0);    
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sl1*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sl2*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);  
    }
    else if(this->isIrreversible[i]==1)
    {
        Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
        dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
        invKc = Kr/Kf;
        const double dCrdC0 = C[sr1];
        const double dCrdC1 = C[sr0];
        ddNdtByVdcTp[sl0*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl1*(this->alignN)+sr0] -= (-Kr*dCrdC0);
        ddNdtByVdcTp[sl2*(this->alignN)+sr0] -= (-Kr*dCrdC0);    
        ddNdtByVdcTp[sr0*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sr1*(this->alignN)+sr0] += (-Kr*dCrdC0);
        ddNdtByVdcTp[sl0*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sl1*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sl2*(this->alignN)+sr1] -= (-Kr*dCrdC1);
        ddNdtByVdcTp[sr0*(this->alignN)+sr1] += (-Kr*dCrdC1);
        ddNdtByVdcTp[sr1*(this->alignN)+sr1] += (-Kr*dCrdC1);  
    }
    const double CF = C[sl0]*C[sl1]*C[sl2];
    const double CR = C[sr0]*C[sr1];

    
    const double q = (Kf*CF) - (Kr*CR);
    dNdtByV[sl0] = dNdtByV[sl0] - q; 
    dNdtByV[sl1] = dNdtByV[sl1] - q;
    dNdtByV[sl2] = dNdtByV[sl2] - q;
    dNdtByV[sr0] = dNdtByV[sr0] + q;
    dNdtByV[sr1] = dNdtByV[sr1] + q;

    const double dqdT = (dKfdT*CF)-(dKrdT*CR);
    ddNdtByVdcTp[sl0*(this->alignN) + (this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sl1*(this->alignN) + (this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sl2*(this->alignN) + (this->nSpecies)] -= dqdT;
    ddNdtByVdcTp[sr0*(this->alignN) + (this->nSpecies)] += dqdT;
    ddNdtByVdcTp[sr1*(this->alignN) + (this->nSpecies)] += dqdT;
    
    const double dCfdC0 = C[sl1]*C[sl2];
    const double dCfdC1 = C[sl0]*C[sl2];
    const double dCfdC2 = C[sl0]*C[sl1];
    ddNdtByVdcTp[sl0*(this->alignN)+sl0] -= Kf*dCfdC0;
    ddNdtByVdcTp[sl1*(this->alignN)+sl0] -= Kf*dCfdC0;
    ddNdtByVdcTp[sl2*(this->alignN)+sl0] -= Kf*dCfdC0;    
    ddNdtByVdcTp[sr0*(this->alignN)+sl0] += Kf*dCfdC0;
    ddNdtByVdcTp[sr1*(this->alignN)+sl0] += Kf*dCfdC0;
    ddNdtByVdcTp[sl0*(this->alignN)+sl1] -= Kf*dCfdC1;
    ddNdtByVdcTp[sl1*(this->alignN)+sl1] -= Kf*dCfdC1;
    ddNdtByVdcTp[sl2*(this->alignN)+sl1] -= Kf*dCfdC1;
    ddNdtByVdcTp[sr0*(this->alignN)+sl1] += Kf*dCfdC1;
    ddNdtByVdcTp[sr1*(this->alignN)+sl1] += Kf*dCfdC1;
    ddNdtByVdcTp[sl0*(this->alignN)+sl2] -= Kf*dCfdC2;
    ddNdtByVdcTp[sl1*(this->alignN)+sl2] -= Kf*dCfdC2;
    ddNdtByVdcTp[sl2*(this->alignN)+sl2] -= Kf*dCfdC2;
    ddNdtByVdcTp[sr0*(this->alignN)+sl2] += Kf*dCfdC2;
    ddNdtByVdcTp[sr1*(this->alignN)+sl2] += Kf*dCfdC2;

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
        } 
    }  
}    

void inline OptReaction::JF33
(
    unsigned int i,
    double Kf,
    double dKfdT,
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    const unsigned int sl0 = lhsSpeciesIndex1D[lhsOffset[i]+0];
    const unsigned int sl1 = lhsSpeciesIndex1D[lhsOffset[i]+1];
    const unsigned int sl2 = lhsSpeciesIndex1D[lhsOffset[i]+2];
    const unsigned int sr0 = rhsSpeciesIndex1D[rhsOffset[i]+0];
    const unsigned int sr1 = rhsSpeciesIndex1D[rhsOffset[i]+1];
    const unsigned int sr2 = rhsSpeciesIndex1D[rhsOffset[i]+2];

    double Kr = 0;
    double dKrdT = 0;
    double invKc = 0;

    if(this->isIrreversible[i]==0)
    {
        const double Kp = (ExpNegGbyRT[sr0]*ExpNegGbyRT[sr1]*ExpNegGbyRT[sr2])/(ExpNegGbyRT[sl0]*ExpNegGbyRT[sl1]*ExpNegGbyRT[sl2]);
        double Kc = Kp;
        Kc = Kc > 1.4901171103413047e-8?Kc:1.4901171103413047e-8;  
        invKc = 1.0/Kc;
        Kr = Kf*invKc;
        const double sumVdBdT = (dBdT[sr0] + dBdT[sr1] + dBdT[sr2]) - (dBdT[sl0] + dBdT[sl1] + dBdT[sl2]);
        const double dKcdTByKc = sumVdBdT;    
        dKrdT = (dKfdT*invKc - (Kc > 1.4901171103413047e-8 ? Kr*dKcdTByKc : 0));  
        
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
}

void inline OptReaction::JFGI
(
    unsigned int i,
    double Kf,
    double dKfdT,
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    double Kr = 0;
    double dKrdT = 0;
    double CF = 1.0;
    double CR = 1.0;
    double invKc = 0;
    if(this->isIrreversible[i]==0)
    {
        int sumVki = 0;
        double Kp = 1.0;
        double sumVdBdT = 0.0;
        for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->lhsSpeciesIndex[i][j];
            Kp = Kp / this->tmp_Exp[si];
            sumVki = sumVki - 1;
            sumVdBdT = sumVdBdT - dBdT[si];
            CF = CF * C[si];
        }

        for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[i][j];
            Kp = Kp * this->tmp_Exp[si];
            sumVki = sumVki + 1;
            sumVdBdT = sumVdBdT + dBdT[si];
            CR = CR * C[si];
        }
        //this->update_Pow_pByRT_SumVki(T);
        const double Kc = std::max(Kp*this->Pow_pByRT_SumVki_I[sumVki],1.4901171103413047e-8);
        invKc = 1.0/Kc;            

        const double dKcdTByKc = sumVdBdT - sumVki*invT;
        Kr = Kf*invKc;
        dKrdT = dKfdT*invKc - (Kc > 1.4901171103413047e-8 ? Kr*dKcdTByKc : 0);
    }
    else if(this->isIrreversible[i]==2)
    {
        Kr = this->Kf_[i - Ikf[1] + Ikf[9]];
        dKrdT = this->dKfdT_[i - Ikf[1] + Ikf[9]];
        invKc = Kr/Kf;        
        for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[i][j];
            CR = CR * C[si];
        }
        for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int sl = this->lhsSpeciesIndex[i][j];
            CF = CF * C[sl];
        }
    }
    else
    {
        for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[i][j];
            CR = CR * C[si];
        }
        for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
        {
            const unsigned int sl = this->lhsSpeciesIndex[i][j];
            CF = CF * C[sl];
        }
    }
    const double q = Kf*CF - Kr*CR;
    const double dqdT = (dKfdT*CF-dKrdT*CR);

    for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
    {
        const unsigned int si = this->lhsSpeciesIndex[i][j];
        ddNdtByVdcTp[si*(this->alignN) + (this->nSpecies)] -= dqdT;
        dNdtByV[si] = dNdtByV[si] - q;                        
    }

    for(unsigned int j = 0; j < this->rhsSpeciesIndex[i].size();j++)
    {
        const unsigned int si = this->rhsSpeciesIndex[i][j];
        ddNdtByVdcTp[si*(this->alignN) + (this->nSpecies)] += dqdT;
        dNdtByV[si] = dNdtByV[si] + q;                        
    }            
                

    for(unsigned int j = 0; j < this->lhsSpeciesIndex[i].size();j++)
    {
        const unsigned int sj = this->lhsSpeciesIndex[i][j];

        double dCfdCj = 1.0;
        for(unsigned int k = 0; k < this->lhsSpeciesIndex[i].size();k++)
        {
            const unsigned int sk = this->lhsSpeciesIndex[i][k];

            dCfdCj = dCfdCj*
            (
                (j==k)? 1:C[sk]
            );
        }

        const double KfdCfdCj = Kf*dCfdCj;

        for(unsigned int k = 0; k < this->lhsSpeciesIndex[i].size();k++)
        {
            const unsigned int sk = this->lhsSpeciesIndex[i][k];
            ddNdtByVdcTp[sk*(this->alignN)+sj] -= KfdCfdCj;
        }

        for(unsigned int k = 0; k < rhsSpeciesIndex[i].size();k++)
        {
            const unsigned int sk = rhsSpeciesIndex[i][k];
            ddNdtByVdcTp[sk*(this->alignN)+sj] += KfdCfdCj;
        }
    }
    if(this->isIrreversible[i]==0||this->isIrreversible[i]==2)
    {
        for(unsigned int j = 0; j < rhsSpeciesIndex[i].size();j++)
        {
            const unsigned int sj = rhsSpeciesIndex[i][j];
    
            double dCrdCj = 1.0;
            for(unsigned int k = 0; k < rhsSpeciesIndex[i].size();k++)
            {
                const unsigned int sk = rhsSpeciesIndex[i][k];
    
                dCrdCj = dCrdCj*
                (
                    (j==k)? 1:C[sk]
                );
            }
    
            const double negKrdCrdCj = -Kr*dCrdCj;

            for(unsigned int k = 0; k < this->lhsSpeciesIndex[i].size();k++)
            {
                const unsigned int sk = this->lhsSpeciesIndex[i][k];
                ddNdtByVdcTp[sk*(this->alignN)+sj] -= negKrdCrdCj;
            }

            for(unsigned int k = 0; k < rhsSpeciesIndex[i].size();k++)
            {
                const unsigned int sk = rhsSpeciesIndex[i][k];
                ddNdtByVdcTp[sk*(this->alignN)+sj] += negKrdCrdCj;                        
            }
        }          
    }

    if(i>=this->Ikf[2] && i <this->Ikf[6])
    {
        const unsigned int kk = i - this->Ikf[2];

        __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[kk]);
        __m256d CF_ = _mm256_set1_pd(CF);            
        __m256d CR_ = _mm256_set1_pd(CR);  
        double dKrdC = 0;
        dKrdC = (this->dKfdC_[kk]*invKc);
        __m256d dKrdC_ = _mm256_set1_pd(dKrdC);

        for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
        {
            __m256d dMdC = load256d(&this->ThirdBodyFactor1D[kk*this->AlignSpecies+j]);
            __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
            WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

            for(unsigned int k = 0; k < this->lhsSpeciesIndex[i].size();k++)
            {
                const unsigned int sk = this->lhsSpeciesIndex[i][k];
                __m256d skv = load256d(&ddNdtByVdcTp[sk*(this->alignN)+j]);
                skv = _mm256_sub_pd(skv,WdMdC);
                store256d(&ddNdtByVdcTp[sk*(this->alignN)+j],skv);
            }

            for(unsigned int k = 0; k < rhsSpeciesIndex[i].size();k++)
            {

                const unsigned int sk = rhsSpeciesIndex[i][k];
                __m256d skv = load256d(&ddNdtByVdcTp[sk*(this->alignN)+j]);
                skv = _mm256_add_pd(skv,WdMdC);
                store256d(&ddNdtByVdcTp[sk*(this->alignN)+j],skv);                                     
            }
        } 
    }   
}

void inline OptReaction::JFGNI
(
    unsigned int iii,
    double Kf,
    double dKfdT,
    const double* __restrict__ C,
    double* __restrict__ dNdtByV,
    double* __restrict__ ddNdtByVdcTp,
    const double* __restrict__ ExpNegGbyRT,
    const double* __restrict__ dBdT
)const noexcept
{

    double Kr = 0;
    double dKrdT = 0;
    double CF = 1.0;
    double CR = 1.0;
    double invKc = 0;
    if(this->isIrreversible[iii]==0)
    {
        double sumVki = 0;
        double Kp = 1.0;
        double sumVdBdT = 0.0;
        for(unsigned int j = 0; j < this->lhsSpeciesIndex[iii].size();j++)
        {
            const unsigned int si = this->lhsSpeciesIndex[iii][j];
            const double sl = this->lhsStoichCoeff[iii][j];
            const double el = this->lhsReactionOrder[iii][j];
            Kp += sl*this->negGstdByRT[si];
            sumVdBdT = sumVdBdT - sl*dBdT[si];
            sumVki = sumVki - sl;
            CF = CF * (C[si] >= small || el >= 1 ? std::pow(std::max(C[si], 0.0), el) : 0.0);
        }

        for(unsigned int j = 0; j < this->rhsSpeciesIndex[iii].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[iii][j];
            const double sr = this->rhsStoichCoeff[iii][j];
            const double er = this->rhsReactionOrder[iii][j];            
            Kp -= sr*this->negGstdByRT[si];
            sumVdBdT = sumVdBdT + sr*dBdT[si];
            sumVki = sumVki + sr;
            CR = CR * (C[si] >= small || er >= 1.0 ? std::pow(std::max(C[si], 0.0), er) : 0.0);    
        }
        Kp = std::exp(Kp);        
        double Kc = Kp*std::pow(this->Pstd/(this->Ru*this->T),sumVki);
        Kc = std::max(Kc,1.4901171103413047e-8);
        invKc = 1.0/Kc;

        const double dKcdTByKc = sumVdBdT - sumVki*invT;

        Kr = Kf*invKc;
        dKrdT = dKfdT*invKc - (Kc > 1.4901171103413047e-8 ? Kr*dKcdTByKc : 0.0);
    }
    else if(this->isIrreversible[iii]==2)
    {
        Kr = this->Kf_[iii - Ikf[1] + Ikf[9]];
        dKrdT = this->dKfdT_[iii - Ikf[1] + Ikf[9]];
        invKc = Kr/Kf;        
        for(unsigned int j = 0; j < this->rhsSpeciesIndex[iii].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[iii][j];
            const double er = lhsReactionOrder[iii][j];
            CR = CR * (C[si] >= small || er >= 1 ? std::pow(std::max(C[si], 0.0), er) : 0.0);   
        }
        for(unsigned int j = 0; j < lhsSpeciesIndex[iii].size();j++)
        {
            const unsigned int si = lhsSpeciesIndex[iii][j];
            const double el = lhsReactionOrder[iii][j];
            CF = CF * (C[si] >= small || el >= 1 ? std::pow(std::max(C[si], 0.0), el) : 0.0);            
        }
    }
    else
    {
        for(unsigned int j = 0; j < this->rhsSpeciesIndex[iii].size();j++)
        {
            const unsigned int si = this->rhsSpeciesIndex[iii][j];
            const double er = lhsReactionOrder[iii][j];
            CR = CR * (C[si] >= small || er >= 1 ? std::pow(std::max(C[si], 0.0), er) : 0.0);   
        }
        for(unsigned int j = 0; j < lhsSpeciesIndex[iii].size();j++)
        {
            const unsigned int si = lhsSpeciesIndex[iii][j];
            const double el = lhsReactionOrder[iii][j];
            CF = CF * (C[si] >= small || el >= 1 ? std::pow(std::max(C[si], 0.0), el) : 0.0);            
        }        
    }

    const double q = Kf*CF - Kr*CR;
    const double dwdT = (dKfdT*CF-dKrdT*CR);

    for(unsigned int j = 0; j < this->lhsSpeciesIndex[iii].size();j++)
    {
        const unsigned int si = this->lhsSpeciesIndex[iii][j];
        const double sl = this->lhsStoichCoeff[iii][j];
        ddNdtByVdcTp[si*(this->alignN) + (this->nSpecies)] -= sl*dwdT;
        dNdtByV[si] = dNdtByV[si] - sl*q;
    }

    for(unsigned int j = 0; j < this->rhsSpeciesIndex[iii].size();j++)
    {
        const unsigned int si = this->rhsSpeciesIndex[iii][j];
        const double sr = this->rhsStoichCoeff[iii][j];
        ddNdtByVdcTp[si*(this->alignN) + (this->nSpecies)] += sr*dwdT;
        dNdtByV[si] = dNdtByV[si] + sr*q;
    }            
                

    for(unsigned int j = 0; j < lhsSpeciesIndex[iii].size();j++)
    {
        const unsigned int sj = lhsSpeciesIndex[iii][j];

        double dCfdCj = 1.0;
        for(unsigned int i = 0; i < lhsSpeciesIndex[iii].size();i++)
        {
            const unsigned int si = lhsSpeciesIndex[iii][i];
            const double el = lhsReactionOrder[iii][i];

            if (i == j)
            {
                dCfdCj *=
                C[si] >= small || el >= 1
                    ? el*std::pow(std::max(C[si], 0.0), el - 1.0)
                    : 0;
            }
            else
            {
                dCfdCj *=
                    C[si] >= small || el >= 1
                    ? std::pow(std::max(C[si], 0.0), el)
                    : 0;
            }            
        }

        const double KfdCfdCj = Kf*dCfdCj;

        for(unsigned int i = 0; i < lhsSpeciesIndex[iii].size();i++)
        {
            const unsigned int si = lhsSpeciesIndex[iii][i];
            const double sl = lhsStoichCoeff[iii][i];
            ddNdtByVdcTp[si*(this->alignN)+sj] -= sl*KfdCfdCj;
        }

        for(unsigned int i = 0; i < rhsSpeciesIndex[iii].size();i++)
        {
            const unsigned int si = rhsSpeciesIndex[iii][i];
            const double sr = rhsStoichCoeff[iii][i];
            ddNdtByVdcTp[si*(this->alignN)+sj] += sr*KfdCfdCj;
        }
    }

    if(this->isIrreversible[iii]==0||this->isIrreversible[iii]==2)
    {
        for(unsigned int j = 0; j < rhsSpeciesIndex[iii].size();j++)
        {
            const unsigned int sj = rhsSpeciesIndex[iii][j];
    
            double dCrdCj = 1.0;
            for(unsigned int i = 0; i < rhsSpeciesIndex[iii].size();i++)
            {

                const unsigned int si = rhsSpeciesIndex[iii][i];
                const double er = rhsReactionOrder[iii][i];
                dCrdCj = dCrdCj*((j==i)? 1:C[si]);

                if (i == j)
                {
                    dCrdCj *=
                        C[si] >= small || er >= 1.0
                      ? er*std::pow(std::max(C[si], 0.0), er - 1.0)
                      : 0;
                }
                else
                {
                    dCrdCj *=
                        C[si] >= small || er >= 1.0
                      ? std::pow(std::max(C[si], 0.0), er)
                      : 0;
                }
            }
    
            const double negKrdCrdCj = -Kr*dCrdCj;

            for(unsigned int i = 0; i < lhsSpeciesIndex[iii].size();i++)
            {
                const unsigned int si = lhsSpeciesIndex[iii][i];
                const double sl = lhsReactionOrder[iii][j];
                ddNdtByVdcTp[si*(this->alignN)+sj] -= sl*negKrdCrdCj;
            }

            for(unsigned int i = 0; i < rhsSpeciesIndex[iii].size();i++)
            {
                const unsigned int si = rhsSpeciesIndex[iii][i];
                const double sr = rhsReactionOrder[iii][j];
                ddNdtByVdcTp[si*(this->alignN)+sj] += sr*negKrdCrdCj;                        
            }
        }
    }

    if(iii>=this->Ikf[2] && iii <this->Ikf[6])
    {
        const unsigned int kk = iii - this->Ikf[2];

        __m256d dKfdC = _mm256_set1_pd(this->dKfdC_[kk]);
        __m256d CF_ = _mm256_set1_pd(CF);            
        __m256d CR_ = _mm256_set1_pd(CR);  
        double dKrdC = 0;
        dKrdC = (this->dKfdC_[kk]*invKc);
        __m256d dKrdC_ = _mm256_set1_pd(dKrdC);

        for(unsigned int j = 0; j < this->AlignSpecies;j=j+4)
        {
            __m256d dMdC = load256d(&this->ThirdBodyFactor1D[kk*this->AlignSpecies+j]);
            __m256d WdMdC = _mm256_mul_pd(_mm256_mul_pd(dKrdC_,dMdC),CR_);
            WdMdC = _mm256_fmsub_pd(_mm256_mul_pd(dKfdC,dMdC),CF_,WdMdC);

            for(unsigned int i = 0; i < lhsSpeciesIndex[iii].size();i++)
            {
                const unsigned int si = lhsSpeciesIndex[iii][i];
                const double sl = lhsStoichCoeff[iii][i];
                __m256d Jac = load256d(&ddNdtByVdcTp[si*(this->alignN)+j]);
                __m256d slv = _mm256_set1_pd(sl);

                Jac = _mm256_sub_pd(Jac,_mm256_mul_pd(WdMdC,slv));
                store256d(&ddNdtByVdcTp[si*(this->alignN)+j],Jac);
            }

            for(unsigned int i = 0; i < rhsSpeciesIndex[iii].size();i++)
            {
                const unsigned int si = rhsSpeciesIndex[iii][i];
                const double sr = rhsStoichCoeff[iii][i];
                __m256d Jac = load256d(&ddNdtByVdcTp[si*(this->alignN)+j]);
                __m256d srv = _mm256_set1_pd(sr);

                Jac = _mm256_add_pd(Jac,_mm256_mul_pd(WdMdC,srv));
                store256d(&ddNdtByVdcTp[si*(this->alignN)+j],Jac);                                     
            }
        } 
    }   
}
