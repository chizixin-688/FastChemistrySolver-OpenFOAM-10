void OptReaction::evalSRIPartialDerivative(double T)const noexcept
{
    for (unsigned int i = 0;i<this->n_SRI;i++)
    {
        const unsigned int j = this->SRIFO[i];
        const unsigned int k = j - this->Ikf[4];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];
        const double Kinf = this->Kf_[j+this->offset_kinf];
        const double invKinf = 1.0/Kinf;
        const double K0 = this->Kf_[j];
        const double dKinfdT = this->dKfdT_[j+this->offset_kinf];
        double F ;
        double dFdT;
        double dFdPr;
        double M = tmp_M[m];
        const double Pr = K0*M*invKinf; 
        const double expbT = this->tmp_Exp[i+this->nSpecies+this->TroeFO.size()*3];
        const double expTc = this->tmp_Exp[i+this->nSpecies+this->TroeFO.size()*3+this->SRIFO.size()];
        const double logPr = std::log10(max(Pr, small));
        const double X = 1/(1 + (logPr*logPr));
        const double psi = this->a_[i]*expbT + expTc;

        const double d = this->d_[i];
        const double e = this->e_[i];

        F = d*std::pow(psi, X)*std::pow(T, e);

        const double dpsidT = this->a_[i]*this->b_[i]/(T*T)*expbT - 1/this->c_[i]*expTc;
        dFdT = F*(X/psi*dpsidT + e/T);
        static const double logTen = std::log(10);
        const double dlogPrdPr = Pr >= small ? 1/(logTen*Pr) : 0;
        const double dXdPr = -(X*X)*2*logPr*dlogPrdPr;
        dFdPr = F*std::log(psi)*dXdPr;
        const double dK0dT =  this->dKfdT_[j]; 
        const double invOnePlusPr = 1.0/(1.0+Pr);
        const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;

        const double dKdT   = Pr*dKinfdT;
        const double N1     = F*invOnePlusPr;
        const double N2     = Pr*dFdPr;
        const double N  = invOnePlusPr*F*K0;
        this->dKfdT_[j] = F*invOnePlusPr*dKdT + F*invOnePlusPr*invOnePlusPr*dPrdT*Kinf + K0*invOnePlusPr*dFdT*M;
        this->dKfdC_[m] =  K0*invOnePlusPr*(N1 + N2); 
        this->Kf_[j] = M*N; 

    }

}