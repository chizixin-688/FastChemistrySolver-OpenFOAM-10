
inline double OptReaction::SRI_F
(
    const double T_,
    const double Pr,
    const unsigned i
) const
{
    const double logPr = std::log10(max(Pr, small));

    const double X = 1/(1 + (logPr*logPr));

    const double psi = this->a_[i]*std::exp(-this->b_[i]/T_) + std::exp(-T_/this->c_[i]);

    return this->d_[i]*std::pow(psi, X)*std::pow(T_, this->e_[i]);
}

inline void OptReaction::SRI_F
(
    const double T_, 
    const double* __restrict__ Pr, 
    double* __restrict__ F,
    const unsigned int i
) const
{
    __m256d small = _mm256_set1_pd(2.2e-16);
    __m256d Pr_ = _mm256_loadu_pd(Pr);
    Pr_ = _mm256_max_pd(small,Pr_);
    const double invLog10 = 1.0/std::log(10);
    __m256d logPr_ = _mm256_mul_pd(vec256_logd(Pr_),_mm256_set1_pd(invLog10));
    __m256d X = _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_fmadd_pd(logPr_,logPr_,_mm256_set1_pd(1.0)));
    __m256d A_ = _mm256_loadu_pd(&this->a_[i]);
    __m256d B_ = _mm256_loadu_pd(&this->b_[i]);
    __m256d C_ = _mm256_loadu_pd(&this->c_[i]);
    __m256d D_ = _mm256_loadu_pd(&this->d_[i]);
    __m256d E_ = _mm256_loadu_pd(&this->e_[i]);
    __m256d R1 = _mm256_mul_pd(B_,_mm256_set1_pd(-this->invT));
    __m256d R2 = _mm256_div_pd(_mm256_set1_pd(-T_),C_);
    __m256d expbT = vec256_expd(R1);
    __m256d expTc = vec256_expd(R2); 
    __m256d psi = _mm256_fmadd_pd(A_,expbT,expTc); 
    __m256d powpX = vec256_powd(psi,X);
    __m256d powTe = vec256_powd(_mm256_set1_pd(T_),E_);
    __m256d F_ = _mm256_mul_pd(_mm256_mul_pd(D_,powpX),powTe);
    _mm256_storeu_pd(&F[0],F_);
}

inline void OptReaction::SRI_F_dFdT_dFdPr
(
    const double T_,
    const double Pr,
    const unsigned int i,
    double& __restrict__ F,
    double& __restrict__ dFdT,
    double& __restrict__ dFdPr
) const
{
    const double expbT = this->tmp_Exp[i+this->nSpecies+this->Troe.size()*3];
    const double expTc = this->tmp_Exp[i+this->nSpecies+this->Troe.size()*3+this->SRI.size()];
    const double logPr = std::log10(max(Pr, small));
    const double X = 1/(1 + (logPr*logPr));
    const double psi = this->a_[i]*expbT + expTc;
    F = this->d_[i]*std::pow(psi, X)*std::pow(T_, this->e_[i]);
    const double dpsidT = this->a_[i]*this->b_[i]/(T_*T_)*expbT - 1/this->c_[i]*expTc;
    dFdT = F*(X/psi*dpsidT + this->e_[i]/T_);
    static const double logTen = std::log(10);
    const double dlogPrdPr = Pr >= small ? 1/(logTen*Pr) : 0;
    const double dXdPr = -(X*X)*2*logPr*dlogPrdPr;
    dFdPr = F*std::log(psi)*dXdPr;
}


inline void OptReaction::SRI_F_dFdT_dFdPr
(
    const double T_,
    __m256d& __restrict__ Pr,
    __m256d& __restrict__ F,
    __m256d& __restrict__ dFdT,
    __m256d& __restrict__ dFdPr,
    const unsigned int i
) const
{
    __m256d small = _mm256_set1_pd(2.2e-16);
    __m256d Pr_ = Pr;
    __m256d cmp = _mm256_cmp_pd(Pr_,small,_CMP_GE_OQ);
    Pr = _mm256_add_pd(Pr,_mm256_set1_pd(1e-100));
    const double invLog10 = 1.0/std::log(10);
    __m256d logPr_ = _mm256_mul_pd(vec256_logd(Pr_),_mm256_set1_pd(invLog10));
    __m256d X = _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_fmadd_pd(logPr_,logPr_,_mm256_set1_pd(1.0)));
    __m256d A_ = _mm256_loadu_pd(&this->a_[i]);
    __m256d B_ = _mm256_loadu_pd(&this->b_[i]);
    __m256d C_ = _mm256_loadu_pd(&this->c_[i]);
    __m256d D_ = _mm256_loadu_pd(&this->d_[i]);
    __m256d E_ = _mm256_loadu_pd(&this->e_[i]);
    __m256d R1 = _mm256_mul_pd(B_,_mm256_set1_pd(-this->invT));
    __m256d R2 = _mm256_div_pd(_mm256_set1_pd(-T_),C_);
    __m256d expbT = vec256_expd(R1);
    __m256d expTc = vec256_expd(R2); 
    __m256d psi = _mm256_fmadd_pd(A_,expbT,expTc); 
    __m256d logpsi = vec256_logd(psi);
    __m256d powpX = vec256_powd(psi,X);
    __m256d powTe = vec256_powd(_mm256_set1_pd(T_),E_);    
    F = _mm256_mul_pd(_mm256_mul_pd(D_,powpX),powTe);
    __m256d dpsidT = _mm256_mul_pd(A_,B_);
    dpsidT = _mm256_mul_pd(dpsidT,expbT);
    dpsidT = _mm256_div_pd(dpsidT,_mm256_set1_pd(T_*T_));
    dpsidT = _mm256_sub_pd(dpsidT,_mm256_div_pd(expTc,C_));
    dFdT = _mm256_fmadd_pd(_mm256_div_pd(X,psi),dpsidT,_mm256_mul_pd(E_,_mm256_set1_pd(this->invT)));
    __m256d logTen = _mm256_set1_pd(std::log(10));
    __m256d invlogTenPr = _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(logTen,Pr));
    __m256d dlogPrdPr = _mm256_blendv_pd(_mm256_setzero_pd(),invlogTenPr,cmp);
    __m256d dXdPr = _mm256_mul_pd(_mm256_mul_pd(X,X),_mm256_set1_pd(-2));
    dXdPr = _mm256_mul_pd(dXdPr,logPr_);
    dXdPr = _mm256_mul_pd(dXdPr,dlogPrdPr);
    dFdPr = _mm256_mul_pd(_mm256_mul_pd(F,dXdPr),logpsi);
}






void OptReaction::Troe_F_3() const noexcept
{
    const size_t i = this->Troe.size()-3;

    const unsigned int j0 = this->Troe[i+0];
    const unsigned int j1 = this->Troe[i+1];
    const unsigned int j2 = this->Troe[i+2];

    const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
    const unsigned int m1 = j1 - this->Ikf[4] + this->Itbr[2];
    const unsigned int m2 = j2 - this->Ikf[4] + this->Itbr[2];

    const unsigned int k0 = j0 - this->Ikf[4];
    const unsigned int k1 = j1 - this->Ikf[4];
    const unsigned int k2 = j2 - this->Ikf[4];

    const double Kinf0 = this->Kf_[j0+this->offset_kinf];
    const double Kinf1 = this->Kf_[j1+this->offset_kinf];
    const double Kinf2 = this->Kf_[j2+this->offset_kinf];

    double M0 = this->tmp_M[m0];   
    double M1 = this->tmp_M[m1]; 
    double M2 = this->tmp_M[m2]; 

    const double K00 = this->Kf_[j0];
    const double K01 = this->Kf_[j1];
    const double K02 = this->Kf_[j2];

    const double Pr0 = K00*M0/Kinf0;  
    const double Pr1 = K01*M1/Kinf1;  
    const double Pr2 = K02*M2/Kinf2;

    const double logPr0 = std::log10(std::max(Pr0, small));
    const double logPr1 = std::log10(std::max(Pr1, small));
    const double logPr2 = std::log10(std::max(Pr2, small));

    const double expTTsss0 = this->tmp_Exp[i + 0 + this->nSpecies];
    const double expTTsss1 = this->tmp_Exp[i + 1 + this->nSpecies];
    const double expTTsss2 = this->tmp_Exp[i + 2 + this->nSpecies];

    const double expTTss0  = this->tmp_Exp[i + 0 + this->nSpecies + this->Troe.size()];
    const double expTTss1  = this->tmp_Exp[i + 1 + this->nSpecies + this->Troe.size()];
    const double expTTss2  = this->tmp_Exp[i + 2 + this->nSpecies + this->Troe.size()];

    const double expTTs0   = this->tmp_Exp[i + 0 + this->nSpecies + this->Troe.size()*2];
    const double expTTs1   = this->tmp_Exp[i + 1 + this->nSpecies + this->Troe.size()*2];
    const double expTTs2   = this->tmp_Exp[i + 2 + this->nSpecies + this->Troe.size()*2];

    const double Fcent0 =(1 - this->alpha_[i + 0])*expTTsss0 + this->alpha_[i + 0]*expTTs0 + expTTss0;
    const double Fcent1 =(1 - this->alpha_[i + 1])*expTTsss1 + this->alpha_[i + 1]*expTTs1 + expTTss1;
    const double Fcent2 =(1 - this->alpha_[i + 2])*expTTsss2 + this->alpha_[i + 2]*expTTs2 + expTTss2;

    const double logFcent0 = std::log10(std::max(Fcent0, small));
    const double logFcent1 = std::log10(std::max(Fcent1, small));
    const double logFcent2 = std::log10(std::max(Fcent2, small));

    const double c0 = -0.4 - 0.67*logFcent0;
    const double c1 = -0.4 - 0.67*logFcent1;
    const double c2 = -0.4 - 0.67*logFcent2;

    const double n0 = 0.75 - 1.27*logFcent0;
    const double n1 = 0.75 - 1.27*logFcent1;
    const double n2 = 0.75 - 1.27*logFcent2;

    const double x1_0 = n0 - 0.14*(logPr0 + c0);
    const double x1_1 = n1 - 0.14*(logPr1 + c1);
    const double x1_2 = n2 - 0.14*(logPr2 + c2);

    const double x2_0 = (logPr0 + c0)/x1_0;
    const double x2_1 = (logPr1 + c1)/x1_1;
    const double x2_2 = (logPr2 + c2)/x1_2;

    const double x3_0 = 1 + (x2_0*x2_0);
    const double x3_1 = 1 + (x2_1*x2_1);
    const double x3_2 = 1 + (x2_2*x2_2);

    const double x4_0 = logFcent0/x3_0;
    const double x4_1 = logFcent1/x3_1;
    const double x4_2 = logFcent2/x3_2;

    const double F0 = std::pow(10, x4_0);
    const double F1 = std::pow(10, x4_1);
    const double F2 = std::pow(10, x4_2);

    const double N0  = 1/(1+Pr0)*F0*K00;
    const double N1  = 1/(1+Pr1)*F1*K01;
    const double N2  = 1/(1+Pr2)*F2*K02;

    this->Kf_[j0] = k0<this->n_Fall_Off_Reaction ? M0*N0 : N0;   
    this->Kf_[j1] = k1<this->n_Fall_Off_Reaction ? M1*N1 : N1;
    this->Kf_[j2] = k2<this->n_Fall_Off_Reaction ? M2*N2 : N2;   
}
void OptReaction::Troe_F_2() const noexcept
{
    const size_t i = this->Troe.size()-2;
    
    const unsigned int j0 = this->Troe[i+0];
    const unsigned int j1 = this->Troe[i+1];

    const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
    const unsigned int m1 = j1 - this->Ikf[4] + this->Itbr[2];

    const unsigned int k0 = j0 - this->Ikf[4];
    const unsigned int k1 = j1 - this->Ikf[4];

    const double Kinf0 = this->Kf_[j0+this->offset_kinf];
    const double Kinf1 = this->Kf_[j1+this->offset_kinf];

    double M0 = this->tmp_M[m0];   
    double M1 = this->tmp_M[m1]; 

    const double K00 = this->Kf_[j0];
    const double K01 = this->Kf_[j1];

    const double Pr0 = K00*M0/Kinf0;  
    const double Pr1 = K01*M1/Kinf1;  

    const double logPr0 = std::log10(std::max(Pr0, small));
    const double logPr1 = std::log10(std::max(Pr1, small));

    const double expTTsss0 = this->tmp_Exp[i + 0 + this->nSpecies];
    const double expTTsss1 = this->tmp_Exp[i + 1 + this->nSpecies];

    const double expTTss0  = this->tmp_Exp[i + 0 + this->nSpecies + this->Troe.size()];
    const double expTTss1  = this->tmp_Exp[i + 1 + this->nSpecies + this->Troe.size()];

    const double expTTs0   = this->tmp_Exp[i + 0 + this->nSpecies + this->Troe.size()*2];
    const double expTTs1   = this->tmp_Exp[i + 1 + this->nSpecies + this->Troe.size()*2];

    const double Fcent0 =(1 - this->alpha_[i + 0])*expTTsss0 + this->alpha_[i + 0]*expTTs0 + expTTss0;
    const double Fcent1 =(1 - this->alpha_[i + 1])*expTTsss1 + this->alpha_[i + 1]*expTTs1 + expTTss1;

    const double logFcent0 = std::log10(std::max(Fcent0, small));
    const double logFcent1 = std::log10(std::max(Fcent1, small));

    const double c0 = -0.4 - 0.67*logFcent0;
    const double c1 = -0.4 - 0.67*logFcent1;

    const double n0 = 0.75 - 1.27*logFcent0;
    const double n1 = 0.75 - 1.27*logFcent1;

    const double x1_0 = n0 - 0.14*(logPr0 + c0);
    const double x1_1 = n1 - 0.14*(logPr1 + c1);

    const double x2_0 = (logPr0 + c0)/x1_0;
    const double x2_1 = (logPr1 + c1)/x1_1;

    const double x3_0 = 1 + (x2_0*x2_0);
    const double x3_1 = 1 + (x2_1*x2_1);

    const double x4_0 = logFcent0/x3_0;
    const double x4_1 = logFcent1/x3_1;

    const double F0 = std::pow(10, x4_0);
    const double F1 = std::pow(10, x4_1);

    const double N0  = 1/(1+Pr0)*F0*K00;
    const double N1  = 1/(1+Pr1)*F1*K01;

    this->Kf_[j0] = k0<this->n_Fall_Off_Reaction ? M0*N0 : N0;
    this->Kf_[j1] = k1<this->n_Fall_Off_Reaction ? M1*N1 : N1;
}

void OptReaction::Troe_F_1() const noexcept
{
    const size_t i = this->Troe.size()-1;
    const unsigned int j0 = this->Troe[i+0];
    const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
    const unsigned int k0 = j0 - this->Ikf[4];
    const double Kinf0 = this->Kf_[j0+this->offset_kinf];
    double M0 = this->tmp_M[m0];   
    const double K00 = this->Kf_[j0];
    const double Pr0 = K00*M0/Kinf0;  
    const double logPr0 = std::log10(std::max(Pr0, small));
    const double expTTsss0 = this->tmp_Exp[i + 0 + this->nSpecies];
    const double expTTss0  = this->tmp_Exp[i + 0 + this->nSpecies + this->Troe.size()];
    const double expTTs0   = this->tmp_Exp[i + 0 + this->nSpecies + this->Troe.size()*2];
    const double Fcent0 =(1 - this->alpha_[i + 0])*expTTsss0 + this->alpha_[i + 0]*expTTs0 + expTTss0;
    const double logFcent0 = std::log10(std::max(Fcent0, small));
    const double c0 = -0.4 - 0.67*logFcent0;
    const double n0 = 0.75 - 1.27*logFcent0;
    const double x1_0 = n0 - 0.14*(logPr0 + c0);
    const double x2_0 = (logPr0 + c0)/x1_0;
    const double x3_0 = 1 + (x2_0*x2_0);
    const double x4_0 = logFcent0/x3_0;
    const double F0 = std::pow(10, x4_0);
    const double N0  = 1/(1+Pr0)*F0*K00;
    this->Kf_[j0] = k0<this->n_Fall_Off_Reaction ? M0*N0 : N0;   

}
void OptReaction::Troe_Jac_1() const noexcept
{
    size_t i = this->Troe.size()-1;
    const unsigned int j = this->Troe[i];
    const unsigned int k = j - this->Ikf[4];
    const unsigned int m = j - this->Ikf[4] + this->Itbr[2];
    const double Kinf = this->Kf_[j+this->offset_kinf];
    const double dKinfdT = this->dKfdT_[j+this->offset_kinf];
    const double K0 = this->Kf_[j];
    const double M = this->tmp_M[m];
    const double Pr = K0*M/Kinf; 
    const double invKinf = 1.0/Kinf;
    const double logPr = std::log10(std::max(Pr, small));
    const double expTTsss = this->tmp_Exp[i+this->nSpecies];
    const double expTTss  = this->tmp_Exp[i+this->nSpecies+this->Troe.size()];
    const double expTTs   = this->tmp_Exp[i+this->nSpecies+this->Troe.size()*2];
    const double Fcent = (1 - this->alpha_[i])*expTTsss + this->alpha_[i]*expTTs + expTTss;
    const double logFcent = std::log10(std::max(Fcent, small));
    const double c = -0.4 - 0.67*logFcent;
    const double n = 0.75 - 1.27*logFcent;
    const double x1 = n - 0.14*(logPr + c);
    const double invx1 = 1.0/x1;
    const double x2 = (logPr + c)*invx1;
    const double x3 = 1 + (x2*x2);
    const double invx3 = 1.0/x3;
    const double x4 = logFcent*invx3;
    const double F = std::pow(10, x4);
    static const double logTen = std::log(10);
    const double dFcentdT = - (1 - this->alpha_[i])*this->invTsss_[i]*expTTsss
      - this->alpha_[i]*this->invTs_[i]*expTTs
      + this->Tss_[i]*invT*invT*expTTss;
    const double dlogFcentdT = Fcent >= small ? dFcentdT/Fcent/logTen : 0;
    const double dcdT = -0.67*dlogFcentdT;
    const double dndT = - 1.27*dlogFcentdT;
    const double dx1dT = dndT - 0.14*dcdT;
    const double dx2dT = (dcdT - x2*dx1dT)*invx1;
    const double dx3dT = 2*x2*dx2dT;
    const double dx4dT = (dlogFcentdT - x4*dx3dT)*invx3;
    const double dlogPrdPr = Pr >= small ? 1/(logTen*Pr) : 0;
    const double dx1dPr = -0.14*dlogPrdPr;
    const double dx2dPr = (dlogPrdPr - x2*dx1dPr)*invx1;
    const double dx3dPr = 2*x2*dx2dPr;
    const double dx4dPr = -x4*dx3dPr*invx3;
    const double dFdPr = logTen*F*dx4dPr;
    const double dK0dT =  this->dKfdT_[j]; 
    const double dPrdT = (M*dK0dT-Pr*dKinfdT)*invKinf;
    const double dFdT = logTen*F*dx4dT + dFdPr*dPrdT;
    const double dKdT   = k<this->n_Fall_Off_Reaction?Pr*dKinfdT:dK0dT;
    const double K      = k<this->n_Fall_Off_Reaction?Kinf      :K0;
    const double MM     = k<this->n_Fall_Off_Reaction?M         :1;
    const double invOneplusPr = 1.0/(1+Pr);
    const double KK     = k<this->n_Fall_Off_Reaction?1         :K0*invKinf;
    const double N1     = k<this->n_Fall_Off_Reaction?F*invOneplusPr  :-F*invOneplusPr;
    const double N2     = k<this->n_Fall_Off_Reaction?Pr*dFdPr  :dFdPr;
    const double N  = invOneplusPr*F*K0;
    this->dKfdT_[j] = F*invOneplusPr*dKdT + F*invOneplusPr*invOneplusPr*dPrdT*K + K0*invOneplusPr*dFdT*MM;
    this->dKfdC_[m] =  K0*invOneplusPr*KK*(N1 + N2); 
    this->Kf_[j] = k<this->n_Fall_Off_Reaction ? M*N : N;   
}

void OptReaction::Troe_Jac_2() const noexcept
{
    size_t i = this->Troe.size()-2;
    const unsigned int j0 = this->Troe[i+0];
    const unsigned int j1 = j0 + 1;
    const unsigned int k0 = j0 - this->Ikf[4];
    const unsigned int k1 = j1 - this->Ikf[4];
    const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
    const unsigned int m1 = j1 - this->Ikf[4] + this->Itbr[2];
    __m256d Kinf = _mm256_setr_pd(Kf_[j0+this->offset_kinf],Kf_[j1+this->offset_kinf],1,1);
    __m256d invKinf = _mm256_div_pd(_mm256_set1_pd(1.0),Kinf);
    __m256d dKinfdT = _mm256_setr_pd(dKfdT_[j0+this->offset_kinf],dKfdT_[j1+this->offset_kinf],1,1);
    __m256d K0 = _mm256_setr_pd(Kf_[j0],Kf_[j1],1,1);           
    __m256d M = _mm256_setr_pd(tmp_M[m0],tmp_M[m1],1,1);
    __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(M,K0),invKinf);
    __m256d small = _mm256_set1_pd(2.2e-16);
    __m256d cmp_result_Pr = _mm256_cmp_pd(Pr,small,_CMP_GE_OQ);
    Pr = _mm256_add_pd(Pr,_mm256_set1_pd(1e-100));
    const double invLog10 = 1.0/std::log(10);
    __m256d logPr_ = _mm256_mul_pd(vec256_logd(_mm256_max_pd(small,Pr)),_mm256_set1_pd(invLog10));
    __m256d InvTsss = _mm256_setr_pd(invTsss_[i+0],invTsss_[i+1],1,1);
    __m256d InvTs = _mm256_setr_pd(invTs_[i+0],invTs_[i+1],1,1);
    __m256d Tss = _mm256_setr_pd(Tss_[i+0],Tss_[i+1],1,1);
    __m256d expTTsss = _mm256_setr_pd(tmp_Exp[i+0+this->nSpecies],tmp_Exp[i+1+this->nSpecies],1,1);
    __m256d expTTss = _mm256_setr_pd(tmp_Exp[i+0+this->nSpecies+this->Troe.size()],tmp_Exp[i+1+this->nSpecies+this->Troe.size()],1,1);
    __m256d expTTs = _mm256_setr_pd(tmp_Exp[i+0+this->nSpecies+this->Troe.size()*2],tmp_Exp[i+1+this->nSpecies+this->Troe.size()*2],1,1);
    __m256d one = _mm256_set1_pd(1.0);
    __m256d alpha = _mm256_setr_pd(alpha_[i+0],alpha_[i+1],1,1);
    __m256d Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha),expTTsss);
    Fcent = _mm256_fmadd_pd(alpha,expTTs,Fcent);
    Fcent = _mm256_add_pd(expTTss,Fcent);
    __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),_mm256_set1_pd(invLog10));
    __m256d c = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(0.67),_mm256_set1_pd(0.4));
    __m256d n = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(-1.27),_mm256_set1_pd(0.75));
    __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(c,logPr_),_mm256_set1_pd(0.14),n);
    __m256d invx1 = _mm256_div_pd(one,x1);
    __m256d x2 = _mm256_mul_pd(_mm256_sub_pd(logPr_,c),invx1);
    __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
    __m256d invx3 = _mm256_div_pd(one,x3);
    __m256d x4 = _mm256_mul_pd(logFcent,invx3);
    __m256d  F = vec256_powd(_mm256_set1_pd(10),x4);
    __m256d logTen = _mm256_set1_pd(std::log(10));
    __m256d dFcentdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_sub_pd(alpha,one),InvTsss),expTTsss);
    dFcentdT = _mm256_sub_pd(dFcentdT,_mm256_mul_pd(_mm256_mul_pd(alpha,InvTs),expTTs));
    __m256d invT2 = _mm256_set1_pd(invT*invT);
    dFcentdT = _mm256_fmadd_pd(expTTss,_mm256_mul_pd(Tss,invT2),dFcentdT);
    __m256d cmp2 = _mm256_cmp_pd(Fcent,small,_CMP_GE_OQ);
    __m256d dlogFcentdT = _mm256_div_pd(_mm256_div_pd(dFcentdT,_mm256_max_pd(Fcent,small)),logTen);
    dlogFcentdT = _mm256_blendv_pd(_mm256_setzero_pd(), dlogFcentdT, cmp2);
    __m256d dcdT = _mm256_mul_pd(dlogFcentdT,_mm256_set1_pd(-0.67));
    __m256d dndT = _mm256_mul_pd(dlogFcentdT,_mm256_set1_pd(-1.27));
    __m256d dx1dT = _mm256_fmadd_pd(dcdT,_mm256_set1_pd(-0.14),dndT);
    __m256d dx2dT = _mm256_mul_pd(_mm256_sub_pd(dcdT,_mm256_mul_pd(x2,dx1dT)),invx1);
    __m256d dx3dT = _mm256_mul_pd(_mm256_mul_pd(x2,dx2dT),_mm256_set1_pd(2.0));
    __m256d dx4dT = _mm256_mul_pd(_mm256_sub_pd(dlogFcentdT,_mm256_mul_pd(x4,dx3dT)),invx3);
    __m256d dFdT = _mm256_mul_pd(logTen,_mm256_mul_pd(F,dx4dT));
    __m256d dlogPrdPr = _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(Pr,logTen));
    dlogPrdPr = _mm256_blendv_pd(_mm256_setzero_pd(), dlogPrdPr, cmp_result_Pr);
    __m256d dx1dPr = _mm256_mul_pd(dlogPrdPr,_mm256_set1_pd(-0.14));
    __m256d dx2dPr = _mm256_mul_pd(_mm256_sub_pd(dlogPrdPr,_mm256_mul_pd(x2,dx1dPr)),invx1);
    __m256d dx3dPr = _mm256_mul_pd(_mm256_mul_pd(x2,dx2dPr),_mm256_set1_pd(2.0));
    __m256d dx4dPr = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(-1.0),x4),dx3dPr),invx3);
    __m256d dFdPr  = _mm256_mul_pd(_mm256_mul_pd(logTen,F),dx4dPr);    
    __m256d dK0dT = _mm256_setr_pd(dKfdT_[j0],dKfdT_[j1],1,1);            
    __m256d dPrdT = _mm256_mul_pd(_mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT)),invKinf);
    dFdT = _mm256_fmadd_pd(dFdPr,dPrdT,dFdT);
    __m256d k = _mm256_setr_pd(k0,k1,1,1);
    __m256d cmp = _mm256_cmp_pd(k,_mm256_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
    __m256d dKdT = _mm256_blendv_pd(dK0dT, _mm256_mul_pd(Pr,dKinfdT), cmp);
    __m256d K = _mm256_blendv_pd(K0, Kinf, cmp);
    __m256d MM = _mm256_blendv_pd(_mm256_set1_pd(1), M, cmp);
    __m256d KK = _mm256_blendv_pd(_mm256_mul_pd(K0,invKinf), _mm256_set1_pd(1), cmp);
    __m256d invOnePlusPr = _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),Pr));
    __m256d N1 = _mm256_mul_pd(F,invOnePlusPr);
    N1 = _mm256_blendv_pd(_mm256_sub_pd(_mm256_setzero_pd(),N1),N1,cmp);            
    __m256d N2 = _mm256_blendv_pd(dFdPr,_mm256_mul_pd(Pr,dFdPr),cmp);
    __m256d N = _mm256_mul_pd(_mm256_mul_pd(F,K0),invOnePlusPr);
    __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(F,invOnePlusPr),dKdT);
    dKfdT = _mm256_fmadd_pd(K,_mm256_mul_pd(_mm256_mul_pd(F,_mm256_mul_pd(invOnePlusPr,invOnePlusPr)),dPrdT),dKfdT);
    dKfdT = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),dFdT),MM,dKfdT);
    this->dKfdT_[j0] = get_elem0(dKfdT);
    this->dKfdT_[j1] = get_elem1(dKfdT);
    __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),KK),_mm256_add_pd(N1,N2));
    this->dKfdC_[m0] = get_elem0(dKfdC);
    this->dKfdC_[m1] = get_elem1(dKfdC);
    __m256d KF = _mm256_blendv_pd(N,_mm256_mul_pd(N,M),cmp); 
    this->Kf_[j0] = get_elem0(KF);
    this->Kf_[j1] = get_elem1(KF);
}

void OptReaction::Troe_Jac_3() const noexcept
{

    size_t i = this->Troe.size() - 3;
    const unsigned int j0 = this->Troe[i+0];
    const unsigned int j1 = j0 + 1;
    const unsigned int j2 = j1 + 1;
    const unsigned int k0 = j0 - this->Ikf[4];
    const unsigned int k1 = j1 - this->Ikf[4];
    const unsigned int k2 = j2 - this->Ikf[4];
    const unsigned int m0 = j0 - this->Ikf[4] + this->Itbr[2];
    const unsigned int m1 = j1 - this->Ikf[4] + this->Itbr[2];
    const unsigned int m2 = j2 - this->Ikf[4] + this->Itbr[2];


    __m256d Kinf = _mm256_setr_pd(Kf_[j0+this->offset_kinf],Kf_[j1+this->offset_kinf],Kf_[j2+this->offset_kinf],1);
    __m256d invKinf = _mm256_div_pd(_mm256_set1_pd(1.0),Kinf);
    __m256d dKinfdT = _mm256_setr_pd(dKfdT_[j0+this->offset_kinf],dKfdT_[j1+this->offset_kinf],dKfdT_[j2+this->offset_kinf],1);
    __m256d K0 = _mm256_setr_pd(Kf_[j0],Kf_[j1],Kf_[j2],1);           
    __m256d M = _mm256_setr_pd(tmp_M[m0],tmp_M[m1],tmp_M[m2],1);
    __m256d Pr = _mm256_mul_pd(_mm256_mul_pd(M,K0),invKinf);
    __m256d small = _mm256_set1_pd(2.2e-16);
    __m256d cmp_result_Pr = _mm256_cmp_pd(Pr,small,_CMP_GE_OQ);
    Pr = _mm256_add_pd(Pr,_mm256_set1_pd(1e-100));
    const double invLog10 = 1.0/std::log(10);
    __m256d logPr_ = _mm256_mul_pd(vec256_logd(_mm256_max_pd(small,Pr)),_mm256_set1_pd(invLog10));
    __m256d InvTsss = _mm256_setr_pd(invTsss_[i+0],invTsss_[i+1],invTsss_[i+2],1);
    __m256d InvTs = _mm256_setr_pd(invTs_[i+0],invTs_[i+1],invTs_[i+2],1);
    __m256d Tss = _mm256_setr_pd(Tss_[i+0],Tss_[i+1],Tss_[i+2],1);
    __m256d expTTsss = _mm256_setr_pd(  tmp_Exp[i+0+this->nSpecies],tmp_Exp[i+1+this->nSpecies],tmp_Exp[i+2+this->nSpecies],1);
    __m256d expTTss = _mm256_setr_pd(   tmp_Exp[i+0+this->nSpecies+this->Troe.size()],tmp_Exp[i+1+this->nSpecies+this->Troe.size()],tmp_Exp[i+2+this->nSpecies+this->Troe.size()],1);
    __m256d expTTs = _mm256_setr_pd(    tmp_Exp[i+0+this->nSpecies+this->Troe.size()*2],tmp_Exp[i+1+this->nSpecies+this->Troe.size()*2],tmp_Exp[i+2+this->nSpecies+this->Troe.size()*2],1);
    __m256d one = _mm256_set1_pd(1.0);
    __m256d alpha = _mm256_setr_pd(alpha_[i+0],alpha_[i+1],alpha_[i+2],1);
    __m256d Fcent  = _mm256_mul_pd(_mm256_sub_pd(one,alpha),expTTsss);
    Fcent = _mm256_fmadd_pd(alpha,expTTs,Fcent);
    Fcent = _mm256_add_pd(expTTss,Fcent);
    __m256d logFcent = _mm256_mul_pd(vec256_logd(_mm256_max_pd(Fcent,small)),_mm256_set1_pd(invLog10));
    __m256d c = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(0.67),_mm256_set1_pd(0.4));
    __m256d n = _mm256_fmadd_pd(logFcent,_mm256_set1_pd(-1.27),_mm256_set1_pd(0.75));
    __m256d x1 = _mm256_fmadd_pd(_mm256_sub_pd(c,logPr_),_mm256_set1_pd(0.14),n);
    __m256d invx1 = _mm256_div_pd(one,x1);
    __m256d x2 = _mm256_mul_pd(_mm256_sub_pd(logPr_,c),invx1);
    __m256d x3 = _mm256_fmadd_pd(x2,x2,one);
    __m256d invx3 = _mm256_div_pd(one,x3);
    __m256d x4 = _mm256_mul_pd(logFcent,invx3);
    __m256d  F = vec256_powd(_mm256_set1_pd(10),x4);
    __m256d logTen = _mm256_set1_pd(std::log(10));
    __m256d dFcentdT = _mm256_mul_pd(_mm256_mul_pd(_mm256_sub_pd(alpha,one),InvTsss),expTTsss);
    dFcentdT = _mm256_sub_pd(dFcentdT,_mm256_mul_pd(_mm256_mul_pd(alpha,InvTs),expTTs));
    __m256d invT2 = _mm256_set1_pd(invT*invT);
    dFcentdT = _mm256_fmadd_pd(expTTss,_mm256_mul_pd(Tss,invT2),dFcentdT);
    __m256d cmp2 = _mm256_cmp_pd(Fcent,small,_CMP_GE_OQ);
    __m256d dlogFcentdT = _mm256_div_pd(_mm256_div_pd(dFcentdT,_mm256_max_pd(Fcent,small)),logTen);
    dlogFcentdT = _mm256_blendv_pd(_mm256_setzero_pd(), dlogFcentdT, cmp2);
    __m256d dcdT = _mm256_mul_pd(dlogFcentdT,_mm256_set1_pd(-0.67));
    __m256d dndT = _mm256_mul_pd(dlogFcentdT,_mm256_set1_pd(-1.27));
    __m256d dx1dT = _mm256_fmadd_pd(dcdT,_mm256_set1_pd(-0.14),dndT);
    __m256d dx2dT = _mm256_mul_pd(_mm256_sub_pd(dcdT,_mm256_mul_pd(x2,dx1dT)),invx1);
    __m256d dx3dT = _mm256_mul_pd(_mm256_mul_pd(x2,dx2dT),_mm256_set1_pd(2.0));
    __m256d dx4dT = _mm256_mul_pd(_mm256_sub_pd(dlogFcentdT,_mm256_mul_pd(x4,dx3dT)),invx3);
    __m256d dFdT = _mm256_mul_pd(logTen,_mm256_mul_pd(F,dx4dT));
    __m256d dlogPrdPr = _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_mul_pd(Pr,logTen));
    dlogPrdPr = _mm256_blendv_pd(_mm256_setzero_pd(), dlogPrdPr, cmp_result_Pr);
    __m256d dx1dPr = _mm256_mul_pd(dlogPrdPr,_mm256_set1_pd(-0.14));
    __m256d dx2dPr = _mm256_mul_pd(_mm256_sub_pd(dlogPrdPr,_mm256_mul_pd(x2,dx1dPr)),invx1);
    __m256d dx3dPr = _mm256_mul_pd(_mm256_mul_pd(x2,dx2dPr),_mm256_set1_pd(2.0));
    __m256d dx4dPr = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(-1.0),x4),dx3dPr),invx3);
    __m256d dFdPr  = _mm256_mul_pd(_mm256_mul_pd(logTen,F),dx4dPr);    
    __m256d dK0dT = _mm256_setr_pd(dKfdT_[j0],dKfdT_[j1],dKfdT_[j2],1);            
    __m256d dPrdT = _mm256_mul_pd(_mm256_fmsub_pd(M,dK0dT,_mm256_mul_pd(Pr,dKinfdT)),invKinf);
    dFdT = _mm256_fmadd_pd(dFdPr,dPrdT,dFdT);
    __m256d k = _mm256_setr_pd(k0,k1,k2,1);
    __m256d cmp = _mm256_cmp_pd(k,_mm256_set1_pd(this->n_Fall_Off_Reaction),_CMP_LT_OQ);
    __m256d dKdT = _mm256_blendv_pd(dK0dT, _mm256_mul_pd(Pr,dKinfdT), cmp);
    __m256d K = _mm256_blendv_pd(K0, Kinf, cmp);
    __m256d MM = _mm256_blendv_pd(_mm256_set1_pd(1), M, cmp);
    __m256d KK = _mm256_blendv_pd(_mm256_mul_pd(K0,invKinf), _mm256_set1_pd(1), cmp);
    __m256d invOnePlusPr = _mm256_div_pd(_mm256_set1_pd(1.0),_mm256_add_pd(_mm256_set1_pd(1.0),Pr));
    __m256d N1 = _mm256_mul_pd(F,invOnePlusPr);
    N1 = _mm256_blendv_pd(_mm256_sub_pd(_mm256_setzero_pd(),N1),N1,cmp);            
    __m256d N2 = _mm256_blendv_pd(dFdPr,_mm256_mul_pd(Pr,dFdPr),cmp);
    __m256d N = _mm256_mul_pd(_mm256_mul_pd(F,K0),invOnePlusPr);
    __m256d dKfdT = _mm256_mul_pd(_mm256_mul_pd(F,invOnePlusPr),dKdT);
    dKfdT = _mm256_fmadd_pd(K,_mm256_mul_pd(_mm256_mul_pd(F,_mm256_mul_pd(invOnePlusPr,invOnePlusPr)),dPrdT),dKfdT);
    dKfdT = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),dFdT),MM,dKfdT);
    this->dKfdT_[j0] = get_elem0(dKfdT);
    this->dKfdT_[j1] = get_elem1(dKfdT);
    this->dKfdT_[j2] = get_elem2(dKfdT);
    __m256d dKfdC = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(K0,invOnePlusPr),KK),_mm256_add_pd(N1,N2));
    this->dKfdC_[m0] = get_elem0(dKfdC);
    this->dKfdC_[m1] = get_elem1(dKfdC);
    this->dKfdC_[m2] = get_elem2(dKfdC);
    __m256d KF = _mm256_blendv_pd(N,_mm256_mul_pd(N,M),cmp);   
    this->Kf_[j0] = get_elem0(KF);
    this->Kf_[j1] = get_elem1(KF);
    this->Kf_[j2] = get_elem2(KF);
}
