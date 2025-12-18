void OptReaction::evalSRIRateConstant(double T)const noexcept
{

    unsigned int remain = this->n_SRI%4;
    __m256d invLog10v = _mm256_set1_pd(0.43429448190325182765112891891661);
    __m256d SRILimiterv = _mm256_set1_pd(SRILimiter);
    __m256d onev = _mm256_set1_pd(1);
    __m256d invTv = _mm256_set1_pd(this->invT);
    __m256d Tv = _mm256_set1_pd(T);
    __m256d logTv = _mm256_set1_pd(this->logT);

    for (unsigned int i = 0; i<this->n_SRI-remain; i=i+4)
    {
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];


        // Compute psi
        __m256d psiv = _mm256_setzero_pd();
        {
            __m256d av = _mm256_loadu_pd(&this->a_[i+0]);
            //const double a0 = this->a_[i+0];
            //const double a1 = this->a_[i+1];
            //const double a2 = this->a_[i+2];
            //const double a3 = this->a_[i+3];

            __m256d bv = _mm256_loadu_pd(&this->b_[i+0]);
            //const double b0 = -this->b_[i+0];
            //const double b1 = -this->b_[i+1];
            //const double b2 = -this->b_[i+2];
            //const double b3 = -this->b_[i+3];

            __m256d invcv = -_mm256_loadu_pd(&this->invc_[i+0]);
            //const double invc0 = -this->invc_[i+0];
            //const double invc1 = -this->invc_[i+1];
            //const double invc2 = -this->invc_[i+2];
            //const double invc3 = -this->invc_[i+3];

            __m256d r0 = _mm256_mul_pd(invcv,Tv);
            r0 = vec256_expd(r0);
            __m256d r1 = _mm256_mul_pd(bv,invTv);
            r1 = vec256_expd(r1);
            psiv = _mm256_fmadd_pd(av,r1,r0);
        }


        //const double Kinf0 = this->Kf_[j+0+this->offset_kinf];
        //const double Kinf1 = this->Kf_[j+1+this->offset_kinf];
        //const double Kinf2 = this->Kf_[j+2+this->offset_kinf];
        //const double Kinf3 = this->Kf_[j+3+this->offset_kinf];
        __m256d Kinf03v = _mm256_loadu_pd(&this->Kf_[j+0+this->offset_kinf]);

        //double M0 = this->tmp_M[m+0];
        //double M1 = this->tmp_M[m+1];
        //double M2 = this->tmp_M[m+2];
        //double M3 = this->tmp_M[m+3];
        __m256d M03v = _mm256_loadu_pd(&this->tmp_M[m+0]);

        //const double K0_0 = this->Kf_[j+0];
        //const double K0_1 = this->Kf_[j+1];
        //const double K0_2 = this->Kf_[j+2];
        //const double K0_3 = this->Kf_[j+3];
        __m256d K0v = _mm256_loadu_pd(&this->Kf_[j+0]);

        //const double Pr0 = K0_0*M0/Kinf0;
        //const double Pr1 = K0_1*M1/Kinf1;
        //const double Pr2 = K0_2*M2/Kinf2;
        //const double Pr3 = K0_3*M3/Kinf3;
        __m256d Prv = _mm256_mul_pd(K0v,M03v);
        Prv = _mm256_div_pd(Prv,Kinf03v);

        //double F0  = 0;
        //double F1  = 0;
        //double F2  = 0;
        //double F3  = 0;
        __m256d Fv = _mm256_setzero_pd();
        {
            __m256d limitedPrv = _mm256_max_pd(Prv,SRILimiterv);
            __m256d logPrv = vec256_logd(limitedPrv);
            logPrv = _mm256_mul_pd(logPrv,invLog10v);
            //const double logPr0 = std::log10(max(Pr0, SRILimiter));
            //const double logPr1 = std::log10(max(Pr1, SRILimiter));
            //const double logPr2 = std::log10(max(Pr2, SRILimiter));
            //const double logPr3 = std::log10(max(Pr3, SRILimiter));

            __m256d Xv = _mm256_mul_pd(logPrv,logPrv);
            Xv = _mm256_add_pd(onev,Xv);
            Xv = _mm256_div_pd(onev,Xv);
            //const double X0 = 1/(1 + (logPr0*logPr0));
            //const double X1 = 1/(1 + (logPr1*logPr1));
            //const double X2 = 1/(1 + (logPr2*logPr2));
            //const double X3 = 1/(1 + (logPr3*logPr3));

            __m256d r0 = vec256_powd(psiv,Xv);


            //const double psi0 = a0*std::exp(b0*invT) + std::exp(-T*invc0);
            //const double psi1 = a1*std::exp(b1*invT) + std::exp(-T*invc1);
            //const double psi2 = a2*std::exp(b2*invT) + std::exp(-T*invc2);
            //const double psi3 = a3*std::exp(b3*invT) + std::exp(-T*invc3);


            __m256d ev = _mm256_loadu_pd(&this->e_[i+0]);
            __m256d r1 = _mm256_mul_pd(ev,logTv);
            //const double e0 = this->e_[i+0];
            //const double e1 = this->e_[i+1];
            //const double e2 = this->e_[i+2];
            //const double e3 = this->e_[i+3];

            //const double logT = this->logT;


            //const double d0 = this->d_[i+0];
            //const double d1 = this->d_[i+1];
            //const double d2 = this->d_[i+2];
            //const double d3 = this->d_[i+3];

            __m256d dv = _mm256_loadu_pd(&this->d_[i+0]);
            Fv = _mm256_mul_pd(r0,vec256_expd(r1));
            Fv = _mm256_mul_pd(Fv,dv);
            //F0 = d0*std::pow(psi0, X0)*std::exp(e0*logT);
            //F1 = d1*std::pow(psi1, X1)*std::exp(e1*logT);
            //F2 = d2*std::pow(psi2, X2)*std::exp(e2*logT);
            //F3 = d3*std::pow(psi3, X3)*std::exp(e3*logT);
        }

        __m256d Nv = _mm256_mul_pd(Fv,K0v);
        Nv = _mm256_div_pd(onev,Nv);
        //const double N0 = 1/(1+Pr0)*F0*K0_0;
        //const double N1 = 1/(1+Pr1)*F1*K0_1;
        //const double N2 = 1/(1+Pr2)*F2*K0_2;
        //const double N3 = 1/(1+Pr3)*F3*K0_3;

        __m256d Kf = _mm256_mul_pd(M03v,Nv);

        _mm256_storeu_pd(&this->Kf_[j],Kf);
        //this->Kf_[j+0] = k+0<this->n_Fall_Off_Reaction ? M0*N0:N0;
        //this->Kf_[j+1] = k+1<this->n_Fall_Off_Reaction ? M1*N1:N1;
        //this->Kf_[j+2] = k+2<this->n_Fall_Off_Reaction ? M2*N2:N2;
        //this->Kf_[j+3] = k+3<this->n_Fall_Off_Reaction ? M3*N3:N3;
        std::cout<<this->Kf_[j+0]<<" "<<this->Kf_[j+1]<<" "<<this->Kf_[j+2]<<" "<<this->Kf_[j+3]<<std::endl;
        std::exit(0);
    }
    if(remain==1)
    {
        unsigned int i = this->n_SRI-1;
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];
        const unsigned int k = j - this->Ikf[4];
        
        const double Kinf0 = this->Kf_[j+0+this->offset_kinf];

        double M0 = this->tmp_M[m+0];

        const double K0_0 = this->Kf_[j+0];

        const double Pr0 = K0_0*M0/Kinf0;

        double F0  = 0;

        {
            const double logPr0 = std::log10(max(Pr0, small));

            const double X0 = 1/(1 + (logPr0*logPr0));

            const double invT = this->invT;

            const double a0 = this->a_[i+0];

            const double b0 = -this->b_[i+0];

            const double invc0 = -this->b_[i+0];


            const double psi0 = a0*std::exp(b0*invT) + std::exp(-T*invc0);


            const double d0 = this->d_[i+0];


            const double e0 = this->e_[i+0];


            const double logT = this->logT;
            F0 = d0*std::pow(psi0, X0)*std::exp(e0*logT);

        }

        const double N0 = 1/(1+Pr0)*F0*K0_0;


        this->Kf_[j+0] = k+0<this->n_Fall_Off_Reaction ? M0*N0:N0;
    }
    else if(remain==2)
    {
        unsigned int i = this->n_SRI-2;
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];
        const unsigned int k = j - this->Ikf[4];
        
        const double Kinf0 = this->Kf_[j+0+this->offset_kinf];
        const double Kinf1 = this->Kf_[j+1+this->offset_kinf];

        double M0 = this->tmp_M[m+0];
        double M1 = this->tmp_M[m+1];

        const double K0_0 = this->Kf_[j+0];
        const double K0_1 = this->Kf_[j+1];

        const double Pr0 = K0_0*M0/Kinf0;
        const double Pr1 = K0_1*M1/Kinf1;

        double F0  = 0;
        double F1  = 0;

        {
            const double logPr0 = std::log10(max(Pr0, small));
            const double logPr1 = std::log10(max(Pr1, small));

            const double X0 = 1/(1 + (logPr0*logPr0));
            const double X1 = 1/(1 + (logPr1*logPr1));

            const double invT = this->invT;


            const double a0 = this->a_[i+0];
            const double a1 = this->a_[i+1];

            const double b0 = -this->b_[i+0];
            const double b1 = -this->b_[i+1];

            const double invc0 = -this->b_[i+0];
            const double invc1 = -this->b_[i+1];

            const double psi0 = a0*std::exp(b0*invT) + std::exp(-T*invc0);
            const double psi1 = a1*std::exp(b1*invT) + std::exp(-T*invc1);

            const double d0 = this->d_[i+0];
            const double d1 = this->d_[i+1];

            const double e0 = this->e_[i+0];
            const double e1 = this->e_[i+1];

            const double logT = this->logT;
            F0 = d0*std::pow(psi0, X0)*std::exp(e0*logT);
            F1 = d1*std::pow(psi1, X1)*std::exp(e1*logT);

        }

        const double N0 = 1/(1+Pr0)*F0*K0_0;
        const double N1 = 1/(1+Pr1)*F1*K0_1;

        this->Kf_[j+0] = k+0<this->n_Fall_Off_Reaction ? M0*N0:N0;
        this->Kf_[j+1] = k+1<this->n_Fall_Off_Reaction ? M1*N1:N1;
    }
    else if(remain==3)
    {
        unsigned int i = this->n_SRI-3;
        const unsigned int j = this->SRIFO[i];
        const unsigned int m = j - this->Ikf[4] + this->Itbr[2];
        const unsigned int k = j - this->Ikf[4];
        
        const double Kinf0 = this->Kf_[j+0+this->offset_kinf];
        const double Kinf1 = this->Kf_[j+1+this->offset_kinf];
        const double Kinf2 = this->Kf_[j+2+this->offset_kinf];

        double M0 = this->tmp_M[m+0];
        double M1 = this->tmp_M[m+1];
        double M2 = this->tmp_M[m+2];

        const double K0_0 = this->Kf_[j+0];
        const double K0_1 = this->Kf_[j+1];
        const double K0_2 = this->Kf_[j+2];

        const double Pr0 = K0_0*M0/Kinf0;
        const double Pr1 = K0_1*M1/Kinf1;
        const double Pr2 = K0_2*M2/Kinf2;

        double F0  = 0;
        double F1  = 0;
        double F2  = 0;

        {
            const double logPr0 = std::log10(max(Pr0, small));
            const double logPr1 = std::log10(max(Pr1, small));
            const double logPr2 = std::log10(max(Pr2, small));


            const double X0 = 1/(1 + (logPr0*logPr0));
            const double X1 = 1/(1 + (logPr1*logPr1));
            const double X2 = 1/(1 + (logPr2*logPr2));
            const double invT = this->invT;


            const double a0 = this->a_[i+0];
            const double a1 = this->a_[i+1];
            const double a2 = this->a_[i+2];

            const double b0 = -this->b_[i+0];
            const double b1 = -this->b_[i+1];
            const double b2 = -this->b_[i+2];

            const double invc0 = -this->b_[i+0];
            const double invc1 = -this->b_[i+1];
            const double invc2 = -this->b_[i+2];

            const double psi0 = a0*std::exp(b0*invT) + std::exp(-T*invc0);
            const double psi1 = a1*std::exp(b1*invT) + std::exp(-T*invc1);
            const double psi2 = a2*std::exp(b2*invT) + std::exp(-T*invc2);

            const double d0 = this->d_[i+0];
            const double d1 = this->d_[i+1];
            const double d2 = this->d_[i+2];

            const double e0 = this->e_[i+0];
            const double e1 = this->e_[i+1];
            const double e2 = this->e_[i+2];


            const double logT = this->logT;
            F0 = d0*std::pow(psi0, X0)*std::exp(e0*logT);
            F1 = d1*std::pow(psi1, X1)*std::exp(e1*logT);
            F2 = d2*std::pow(psi2, X2)*std::exp(e2*logT);

        }

        const double N0 = 1/(1+Pr0)*F0*K0_0;
        const double N1 = 1/(1+Pr1)*F1*K0_1;
        const double N2 = 1/(1+Pr2)*F2*K0_2;


        this->Kf_[j+0] = k+0<this->n_Fall_Off_Reaction ? M0*N0:N0;
        this->Kf_[j+1] = k+1<this->n_Fall_Off_Reaction ? M1*N1:N1;
        this->Kf_[j+2] = k+2<this->n_Fall_Off_Reaction ? M2*N2:N2;
    }

}