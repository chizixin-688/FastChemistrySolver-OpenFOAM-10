/*---------------------------------------------------------------------------*\
  Description
      Computing the forward rate constant and the partial derivatives of Plog
      reactions. Including Kf and dKfdT

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"

//=============================================================================//

void FastChemistry::OptReaction::evalPlogPartialDerivative()const noexcept
{
    const double logT_ = this->logT;
    const double invT_ = this->invT;
    __m256d logTv = _mm256_set1_pd(this->logT);
    __m256d invTv = _mm256_set1_pd(this->invT);
    const unsigned int remain = this->ActivePlogReactionNumber%4;
    for(unsigned int i=0; i<this->ActivePlogReactionNumber-remain; i=i+4)
    {
        unsigned j0 = ActivePlogReactionIndex[i+0];
        unsigned index0 = this->Pindex[j0];
        double weight0 = (this->logP - this->logPi[j0][index0])*this->rDeltaP_[j0][index0];
        double Kf00 = this->Kf_[j0+this->Ikf[6]];



        unsigned j1 = ActivePlogReactionIndex[i+1];
        unsigned index1 = this->Pindex[j1];
        double weight1 = (this->logP - this->logPi[j1][index1])*this->rDeltaP_[j1][index1];
        double Kf01 = this->Kf_[j1+this->Ikf[6]];


        unsigned j2 = ActivePlogReactionIndex[i+2];
        unsigned index2 = this->Pindex[j2];
        double weight2 = (this->logP - this->logPi[j2][index2])*this->rDeltaP_[j2][index2];
        double Kf02 = this->Kf_[j2+this->Ikf[6]];


        unsigned j3 = ActivePlogReactionIndex[i+3];
        unsigned index3 = this->Pindex[j3];
        double weight3 = (this->logP - this->logPi[j3][index3])*this->rDeltaP_[j3][index3];
        double Kf03 = this->Kf_[j3+this->Ikf[6]];


        __m256d weightv = _mm256_setr_pd(weight0,weight2,weight1,weight3);
        __m256d kf0v = _mm256_setr_pd(Kf00,Kf02,Kf01,Kf03);
        
        __m256d tmp0v = _mm256_setzero_pd();
        __m128d a0 = _mm_loadu_pd(&this->logAPlog[j0][index0]);
        __m128d a1 = _mm_loadu_pd(&this->logAPlog[j1][index1]);
        tmp0v = _mm256_set_m128d(a1, a0);//a_logA0 a_logA1 b_logA0 b_logA1 
            
        __m256d tmp1v = _mm256_setzero_pd();
        __m128d a2 = _mm_loadu_pd(&this->logAPlog[j2][index2]);
        __m128d a3 = _mm_loadu_pd(&this->logAPlog[j3][index3]);
        tmp1v = _mm256_set_m128d(a3, a2);//c_logA0 c_logA1 d_logA0 d_logA1 
            

        __m256d tmp2v = _mm256_setzero_pd();
        __m128d a4 = _mm_loadu_pd(&this->betaPlog[j0][index0]);
        __m128d a5 = _mm_loadu_pd(&this->betaPlog[j1][index1]);
        tmp2v = _mm256_set_m128d(a5, a4);//a_logA0 a_logA1 b_logA0 b_logA1 
            
        __m256d tmp3v = _mm256_setzero_pd();
        __m128d a6 = _mm_loadu_pd(&this->betaPlog[j2][index2]);
        __m128d a7 = _mm_loadu_pd(&this->betaPlog[j3][index3]);
        tmp3v = _mm256_set_m128d(a7, a6);//a_logA0 a_logA1 b_logA0 b_logA1 
            

        __m256d tmp4v = _mm256_setzero_pd();
        __m128d a8 = _mm_loadu_pd(&this->TaPlog[j0][index0]);
        __m128d a9 = _mm_loadu_pd(&this->TaPlog[j1][index1]);
        tmp4v = _mm256_set_m128d(a9, a8);//a_logA0 a_logA1 b_logA0 b_logA1 
            
        __m256d tmp5v = _mm256_setzero_pd();
        __m128d a10 = _mm_loadu_pd(&this->TaPlog[j2][index2]);
        __m128d a11 = _mm_loadu_pd(&this->TaPlog[j3][index3]);
        tmp5v = _mm256_set_m128d(a11, a10);//a_logA0 a_logA1 b_logA0 b_logA1 
            
            //a_logA0 a_logA1 b_logA0 b_logA1 
            //c_logA0 c_logA1 d_logA0 d_logA1 

            // a_logA0 c_logA0 b_logA0 d_logA0
        __m256d logA0v = _mm256_shuffle_pd (tmp0v, tmp1v, 0b0000);

            // a_logA1 c_logA1 b_logA1 d_logA1
        __m256d logA1v = _mm256_shuffle_pd (tmp0v, tmp1v, 0b1111);

            // a_beta0 c_beta0 b_beta0 d_beta0
        __m256d beta0v = _mm256_shuffle_pd (tmp2v, tmp3v, 0b0000);

            // a_beta1 c_beta1 b_beta1 d_beta1
        __m256d beta1v = _mm256_shuffle_pd (tmp2v, tmp3v, 0b1111);

            // a_Ta0 c_Ta0 b_Ta0 d_Ta0
        __m256d Ta0v = _mm256_shuffle_pd (tmp4v, tmp5v, 0b0000);

            // a_Ta1 c_Ta1 b_Ta1 d_Ta1
        __m256d Ta1v = _mm256_shuffle_pd (tmp4v, tmp5v, 0b1111);


        __m256d logA1A0v = _mm256_sub_pd(logA1v,logA0v);//logA1 - logA0
        __m256d beta1beta0v = _mm256_sub_pd(beta1v,beta0v);//beta1 - beta0
        __m256d Ta0Ta1v = _mm256_sub_pd(Ta0v,Ta1v);  //Ta0 - Ta1

            //double logk1k00 = logA1A00 + beta1beta00*logT_ + Ta0Ta10*invT_;
            //double logk1k03 = logA1A03 + beta1beta03*logT_ + Ta0Ta13*invT_;
        __m256d logk1k0v = _mm256_fmadd_pd(beta1beta0v,logTv,logA1A0v);
        logk1k0v = _mm256_fmadd_pd(Ta0Ta1v,invTv,logk1k0v);
        logk1k0v = _mm256_mul_pd(logk1k0v,weightv);
        logk1k0v = vec256_expd(logk1k0v);
        __m256d kfv = _mm256_mul_pd(kf0v,logk1k0v);
            //double k3 = Kf03*std::exp(weight3*logk1k03);

        double kfa = this->get_elem0(kfv);
        double kfb = this->get_elem2(kfv);
        double kfc = this->get_elem1(kfv);
        double kfd = this->get_elem3(kfv);
        this->Kf_[j0+this->Ikf[6]] = kfa;
        this->Kf_[j1+this->Ikf[6]] = kfb;
        this->Kf_[j2+this->Ikf[6]] = kfc;
        this->Kf_[j3+this->Ikf[6]] = kfd;

                //double beta0 = this->beta[i+this->Ikf[6]];
                //double beta1 = this->beta[i+this->Ikf[11]];
                //double Ta0 = this->Ta[i+this->Ikf[6]];
                //double Ta1 = this->Ta[i+this->Ikf[11]];
                //double invt = this->invT;
                //double dKfdT = Kf*invt*(beta0 + Ta0*invt + (beta1-beta0+(Ta1-Ta0)*invt)*weight);
                //this->dKfdT_[i+this->Ikf[6]] = dKfdT;

        __m256d dKfdTv = _mm256_fmadd_pd(-Ta0Ta1v,invTv,beta1beta0v);
        __m256d beta0Ta0invTv = _mm256_fmadd_pd(invTv,Ta0v,beta0v);
        dKfdTv = _mm256_fmadd_pd(dKfdTv,weightv,beta0Ta0invTv);
        __m256d KfinvTv = _mm256_mul_pd(kfv,invTv);
        dKfdTv = _mm256_mul_pd(dKfdTv,KfinvTv);
        double dkfdTa = this->get_elem0(dKfdTv);
        double dkfdTb = this->get_elem2(dKfdTv);
        double dkfdTc = this->get_elem1(dKfdTv);
        double dkfdTd = this->get_elem3(dKfdTv);
        this->dKfdT_[j0+this->Ikf[6]] = dkfdTa;
        this->dKfdT_[j1+this->Ikf[6]] = dkfdTb;
        this->dKfdT_[j2+this->Ikf[6]] = dkfdTc;
        this->dKfdT_[j3+this->Ikf[6]] = dkfdTd;

    }
    if(remain==1)
    {
        unsigned int i=this->ActivePlogReactionNumber-1;
        unsigned j0 = ActivePlogReactionIndex[i+0];
        unsigned index0 = this->Pindex[j0];
        double weight0 = (this->logP - this->logPi[j0][index0])*this->rDeltaP_[j0][index0];
        double Kf00 = this->Kf_[j0+this->Ikf[6]];
        double logA1A00 = this->logAPlog[j0][index0+1] - this->logAPlog[j0][index0];
        double beta1beta00 = this->betaPlog[j0][index0+1] - this->betaPlog[j0][index0];
        double Ta0Ta10 = this->TaPlog[j0][index0] - this->TaPlog[j0][index0+1];
        double logk1k00 = logA1A00 + beta1beta00*logT_ + Ta0Ta10*invT_;
        double k0 = Kf00*std::exp(weight0*logk1k00);
        this->Kf_[j0+this->Ikf[6]] = k0;

                //double beta0 = this->beta[i+this->Ikf[6]];
                //double beta1 = this->beta[i+this->Ikf[11]];
                //double Ta0 = this->Ta[i+this->Ikf[6]];
                //double Ta1 = this->Ta[i+this->Ikf[11]];
                //double invt = this->invT;
                //double dKfdT = Kf*invt*(beta0 + Ta0*invt + (beta1-beta0+(Ta1-Ta0)*invt)*weight);
                //this->dKfdT_[i+this->Ikf[6]] = dKfdT;
        double beta0 = this->betaPlog[j0][index0];
        double Ta0 = this->TaPlog[j0][index0];
        double dKfdT = k0*invT_*(beta0 + Ta0*invT_ + (beta1beta00-(Ta0Ta10)*invT_)*weight0);
        this->dKfdT_[j0+this->Ikf[6]] = dKfdT;
    }
    else if(remain==2)
    {
        unsigned int i=this->ActivePlogReactionNumber-2;
        unsigned j0 = ActivePlogReactionIndex[i+0];
        unsigned j1 = ActivePlogReactionIndex[i+1];
        unsigned index0 = this->Pindex[j0];
        unsigned index1 = this->Pindex[j1];




        __m128d a0 = _mm_loadu_pd(&this->logAPlog[j0][index0]);
        __m128d a1 = _mm_loadu_pd(&this->logAPlog[j1][index1]);
        __m256d  tmp0v = _mm256_set_m128d(a1, a0);//a_logA0 a_logA1 b_logA0 b_logA1 
            // a_logA0 a_logA0 b_logA0 b_logA0
        __m256d logA0v = _mm256_shuffle_pd (tmp0v, tmp0v, 0b0000);
            // a_logA1 a_logA1 b_logA1 b_logA1
        __m256d logA1v = _mm256_shuffle_pd (tmp0v, tmp0v, 0b1111);


        __m128d a2 = _mm_loadu_pd(&this->betaPlog[j0][index0]);
        __m128d a3 = _mm_loadu_pd(&this->betaPlog[j1][index1]);
        __m256d  tmp1v = _mm256_set_m128d(a3, a2);//a_logA0 a_logA1 b_logA0 b_logA1 
            // a_logA0 a_logA0 b_logA0 b_logA0
        __m256d beta0v = _mm256_shuffle_pd (tmp1v, tmp1v, 0b0000);
            // a_logA1 a_logA1 b_logA1 b_logA1
        __m256d beta1v = _mm256_shuffle_pd (tmp1v, tmp1v, 0b1111);



        __m128d a4 = _mm_loadu_pd(&this->TaPlog[j0][index0]);
        __m128d a5 = _mm_loadu_pd(&this->TaPlog[j1][index1]);
        __m256d  tmp2v = _mm256_set_m128d(a5, a4);//a_logA0 a_logA1 b_logA0 b_logA1 
            // a_logA0 a_logA0 b_logA0 b_logA0
        __m256d Ta0v = _mm256_shuffle_pd (tmp2v, tmp2v, 0b0000);
            // a_logA1 a_logA1 b_logA1 b_logA1
        __m256d Ta1v = _mm256_shuffle_pd (tmp2v, tmp2v, 0b1111);

        double weight0 = (this->logP - this->logPi[j0][index0])*this->rDeltaP_[j0][index0];
        double weight1 = (this->logP - this->logPi[j1][index1])*this->rDeltaP_[j1][index1];
        double Kf00 = this->Kf_[j0+this->Ikf[6]];
        double Kf01 = this->Kf_[j1+this->Ikf[6]];
        __m256d weightv = _mm256_setr_pd(weight0,weight0,weight1,weight1);
        __m256d kf0v = _mm256_setr_pd(Kf00,Kf00,Kf01,Kf01);

        __m256d logA1A0v = _mm256_sub_pd(logA1v,logA0v);
        __m256d beta1beta0v = _mm256_sub_pd(beta1v,beta0v);
        __m256d Ta0Ta1v = _mm256_sub_pd(Ta0v,Ta1v);

        __m256d logk1k0v = _mm256_fmadd_pd(beta1beta0v,logTv,logA1A0v);
        logk1k0v = _mm256_fmadd_pd(Ta0Ta1v,invTv,logk1k0v);
        logk1k0v = _mm256_mul_pd(logk1k0v,weightv);
        logk1k0v = vec256_expd(logk1k0v);
        __m256d kfv = _mm256_mul_pd(kf0v,logk1k0v);

        double kfa = this->get_elem0(kfv);
        double kfb = this->get_elem2(kfv);

        this->Kf_[j0+this->Ikf[6]] = kfa;
        this->Kf_[j1+this->Ikf[6]] = kfb;

        __m256d dKfdTv = _mm256_fmadd_pd(-Ta0Ta1v,invTv,beta1beta0v);
        __m256d beta0Ta0invTv = _mm256_fmadd_pd(invTv,Ta0v,beta0v);
        dKfdTv = _mm256_fmadd_pd(dKfdTv,weightv,beta0Ta0invTv);
        __m256d KfinvTv = _mm256_mul_pd(kfv,invTv);
        dKfdTv = _mm256_mul_pd(dKfdTv,KfinvTv);
        double dkfdTa = this->get_elem0(dKfdTv);
        double dkfdTb = this->get_elem2(dKfdTv);

        this->dKfdT_[j0+this->Ikf[6]] = dkfdTa;
        this->dKfdT_[j1+this->Ikf[6]] = dkfdTb;

    }
    else if(remain==3)
    {
        unsigned int i=this->ActivePlogReactionNumber-3;


        unsigned j0 = ActivePlogReactionIndex[i+0];
        unsigned index0 = this->Pindex[j0];
        double weight0 = (this->logP - this->logPi[j0][index0])*this->rDeltaP_[j0][index0];
        double Kf00 = this->Kf_[j0+this->Ikf[6]];



        unsigned j1 = ActivePlogReactionIndex[i+1];
        unsigned index1 = this->Pindex[j1];
        double weight1 = (this->logP - this->logPi[j1][index1])*this->rDeltaP_[j1][index1];
        double Kf01 = this->Kf_[j1+this->Ikf[6]];


        unsigned j2 = ActivePlogReactionIndex[i+2];
        unsigned index2 = this->Pindex[j2];
        double weight2 = (this->logP - this->logPi[j2][index2])*this->rDeltaP_[j2][index2];
        double Kf02 = this->Kf_[j2+this->Ikf[6]];

        __m256d weightv = _mm256_setr_pd(weight0,weight2,weight1,weight2);
        __m256d kf0v = _mm256_setr_pd(Kf00,Kf02,Kf01,Kf02);
        
        __m256d tmp0v = _mm256_setzero_pd();
        __m128d a0 = _mm_loadu_pd(&this->logAPlog[j0][index0]);
        __m128d a1 = _mm_loadu_pd(&this->logAPlog[j1][index1]);
        tmp0v = _mm256_set_m128d(a1, a0);//a_logA0 a_logA1 b_logA0 b_logA1 
            
        __m256d tmp1v = _mm256_setzero_pd();
        __m128d a2 = _mm_loadu_pd(&this->logAPlog[j2][index2]);
        //__m128d a3 = _mm_loadu_pd(&this->logAPlog[j2][index2]);
        tmp1v = _mm256_set_m128d(a2, a2);//c_logA0 c_logA1 d_logA0 d_logA1 
            

        __m256d tmp2v = _mm256_setzero_pd();
        __m128d a4 = _mm_loadu_pd(&this->betaPlog[j0][index0]);
        __m128d a5 = _mm_loadu_pd(&this->betaPlog[j1][index1]);
        tmp2v = _mm256_set_m128d(a5, a4);//a_logA0 a_logA1 b_logA0 b_logA1 
            
        __m256d tmp3v = _mm256_setzero_pd();
        __m128d a6 = _mm_loadu_pd(&this->betaPlog[j2][index2]);
        //__m128d a7 = _mm_loadu_pd(&this->betaPlog[j2][index2]);
        tmp3v = _mm256_set_m128d(a6, a6);//a_logA0 a_logA1 b_logA0 b_logA1 
            

        __m256d tmp4v = _mm256_setzero_pd();
        __m128d a8 = _mm_loadu_pd(&this->TaPlog[j0][index0]);
        __m128d a9 = _mm_loadu_pd(&this->TaPlog[j1][index1]);
        tmp4v = _mm256_set_m128d(a9, a8);//a_logA0 a_logA1 b_logA0 b_logA1 
            
        __m256d tmp5v = _mm256_setzero_pd();
        __m128d a10 = _mm_loadu_pd(&this->TaPlog[j2][index2]);
        //__m128d a11 = _mm_loadu_pd(&this->TaPlog[j3][index3]);
        tmp5v = _mm256_set_m128d(a10, a10);//a_logA0 a_logA1 b_logA0 b_logA1 
            
            //a_logA0 a_logA1 b_logA0 b_logA1 
            //c_logA0 c_logA1 d_logA0 d_logA1 

            // a_logA0 c_logA0 b_logA0 d_logA0
        __m256d logA0v = _mm256_shuffle_pd (tmp0v, tmp1v, 0b0000);

            // a_logA1 c_logA1 b_logA1 d_logA1
        __m256d logA1v = _mm256_shuffle_pd (tmp0v, tmp1v, 0b1111);

            // a_beta0 c_beta0 b_beta0 d_beta0
        __m256d beta0v = _mm256_shuffle_pd (tmp2v, tmp3v, 0b0000);

            // a_beta1 c_beta1 b_beta1 d_beta1
        __m256d beta1v = _mm256_shuffle_pd (tmp2v, tmp3v, 0b1111);

            // a_Ta0 c_Ta0 b_Ta0 d_Ta0
        __m256d Ta0v = _mm256_shuffle_pd (tmp4v, tmp5v, 0b0000);

            // a_Ta1 c_Ta1 b_Ta1 d_Ta1
        __m256d Ta1v = _mm256_shuffle_pd (tmp4v, tmp5v, 0b1111);


        __m256d logA1A0v = _mm256_sub_pd(logA1v,logA0v);
        __m256d beta1beta0v = _mm256_sub_pd(beta1v,beta0v);
        __m256d Ta0Ta1v = _mm256_sub_pd(Ta0v,Ta1v);

            //double logk1k00 = logA1A00 + beta1beta00*logT_ + Ta0Ta10*invT_;
            //double logk1k03 = logA1A03 + beta1beta03*logT_ + Ta0Ta13*invT_;
        __m256d logk1k0v = _mm256_fmadd_pd(beta1beta0v,logTv,logA1A0v);
        logk1k0v = _mm256_fmadd_pd(Ta0Ta1v,invTv,logk1k0v);
        logk1k0v = _mm256_mul_pd(logk1k0v,weightv);
        logk1k0v = vec256_expd(logk1k0v);
        __m256d kfv = _mm256_mul_pd(kf0v,logk1k0v);
            //double k3 = Kf03*std::exp(weight3*logk1k03);

        double kfa = this->get_elem0(kfv);
        double kfb = this->get_elem2(kfv);
        double kfc = this->get_elem1(kfv);
        this->Kf_[j0+this->Ikf[6]] = kfa;
        this->Kf_[j1+this->Ikf[6]] = kfb;
        this->Kf_[j2+this->Ikf[6]] = kfc;

        __m256d dKfdTv = _mm256_fmadd_pd(-Ta0Ta1v,invTv,beta1beta0v);
        __m256d beta0Ta0invTv = _mm256_fmadd_pd(invTv,Ta0v,beta0v);
        dKfdTv = _mm256_fmadd_pd(dKfdTv,weightv,beta0Ta0invTv);
        __m256d KfinvTv = _mm256_mul_pd(kfv,invTv);
        dKfdTv = _mm256_mul_pd(dKfdTv,KfinvTv);
        double dkfdTa = this->get_elem0(dKfdTv);
        double dkfdTb = this->get_elem2(dKfdTv);
        double dkfdTc = this->get_elem1(dKfdTv);
        this->dKfdT_[j0+this->Ikf[6]] = dkfdTa;
        this->dKfdT_[j1+this->Ikf[6]] = dkfdTb;
        this->dKfdT_[j2+this->Ikf[6]] = dkfdTc;
    }
}
