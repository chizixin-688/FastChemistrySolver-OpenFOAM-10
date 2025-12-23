#include "OptReaction.H"
#include <immintrin.h>  

void FastChemistry::OptReaction::findPlogPressureRange(double p)const noexcept
{
    int k = 0;
    for(unsigned int i = 0; i< this->n_PlogReaction; i++)
    {
        const size_t length = this->Prange[i].size();
        if(p<=this->Prange[i][0])
        {
            double A0 = this->APlog[i][0];
            double beta0 = this->betaPlog[i][0];
            double Ta0 = this->TaPlog[i][0];

            this->A[i+this->Ikf[6]] = A0;
            this->A[i+this->Ikf[11]] = A0;
            this->beta[i+this->Ikf[6]] = beta0;
            this->beta[i+this->Ikf[11]] = beta0;
            this->Ta[i+this->Ikf[6]] = Ta0;
            this->Ta[i+this->Ikf[11]] = Ta0;
            this->Pindex[i] = 0;
        }
        else if(p>=this->Prange[i][length-1])
        {
            double A1 = this->APlog[i][length-1];
            double beta1 = this->betaPlog[i][length-1];
            double Ta1 = this->TaPlog[i][length-1];
                
            this->A[i+this->Ikf[6]] = A1;
            this->A[i+this->Ikf[11]] = A1;
            this->beta[i+this->Ikf[6]] = beta1;
            this->beta[i+this->Ikf[11]] = beta1;
            this->Ta[i+this->Ikf[6]] = Ta1;
            this->Ta[i+this->Ikf[11]] = Ta1;
            this->Pindex[i] = static_cast<unsigned int>(length-1);
        }
        else
        {
            this->ActivePlogReactionIndex[k] = i;
            k++;
            unsigned int index = 0;
            for(unsigned int j = 0; j < length-1;j++)
            {
                if(this->Prange[i][j]<=p && p<this->Prange[i][j+1])
                {
                    index = j;
                    break;
                }
            }
            this->Pindex[i] = index;
            this->A[i+this->Ikf[6]] = this->APlog[i][index+0];
            this->A[i+this->Ikf[11]] = this->APlog[i][index+1];
            this->beta[i+this->Ikf[6]] = this->betaPlog[i][index+0];
            this->beta[i+this->Ikf[11]] = this->betaPlog[i][index+1];
            this->Ta[i+this->Ikf[6]] = this->TaPlog[i][index+0];
            this->Ta[i+this->Ikf[11]] = this->TaPlog[i][index+1];
        }
    }
    this->ActivePlogReactionNumber = k;
}