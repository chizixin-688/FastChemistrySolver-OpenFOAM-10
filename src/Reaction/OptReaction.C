/*---------------------------------------------------------------------------*\
  Description
      Implementation file for OptReaction class. This file contains the 
      constructor and destructor.

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "OptReaction.H"

//=============================================================================//

FastChemistry::OptReaction::OptReaction()
{}


FastChemistry::OptReaction::~OptReaction()
{
    if(this->tmp_Exp!=nullptr)
    {
        free(this->tmp_Exp);
    }
    if(this->negGstdByRT!=nullptr)
    {
        free(this->negGstdByRT);
    }
    if(this->invNegGstdByRT!=nullptr)
    {
        free(this->invNegGstdByRT);
    }
    if(this->ThirdBodyFactor1D!=nullptr)
    {
        free(this->ThirdBodyFactor1D);
    }
}
