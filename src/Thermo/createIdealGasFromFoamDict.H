/*---------------------------------------------------------------------------*\
  Description
      Header file for ideal gas objects. 
      Initializing idealGas objects using OpenFOAM dictionary

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. Standard C++ library headers
//---------------------------------
#include <memory>

//---------------------------------
// 2. OpenFOAM library headers
//---------------------------------
#include "fvCFD.H"
#include "hashedWordList.H"

//---------------------------------
// 3. FastChemistry headers
//---------------------------------
#include "idealGas.H"

//=============================================================================//

#ifndef CreateIdealGasFromFoamDict_H
#define CreateIdealGasFromFoamDict_H


namespace FastChemistry
{

    void createIdealGasFromFoamDict
    (
        const Foam::dictionary& physicalDict,
        std::unique_ptr<idealGas>& gas
    );

}

#endif