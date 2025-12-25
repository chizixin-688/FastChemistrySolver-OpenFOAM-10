/*---------------------------------------------------------------------------*\
  Description
      Header file for gas phase reaction objects. 
      Initializing OptReaction objects using OpenFOAM dictionary

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

#ifndef CreateGasPhaseReactionFromFoamDict_H
#define CreateGasPhaseReactionFromFoamDict_H

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
#include "OptReaction.H"

//=============================================================================//

namespace FastChemistry
{



    // Read reaction information from OpenFOAM, read reaction with integer number
    // e.g. A+B=C
    void readReactionInfo
    (
        std::vector<unsigned int>& lhsIndex,
        std::vector<unsigned int>& rhsIndex,
        const dictionary& nthreaction,
        const hashedWordList& speciesTable,
        std::unique_ptr<OptReaction>& reaction
    );

    // Read reaction information from OpenFOAM, read reaction with non-integer number
    // e.g. A+1.5B=C, A+B^1.5=B
    void readReactionInfo
    (
        std::vector<unsigned int>& lhsIndex,
        std::vector<double>& lhsStoichCoeff,
        std::vector<double>& lhsReactionOrder,
        std::vector<unsigned int>& rhsIndex,
        std::vector<double>& rhsStoichCoeff,
        std::vector<double>& rhsRectionOrder,
        const dictionary& nthreaction,
        const hashedWordList& speciesTable,
        std::unique_ptr<OptReaction>& reaction
    );

    // Check whether the stoichiometric number is an integer, this class only treat 
    // integer.
    bool checkInteger
    (
        const dictionary& nthreaction
    );


    void createGasPhaseReactionFromFoamDict
    (
        const Foam::dictionary& chemistryDict,
        const Foam::dictionary& thermoDict,
        std::unique_ptr<OptReaction>& reaction
    );

}

#endif