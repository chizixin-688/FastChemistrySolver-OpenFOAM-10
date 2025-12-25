/*---------------------------------------------------------------------------*\
  Description
      Template instantiation file for OpenFOAM-10

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. OpenFOAM library headers
//---------------------------------
#include "basicChemistryModel.H"
#include "forGases.H"
#include "makeChemistrySolver.H"

//---------------------------------
// 2. FastChemistry headers
//---------------------------------
#include "OptRosenbrock34.H"
#include "OptRodas34.H"
#include "OptSeulex.H"
#include "FastChemistryModel.H"

//=============================================================================//

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedefThermo(sutherlandTransport,   sensibleEnthalpy,   janafThermo,    perfectGas, specie);
    defineChemistrySolver (FastChemistryModel,  sutherlandTransportsensibleEnthalpyjanafThermoperfectGasspecie)
    makeChemistrySolver(OptRosenbrock34,FastChemistryModel,sutherlandTransportsensibleEnthalpyjanafThermoperfectGasspecie);
    makeChemistrySolver(OptRodas34,     FastChemistryModel,sutherlandTransportsensibleEnthalpyjanafThermoperfectGasspecie);
    makeChemistrySolver(OptSeulex,     FastChemistryModel,sutherlandTransportsensibleEnthalpyjanafThermoperfectGasspecie);
    
}


// ************************************************************************* //
