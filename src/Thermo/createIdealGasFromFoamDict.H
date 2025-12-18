#ifndef CreateIdealGasFromFoamDict_H
#define CreateIdealGasFromFoamDict_H
#include "fvCFD.H"
#include "idealGas.H"
#include <memory>
namespace FastChemistry
{

    void createIdealGasFromFoamDict
    (
        const Foam::dictionary& physicalDict,
        std::unique_ptr<idealGas>& gas
    );

}

#endif