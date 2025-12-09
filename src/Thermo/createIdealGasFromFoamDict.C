
#include "createIdealGasFromFoamDict.H"
#include "fvCFD.H"
#include "hashedWordList.H"

namespace FastThermo
{

void createIdealGasFromFoamDict
(
    const Foam::dictionary& physicalDict,
    std::unique_ptr<idealGas>& gas
)
{

    hashedWordList speciesTable(physicalDict.lookup("species"));
    int nSpecies = speciesTable.size();

    // Create temporary idealgas object
    gas = std::make_unique<idealGas>(nSpecies);


    double TcommonMax_ = 0;
    double TcommonMin_ = 1e100;
    for(int i = 0; i < speciesTable.size();i++)
    {
        const Foam::dictionary specieDict(physicalDict.subDict(speciesTable[i]));
        const Foam::dictionary thermodynamicsDict(specieDict.subDict("thermodynamics"));
        gas->Tcommon[i] = thermodynamicsDict.lookup<scalar>("Tcommon");
        gas->Tlow[i] = thermodynamicsDict.lookup<scalar>("Tlow");
        gas->Thigh[i] = thermodynamicsDict.lookup<scalar>("Thigh");
        Foam::FixedList<scalar,7> temp1(thermodynamicsDict.lookup("highCpCoeffs"));

        const Foam::dictionary species(specieDict.subDict("specie"));
        gas->W[i] = species.lookup<scalar>("molWeight");
        gas->invW[i] = 1.0/gas->W[i];
        for(unsigned int j = 0; j < 7; j ++)
        {
            gas->HCoeffs[i][j] = temp1[j];
        }
        Foam::FixedList<scalar,7> temp2(thermodynamicsDict.lookup("lowCpCoeffs")); 
        for(unsigned int j = 0; j < 7; j ++)
        {
            gas->LCoeffs[i][j] = temp2[j];
        }               
        TcommonMax_ = (gas->Tcommon[i]>TcommonMax_)?gas->Tcommon[i]:TcommonMax_;
        TcommonMin_ = (gas->Tcommon[i]<TcommonMin_)?gas->Tcommon[i]:TcommonMin_;
    }
    gas->TcommonMin = TcommonMin_;
    gas->TcommonMax = TcommonMax_;

    for(int i = 0; i < gas->nSpecies;i++)
    {
        if (gas->TlowMin > gas->Tlow[i])
        {
            gas->TlowMin = gas->Tlow[i];
        }
        if(gas->ThighMax < gas->Thigh[i])
        {
            gas->ThighMax = gas->Thigh[i];
        }
    }
}
}