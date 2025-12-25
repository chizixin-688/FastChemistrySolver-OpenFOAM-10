/*---------------------------------------------------------------------------*\
  Description
      Initializing OptReaction objects using OpenFOAM dictionary

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. FastChemistry headers
//---------------------------------
#include "createGasPhaseReactionFromFoamDict.H"

//=============================================================================//


namespace FastChemistry
{

void readReactionInfo
(
    std::vector<unsigned int>& inputLhsIndex,
    std::vector<unsigned int>& inputRhsIndex,
    const dictionary& nthreaction,
    const hashedWordList& speciesTable,
    std::unique_ptr<OptReaction>& reaction
)
{
    inputLhsIndex.clear();
    inputRhsIndex.clear();

    Foam::string reactionName = nthreaction.lookup("reaction");
    reaction->reactionTable_.push_back(reactionName);    
    std::string stdReactionName(reactionName);

    std::istringstream iss(stdReactionName);
    std::vector<std::string> words;

    std::vector<std::string> ReactantStr;
    std::vector<std::string> ProductStr;

    std::string Word;
    while (iss >> Word) 
    {
        words.push_back(Word);
    }

    size_t index = 0;
    for (size_t i = 0; i < words.size();i++)
    {
        if(words[i]=="=")
        {
            index =i;
        }
    }
    for (size_t i = 0; i < index;i++)
    {
        if(words[i]!="+")
        {
            ReactantStr.push_back(words[i]);
        }
    }
    for (size_t i = index+1; i < words.size();i++)
    {
        if(words[i]!="+")
        {
            ProductStr.push_back(words[i]);
        }
    }

    // Reactant
    for(size_t  i = 0; i < ReactantStr.size();i++)
    {
        size_t first = 0;
        size_t second = ReactantStr[i].size();
        for(unsigned int  j = 0; j < ReactantStr[i].size();j++)
        {
            if(!std::isdigit(ReactantStr[i][j]) && ReactantStr[i][j]!='.')
            {
                first = j;
                break;
            }
        }
        for(unsigned int  j = 0; j < ReactantStr[i].size();j++)
        {
            if(ReactantStr[i][j]=='^')
            {
                second = j;
            }
        }
        std::string coeffStr = ReactantStr[i].substr(0, first);
        std::string speciesStr = ReactantStr[i].substr(first,second-first);

        unsigned int sl = 0;
        if(first==0)
        {
            sl=1;
        }
        else
        {
            double val = std::round(std::stod(coeffStr));
            sl = static_cast<unsigned int>(val);
        }
        const int newSpecIndex = speciesTable[speciesStr];

        while(sl!=0)
        {
            inputLhsIndex.push_back(newSpecIndex);  
            sl--;
        }
    }

    // Product
    for(size_t  i = 0; i < ProductStr.size();i++)
    {

        size_t first = 0;
        size_t second = ProductStr[i].size();
        for(size_t  j = 0; j < ProductStr[i].size();j++)
        {
            if(!std::isdigit(ProductStr[i][j]) && ProductStr[i][j]!='.')
            {
                first = j;
                break;
            }
        }
        for(size_t  j = 0; j < ProductStr[i].size();j++)
        {
            if(ProductStr[i][j]=='^')
            {
                second = j;
            }
        }        
        std::string coeffStr = ProductStr[i].substr(0, first);
        std::string speciesStr = ProductStr[i].substr(first,second-first);


        unsigned int sr = 0;
        if(first==0)
        {
            sr=1;
        }
        else
        {
            double val = std::round(std::stod(coeffStr));
            sr = static_cast<unsigned int>(val);
        }

        const int newSpecIndex = speciesTable[speciesStr];

        while(sr!=0)
        {
            inputRhsIndex.push_back(newSpecIndex);            
            sr--;
        }
    }
}
void readReactionInfo
(
    std::vector<unsigned int>& inputLhsIndex,
    std::vector<double>& inputLhsstoichCoeff,
    std::vector<double>& inputLhsReactionOrder,
    std::vector<unsigned int>& inputRhsIndex,
    std::vector<double>& inputRhsstoichCoeff,
    std::vector<double>& inputRhsReactionOrder,
    const dictionary& nthreaction,
    const hashedWordList& speciesTable,
    std::unique_ptr<OptReaction>& reaction
)
{
    inputLhsIndex.clear();
    inputLhsstoichCoeff.clear();
    inputLhsReactionOrder.clear();
    inputRhsIndex.clear();
    inputRhsstoichCoeff.clear();
    inputRhsReactionOrder.clear();

    string reactionName = nthreaction.lookup("reaction");
    reaction->reactionTable_.push_back(reactionName);    
    std::string stdReactionName(reactionName);

    std::istringstream iss(stdReactionName);
    std::vector<std::string> words;

    std::vector<std::string> ReactantStr;
    std::vector<std::string> ProductStr;

    std::string Word;
    while (iss >> Word) 
    {
        words.push_back(Word);
    }

    size_t index = 0;
    for (size_t i = 0; i < words.size();i++)
    {
        if(words[i]=="=")
        {
            index =i;
        }
    }
    for (size_t i = 0; i < index;i++)
    {
        if(words[i]!="+")
        {
            ReactantStr.push_back(words[i]);
        }
    }
    for (size_t i = index+1; i < words.size();i++)
    {
        if(words[i]!="+")
        {
            ProductStr.push_back(words[i]);
        }
    }

    // Reactant
    for(size_t  i = 0; i < ReactantStr.size();i++)
    {
        size_t first = 0;
        size_t second = ReactantStr[i].size();
        for(unsigned int  j = 0; j < ReactantStr[i].size();j++)
        {
            if(!std::isdigit(ReactantStr[i][j]) && ReactantStr[i][j]!='.')
            {
                first = j;
                break;
            }
        }
        for(unsigned int  j = 0; j < ReactantStr[i].size();j++)
        {
            if(ReactantStr[i][j]=='^')
            {
                second = j;
                break;
            }
        }
        std::string coeffStr = ReactantStr[i].substr(0, first);
        std::string speciesStr = ReactantStr[i].substr(first,second-first);
        std::string orderStr = ReactantStr[i].substr(second);

        double sl = 0;
        if(first==0)
        {
            sl=1.0;
        }
        else
        {
            sl = (std::stod(coeffStr));
        }
        const unsigned int newSpecIndex = speciesTable[speciesStr];
        double el = 0;
        if(orderStr.empty())
        {
            el = sl;
        }
        else
        {
            el = std::stod(orderStr.substr(1));
        }

        inputLhsIndex.push_back(newSpecIndex);
        inputLhsstoichCoeff.push_back(sl);
        inputLhsReactionOrder.push_back(el);
    }

    // Product
    for(size_t  i = 0; i < ProductStr.size();i++)
    {

        size_t first = 0;
        size_t second = ProductStr[i].size();
        for(size_t  j = 0; j < ProductStr[i].size();j++)
        {
            if(!std::isdigit(ProductStr[i][j]) && ProductStr[i][j]!='.')
            {
                first = j;
                break;
            }
        }
        for(size_t  j = 0; j < ProductStr[i].size();j++)
        {
            if(ProductStr[i][j]=='^')
            {
                second = j;
                break;
            }
        }        
        std::string coeffStr = ProductStr[i].substr(0, first);
        std::string speciesStr = ProductStr[i].substr(first,second-first);
        std::string orderStr = ProductStr[i].substr(second);


        double sr = 0;
        if(first==0)
        {
            sr=1.0;
        }
        else
        {
            sr = std::stod(coeffStr);
        }
        const unsigned int newSpecIndex = speciesTable[speciesStr];
        double er = 0;
        if(orderStr.empty())
        {
            er = sr;
        }
        else
        {
            er = std::stod(orderStr.substr(1));
        }

        inputRhsIndex.push_back(newSpecIndex);
        inputRhsstoichCoeff.push_back(sr);
        inputRhsReactionOrder.push_back(er);
    }
}



bool checkInteger
(
    const dictionary& nthreaction
)
{
    bool isInteger = true;
    string reactionName = nthreaction.lookup("reaction");
    std::string stdReactionName(reactionName);

    std::istringstream iss(stdReactionName);
    std::vector<std::string> words;

    std::vector<std::string> ReactantStr;
    std::vector<std::string> ProductStr;

    std::string Word;
    while (iss >> Word) 
    {
        words.push_back(Word);
    }

    size_t index = 0;
    for (size_t i = 0; i < words.size();i++)
    {
        if(words[i]=="=")
        {
            index =i;
        }
    }
    for (size_t i = 0; i < index;i++)
    {
        if(words[i]!="+")
        {
            ReactantStr.push_back(words[i]);
        }
    }
    for (size_t i = index+1; i < words.size();i++)
    {
        if(words[i]!="+")
        {
            ProductStr.push_back(words[i]);
        }
    }

    //ReactantStr example: ["CH4", "2O2^1.0", "0.5O2^1.0", "0.5O2^1.5", "O2^1.0"]
    for(size_t i = 0; i < ReactantStr.size();i++)
    {
        size_t first = 0;
        size_t second = ReactantStr[i].size();
        for(size_t  j = 0; j < ReactantStr[i].size();j++)
        {
            if(!std::isdigit(ReactantStr[i][j]) && ReactantStr[i][j]!='.')
            {
                first = j;
                break;
            }
        }
        for(size_t  j = 0; j < ReactantStr[i].size();j++)
        {
            if(ReactantStr[i][j]=='^')
            {
                second = j;
                break;
            }
        }
        std::string coeffStr = ReactantStr[i].substr(0, first);
        std::string speciesStr = ReactantStr[i].substr(first,second-first);
        std::string reactionOrderStr = ReactantStr[i].substr(second);

        // coeffStr e.g. "1", "1.0", "1.2", ""
        if(!coeffStr.empty())
        {
            double val = std::stod(coeffStr); 
            if (fabs(val - round(val)) > 2.22e-16)
            {
                isInteger = false;// Stoichiometric number is not an integer
            }
        }

        if(!reactionOrderStr.empty())
        {
            reactionOrderStr = reactionOrderStr.substr(1);
            double val = std::stod(reactionOrderStr);
            if (fabs(val - round(1.0)) > 2.22e-16)
            {
                isInteger = false;// Reaction order is not 1.0
            }
        }
    }

    for(size_t i = 0; i < ProductStr.size();i++)
    {
        size_t first = 0;
        size_t second = ProductStr[i].size();
        for(size_t  j = 0; j < ProductStr[i].size();j++)
        {
            if(!std::isdigit(ProductStr[i][j]) && ProductStr[i][j]!='.')
            {
                first = j;
                break;
            }

        }

        for(size_t  j = 0; j < ProductStr[i].size();j++)
        {
            if(ProductStr[i][j]=='^')
            {
                second = j;
                break;
            }
        }

        std::string coeffStr = ProductStr[i].substr(0, first);
        std::string speciesStr = ProductStr[i].substr(first,second-first);
        std::string reactionOrderStr = ProductStr[i].substr(second);
        if(!coeffStr.empty())
        {
            double val = std::stod(coeffStr);
            if (fabs(val - round(val)) > 2.22e-16)
            {
                isInteger = false;// Stoichiometric number is not an integer
            }
        }
        if(!reactionOrderStr.empty())
        {
            reactionOrderStr = reactionOrderStr.substr(1);
            double val = std::stod(reactionOrderStr);
            if (fabs(val - round(val) > 2.22e-16))
            {
                isInteger = false;// Reaction order is not 1.0
            }
        }
    }
    return isInteger;
}

void createGasPhaseReactionFromFoamDict
(
    const Foam::dictionary& chemistryDict,
    const Foam::dictionary& thermoDict,
    std::unique_ptr<OptReaction>& gasPhaseReaction
)
{
    gasPhaseReaction = std::make_unique<OptReaction>();

    hashedWordList speciesTable(thermoDict.lookup("species"));

    gasPhaseReaction->speciesTable_.resize(speciesTable.size());
    for(unsigned int i=0; i<gasPhaseReaction->speciesTable_.size();i++)
    {
        gasPhaseReaction->speciesTable_[i] = speciesTable[i];
    }

    const dictionary& reactions(chemistryDict.subDict("reactions"));

    unsigned int nLindemann = 0;
    unsigned int nTroe = 0;
    unsigned int nSRI = 0;
    unsigned int nLindemannFO = 0;
    unsigned int nTroeFO = 0;
    unsigned int nSRIFO = 0;
    unsigned int nLindemannCA = 0;
    unsigned int nTroeCA = 0;
    unsigned int nSRICA = 0;

    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& nthreaction = reactions.subDict(key);
        const word reactionTypeName = nthreaction.lookup("type");
        Foam::string reactionName = nthreaction.lookup("reaction");

        gasPhaseReaction->n_Reactions++;

        if(reactionTypeName == "irreversibleArrhenius")
        {
            gasPhaseReaction->n_Arrhenius++;
        }
        else if(reactionTypeName == "reversibleArrhenius")
        {
            gasPhaseReaction->n_Arrhenius++;
        }
        else if(reactionTypeName == "nonEquilibriumReversibleArrhenius")
        {
            gasPhaseReaction->n_NonEquilibriumReversibleArrhenius++;
        }
        else if(reactionTypeName == "nonEquilibriumReversibleThirdBodyArrhenius")
        {
            gasPhaseReaction->n_NonEquilibriumThirdBodyReaction++;
        }
        else if
        (
            reactionTypeName == "reversibleThirdBodyArrhenius"||
            reactionTypeName == "irreversibleThirdBodyArrhenius"
        )
        {
            gasPhaseReaction->n_ThirdBodyReaction++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusLindemannFallOff"||
            reactionTypeName == "irreversibleArrheniusLindemannFallOff"
        )
        {
            gasPhaseReaction->n_Fall_Off_Reaction++;nLindemann++;nLindemannFO++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusTroeFallOff"||
            reactionTypeName == "irreversibleArrheniusTroeFallOff"
        )
        {
            gasPhaseReaction->n_Fall_Off_Reaction++;nTroe++;nTroeFO++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusSRIFallOff"||
            reactionTypeName == "irreversibleArrheniusSRIFallOff"
        )
        {
            gasPhaseReaction->n_Fall_Off_Reaction++;nSRI++;nSRIFO++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusLindemannChemicallyActivated"||
            reactionTypeName == "irreversibleArrheniusLindemannChemicallyActivated"
        )
        {
            gasPhaseReaction->n_ChemicallyActivated_Reaction++;nLindemann++;nLindemannCA++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusTroeChemicallyActivated"||
            reactionTypeName == "irreversibleArrheniusTroeChemicallyActivated"
        )
        {
            gasPhaseReaction->n_ChemicallyActivated_Reaction++;nTroe++;nTroeCA++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusSRIChemicallyActivated"||
            reactionTypeName == "irreversibleArrheniusSRIChemicallyActivated"
        )
        {
            gasPhaseReaction->n_ChemicallyActivated_Reaction++;nSRI++;nSRICA++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusPLOG"||
            reactionTypeName == "irreversibleArrheniusPLOG"
        )
        {
            gasPhaseReaction->n_PlogReaction++;
        }
        else
        {
            FatalErrorInFunction<< "unknown reaction type:"
                << reactionTypeName << exit(FatalError);
        }
    }
    gasPhaseReaction->n_Lindemann = nLindemann;
    gasPhaseReaction->n_Troe = nTroe;
    gasPhaseReaction->n_SRI = nSRI;

    gasPhaseReaction->n_LindemannFO = nLindemannFO;
    gasPhaseReaction->n_TroeFO = nTroeFO;
    gasPhaseReaction->n_SRIFO = nSRIFO;
    
    gasPhaseReaction->n_LindemannCA = nLindemannCA;
    gasPhaseReaction->n_TroeCA = nTroeCA;
    gasPhaseReaction->n_SRICA = nSRICA;
    {
        gasPhaseReaction->Itbr[0] = 0;
        gasPhaseReaction->Itbr[1] = gasPhaseReaction->n_NonEquilibriumThirdBodyReaction;
        gasPhaseReaction->Itbr[2] = gasPhaseReaction->Itbr[1] + gasPhaseReaction->n_ThirdBodyReaction;
        gasPhaseReaction->Itbr[3] = gasPhaseReaction->Itbr[2] + gasPhaseReaction->n_Fall_Off_Reaction;   
        gasPhaseReaction->Itbr[4] = gasPhaseReaction->Itbr[3] + gasPhaseReaction->n_ChemicallyActivated_Reaction;       
        gasPhaseReaction->Itbr[5] = gasPhaseReaction->Itbr[4] + gasPhaseReaction->n_NonEquilibriumThirdBodyReaction; 
    }

    {
        auto& IKfref = gasPhaseReaction->Ikf;
        IKfref[0] = 0;
        IKfref[1] = gasPhaseReaction->n_Arrhenius;
        IKfref[2] = IKfref[1] + gasPhaseReaction->n_NonEquilibriumReversibleArrhenius;
        IKfref[3] = IKfref[2] + gasPhaseReaction->n_NonEquilibriumThirdBodyReaction;   
        IKfref[4] = IKfref[3] + gasPhaseReaction->n_ThirdBodyReaction;       
        IKfref[5] = IKfref[4] + gasPhaseReaction->n_Fall_Off_Reaction; 
        IKfref[6] = IKfref[5] + gasPhaseReaction->n_ChemicallyActivated_Reaction; 
        IKfref[7] = IKfref[6] + gasPhaseReaction->n_PlogReaction;
        IKfref[8] = IKfref[7] + gasPhaseReaction->n_Fall_Off_Reaction;   
        IKfref[9] = IKfref[8] + gasPhaseReaction->n_ChemicallyActivated_Reaction;   
        IKfref[10] = IKfref[9] + gasPhaseReaction->n_NonEquilibriumReversibleArrhenius;        
        IKfref[11] = IKfref[10] + gasPhaseReaction->n_NonEquilibriumThirdBodyReaction;
        IKfref[12] = IKfref[11] + gasPhaseReaction->n_PlogReaction;

        gasPhaseReaction->offset_kinf = - IKfref[4] + IKfref[7];
    }


    //gasPhaseReaction->n_Reactions                        = n_Reactions;
    gasPhaseReaction->nSpecies = speciesTable.size();

    gasPhaseReaction->A.resize(gasPhaseReaction->Ikf[12]);
    gasPhaseReaction->beta.resize(gasPhaseReaction->Ikf[12]);
    gasPhaseReaction->Ta.resize(gasPhaseReaction->Ikf[12]);
    gasPhaseReaction->lhsSpeciesIndex.resize(gasPhaseReaction->n_Reactions);
    gasPhaseReaction->rhsSpeciesIndex.resize(gasPhaseReaction->n_Reactions);
    gasPhaseReaction->lhsStoichCoeff.resize(gasPhaseReaction->n_Reactions);    
    gasPhaseReaction->rhsStoichCoeff.resize(gasPhaseReaction->n_Reactions);
    gasPhaseReaction->lhsReactionOrder.resize(gasPhaseReaction->n_Reactions);
    gasPhaseReaction->rhsReactionOrder.resize(gasPhaseReaction->n_Reactions);
    std::vector<std::vector<double>> ThirdBodyFactor(gasPhaseReaction->Itbr[5]);

    gasPhaseReaction->alpha_.resize(0);
    gasPhaseReaction->alpha_.reserve(nTroe);
    gasPhaseReaction->Ts_.resize(0);
    gasPhaseReaction->Ts_.reserve(nTroe);    
    gasPhaseReaction->Tss_.resize(0);
    gasPhaseReaction->Tss_.reserve(nTroe);    
    gasPhaseReaction->Tsss_.resize(0);
    gasPhaseReaction->Tsss_.reserve(nTroe);

    gasPhaseReaction->a_.resize(0);
    gasPhaseReaction->b_.resize(0);
    gasPhaseReaction->c_.resize(0);
    gasPhaseReaction->d_.resize(0);
    gasPhaseReaction->e_.resize(0);
    gasPhaseReaction->a_.reserve(nSRI);
    gasPhaseReaction->b_.reserve(nSRI);
    gasPhaseReaction->c_.reserve(nSRI);
    gasPhaseReaction->d_.reserve(nSRI);
    gasPhaseReaction->e_.reserve(nSRI);

    unsigned int alignSpecies = ((gasPhaseReaction->nSpecies+3)/4)*4;


    if (posix_memalign(reinterpret_cast<void**>(&gasPhaseReaction->negGstdByRT), 32, alignSpecies * sizeof(double)))
    {
        throw std::bad_alloc();
    }
    memset(gasPhaseReaction->negGstdByRT, 0, alignSpecies * sizeof(double));

    if (posix_memalign(reinterpret_cast<void**>(&gasPhaseReaction->invNegGstdByRT), 32, alignSpecies * sizeof(double)))
    {
        throw std::bad_alloc();
    }
    memset(gasPhaseReaction->invNegGstdByRT, 1, alignSpecies * sizeof(double));


    gasPhaseReaction->isIrreversible.resize(gasPhaseReaction->n_Reactions,0); 
    gasPhaseReaction->isGlobal.resize(gasPhaseReaction->n_Reactions,0); 

    // Find temperature independent Arrhenius reaction
    int iArrhenius = 0;
    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");

        bool isInteger = checkInteger(reactDict);

        
        if
        (
            reactionTypeName=="irreversibleArrhenius"||
            reactionTypeName=="reversibleArrhenius"
        )
        {
            auto a = reactDict.lookup<scalar>("beta");
            auto b = reactDict.lookup<scalar>("Ta");
            if(a==0&&b==0)
            {
                if(reactionTypeName.find("irreversible",0)!=std::string::npos)
                {gasPhaseReaction->isIrreversible[iArrhenius]=1;}
                gasPhaseReaction->reactionType_.push_back(reactionTypeName);
                gasPhaseReaction->reactionName_.push_back(key);
                gasPhaseReaction->A[iArrhenius] = reactDict.lookup<scalar>("A");
                gasPhaseReaction->beta[iArrhenius] = reactDict.lookup<scalar>("beta");
                gasPhaseReaction->Ta[iArrhenius] = reactDict.lookup<scalar>("Ta");

                if(isInteger==true)
                {
                    FastChemistry::readReactionInfo
                    (
                        gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                        gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                        reactDict,
                        speciesTable,
                        gasPhaseReaction
                    );
                }
                else
                {
                    FastChemistry::readReactionInfo
                    (
                        gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                        gasPhaseReaction->lhsStoichCoeff[iArrhenius],
                        gasPhaseReaction->lhsReactionOrder[iArrhenius],
                        gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                        gasPhaseReaction->rhsStoichCoeff[iArrhenius],
                        gasPhaseReaction->rhsReactionOrder[iArrhenius],
                        reactDict,
                        speciesTable,
                        gasPhaseReaction
                    );
                    gasPhaseReaction->isGlobal[iArrhenius] = 1;
                }

                iArrhenius++;
            }
        }
    }


    // Find temperature related reaction
    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = checkInteger(reactDict);

        if
        (
            reactionTypeName=="irreversibleArrhenius"||
            reactionTypeName=="reversibleArrhenius"
        )
        {
            auto a = reactDict.lookup<scalar>("beta");
            auto b = reactDict.lookup<scalar>("Ta");

            if(!(a==0&&b==0))
            {
                if(reactionTypeName.find("irreversible",0)!=std::string::npos)
                {gasPhaseReaction->isIrreversible[iArrhenius]=1;}
                gasPhaseReaction->reactionType_.push_back(reactionTypeName);
                gasPhaseReaction->reactionName_.push_back(key);
                gasPhaseReaction->A[iArrhenius] = reactDict.lookup<scalar>("A");
                gasPhaseReaction->beta[iArrhenius] = reactDict.lookup<scalar>("beta");
                gasPhaseReaction->Ta[iArrhenius] = reactDict.lookup<scalar>("Ta");

                if(isInteger==true)
                {
                    FastChemistry::readReactionInfo
                    (
                        gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                        gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                        reactDict,
                        speciesTable,
                        gasPhaseReaction
                    );
                }
                else
                {
                    FastChemistry::readReactionInfo
                    (
                        gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                        gasPhaseReaction->lhsStoichCoeff[iArrhenius],
                        gasPhaseReaction->lhsReactionOrder[iArrhenius],
                        gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                        gasPhaseReaction->rhsStoichCoeff[iArrhenius],
                        gasPhaseReaction->rhsReactionOrder[iArrhenius],
                        reactDict,
                        speciesTable,
                        gasPhaseReaction
                    );
                    gasPhaseReaction->isGlobal[iArrhenius] = 1;                    
                }
                iArrhenius++;
            }
        }
    }


    auto j = gasPhaseReaction->Ikf[9];
    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = checkInteger(reactDict);

        if(reactionTypeName=="nonEquilibriumReversibleArrhenius")
        {
            gasPhaseReaction->isIrreversible[iArrhenius]=2;
            const dictionary& forwardDict = reactDict.subDict("forward");
            const dictionary& reverseDict = reactDict.subDict("reverse");

            gasPhaseReaction->reactionType_.push_back(reactionTypeName);            
            gasPhaseReaction->reactionName_.push_back(key);
            gasPhaseReaction->A[iArrhenius] = forwardDict.lookup<scalar>("A");       
            gasPhaseReaction->beta[iArrhenius] = forwardDict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[iArrhenius] = forwardDict.lookup<scalar>("Ta");

            gasPhaseReaction->A[j] = reverseDict.lookup<scalar>("A");       
            gasPhaseReaction->beta[j] = reverseDict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[j] = reverseDict.lookup<scalar>("Ta");
            if(isInteger==true)
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
            }
            else
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->lhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->lhsReactionOrder[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );                
                    gasPhaseReaction->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;j++;
        }
    }

    unsigned int k = 0;
    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = checkInteger(reactDict);

        if(reactionTypeName=="nonEquilibriumReversibleThirdBodyArrhenius")
        {
            gasPhaseReaction->isIrreversible[iArrhenius]=2;
            const dictionary& forwardDict = reactDict.subDict("forward");
            const dictionary& reverseDict = reactDict.subDict("reverse");

            gasPhaseReaction->reactionType_.push_back(reactionTypeName); 
            gasPhaseReaction->reactionName_.push_back(key);

            gasPhaseReaction->A[iArrhenius] = forwardDict.lookup<scalar>("A");     
            gasPhaseReaction->beta[iArrhenius] = forwardDict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[iArrhenius] = forwardDict.lookup<scalar>("Ta");

            gasPhaseReaction->A[j] = reverseDict.lookup<scalar>("A");        
            gasPhaseReaction->beta[j] = reverseDict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[j] = reverseDict.lookup<scalar>("Ta");     

            List<Tuple2<word, scalar>> forwardCoeffs(forwardDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(forwardCoeffs.size());
            forAll(forwardCoeffs, n)
            {
                const int l = speciesTable[(forwardCoeffs[n].first())];
                const scalar ThirdBodyFactor_n = forwardCoeffs[n].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_n;
            }
            
            List<Tuple2<word, scalar>> reverseCoeffs(reverseDict.lookup("coeffs"));

            auto begin = k + gasPhaseReaction->Itbr[4];
            

            ThirdBodyFactor[begin].resize(reverseCoeffs.size());
            forAll(reverseCoeffs, n)
            {
                const int l = speciesTable[(reverseCoeffs[n].first())];
                const scalar ThirdBodyFactor_n = reverseCoeffs[n].second();
                ThirdBodyFactor[begin][l] = ThirdBodyFactor_n;
            }
            if(isInteger==true)
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
            }
            else
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->lhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->lhsReactionOrder[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                ); 
                gasPhaseReaction->isGlobal[iArrhenius] = 1;
            }
            iArrhenius++;j++;k++;
        }
    }

    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = checkInteger(reactDict);

        if
        (
            reactionTypeName=="reversibleThirdBodyArrhenius"||
            reactionTypeName=="irreversibleThirdBodyArrhenius"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {gasPhaseReaction->isIrreversible[iArrhenius]=1;}
            gasPhaseReaction->reactionType_.push_back(reactionTypeName); 
            gasPhaseReaction->reactionName_.push_back(key);

            gasPhaseReaction->A[iArrhenius] = reactDict.lookup<scalar>("A");        
            gasPhaseReaction->beta[iArrhenius] = reactDict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[iArrhenius] = reactDict.lookup<scalar>("Ta");

            List<Tuple2<word, scalar>> coeffs(reactDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }
            if(isInteger==true)
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
            }
            else
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->lhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->lhsReactionOrder[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
                gasPhaseReaction->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;k++;
        }
    }


    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = checkInteger(reactDict);

        if
        (
            reactionTypeName=="reversibleArrheniusLindemannFallOff"||
            reactionTypeName=="irreversibleArrheniusLindemannFallOff"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {gasPhaseReaction->isIrreversible[iArrhenius]=1;}
            gasPhaseReaction->reactionType_.push_back(reactionTypeName);
            gasPhaseReaction->reactionName_.push_back(key);
            const dictionary& k0Dict = reactDict.subDict("k0");
            const dictionary& kInfDict = reactDict.subDict("kInf");

            const dictionary& thirdBodyEfficienciesDict =
            reactDict.subDict("thirdBodyEfficiencies");
    
            gasPhaseReaction->A[iArrhenius] = k0Dict.lookup<scalar>("A");        
            gasPhaseReaction->beta[iArrhenius] = k0Dict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[iArrhenius] = k0Dict.lookup<scalar>("Ta");
    
            auto begin = iArrhenius - gasPhaseReaction->Ikf[4] + gasPhaseReaction->Ikf[7];
            gasPhaseReaction->A[begin] = (kInfDict.lookup<scalar>("A")) ;         
            gasPhaseReaction->beta[begin] = (kInfDict.lookup<scalar>("beta")) ;
            gasPhaseReaction->Ta[begin] = (kInfDict.lookup<scalar>("Ta")) ;            
    
            List<Tuple2<word, scalar>> coeffs(thirdBodyEfficienciesDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
    
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }
            gasPhaseReaction->LindemannFO.push_back(iArrhenius);

            if(isInteger==true)
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
            }
            else
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->lhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->lhsReactionOrder[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
                    gasPhaseReaction->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;k++;
       }
    }

    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = checkInteger(reactDict);
        if
        (
            reactionTypeName=="reversibleArrheniusTroeFallOff"||
            reactionTypeName=="irreversibleArrheniusTroeFallOff"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {gasPhaseReaction->isIrreversible[iArrhenius]=1;}
            gasPhaseReaction->reactionType_.push_back(reactionTypeName);
            gasPhaseReaction->reactionName_.push_back(key);
            const dictionary& k0Dict = reactDict.subDict("k0");
            const dictionary& kInfDict = reactDict.subDict("kInf");
            const dictionary& FDict = reactDict.subDict("F");
            const dictionary& thirdBodyEfficienciesDict =
            reactDict.subDict("thirdBodyEfficiencies");
    
            gasPhaseReaction->A[iArrhenius] = k0Dict.lookup<scalar>("A");        
            gasPhaseReaction->beta[iArrhenius] = k0Dict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[iArrhenius] = k0Dict.lookup<scalar>("Ta");
    
            auto begin = iArrhenius - gasPhaseReaction->Ikf[4] + gasPhaseReaction->Ikf[7];
            gasPhaseReaction->A[begin] = (kInfDict.lookup<scalar>("A")) ;         
            gasPhaseReaction->beta[begin] = (kInfDict.lookup<scalar>("beta")) ;
            gasPhaseReaction->Ta[begin] = (kInfDict.lookup<scalar>("Ta")) ;            
    
            List<Tuple2<word, scalar>> coeffs(thirdBodyEfficienciesDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
    
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }

            gasPhaseReaction->TroeFO.push_back(iArrhenius);
            gasPhaseReaction->alpha_.push_back(FDict.lookup<scalar>("alpha"));    
            gasPhaseReaction->Ts_.push_back(FDict.lookup<scalar>("Ts"));    
            gasPhaseReaction->Tss_.push_back(FDict.lookup<scalar>("Tss"));    
            gasPhaseReaction->Tsss_.push_back(FDict.lookup<scalar>("Tsss"));
            if(isInteger==true)
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
            }
            else
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->lhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->lhsReactionOrder[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
                gasPhaseReaction->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;k++;
       }
    }

    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = checkInteger(reactDict); 
        if
        (
            reactionTypeName=="reversibleArrheniusSRIFallOff"||
            reactionTypeName=="irreversibleArrheniusSRIFallOff"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {gasPhaseReaction->isIrreversible[iArrhenius]=1;}
            gasPhaseReaction->reactionType_.push_back(reactionTypeName);
            gasPhaseReaction->reactionName_.push_back(key);
            const dictionary& k0Dict = reactDict.subDict("k0");
            const dictionary& kInfDict = reactDict.subDict("kInf");
            const dictionary& FDict = reactDict.subDict("F");
            const dictionary& thirdBodyEfficienciesDict =
            reactDict.subDict("thirdBodyEfficiencies");
    
            gasPhaseReaction->A[iArrhenius] = k0Dict.lookup<scalar>("A");        
            gasPhaseReaction->beta[iArrhenius] = k0Dict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[iArrhenius] = k0Dict.lookup<scalar>("Ta");
    
            auto begin = iArrhenius - gasPhaseReaction->Ikf[4] + gasPhaseReaction->Ikf[7];
            gasPhaseReaction->A[begin] = (kInfDict.lookup<scalar>("A")) ;         
            gasPhaseReaction->beta[begin] = (kInfDict.lookup<scalar>("beta")) ;
            gasPhaseReaction->Ta[begin] = (kInfDict.lookup<scalar>("Ta")) ;            
    
            List<Tuple2<word, scalar>> coeffs(thirdBodyEfficienciesDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
    
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }
            
            gasPhaseReaction->SRIFO.push_back(iArrhenius);
            gasPhaseReaction->a_.push_back(FDict.lookup<scalar>("a"));    
            gasPhaseReaction->b_.push_back(FDict.lookup<scalar>("b"));    
            gasPhaseReaction->c_.push_back(FDict.lookup<scalar>("c"));    
            gasPhaseReaction->d_.push_back(FDict.lookup<scalar>("d"));  
            gasPhaseReaction->e_.push_back(FDict.lookup<scalar>("e"));  
            
            if(isInteger==true)
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
            }
            else
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->lhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->lhsReactionOrder[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
                gasPhaseReaction->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;k++;
       }
    }

    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = checkInteger(reactDict);

        if(
            reactionTypeName=="reversibleArrheniusLindemannChemicallyActivated"||
            reactionTypeName=="irreversibleArrheniusLindemannChemicallyActivated"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {gasPhaseReaction->isIrreversible[iArrhenius]=1;}
            gasPhaseReaction->reactionType_.push_back(reactionTypeName); 
            gasPhaseReaction->reactionName_.push_back(key);   

            const dictionary& k0Dict = reactDict.subDict("k0");
            const dictionary& kInfDict = reactDict.subDict("kInf");

            const dictionary& thirdBodyEfficienciesDict = 
            reactDict.subDict("thirdBodyEfficiencies");

            gasPhaseReaction->A[iArrhenius] = k0Dict.lookup<scalar>("A");
            gasPhaseReaction->beta[iArrhenius] = k0Dict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[iArrhenius] = k0Dict.lookup<scalar>("Ta");

            auto begin = iArrhenius - gasPhaseReaction->Ikf[5] + gasPhaseReaction->Ikf[8];
            gasPhaseReaction->A[begin] = kInfDict.lookup<scalar>("A");
            gasPhaseReaction->beta[begin] = kInfDict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[begin] = kInfDict.lookup<scalar>("Ta");

            List<Tuple2<word, scalar>> coeffs(thirdBodyEfficienciesDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }
            gasPhaseReaction->LindemannCA.push_back(iArrhenius);

            if(isInteger==true)
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
            }
            else
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->lhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->lhsReactionOrder[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
                gasPhaseReaction->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;
            k++;
        }

    }

    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = checkInteger(reactDict);

        if(
            reactionTypeName=="reversibleArrheniusTroeChemicallyActivated"||
            reactionTypeName=="irreversibleArrheniusTroeChemicallyActivated"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {gasPhaseReaction->isIrreversible[iArrhenius]=1;}

            gasPhaseReaction->reactionType_.push_back(reactionTypeName); 
            gasPhaseReaction->reactionName_.push_back(key);   

            const dictionary& k0Dict = reactDict.subDict("k0");
            const dictionary& kInfDict = reactDict.subDict("kInf");
            const dictionary& FDict = reactDict.subDict("F");
            const dictionary& thirdBodyEfficienciesDict = 
            reactDict.subDict("thirdBodyEfficiencies");

            gasPhaseReaction->A[iArrhenius] = k0Dict.lookup<scalar>("A");
            gasPhaseReaction->beta[iArrhenius] = k0Dict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[iArrhenius] = k0Dict.lookup<scalar>("Ta");

            auto begin = iArrhenius - gasPhaseReaction->Ikf[5] + gasPhaseReaction->Ikf[8];
            gasPhaseReaction->A[begin] = kInfDict.lookup<scalar>("A");
            gasPhaseReaction->beta[begin] = kInfDict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[begin] = kInfDict.lookup<scalar>("Ta");

            List<Tuple2<word, scalar>> coeffs(thirdBodyEfficienciesDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }

            gasPhaseReaction->TroeCA.push_back(iArrhenius);
            gasPhaseReaction->alpha_.push_back(FDict.lookup<scalar>("alpha"));    
            gasPhaseReaction->Ts_.push_back(FDict.lookup<scalar>("Ts"));    
            gasPhaseReaction->Tss_.push_back(FDict.lookup<scalar>("Tss"));    
            gasPhaseReaction->Tsss_.push_back(FDict.lookup<scalar>("Tsss"));                
            if(isInteger==true)
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
            }
            else
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->lhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->lhsReactionOrder[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
                gasPhaseReaction->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;
            k++;
        }
    }    

    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = checkInteger(reactDict);

        if(
            reactionTypeName=="reversibleArrheniusSRIChemicallyActivated" ||
            reactionTypeName=="irreversibleArrheniusSRIChemicallyActivated" 
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {gasPhaseReaction->isIrreversible[iArrhenius]=1;}

            gasPhaseReaction->reactionType_.push_back(reactionTypeName); 
            gasPhaseReaction->reactionName_.push_back(key);   

            const dictionary& k0Dict = reactDict.subDict("k0");
            const dictionary& kInfDict = reactDict.subDict("kInf");
            const dictionary& FDict = reactDict.subDict("F");
            const dictionary& thirdBodyEfficienciesDict = 
            reactDict.subDict("thirdBodyEfficiencies");

            gasPhaseReaction->A[iArrhenius] = k0Dict.lookup<scalar>("A");
            gasPhaseReaction->beta[iArrhenius] = k0Dict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[iArrhenius] = k0Dict.lookup<scalar>("Ta");

            auto begin = iArrhenius - gasPhaseReaction->Ikf[5] + gasPhaseReaction->Ikf[8];
            gasPhaseReaction->A[begin] = kInfDict.lookup<scalar>("A");
            gasPhaseReaction->beta[begin] = kInfDict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[begin] = kInfDict.lookup<scalar>("Ta");

            List<Tuple2<word, scalar>> coeffs(thirdBodyEfficienciesDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }
            gasPhaseReaction->SRICA.push_back(iArrhenius);
            gasPhaseReaction->a_.push_back(FDict.lookup<scalar>("a"));    
            gasPhaseReaction->b_.push_back(FDict.lookup<scalar>("b"));    
            gasPhaseReaction->c_.push_back(FDict.lookup<scalar>("c"));    
            gasPhaseReaction->d_.push_back(FDict.lookup<scalar>("d"));  
            gasPhaseReaction->e_.push_back(FDict.lookup<scalar>("e"));
            if(isInteger==true)
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
            }
            else
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->lhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->lhsReactionOrder[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
                gasPhaseReaction->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;
            k++;
        }
    }

    {
        unsigned int remain = 4 - gasPhaseReaction->nSpecies%4;

        gasPhaseReaction->AlignSpecies = gasPhaseReaction->nSpecies+remain;

        if (
            posix_memalign
            (
                reinterpret_cast<void**>(&gasPhaseReaction->ThirdBodyFactor1D), 
                32, 
                ThirdBodyFactor.size()*gasPhaseReaction->AlignSpecies*sizeof(double)
            )
            )
        {
            throw std::bad_alloc();
        }
        memset(gasPhaseReaction->ThirdBodyFactor1D, 0, ThirdBodyFactor.size()*gasPhaseReaction->AlignSpecies*sizeof(double));

        auto& ThirdBodyFactor1Dref = gasPhaseReaction->ThirdBodyFactor1D;
        unsigned int count = 0;
        for(unsigned int i = 0; i < ThirdBodyFactor.size();i++)
        {
            for(unsigned int J = 0; J < ThirdBodyFactor[i].size();J++)
            {
                ThirdBodyFactor1Dref[count] = ThirdBodyFactor[i][J];
                count++;
            }
            for(unsigned int J = 0; J < remain;J++)
            {
                ThirdBodyFactor1Dref[count] = 0;
                count++;                
            }
        }
   }



    gasPhaseReaction->APlog.resize(gasPhaseReaction->n_PlogReaction);
    gasPhaseReaction->logAPlog.resize(gasPhaseReaction->n_PlogReaction);
    gasPhaseReaction->betaPlog.resize(gasPhaseReaction->n_PlogReaction);
    gasPhaseReaction->TaPlog.resize(gasPhaseReaction->n_PlogReaction);
    gasPhaseReaction->Prange.resize(gasPhaseReaction->n_PlogReaction);
    gasPhaseReaction->rDeltaP_.resize(gasPhaseReaction->n_PlogReaction);
    gasPhaseReaction->logPi.resize(gasPhaseReaction->n_PlogReaction);
    gasPhaseReaction->Pindex.resize(gasPhaseReaction->n_PlogReaction);
    gasPhaseReaction->ActivePlogReactionIndex.resize(gasPhaseReaction->n_PlogReaction);

    unsigned int a = 0;
    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = checkInteger(reactDict);
        if
        (
            reactionTypeName=="reversibleArrheniusPLOG"||
            reactionTypeName=="irreversibleArrheniusPLOG"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {
                gasPhaseReaction->isIrreversible[iArrhenius]=1;
            }
            gasPhaseReaction->reactionType_.push_back(reactionTypeName); 
            gasPhaseReaction->reactionName_.push_back(key);   

            List<List<double>> PlogData(reactDict.lookup("ArrheniusData"));
            unsigned int pSize = PlogData.size();

            gasPhaseReaction->APlog[a].resize(pSize);
            gasPhaseReaction->logAPlog[a].resize(pSize); 
            gasPhaseReaction->betaPlog[a].resize(pSize);
            gasPhaseReaction->TaPlog[a].resize(pSize);
            gasPhaseReaction->Prange[a].resize(pSize);
            gasPhaseReaction->rDeltaP_[a].resize(pSize-1);
            gasPhaseReaction->logPi[a].resize(pSize);

            for(unsigned int i = 0; i < pSize; i ++)
            {
                gasPhaseReaction->Prange[a][i] = PlogData[i][0];  
                gasPhaseReaction->APlog[a][i] = PlogData[i][1];
                if(PlogData[i][1]<=0)
                {
                    std::cout<<"Pre-factor A should be larger "<<
                    "than zero for Plog reaction"<<std::endl;
                    std::exit(0);
                }
                gasPhaseReaction->logAPlog[a][i] = std::log(gasPhaseReaction->APlog[a][i]);
                gasPhaseReaction->betaPlog[a][i] = PlogData[i][2];
                gasPhaseReaction->TaPlog[a][i] = PlogData[i][3];
                gasPhaseReaction->logPi[a][i] = std::log(gasPhaseReaction->Prange[a][i]);
            }

            for(unsigned int i = 0; i < pSize-1; i ++)
            {
                gasPhaseReaction->rDeltaP_[a][i] = 1.0/
                    (gasPhaseReaction->logPi[a][i+1]-gasPhaseReaction->logPi[a][i]);
            }
            if(isInteger==true)
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
            }
            else
            {
                FastChemistry::readReactionInfo
                (
                    gasPhaseReaction->lhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->lhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->lhsReactionOrder[iArrhenius],
                    gasPhaseReaction->rhsSpeciesIndex[iArrhenius],
                    gasPhaseReaction->rhsStoichCoeff[iArrhenius],
                    gasPhaseReaction->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable,
                    gasPhaseReaction
                );
                gasPhaseReaction->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;a++;
        }
    }

    List<int> sumVki(gasPhaseReaction->n_Reactions);

     for(unsigned int i = 0;i<gasPhaseReaction->n_Reactions;i++)
    {
        sumVki[i] = 0;
        for(unsigned int jj = 0; jj<gasPhaseReaction->rhsSpeciesIndex[i].size();jj++)
        {
            sumVki[i] = sumVki[i] + 1;
        }  
        for(unsigned int jj = 0; jj<gasPhaseReaction->lhsSpeciesIndex[i].size();jj++)
        {
            sumVki[i] = sumVki[i] - 1;
        }  
    } 

    //gasPhaseReaction->Pow_pByRT_SumVki_I.insert({sumVki[0],0.0});
    for(int i = 0; i < sumVki.size();i++)
    {
        auto it = gasPhaseReaction->Pow_pByRT_SumVki_I.find(sumVki[i]);
        if(it==gasPhaseReaction->Pow_pByRT_SumVki_I.end())
        {
            if
            (
                sumVki[i] !=0 && 
                sumVki[i] !=1 && 
                sumVki[i] !=2 && 
                sumVki[i] !=-1 && 
                sumVki[i] !=-2 
            )
            {
                gasPhaseReaction->Pow_pByRT_SumVki_I.insert({sumVki[i],0.0});
            }

        }
    }

    int alignNKf = ((gasPhaseReaction->Ikf[12]+3)/4)*4;
    gasPhaseReaction->Kf_.resize(alignNKf);
    gasPhaseReaction->dKfdT_.resize(alignNKf);


    int alignNtmp_M = ((gasPhaseReaction->Itbr[5]+3)/4)*4;
    gasPhaseReaction->dKfdC_.resize(alignNtmp_M);
    gasPhaseReaction->tmp_M.resize(alignNtmp_M);

    {
        gasPhaseReaction->tmp_ExpSize = (gasPhaseReaction->nSpecies + nTroe*3 + nSRI*2);
        const unsigned int bytes = (gasPhaseReaction->nSpecies + nTroe*3 + nSRI*2)  * sizeof(double);
        if (posix_memalign(reinterpret_cast<void**>(&gasPhaseReaction->tmp_Exp), 32, bytes))
        {
            throw std::bad_alloc();
        }
        memset(gasPhaseReaction->tmp_Exp, 0, bytes);
    }



    gasPhaseReaction->invTs_.resize(nTroe);
    gasPhaseReaction->invTsss_.resize(nTroe);
    for(unsigned int i = 0; i < nTroe;i++)
    {
        gasPhaseReaction->invTs_[i] = 1.0/gasPhaseReaction->Ts_[i];
    }
    for(unsigned int i = 0; i < nTroe;i++)
    {
        gasPhaseReaction->invTsss_[i] = 1.0/gasPhaseReaction->Tsss_[i];
    }

    gasPhaseReaction->invc_.resize(nSRI);
    for(unsigned int i = 0; i < nSRI;i++)
    {
        gasPhaseReaction->invc_[i] = 1.0/gasPhaseReaction->c_[i];
    }    

    gasPhaseReaction->n_Temperature_Independent_Reaction =0;
    if(gasPhaseReaction->n_Arrhenius>0)
    {
        for(unsigned int ii = 0; ii < gasPhaseReaction->n_Arrhenius;ii++)
        {
            if(gasPhaseReaction->beta[ii]==0&&gasPhaseReaction->Ta[ii]==0)
            {gasPhaseReaction->n_Temperature_Independent_Reaction++;}
        }
    }

    for(unsigned int i0 = 0; i0 < gasPhaseReaction->n_Temperature_Independent_Reaction;i0++)
    {
        gasPhaseReaction->Kf_[i0]=gasPhaseReaction->A[i0];
        gasPhaseReaction->dKfdT_[i0] = 0;
    }

    gasPhaseReaction->n_ = gasPhaseReaction->nSpecies+1;

    {
        unsigned int lhsAll=0;
        auto& lhsSpeciesIndexRef = gasPhaseReaction->lhsSpeciesIndex;
        auto& lhsOffsetRef = gasPhaseReaction->lhsOffset;
        auto& lhsSpeciesIndex1DRef = gasPhaseReaction->lhsSpeciesIndex1D;
        for(size_t i = 0; i < lhsSpeciesIndexRef.size();i++)
        {
            for(size_t J = 0; J < lhsSpeciesIndexRef[i].size();J++)
            {
                lhsAll++;
            }
        }
        lhsSpeciesIndex1DRef.resize(lhsAll);    
        lhsOffsetRef.resize(lhsSpeciesIndexRef.size()+1);
        lhsAll=0;
        for(size_t i = 0; i < lhsSpeciesIndexRef.size();i++)
        {
            lhsOffsetRef[i+1] = lhsOffsetRef[i] + static_cast<unsigned int>(lhsSpeciesIndexRef[i].size());
            for(size_t J = 0; J < lhsSpeciesIndexRef[i].size();J++)
            {
                lhsSpeciesIndex1DRef[lhsAll] = lhsSpeciesIndexRef[i][J];
                lhsAll++;
            }
        }         
        lhsOffsetRef[lhsSpeciesIndexRef.size()] = 
            static_cast<unsigned int>(lhsSpeciesIndex1DRef.size());

        auto& rhsSpeciesIndexRef = gasPhaseReaction->rhsSpeciesIndex;
        auto& rhsOffsetRef = gasPhaseReaction->rhsOffset;
        auto& rhsSpeciesIndex1DRef = gasPhaseReaction->rhsSpeciesIndex1D;


        unsigned int rhsAll=0;
        for(size_t i = 0; i < rhsSpeciesIndexRef.size();i++)
        {
            for(size_t J = 0; J < rhsSpeciesIndexRef[i].size();J++)
            {
                rhsAll++;
            }
        }
        rhsSpeciesIndex1DRef.resize(rhsAll);    
        rhsOffsetRef.resize(rhsSpeciesIndexRef.size()+1);
        rhsAll=0;
        for(size_t i = 0; i < rhsSpeciesIndexRef.size();i++)
        {
            rhsOffsetRef[i+1] = rhsOffsetRef[i] + static_cast<unsigned int>(rhsSpeciesIndexRef[i].size());
            for(size_t J = 0; J < rhsSpeciesIndexRef[i].size();J++)
            {
                rhsSpeciesIndex1DRef[rhsAll] = rhsSpeciesIndexRef[i][J];
                rhsAll++;
            }
        }       
        rhsOffsetRef[rhsSpeciesIndexRef.size()] = static_cast<unsigned int>(rhsSpeciesIndex1DRef.size());
    }


    for(std::size_t i = 0; i < gasPhaseReaction->lhsSpeciesIndex.size();i++)
    {
        std::size_t lhsNumber = gasPhaseReaction->lhsSpeciesIndex[i].size();
        std::size_t rhsNumber = gasPhaseReaction->rhsSpeciesIndex[i].size();

        if(gasPhaseReaction->isGlobal[i]==1)
        {
            gasPhaseReaction->reactionGNIindex.push_back(static_cast<unsigned int>(i));

            continue;
        }


        if(lhsNumber==1 && rhsNumber==1)
        {
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                gasPhaseReaction->reversibleReaction11index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D11RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D11RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                gasPhaseReaction->irreversibleReaction11index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D11IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D11IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                gasPhaseReaction->nonEquilibriumReaction11index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D11NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D11NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
            }
        }
        else if(lhsNumber==1 && rhsNumber==2)
        {
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                gasPhaseReaction->reversibleReaction12index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D12RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D12RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D12RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                gasPhaseReaction->irreversibleReaction12index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D12IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D12IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D12IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                gasPhaseReaction->nonEquilibriumReaction12index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D12NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D12NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D12NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
            }
        }
        else if(lhsNumber==1 && rhsNumber==3)
        {
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                gasPhaseReaction->reversibleReaction13index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D13RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D13RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D13RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D13RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][2]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                gasPhaseReaction->irreversibleReaction13index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D13IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D13IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D13IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D13IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][2]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                gasPhaseReaction->nonEquilibriumReaction13index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D13NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D13NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D13NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D13NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][2]);
            }
        }
        else if(lhsNumber==2 && rhsNumber==1)
        {
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                gasPhaseReaction->reversibleReaction21index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D21RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D21RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D21RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                gasPhaseReaction->irreversibleReaction21index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D21IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D21IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D21IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                gasPhaseReaction->nonEquilibriumReaction21index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D21NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D21NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D21NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
            }
        }
        else if(lhsNumber==2 && rhsNumber==2)
        {
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                gasPhaseReaction->reversibleReaction22index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D22RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D22RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D22RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D22RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                gasPhaseReaction->irreversibleReaction22index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D22IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D22IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D22IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D22IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                gasPhaseReaction->nonEquilibriumReaction22index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D22NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D22NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D22NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D22NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
            }
        }
        else if(lhsNumber==2 && rhsNumber==3)
        {
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                gasPhaseReaction->reversibleReaction23index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D23RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D23RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D23RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D23RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D23RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][2]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                gasPhaseReaction->irreversibleReaction23index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D23IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D23IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D23IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D23IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D23IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][2]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                gasPhaseReaction->nonEquilibriumReaction23index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D23NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D23NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D23NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D23NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D23NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][2]);
            }
        }
        else if(lhsNumber==3 && rhsNumber==1)
        {
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                gasPhaseReaction->reversibleReaction31index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D31RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D31RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->lhsSpeciesIndex1D31RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][2]);
                gasPhaseReaction->rhsSpeciesIndex1D31RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                gasPhaseReaction->irreversibleReaction31index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D31IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D31IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->lhsSpeciesIndex1D31IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][2]);
                gasPhaseReaction->rhsSpeciesIndex1D31IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                gasPhaseReaction->nonEquilibriumReaction31index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D31NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D31NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->lhsSpeciesIndex1D31NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][2]);
                gasPhaseReaction->rhsSpeciesIndex1D31NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
            }
        }
        else if(lhsNumber==3 && rhsNumber==2)
        {
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                gasPhaseReaction->reversibleReaction32index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D32RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D32RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->lhsSpeciesIndex1D32RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][2]);
                gasPhaseReaction->rhsSpeciesIndex1D32RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D32RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                gasPhaseReaction->irreversibleReaction32index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D32IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D32IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->lhsSpeciesIndex1D32IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][2]);
                gasPhaseReaction->rhsSpeciesIndex1D32IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D32IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                gasPhaseReaction->nonEquilibriumReaction32index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D32NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D32NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->lhsSpeciesIndex1D32NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][2]);
                gasPhaseReaction->rhsSpeciesIndex1D32NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D32NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
            }

        }
        else if(lhsNumber==3 && rhsNumber==3)
        {
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                gasPhaseReaction->reversibleReaction33index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D33RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D33RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->lhsSpeciesIndex1D33RR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][2]);
                gasPhaseReaction->rhsSpeciesIndex1D33RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D33RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D33RR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][2]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                gasPhaseReaction->irreversibleReaction33index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D33IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D33IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->lhsSpeciesIndex1D33IR.push_back(gasPhaseReaction->lhsSpeciesIndex[i][2]);
                gasPhaseReaction->rhsSpeciesIndex1D33IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D33IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D33IR.push_back(gasPhaseReaction->rhsSpeciesIndex[i][2]);
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                gasPhaseReaction->nonEquilibriumReaction33index.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->lhsSpeciesIndex1D33NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][0]);
                gasPhaseReaction->lhsSpeciesIndex1D33NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][1]);
                gasPhaseReaction->lhsSpeciesIndex1D33NER.push_back(gasPhaseReaction->lhsSpeciesIndex[i][2]);
                gasPhaseReaction->rhsSpeciesIndex1D33NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][0]);
                gasPhaseReaction->rhsSpeciesIndex1D33NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][1]);
                gasPhaseReaction->rhsSpeciesIndex1D33NER.push_back(gasPhaseReaction->rhsSpeciesIndex[i][2]);
            }
        }
        if(lhsNumber>3 || rhsNumber>3)
        {
            gasPhaseReaction->reactionGIindex.push_back(static_cast<unsigned int>(i));
        }
    }

    {
        if(gasPhaseReaction->reactionGNIindex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RFGNI);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JFGNI);
        }
        if(gasPhaseReaction->reactionGIindex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RFGI);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JFGI);
        }
        if(gasPhaseReaction->reversibleReaction11index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF11RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF11RR);
        }
        if(gasPhaseReaction->irreversibleReaction11index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF11IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF11IR);
        }
        if(gasPhaseReaction->nonEquilibriumReaction11index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF11NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF11NER);
        }
        if(gasPhaseReaction->reversibleReaction12index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF12RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF12RR);
        }
        if(gasPhaseReaction->irreversibleReaction12index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF12IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF12IR);
        }
        if(gasPhaseReaction->nonEquilibriumReaction12index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF12NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF12NER);
        }
        if(gasPhaseReaction->reversibleReaction13index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF13RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF13RR);
        }
        if(gasPhaseReaction->irreversibleReaction13index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF13IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF13IR);
        }
        if(gasPhaseReaction->nonEquilibriumReaction13index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF13NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF13NER);
        }
        if(gasPhaseReaction->reversibleReaction21index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF21RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF21RR);
        }
        if(gasPhaseReaction->irreversibleReaction21index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF21IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF21IR);
        }
        if(gasPhaseReaction->nonEquilibriumReaction21index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF21NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF21NER);
        }
        if(gasPhaseReaction->reversibleReaction22index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF22RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF22RR);
        }
        if(gasPhaseReaction->irreversibleReaction22index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF22IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF22IR);
        }
        if(gasPhaseReaction->nonEquilibriumReaction22index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF22NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF22NER);
        }
        if(gasPhaseReaction->reversibleReaction23index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF23RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF23RR);
        }
        if(gasPhaseReaction->irreversibleReaction23index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF23IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF23IR);
        }
        if(gasPhaseReaction->nonEquilibriumReaction23index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF23NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF23NER);
        }

        if(gasPhaseReaction->reversibleReaction31index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF31RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF31RR);
        }
        if(gasPhaseReaction->irreversibleReaction31index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF31IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF31IR);
        }
        if(gasPhaseReaction->nonEquilibriumReaction31index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF31NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF31NER);
        }

        if(gasPhaseReaction->reversibleReaction32index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF32RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF32RR);
        }
        if(gasPhaseReaction->irreversibleReaction32index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF32IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF32IR);
        }
        if(gasPhaseReaction->nonEquilibriumReaction32index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF32NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF32NER);
        }
        if(gasPhaseReaction->reversibleReaction33index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF33RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF33RR);
        }
        if(gasPhaseReaction->irreversibleReaction33index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF33IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF33IR);
        }
        if(gasPhaseReaction->nonEquilibriumReaction33index.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF33NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF33NER);
        }
    }

    /*********************************************************************************************/

}


}   // End namespace FastChemistry