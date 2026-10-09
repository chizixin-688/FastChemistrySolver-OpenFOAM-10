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
#include <set>
//=============================================================================//


namespace FastChemistry
{

void readElementaryReactionInfo
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
    reaction->reactionEquation_.push_back(reactionName);    
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
    reaction->reactionEquation_.push_back(reactionName);    
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



void validateThirdBodyEfficiencies
(
    const Foam::List<Foam::Tuple2<Foam::word, Foam::scalar>>& coeffs,
    const Foam::hashedWordList& speciesTable,
    const Foam::word& reactionKey,
    const char* what
)
{
    const Foam::label nS = speciesTable.size();

    if (coeffs.size() != static_cast<Foam::label>(nS))
    {
        FatalErrorInFunction
            << "Reaction \"" << reactionKey << "\" (" << what << "): "
            << "the third-body efficiency list must enumerate every species. "
            << "Number of listed entries = " << coeffs.size()
            << ", number of species = " << nS << "."
            << exit(FatalError);
    }

    std::vector<Foam::label> seen(nS, 0);
    for (const auto& item : coeffs)
    {
        const Foam::label l = speciesTable[item.first()];
        if (l < 0 || l >= nS)
        {
            FatalErrorInFunction
                << "Reaction \"" << reactionKey << "\" (" << what << "): "
                << "unknown species \"" << item.first()
                << "\" in the third-body efficiency list."
                << exit(FatalError);
        }
        seen[l] = 1;
    }
    for (Foam::label i = 0; i < nS; ++i)
    {
        if (seen[i] == 0)
        {
            FatalErrorInFunction
                << "Reaction \"" << reactionKey << "\" (" << what << "): "
                << "missing third-body efficiency entry for species index "
                << i << " (every species must be listed explicitly)."
                << exit(FatalError);
        }
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
            if (fabs(val - round(val)) > 2.22e-16)
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
    gasPhaseReaction->AlignSpecies = alignSpecies;

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
                    FastChemistry::readElementaryReactionInfo
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
                    FastChemistry::readElementaryReactionInfo
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

    // Index of reverse rate constant part of non equilibrium reaction,
    //  these reactions uses Arrhenius model instead of equilibrium
    //  constant for reverse rate constant
    auto NERidx = gasPhaseReaction->Ikf[9];
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

            gasPhaseReaction->A[NERidx] = reverseDict.lookup<scalar>("A");       
            gasPhaseReaction->beta[NERidx] = reverseDict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[NERidx] = reverseDict.lookup<scalar>("Ta");
            if(isInteger==true)
            {
                FastChemistry::readElementaryReactionInfo
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
            iArrhenius++;NERidx++;
        }
    }

    // Index of third body factor
    unsigned int TBFidx = 0;
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

            gasPhaseReaction->A[NERidx] = reverseDict.lookup<scalar>("A");        
            gasPhaseReaction->beta[NERidx] = reverseDict.lookup<scalar>("beta");
            gasPhaseReaction->Ta[NERidx] = reverseDict.lookup<scalar>("Ta");     

            List<Tuple2<word, scalar>> forwardCoeffs(forwardDict.lookup("coeffs"));
            validateThirdBodyEfficiencies(forwardCoeffs, speciesTable, key, "third-body forward efficiencies");
            ThirdBodyFactor[TBFidx].resize(forwardCoeffs.size());
            forAll(forwardCoeffs, n)
            {
                const int l = speciesTable[(forwardCoeffs[n].first())];
                const scalar ThirdBodyFactor_n = forwardCoeffs[n].second();
                ThirdBodyFactor[TBFidx][l] = ThirdBodyFactor_n;
            }
            
            List<Tuple2<word, scalar>> reverseCoeffs(reverseDict.lookup("coeffs"));
            validateThirdBodyEfficiencies(reverseCoeffs, speciesTable, key, "third-body reverse efficiencies");

            auto begin = TBFidx + gasPhaseReaction->Itbr[4];
            

            ThirdBodyFactor[begin].resize(reverseCoeffs.size());
            forAll(reverseCoeffs, n)
            {
                const int l = speciesTable[(reverseCoeffs[n].first())];
                const scalar ThirdBodyFactor_n = reverseCoeffs[n].second();
                ThirdBodyFactor[begin][l] = ThirdBodyFactor_n;
            }
            if(isInteger==true)
            {
                FastChemistry::readElementaryReactionInfo
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

            // raiseError if reverse ThirdBodyFactor and forward ThirdBodyFactor;
            // is not same
            for(size_t i=0; i<ThirdBodyFactor[TBFidx].size();i++)
            {
                if(ThirdBodyFactor[TBFidx][i]!=ThirdBodyFactor[begin][i])
                {
                    Info<<"************FastChemistry Error************ "<<endl;
                    Info<<"different third body factor for forward part "<<endl;
                    Info<<"and reverse part is not allowed, if you need, "<<endl;
                    Info<<"split this reaction into two irreversible "<<endl;
                    Info<<"reactions. The invalid reaction is: "<<endl;
                    Info<<gasPhaseReaction->reactionEquation_[iArrhenius]<<endl;
                    FatalErrorInFunction<<exit(FatalError);
                }
            }
            iArrhenius++;NERidx++;TBFidx++;
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
            validateThirdBodyEfficiencies(coeffs, speciesTable, key, "third-body reaction");
            ThirdBodyFactor[TBFidx].resize(coeffs.size());
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[TBFidx][l] = ThirdBodyFactor_m;
            }
            if(isInteger==true)
            {
                FastChemistry::readElementaryReactionInfo
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
            iArrhenius++;TBFidx++;
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

            if (!reactDict.found("thirdBodyEfficiencies"))
            {
                FatalErrorInFunction
                    << "Reaction \"" << key << "\": required dictionary entry "
                    << "thirdBodyEfficiencies is missing."
                    << exit(FatalError);
            }
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
            validateThirdBodyEfficiencies(coeffs, speciesTable, key, "third-body reaction");
            ThirdBodyFactor[TBFidx].resize(coeffs.size());
    
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[TBFidx][l] = ThirdBodyFactor_m;
            }
            gasPhaseReaction->LindemannFO.push_back(iArrhenius);

            if(isInteger==true)
            {
                FastChemistry::readElementaryReactionInfo
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
            iArrhenius++;TBFidx++;
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
            if (!reactDict.found("thirdBodyEfficiencies"))
            {
                FatalErrorInFunction
                    << "Reaction \"" << key << "\": required dictionary entry "
                    << "thirdBodyEfficiencies is missing."
                    << exit(FatalError);
            }
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
            validateThirdBodyEfficiencies(coeffs, speciesTable, key, "third-body reaction");
            ThirdBodyFactor[TBFidx].resize(coeffs.size());
    
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[TBFidx][l] = ThirdBodyFactor_m;
            }

            gasPhaseReaction->TroeFO.push_back(iArrhenius);
            gasPhaseReaction->alpha_.push_back(FDict.lookup<scalar>("alpha"));    
            gasPhaseReaction->Ts_.push_back(FDict.lookup<scalar>("Ts"));    
            gasPhaseReaction->Tss_.push_back(FDict.lookup<scalar>("Tss"));    
            gasPhaseReaction->Tsss_.push_back(FDict.lookup<scalar>("Tsss"));
            if(isInteger==true)
            {
                FastChemistry::readElementaryReactionInfo
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
            iArrhenius++;TBFidx++;
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
            if (!reactDict.found("thirdBodyEfficiencies"))
            {
                FatalErrorInFunction
                    << "Reaction \"" << key << "\": required dictionary entry "
                    << "thirdBodyEfficiencies is missing."
                    << exit(FatalError);
            }
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
            validateThirdBodyEfficiencies(coeffs, speciesTable, key, "third-body reaction");
            ThirdBodyFactor[TBFidx].resize(coeffs.size());
    
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[TBFidx][l] = ThirdBodyFactor_m;
            }
            
            gasPhaseReaction->SRIFO.push_back(iArrhenius);
            gasPhaseReaction->a_.push_back(FDict.lookup<scalar>("a"));    
            gasPhaseReaction->b_.push_back(FDict.lookup<scalar>("b"));    
            gasPhaseReaction->c_.push_back(FDict.lookup<scalar>("c"));    
            gasPhaseReaction->d_.push_back(FDict.lookup<scalar>("d"));  
            gasPhaseReaction->e_.push_back(FDict.lookup<scalar>("e"));  
            
            if(isInteger==true)
            {
                FastChemistry::readElementaryReactionInfo
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
            iArrhenius++;TBFidx++;
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

            if (!reactDict.found("thirdBodyEfficiencies"))
            {
                FatalErrorInFunction
                    << "Reaction \"" << key << "\": required dictionary entry "
                    << "thirdBodyEfficiencies is missing."
                    << exit(FatalError);
            }
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
            validateThirdBodyEfficiencies(coeffs, speciesTable, key, "third-body reaction");
            ThirdBodyFactor[TBFidx].resize(coeffs.size());
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[TBFidx][l] = ThirdBodyFactor_m;
            }
            gasPhaseReaction->LindemannCA.push_back(iArrhenius);

            if(isInteger==true)
            {
                FastChemistry::readElementaryReactionInfo
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
            TBFidx++;
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
            if (!reactDict.found("thirdBodyEfficiencies"))
            {
                FatalErrorInFunction
                    << "Reaction \"" << key << "\": required dictionary entry "
                    << "thirdBodyEfficiencies is missing."
                    << exit(FatalError);
            }
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
            validateThirdBodyEfficiencies(coeffs, speciesTable, key, "third-body reaction");
            ThirdBodyFactor[TBFidx].resize(coeffs.size());
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[TBFidx][l] = ThirdBodyFactor_m;
            }

            gasPhaseReaction->TroeCA.push_back(iArrhenius);
            gasPhaseReaction->alpha_.push_back(FDict.lookup<scalar>("alpha"));    
            gasPhaseReaction->Ts_.push_back(FDict.lookup<scalar>("Ts"));    
            gasPhaseReaction->Tss_.push_back(FDict.lookup<scalar>("Tss"));    
            gasPhaseReaction->Tsss_.push_back(FDict.lookup<scalar>("Tsss"));                
            if(isInteger==true)
            {
                FastChemistry::readElementaryReactionInfo
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
            TBFidx++;
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
            if (!reactDict.found("thirdBodyEfficiencies"))
            {
                FatalErrorInFunction
                    << "Reaction \"" << key << "\": required dictionary entry "
                    << "thirdBodyEfficiencies is missing."
                    << exit(FatalError);
            }
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
            validateThirdBodyEfficiencies(coeffs, speciesTable, key, "third-body reaction");
            ThirdBodyFactor[TBFidx].resize(coeffs.size());
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[TBFidx][l] = ThirdBodyFactor_m;
            }
            gasPhaseReaction->SRICA.push_back(iArrhenius);
            gasPhaseReaction->a_.push_back(FDict.lookup<scalar>("a"));    
            gasPhaseReaction->b_.push_back(FDict.lookup<scalar>("b"));    
            gasPhaseReaction->c_.push_back(FDict.lookup<scalar>("c"));    
            gasPhaseReaction->d_.push_back(FDict.lookup<scalar>("d"));  
            gasPhaseReaction->e_.push_back(FDict.lookup<scalar>("e"));
            if(isInteger==true)
            {
                FastChemistry::readElementaryReactionInfo
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
            TBFidx++;
        }
    }


    // Set third body factor
    {


        //gasPhaseReaction->AlignSpecies = gasPhaseReaction->nSpecies+remain;

        


        std::size_t totalSize = ThirdBodyFactor.size()*
            gasPhaseReaction->AlignSpecies*sizeof(double);
        if 
        (
            posix_memalign
            (
                reinterpret_cast<void**>(&gasPhaseReaction->ThirdBodyFactor1D), 
                32, 
                totalSize
            )
        )
        {
            throw std::bad_alloc();
        }
        memset(gasPhaseReaction->ThirdBodyFactor1D, 0, totalSize);

        auto& ThirdBodyFactor1Dref = gasPhaseReaction->ThirdBodyFactor1D;
        unsigned int nCols = gasPhaseReaction->AlignSpecies;
        for(unsigned int i = 0; i < ThirdBodyFactor.size();i++)
        {
            for(unsigned int J = 0; J < ThirdBodyFactor[i].size();J++)
            {
                ThirdBodyFactor1Dref[i*nCols+J] = ThirdBodyFactor[i][J];
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
                    FatalErrorInFunction<< "Pre-factor A should be larger "
                        << "than zero for Plog reaction"
                        << exit(FatalError);
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
                FastChemistry::readElementaryReactionInfo
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
            unsigned int sp0 = gasPhaseReaction->lhsSpeciesIndex[i][0];
            unsigned int sp1 = gasPhaseReaction->rhsSpeciesIndex[i][0];
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                    gasPhaseReaction->RR11AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR11AllIndex.push_back(sp0);
                    gasPhaseReaction->RR11AllIndex.push_back(sp1);
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                    gasPhaseReaction->IR11AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR11AllIndex.push_back(sp0);
                    gasPhaseReaction->IR11AllIndex.push_back(sp1);
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                gasPhaseReaction->NER11AllIndex.push_back(static_cast<unsigned int>(i));
                gasPhaseReaction->NER11AllIndex.push_back(sp0);
                gasPhaseReaction->NER11AllIndex.push_back(sp1);
            }
        }
        else if(lhsNumber==1 && rhsNumber==2)
        {
            unsigned int sp0 = gasPhaseReaction->lhsSpeciesIndex[i][0];
            unsigned int sp1 = gasPhaseReaction->rhsSpeciesIndex[i][0];
            unsigned int sp2 = gasPhaseReaction->rhsSpeciesIndex[i][1];
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                if(sp1==sp2)
                {
                    gasPhaseReaction->RR12DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR12DupAllIndex.push_back(sp0);
                    gasPhaseReaction->RR12DupAllIndex.push_back(sp2);
                }
                else
                {
                    gasPhaseReaction->RR12AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR12AllIndex.push_back(sp0);
                    gasPhaseReaction->RR12AllIndex.push_back(sp1);
                    gasPhaseReaction->RR12AllIndex.push_back(sp2);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                if(sp1==sp2)
                {
                    gasPhaseReaction->IR12DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR12DupAllIndex.push_back(sp0);
                    gasPhaseReaction->IR12DupAllIndex.push_back(sp2);
                }
                else
                {
                    gasPhaseReaction->IR12AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR12AllIndex.push_back(sp0);
                    gasPhaseReaction->IR12AllIndex.push_back(sp1);
                    gasPhaseReaction->IR12AllIndex.push_back(sp2);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                if(sp1==sp2)
                {
                    gasPhaseReaction->NER12DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER12DupAllIndex.push_back(sp0);
                    gasPhaseReaction->NER12DupAllIndex.push_back(sp2);
                }
                else
                {
                    gasPhaseReaction->NER12AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER12AllIndex.push_back(sp0);
                    gasPhaseReaction->NER12AllIndex.push_back(sp1);
                    gasPhaseReaction->NER12AllIndex.push_back(sp2);
                }
            }
        }
        else if(lhsNumber==1 && rhsNumber==3)
        {
            unsigned int sp0 = gasPhaseReaction->lhsSpeciesIndex[i][0];
            unsigned int sp1 = gasPhaseReaction->rhsSpeciesIndex[i][0];
            unsigned int sp2 = gasPhaseReaction->rhsSpeciesIndex[i][1];
            unsigned int sp3 = gasPhaseReaction->rhsSpeciesIndex[i][2];
            if(sp1==sp2&&sp2==sp3)
            {
                FatalErrorInFunction<< "unsupport reaction type:"
                    << "A=B+B+B" << exit(FatalError);
            }


            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                if(sp1==sp2)
                {
                    gasPhaseReaction->RR13DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR13DupAllIndex.push_back(sp0);
                    gasPhaseReaction->RR13DupAllIndex.push_back(sp3);
                    gasPhaseReaction->RR13DupAllIndex.push_back(sp1);
                }
                else if(sp1==sp3)
                {
                    gasPhaseReaction->RR13DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR13DupAllIndex.push_back(sp0);
                    gasPhaseReaction->RR13DupAllIndex.push_back(sp2);
                    gasPhaseReaction->RR13DupAllIndex.push_back(sp1);
                }
                else if(sp2==sp3)
                {
                    gasPhaseReaction->RR13DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR13DupAllIndex.push_back(sp0);
                    gasPhaseReaction->RR13DupAllIndex.push_back(sp1);
                    gasPhaseReaction->RR13DupAllIndex.push_back(sp2);
                }
                else
                {
                    gasPhaseReaction->RR13AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR13AllIndex.push_back(sp0);
                    gasPhaseReaction->RR13AllIndex.push_back(sp1);
                    gasPhaseReaction->RR13AllIndex.push_back(sp2);
                    gasPhaseReaction->RR13AllIndex.push_back(sp3);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                if(sp1==sp2)
                {
                    gasPhaseReaction->IR13DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR13DupAllIndex.push_back(sp0);
                    gasPhaseReaction->IR13DupAllIndex.push_back(sp3);
                    gasPhaseReaction->IR13DupAllIndex.push_back(sp1);
                }
                else if(sp1==sp3)
                {
                    gasPhaseReaction->IR13DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR13DupAllIndex.push_back(sp0);
                    gasPhaseReaction->IR13DupAllIndex.push_back(sp2);
                    gasPhaseReaction->IR13DupAllIndex.push_back(sp1);
                }
                else if(sp2==sp3)
                {
                    gasPhaseReaction->IR13DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR13DupAllIndex.push_back(sp0);
                    gasPhaseReaction->IR13DupAllIndex.push_back(sp1);
                    gasPhaseReaction->IR13DupAllIndex.push_back(sp2);
                }
                else
                {
                    gasPhaseReaction->IR13AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR13AllIndex.push_back(sp0);
                    gasPhaseReaction->IR13AllIndex.push_back(sp1);
                    gasPhaseReaction->IR13AllIndex.push_back(sp2);
                    gasPhaseReaction->IR13AllIndex.push_back(sp3);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                if(sp1==sp2)
                {
                    gasPhaseReaction->NER13DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER13DupAllIndex.push_back(sp0);
                    gasPhaseReaction->NER13DupAllIndex.push_back(sp3);
                    gasPhaseReaction->NER13DupAllIndex.push_back(sp1);
                }
                else if(sp1==sp3)
                {
                    gasPhaseReaction->NER13DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER13DupAllIndex.push_back(sp0);
                    gasPhaseReaction->NER13DupAllIndex.push_back(sp2);
                    gasPhaseReaction->NER13DupAllIndex.push_back(sp1);
                }
                else if(sp2==sp3)
                {
                    gasPhaseReaction->NER13DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER13DupAllIndex.push_back(sp0);
                    gasPhaseReaction->NER13DupAllIndex.push_back(sp1);
                    gasPhaseReaction->NER13DupAllIndex.push_back(sp2);
                }
                else
                {
                    gasPhaseReaction->NER13AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER13AllIndex.push_back(sp0);
                    gasPhaseReaction->NER13AllIndex.push_back(sp1);
                    gasPhaseReaction->NER13AllIndex.push_back(sp2);
                    gasPhaseReaction->NER13AllIndex.push_back(sp3);
                }
            }
        }
        else if(lhsNumber==2 && rhsNumber==1)
        {
            unsigned int sp0 = gasPhaseReaction->lhsSpeciesIndex[i][0];
            unsigned int sp1 = gasPhaseReaction->lhsSpeciesIndex[i][1];
            unsigned int sp2 = gasPhaseReaction->rhsSpeciesIndex[i][0];
            std::set<unsigned int> numbers = {sp0, sp1, sp2};
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                if(numbers.size()<3)
                {
                    gasPhaseReaction->RR21DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR21DupAllIndex.push_back(sp0);
                    gasPhaseReaction->RR21DupAllIndex.push_back(sp2);
                }
                else
                {
                    gasPhaseReaction->RR21AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR21AllIndex.push_back(sp0);
                    gasPhaseReaction->RR21AllIndex.push_back(sp1);
                    gasPhaseReaction->RR21AllIndex.push_back(sp2);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                if(numbers.size()<3)
                {
                    gasPhaseReaction->IR21DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR21DupAllIndex.push_back(sp0);
                    gasPhaseReaction->IR21DupAllIndex.push_back(sp2);
                }
                else
                {
                    gasPhaseReaction->IR21AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR21AllIndex.push_back(sp0);
                    gasPhaseReaction->IR21AllIndex.push_back(sp1);
                    gasPhaseReaction->IR21AllIndex.push_back(sp2);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                if(numbers.size()<3)
                {
                    gasPhaseReaction->NER21DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER21DupAllIndex.push_back(sp0);
                    gasPhaseReaction->NER21DupAllIndex.push_back(sp2);
                }
                else
                {
                    gasPhaseReaction->NER21AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER21AllIndex.push_back(sp0);
                    gasPhaseReaction->NER21AllIndex.push_back(sp1);
                    gasPhaseReaction->NER21AllIndex.push_back(sp2);
                }
            }
        }
        else if(lhsNumber==2 && rhsNumber==2)
        {
            unsigned int sp0 = gasPhaseReaction->lhsSpeciesIndex[i][0];
            unsigned int sp1 = gasPhaseReaction->lhsSpeciesIndex[i][1];
            unsigned int sp2 = gasPhaseReaction->rhsSpeciesIndex[i][0];
            unsigned int sp3 = gasPhaseReaction->rhsSpeciesIndex[i][1];
            std::set<unsigned int> numbers = {sp0, sp1, sp2, sp3};
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                if(numbers.size()<4)
                {
                    gasPhaseReaction->RR22DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR22DupAllIndex.push_back(sp0);
                    gasPhaseReaction->RR22DupAllIndex.push_back(sp1);
                    gasPhaseReaction->RR22DupAllIndex.push_back(sp2);
                    gasPhaseReaction->RR22DupAllIndex.push_back(sp3);
                }
                else
                {
                    gasPhaseReaction->RR22AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR22AllIndex.push_back(sp0);
                    gasPhaseReaction->RR22AllIndex.push_back(sp1);
                    gasPhaseReaction->RR22AllIndex.push_back(sp2);
                    gasPhaseReaction->RR22AllIndex.push_back(sp3);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                if(numbers.size()<4)
                {
                    gasPhaseReaction->IR22DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR22DupAllIndex.push_back(sp0);
                    gasPhaseReaction->IR22DupAllIndex.push_back(sp1);
                    gasPhaseReaction->IR22DupAllIndex.push_back(sp2);
                    gasPhaseReaction->IR22DupAllIndex.push_back(sp3);
                }
                else
                {
                    gasPhaseReaction->IR22AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR22AllIndex.push_back(sp0);
                    gasPhaseReaction->IR22AllIndex.push_back(sp1);
                    gasPhaseReaction->IR22AllIndex.push_back(sp2);
                    gasPhaseReaction->IR22AllIndex.push_back(sp3);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                if(numbers.size()<4)
                {
                    gasPhaseReaction->NER22DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER22DupAllIndex.push_back(sp0);
                    gasPhaseReaction->NER22DupAllIndex.push_back(sp1);
                    gasPhaseReaction->NER22DupAllIndex.push_back(sp2);
                    gasPhaseReaction->NER22DupAllIndex.push_back(sp3);
                }
                else
                {
                    gasPhaseReaction->NER22AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER22AllIndex.push_back(sp0);
                    gasPhaseReaction->NER22AllIndex.push_back(sp1);
                    gasPhaseReaction->NER22AllIndex.push_back(sp2);
                    gasPhaseReaction->NER22AllIndex.push_back(sp3);
                }
            }
        }
        else if(lhsNumber==2 && rhsNumber==3)
        {

            unsigned int sp0 = gasPhaseReaction->lhsSpeciesIndex[i][0];
            unsigned int sp1 = gasPhaseReaction->lhsSpeciesIndex[i][1];
            unsigned int sp2 = gasPhaseReaction->rhsSpeciesIndex[i][0];
            unsigned int sp3 = gasPhaseReaction->rhsSpeciesIndex[i][1];
            unsigned int sp4 = gasPhaseReaction->rhsSpeciesIndex[i][2];
            std::set<unsigned int> numbers = {sp0, sp1, sp2, sp3, sp4};

            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                if(numbers.size()<5)
                {
                    gasPhaseReaction->RR23DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR23DupAllIndex.push_back(sp0);
                    gasPhaseReaction->RR23DupAllIndex.push_back(sp1);
                    gasPhaseReaction->RR23DupAllIndex.push_back(sp2);
                    gasPhaseReaction->RR23DupAllIndex.push_back(sp3);
                    gasPhaseReaction->RR23DupAllIndex.push_back(sp4);
                }
                else
                {
                    gasPhaseReaction->RR23AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR23AllIndex.push_back(sp0);
                    gasPhaseReaction->RR23AllIndex.push_back(sp1);
                    gasPhaseReaction->RR23AllIndex.push_back(sp2);
                    gasPhaseReaction->RR23AllIndex.push_back(sp3);
                    gasPhaseReaction->RR23AllIndex.push_back(sp4);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                if(numbers.size()<5)
                {
                    gasPhaseReaction->IR23DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR23DupAllIndex.push_back(sp0);
                    gasPhaseReaction->IR23DupAllIndex.push_back(sp1);
                    gasPhaseReaction->IR23DupAllIndex.push_back(sp2);
                    gasPhaseReaction->IR23DupAllIndex.push_back(sp3);
                    gasPhaseReaction->IR23DupAllIndex.push_back(sp4);
                }
                else
                {
                    gasPhaseReaction->IR23AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR23AllIndex.push_back(sp0);
                    gasPhaseReaction->IR23AllIndex.push_back(sp1);
                    gasPhaseReaction->IR23AllIndex.push_back(sp2);
                    gasPhaseReaction->IR23AllIndex.push_back(sp3);
                    gasPhaseReaction->IR23AllIndex.push_back(sp4);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                if(numbers.size()<5)
                {
                    gasPhaseReaction->NER23DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER23DupAllIndex.push_back(sp0);
                    gasPhaseReaction->NER23DupAllIndex.push_back(sp1);
                    gasPhaseReaction->NER23DupAllIndex.push_back(sp2);
                    gasPhaseReaction->NER23DupAllIndex.push_back(sp3);
                    gasPhaseReaction->NER23DupAllIndex.push_back(sp4);
                }
                else
                {
                    gasPhaseReaction->NER23AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER23AllIndex.push_back(sp0);
                    gasPhaseReaction->NER23AllIndex.push_back(sp1);
                    gasPhaseReaction->NER23AllIndex.push_back(sp2);
                    gasPhaseReaction->NER23AllIndex.push_back(sp3);
                    gasPhaseReaction->NER23AllIndex.push_back(sp4);
                }
            }
        }
        else if(lhsNumber==3 && rhsNumber==1)
        {
            unsigned int sp0 = gasPhaseReaction->lhsSpeciesIndex[i][0];
            unsigned int sp1 = gasPhaseReaction->lhsSpeciesIndex[i][1];
            unsigned int sp2 = gasPhaseReaction->lhsSpeciesIndex[i][2];
            unsigned int sp3 = gasPhaseReaction->rhsSpeciesIndex[i][0];
            std::set<unsigned int> numbers = {sp0, sp1, sp2, sp3};
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                if(numbers.size()<4)
                {
                    gasPhaseReaction->RR31DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR31DupAllIndex.push_back(sp0);
                    gasPhaseReaction->RR31DupAllIndex.push_back(sp1);
                    gasPhaseReaction->RR31DupAllIndex.push_back(sp2);
                    gasPhaseReaction->RR31DupAllIndex.push_back(sp3);
                }
                else
                {
                    gasPhaseReaction->RR31AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR31AllIndex.push_back(sp0);
                    gasPhaseReaction->RR31AllIndex.push_back(sp1);
                    gasPhaseReaction->RR31AllIndex.push_back(sp2);
                    gasPhaseReaction->RR31AllIndex.push_back(sp3);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                if(numbers.size()<4)
                {
                    gasPhaseReaction->IR31DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR31DupAllIndex.push_back(sp0);
                    gasPhaseReaction->IR31DupAllIndex.push_back(sp1);
                    gasPhaseReaction->IR31DupAllIndex.push_back(sp2);
                    gasPhaseReaction->IR31DupAllIndex.push_back(sp3);
                }
                else
                {
                    gasPhaseReaction->IR31AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR31AllIndex.push_back(sp0);
                    gasPhaseReaction->IR31AllIndex.push_back(sp1);
                    gasPhaseReaction->IR31AllIndex.push_back(sp2);
                    gasPhaseReaction->IR31AllIndex.push_back(sp3);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                if(numbers.size()<4)
                {
                    gasPhaseReaction->NER31DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER31DupAllIndex.push_back(sp0);
                    gasPhaseReaction->NER31DupAllIndex.push_back(sp1);
                    gasPhaseReaction->NER31DupAllIndex.push_back(sp2);
                    gasPhaseReaction->NER31DupAllIndex.push_back(sp3);
                }
                else
                {
                    gasPhaseReaction->NER31AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER31AllIndex.push_back(sp0);
                    gasPhaseReaction->NER31AllIndex.push_back(sp1);
                    gasPhaseReaction->NER31AllIndex.push_back(sp2);
                    gasPhaseReaction->NER31AllIndex.push_back(sp3);
                }
            }

        }
        else if(lhsNumber==3 && rhsNumber==2)
        {
            unsigned int sp0 = gasPhaseReaction->lhsSpeciesIndex[i][0];
            unsigned int sp1 = gasPhaseReaction->lhsSpeciesIndex[i][1];
            unsigned int sp2 = gasPhaseReaction->lhsSpeciesIndex[i][2];
            unsigned int sp3 = gasPhaseReaction->rhsSpeciesIndex[i][0];
            unsigned int sp4 = gasPhaseReaction->rhsSpeciesIndex[i][1];
            std::set<unsigned int> numbers = {sp0, sp1, sp2, sp3, sp4};
            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                if(numbers.size()<5)
                {
                    gasPhaseReaction->RR32DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR32DupAllIndex.push_back(sp0);
                    gasPhaseReaction->RR32DupAllIndex.push_back(sp1);
                    gasPhaseReaction->RR32DupAllIndex.push_back(sp2);
                    gasPhaseReaction->RR32DupAllIndex.push_back(sp3);
                    gasPhaseReaction->RR32DupAllIndex.push_back(sp4);
                }
                else
                {
                    gasPhaseReaction->RR32AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR32AllIndex.push_back(sp0);
                    gasPhaseReaction->RR32AllIndex.push_back(sp1);
                    gasPhaseReaction->RR32AllIndex.push_back(sp2);
                    gasPhaseReaction->RR32AllIndex.push_back(sp3);
                    gasPhaseReaction->RR32AllIndex.push_back(sp4);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                if(numbers.size()<5)
                {
                    gasPhaseReaction->IR32DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR32DupAllIndex.push_back(sp0);
                    gasPhaseReaction->IR32DupAllIndex.push_back(sp1);
                    gasPhaseReaction->IR32DupAllIndex.push_back(sp2);
                    gasPhaseReaction->IR32DupAllIndex.push_back(sp3);
                    gasPhaseReaction->IR32DupAllIndex.push_back(sp4);
                }
                else
                {
                    gasPhaseReaction->IR32AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR32AllIndex.push_back(sp0);
                    gasPhaseReaction->IR32AllIndex.push_back(sp1);
                    gasPhaseReaction->IR32AllIndex.push_back(sp2);
                    gasPhaseReaction->IR32AllIndex.push_back(sp3);
                    gasPhaseReaction->IR32AllIndex.push_back(sp4);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                if(numbers.size()<5)
                {
                    gasPhaseReaction->NER32DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER32DupAllIndex.push_back(sp0);
                    gasPhaseReaction->NER32DupAllIndex.push_back(sp1);
                    gasPhaseReaction->NER32DupAllIndex.push_back(sp2);
                    gasPhaseReaction->NER32DupAllIndex.push_back(sp3);
                    gasPhaseReaction->NER32DupAllIndex.push_back(sp4);
                }
                else
                {
                    gasPhaseReaction->NER32AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER32AllIndex.push_back(sp0);
                    gasPhaseReaction->NER32AllIndex.push_back(sp1);
                    gasPhaseReaction->NER32AllIndex.push_back(sp2);
                    gasPhaseReaction->NER32AllIndex.push_back(sp3);
                    gasPhaseReaction->NER32AllIndex.push_back(sp4);
                }
            }

        }
        else if(lhsNumber==3 && rhsNumber==3)
        {
            unsigned int sp0 = gasPhaseReaction->lhsSpeciesIndex[i][0];
            unsigned int sp1 = gasPhaseReaction->lhsSpeciesIndex[i][1];
            unsigned int sp2 = gasPhaseReaction->lhsSpeciesIndex[i][2];
            unsigned int sp3 = gasPhaseReaction->rhsSpeciesIndex[i][0];
            unsigned int sp4 = gasPhaseReaction->rhsSpeciesIndex[i][1];
            unsigned int sp5 = gasPhaseReaction->rhsSpeciesIndex[i][2];
            std::set<unsigned int> numbers = {sp0, sp1, sp2, sp3, sp4, sp5};

            if(gasPhaseReaction->isIrreversible[i]==0)
            {
                if(numbers.size()<6)
                {
                    gasPhaseReaction->RR33DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR33DupAllIndex.push_back(sp0);
                    gasPhaseReaction->RR33DupAllIndex.push_back(sp1);
                    gasPhaseReaction->RR33DupAllIndex.push_back(sp2);
                    gasPhaseReaction->RR33DupAllIndex.push_back(sp3);
                    gasPhaseReaction->RR33DupAllIndex.push_back(sp4);
                    gasPhaseReaction->RR33DupAllIndex.push_back(sp5);
                }
                else
                {
                    gasPhaseReaction->RR33AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->RR33AllIndex.push_back(sp0);
                    gasPhaseReaction->RR33AllIndex.push_back(sp1);
                    gasPhaseReaction->RR33AllIndex.push_back(sp2);
                    gasPhaseReaction->RR33AllIndex.push_back(sp3);
                    gasPhaseReaction->RR33AllIndex.push_back(sp4);
                    gasPhaseReaction->RR33AllIndex.push_back(sp5);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==1)
            {
                if(numbers.size()<6)
                {
                    gasPhaseReaction->IR33DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR33DupAllIndex.push_back(sp0);
                    gasPhaseReaction->IR33DupAllIndex.push_back(sp1);
                    gasPhaseReaction->IR33DupAllIndex.push_back(sp2);
                    gasPhaseReaction->IR33DupAllIndex.push_back(sp3);
                    gasPhaseReaction->IR33DupAllIndex.push_back(sp4);
                    gasPhaseReaction->IR33DupAllIndex.push_back(sp5);
                }
                else
                {
                    gasPhaseReaction->IR33AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->IR33AllIndex.push_back(sp0);
                    gasPhaseReaction->IR33AllIndex.push_back(sp1);
                    gasPhaseReaction->IR33AllIndex.push_back(sp2);
                    gasPhaseReaction->IR33AllIndex.push_back(sp3);
                    gasPhaseReaction->IR33AllIndex.push_back(sp4);
                    gasPhaseReaction->IR33AllIndex.push_back(sp5);
                }
            }
            else if(gasPhaseReaction->isIrreversible[i]==2)
            {
                if(numbers.size()<6)
                {
                    gasPhaseReaction->NER33DupAllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER33DupAllIndex.push_back(sp0);
                    gasPhaseReaction->NER33DupAllIndex.push_back(sp1);
                    gasPhaseReaction->NER33DupAllIndex.push_back(sp2);
                    gasPhaseReaction->NER33DupAllIndex.push_back(sp3);
                    gasPhaseReaction->NER33DupAllIndex.push_back(sp4);
                    gasPhaseReaction->NER33DupAllIndex.push_back(sp5);
                }
                else
                {
                    gasPhaseReaction->NER33AllIndex.push_back(static_cast<unsigned int>(i));
                    gasPhaseReaction->NER33AllIndex.push_back(sp0);
                    gasPhaseReaction->NER33AllIndex.push_back(sp1);
                    gasPhaseReaction->NER33AllIndex.push_back(sp2);
                    gasPhaseReaction->NER33AllIndex.push_back(sp3);
                    gasPhaseReaction->NER33AllIndex.push_back(sp4);
                    gasPhaseReaction->NER33AllIndex.push_back(sp5);
                }
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
        if(gasPhaseReaction->RR11AllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF11RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF11RR);
        }
        if(gasPhaseReaction->IR11AllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF11IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF11IR);
        }
        if(gasPhaseReaction->NER11AllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF11NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF11NER);
        }
        if(gasPhaseReaction->RR12AllIndex.size()>0||gasPhaseReaction->RR12DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF12RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF12RR);
        }
        if(gasPhaseReaction->IR12AllIndex.size()>0||gasPhaseReaction->IR12DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF12IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF12IR);
        }
        if(gasPhaseReaction->NER12AllIndex.size()>0||gasPhaseReaction->NER12DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF12NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF12NER);
        }
        if(gasPhaseReaction->RR13AllIndex.size()>0||gasPhaseReaction->RR13DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF13RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF13RR);
        }
        if(gasPhaseReaction->IR13AllIndex.size()>0||gasPhaseReaction->IR13DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF13IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF13IR);
        }
        if(gasPhaseReaction->NER13AllIndex.size()>0||gasPhaseReaction->NER13DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF13NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF13NER);
        }
        if(gasPhaseReaction->RR21AllIndex.size()>0||gasPhaseReaction->RR21DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF21RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF21RR);
        }
        if(gasPhaseReaction->IR21AllIndex.size()>0||gasPhaseReaction->IR21DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF21IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF21IR);
        }
        if(gasPhaseReaction->NER21AllIndex.size()>0||gasPhaseReaction->NER21DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF21NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF21NER);
        }
        if(gasPhaseReaction->RR22AllIndex.size()>0||gasPhaseReaction->RR22DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF22RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF22RR);
        }
        if(gasPhaseReaction->IR22AllIndex.size()>0||gasPhaseReaction->IR22DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF22IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF22IR);
        }
        if(gasPhaseReaction->NER22AllIndex.size()>0||gasPhaseReaction->NER22DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF22NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF22NER);
        }
        if(gasPhaseReaction->RR23AllIndex.size()>0||gasPhaseReaction->RR23DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF23RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF23RR);
        }
        if(gasPhaseReaction->IR23AllIndex.size()>0||gasPhaseReaction->IR23DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF23IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF23IR);
        }
        if(gasPhaseReaction->NER23AllIndex.size()>0||gasPhaseReaction->NER23DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF23NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF23NER);
        }

        if(gasPhaseReaction->RR31AllIndex.size()>0||gasPhaseReaction->RR31DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF31RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF31RR);
        }
        if(gasPhaseReaction->IR31AllIndex.size()>0||gasPhaseReaction->IR31DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF31IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF31IR);
        }
        if(gasPhaseReaction->NER31AllIndex.size()>0||gasPhaseReaction->NER31DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF31NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF31NER);
        }

        if(gasPhaseReaction->RR32AllIndex.size()>0||gasPhaseReaction->RR32DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF32RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF32RR);
        }
        if(gasPhaseReaction->IR32AllIndex.size()>0||gasPhaseReaction->IR32DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF32IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF32IR);
        }
        if(gasPhaseReaction->NER32AllIndex.size()>0||gasPhaseReaction->NER32DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF32NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF32NER);
        }
        if(gasPhaseReaction->RR33AllIndex.size()>0||gasPhaseReaction->RR33DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF33RR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF33RR);
        }
        if(gasPhaseReaction->IR33AllIndex.size()>0||gasPhaseReaction->IR33DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF33IR);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF33IR);
        }
        if(gasPhaseReaction->NER33AllIndex.size()>0||gasPhaseReaction->NER33DupAllIndex.size()>0)
        {
            gasPhaseReaction->RFptr.push_back(&OptReaction::RF33NER);
            gasPhaseReaction->JFptr.push_back(&OptReaction::JF33NER);
        }
    }


    gasPhaseReaction->RR11Size = static_cast<unsigned int>(gasPhaseReaction->RR11AllIndex.size()/3);
    gasPhaseReaction->IR11Size = static_cast<unsigned int>(gasPhaseReaction->IR11AllIndex.size()/3);
    gasPhaseReaction->NER11Size = static_cast<unsigned int>(gasPhaseReaction->NER11AllIndex.size()/3);

    gasPhaseReaction->RR12Size = static_cast<unsigned int>(gasPhaseReaction->RR12AllIndex.size()/4);
    gasPhaseReaction->RR12DupSize = static_cast<unsigned int>(gasPhaseReaction->RR12DupAllIndex.size()/3);
    gasPhaseReaction->IR12Size = static_cast<unsigned int>(gasPhaseReaction->IR12AllIndex.size()/4);
    gasPhaseReaction->IR12DupSize = static_cast<unsigned int>(gasPhaseReaction->IR12DupAllIndex.size()/3);
    gasPhaseReaction->NER12Size = static_cast<unsigned int>(gasPhaseReaction->NER12AllIndex.size()/4);
    gasPhaseReaction->NER12DupSize = static_cast<unsigned int>(gasPhaseReaction->NER12DupAllIndex.size()/3);

    gasPhaseReaction->RR13Size = static_cast<unsigned int>(gasPhaseReaction->RR13AllIndex.size()/5);
    gasPhaseReaction->RR13DupSize = static_cast<unsigned int>(gasPhaseReaction->RR13DupAllIndex.size()/4);
    gasPhaseReaction->IR13Size = static_cast<unsigned int>(gasPhaseReaction->IR13AllIndex.size()/5);
    gasPhaseReaction->IR13DupSize = static_cast<unsigned int>(gasPhaseReaction->IR13DupAllIndex.size()/4);
    gasPhaseReaction->NER13Size = static_cast<unsigned int>(gasPhaseReaction->NER13AllIndex.size()/5);
    gasPhaseReaction->NER13DupSize = static_cast<unsigned int>(gasPhaseReaction->NER13DupAllIndex.size()/4);

    gasPhaseReaction->RR21Size = static_cast<unsigned int>(gasPhaseReaction->RR21AllIndex.size()/4);
    gasPhaseReaction->RR21DupSize = static_cast<unsigned int>(gasPhaseReaction->RR21DupAllIndex.size()/3);
    gasPhaseReaction->IR21Size = static_cast<unsigned int>(gasPhaseReaction->IR21AllIndex.size()/4);
    gasPhaseReaction->IR21DupSize = static_cast<unsigned int>(gasPhaseReaction->IR21DupAllIndex.size()/3);
    gasPhaseReaction->NER21Size = static_cast<unsigned int>(gasPhaseReaction->NER21AllIndex.size()/4);
    gasPhaseReaction->NER21DupSize = static_cast<unsigned int>(gasPhaseReaction->NER21DupAllIndex.size()/3);

    gasPhaseReaction->RR22Size = static_cast<unsigned int>(gasPhaseReaction->RR22AllIndex.size()/5);
    gasPhaseReaction->RR22DupSize = static_cast<unsigned int>(gasPhaseReaction->RR22DupAllIndex.size()/5);
    gasPhaseReaction->IR22Size = static_cast<unsigned int>(gasPhaseReaction->IR22AllIndex.size()/5);
    gasPhaseReaction->IR22DupSize = static_cast<unsigned int>(gasPhaseReaction->IR22DupAllIndex.size()/5);
    gasPhaseReaction->NER22Size = static_cast<unsigned int>(gasPhaseReaction->NER22AllIndex.size()/5);
    gasPhaseReaction->NER22DupSize = static_cast<unsigned int>(gasPhaseReaction->NER22DupAllIndex.size()/5);

    gasPhaseReaction->RR23Size = static_cast<unsigned int>(gasPhaseReaction->RR23AllIndex.size()/6);
    gasPhaseReaction->RR23DupSize = static_cast<unsigned int>(gasPhaseReaction->RR23DupAllIndex.size()/6);
    gasPhaseReaction->IR23Size = static_cast<unsigned int>(gasPhaseReaction->IR23AllIndex.size()/6);
    gasPhaseReaction->IR23DupSize = static_cast<unsigned int>(gasPhaseReaction->IR23DupAllIndex.size()/6);
    gasPhaseReaction->NER23Size = static_cast<unsigned int>(gasPhaseReaction->NER23AllIndex.size()/6);
    gasPhaseReaction->NER23DupSize = static_cast<unsigned int>(gasPhaseReaction->NER23DupAllIndex.size()/6);

    gasPhaseReaction->RR32Size = static_cast<unsigned int>(gasPhaseReaction->RR32AllIndex.size()/6);
    gasPhaseReaction->RR32DupSize = static_cast<unsigned int>(gasPhaseReaction->RR32DupAllIndex.size()/6);
    gasPhaseReaction->IR32Size = static_cast<unsigned int>(gasPhaseReaction->IR32AllIndex.size()/6);
    gasPhaseReaction->IR32DupSize = static_cast<unsigned int>(gasPhaseReaction->IR32DupAllIndex.size()/6);
    gasPhaseReaction->NER32Size = static_cast<unsigned int>(gasPhaseReaction->NER32AllIndex.size()/6);
    gasPhaseReaction->NER32DupSize = static_cast<unsigned int>(gasPhaseReaction->NER32DupAllIndex.size()/6);

    gasPhaseReaction->RR31Size = static_cast<unsigned int>(gasPhaseReaction->RR31AllIndex.size()/5);
    gasPhaseReaction->RR31DupSize = static_cast<unsigned int>(gasPhaseReaction->RR31DupAllIndex.size()/5);
    gasPhaseReaction->IR31Size = static_cast<unsigned int>(gasPhaseReaction->IR31AllIndex.size()/5);
    gasPhaseReaction->IR31DupSize = static_cast<unsigned int>(gasPhaseReaction->IR31DupAllIndex.size()/5);
    gasPhaseReaction->NER31Size = static_cast<unsigned int>(gasPhaseReaction->NER31AllIndex.size()/5);
    gasPhaseReaction->NER31DupSize = static_cast<unsigned int>(gasPhaseReaction->NER31DupAllIndex.size()/5);

    gasPhaseReaction->RR33Size = static_cast<unsigned int>(gasPhaseReaction->RR33AllIndex.size()/7);
    gasPhaseReaction->RR33DupSize = static_cast<unsigned int>(gasPhaseReaction->RR33DupAllIndex.size()/7);
    gasPhaseReaction->IR33Size = static_cast<unsigned int>(gasPhaseReaction->IR33AllIndex.size()/7);
    gasPhaseReaction->IR33DupSize = static_cast<unsigned int>(gasPhaseReaction->IR33DupAllIndex.size()/7);
    gasPhaseReaction->NER33Size = static_cast<unsigned int>(gasPhaseReaction->NER33AllIndex.size()/7);
    gasPhaseReaction->NER33DupSize = static_cast<unsigned int>(gasPhaseReaction->NER33DupAllIndex.size()/7);



    /*********************************************************************************************/

}


}   // End namespace FastChemistry