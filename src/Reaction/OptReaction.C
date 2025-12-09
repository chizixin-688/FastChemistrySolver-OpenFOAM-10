#include "OptReaction.H"
#include <unordered_set>
#include <numeric>
#include <queue>
bool hasDuplicateFast(const std::vector<unsigned int>& lhs, const std::vector<unsigned int>& rhs)
{

    std::unordered_set<unsigned int> seen;
    seen.reserve(lhs.size() + rhs.size());

    for (auto x : lhs)
    {
        if (!seen.insert(x).second) 
        {
            return true;
        }
    }

    for (auto x : rhs)
    {
        if (!seen.insert(x).second) 
        {
            return true;
        }
    }

    return false;
}
void OptReaction::readReactionInfo
(
    std::vector<unsigned int>& inputLhsIndex,
    std::vector<unsigned int>& inputRhsIndex,
    const dictionary& nthreaction,
    const hashedWordList& speciesTable
)
{
    inputLhsIndex.clear();
    inputRhsIndex.clear();

    string reactionName = nthreaction.lookup("reaction");
    this->reactionTable_.push_back(reactionName);    
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
void OptReaction::readReactionInfo
(
    std::vector<unsigned int>& inputLhsIndex,
    std::vector<double>& inputLhsstoichCoeff,
    std::vector<double>& inputLhsReactionOrder,
    std::vector<unsigned int>& inputRhsIndex,
    std::vector<double>& inputRhsstoichCoeff,
    std::vector<double>& inputRhsReactionOrder,
    const dictionary& nthreaction,
    const hashedWordList& speciesTable
)
{
    inputLhsIndex.clear();
    inputLhsstoichCoeff.clear();
    inputLhsReactionOrder.clear();
    inputRhsIndex.clear();
    inputRhsstoichCoeff.clear();
    inputRhsReactionOrder.clear();

    string reactionName = nthreaction.lookup("reaction");
    this->reactionTable_.push_back(reactionName);    
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
bool OptReaction::checkInteger
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


/*bool OptReaction::findTwoTwoReaction
(
    const dictionary& nthreaction,
    const hashedWordList& speciesTable
)
{
    string reactionName = nthreaction.lookup("reaction");
    std::string stdReactionName(reactionName);
    std::istringstream iss(stdReactionName);
    List<word> words;
    List<word> ReactantStr;
    List<word> ProductStr;
    std::string Word;
    while (iss >> Word) 
    {
        words.append(Word);
    }
    int index = 0;
    for (int i = 0; i < words.size();i++)
    {
        if(words[i]=="=")
        {
            index =i;
        }
    }
    for (int i = 0; i < index;i++)
    {
        if(words[i]!="+")
        {
            ReactantStr.append(words[i]);
        }
    }
    for (int i = index+1; i < words.size();i++)
    {
        if(words[i]!="+")
        {
            ProductStr.append(words[i]);
        }
    }

    std::vector<unsigned int> tmplhsIndex;
    std::vector<unsigned int> tmprhsIndex;

    for(int  i = 0; i < ReactantStr.size();i++)
    {
        int first=0;

        for(unsigned int  j = 0; j < ReactantStr[i].size();j++)
        {
            if(!std::isdigit(ReactantStr[i][j]) && ReactantStr[i][j]!='.')
            {
                first = j;
                break;
            }
        }
        std::string coeffStr = ReactantStr[i].substr(0, first);
        std::string speciesStr = ReactantStr[i].substr(first);    

        unsigned int sl = 0;
        if(first==0)
        {
            sl=1;
        }
        else
        {
            double val = std::stod(coeffStr);
            double intpart;
            if (std::modf(val, &intpart) != 0.0) 
            {
                throw std::runtime_error("Stoichiometric coefficient must be integer-valued: " + coeffStr);
            }

            if (val < 0 || val > static_cast<double>(std::numeric_limits<unsigned int>::max())) 
            {
                throw std::runtime_error("Stoichiometric coefficient out of range: " + coeffStr);
            }
            sl = static_cast<unsigned int>(val);
        }
        const int newSpecIndex = speciesTable[ReactantStr[i].substr(first)];

        while(sl!=0)
        {
            tmplhsIndex.push_back(newSpecIndex);                    
            sl--;
        }
    }
    for(int  i = 0; i < ProductStr.size();i++)
    {

        int first=0;
        for(unsigned int  j = 0; j < ProductStr[i].size();j++)
        {
            if(!std::isdigit(ProductStr[i][j])&& ProductStr[i][j]!='.')
            {
                first = j;
                break;
            }
        }

        std::string coeffStr = ProductStr[i].substr(0, first);
        std::string speciesStr = ProductStr[i].substr(first);    

        unsigned int sr = 0;
        if(first==0)
        {
            sr=1;
        }
        else
        {
            double val = std::stod(coeffStr);
            double intpart;
            if (std::modf(val, &intpart) != 0.0) 
            {
                throw std::runtime_error("Stoichiometric coefficient must be integer-valued: " + coeffStr);
            }

            if (val < 0 || val > static_cast<double>(std::numeric_limits<unsigned int>::max())) 
            {
                throw std::runtime_error("Stoichiometric coefficient out of range: " + coeffStr);
            }
            sr = static_cast<unsigned int>(val);
        }

        const int newSpecIndex = speciesTable[ProductStr[i].substr(first)];

        while(sr!=0)
        {
            tmprhsIndex.push_back(newSpecIndex);                   
            sr--;
        }
    }

    if (tmplhsIndex.size()==2&&tmprhsIndex.size()==2)
    {
        return true;
    }
    else
    {
        return false;
    }
}*/



OptReaction::OptReaction
(
    const dictionary& chemistryDict,
    const dictionary& physicalDict
)
{
    this->readInfo(chemistryDict,physicalDict);
}

OptReaction::OptReaction
()
{}

void OptReaction::readInfo
(
    const dictionary& chemistryDict,
    const dictionary& physicalDict    
)
{

    hashedWordList speciesTable(physicalDict.lookup("species"));

    this->speciesTable_.resize(speciesTable.size());
    for(unsigned int i=0; i<this->speciesTable_.size();i++)
    {
        this->speciesTable_[i] = speciesTable[i];
    }

    const dictionary& reactions(chemistryDict.subDict("reactions"));

    unsigned int nLindemann = 0;
    unsigned int nTroe = 0;
    unsigned int nSRI = 0;   


    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& nthreaction = reactions.subDict(key);
        const word reactionTypeName = nthreaction.lookup("type");
        Foam::string reactionName = nthreaction.lookup("reaction");

        this->n_Reactions++;

        if(reactionTypeName == "irreversibleArrhenius")
        {
            this->n_Arrhenius++;
        }
        else if(reactionTypeName == "reversibleArrhenius")
        {
            this->n_Arrhenius++;
        }
        else if(reactionTypeName == "nonEquilibriumReversibleArrhenius")
        {
            this->n_NonEquilibriumReversibleArrhenius++;
        }
        else if(reactionTypeName == "nonEquilibriumReversibleThirdBodyArrhenius")
        {
            this->n_NonEquilibriumThirdBodyReaction++;
        }
        else if
        (
            reactionTypeName == "reversibleThirdBodyArrhenius"||
            reactionTypeName == "irreversibleThirdBodyArrhenius"
        )
        {
            this->n_ThirdBodyReaction++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusLindemannFallOff"||
            reactionTypeName == "irreversibleArrheniusLindemannFallOff"
        )
        {
            this->n_Fall_Off_Reaction++;nLindemann++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusTroeFallOff"||
            reactionTypeName == "irreversibleArrheniusTroeFallOff"
        )
        {
            this->n_Fall_Off_Reaction++;nTroe++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusSRIFallOff"||
            reactionTypeName == "irreversibleArrheniusSRIFallOff"
        )
        {
            this->n_Fall_Off_Reaction++;nSRI++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusLindemannChemicallyActivated"||
            reactionTypeName == "irreversibleArrheniusLindemannChemicallyActivated"
        )
        {
            this->n_ChemicallyActivated_Reaction++;nLindemann++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusTroeChemicallyActivated"||
            reactionTypeName == "irreversibleArrheniusTroeChemicallyActivated"
        )
        {
            this->n_ChemicallyActivated_Reaction++;nTroe++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusSRIChemicallyActivated"||
            reactionTypeName == "irreversibleArrheniusSRIChemicallyActivated"
        )
        {
            this->n_ChemicallyActivated_Reaction++;nSRI++;
        }
        else if
        (
            reactionTypeName == "reversibleArrheniusPLOG"||
            reactionTypeName == "irreversibleArrheniusPLOG"
        )
        {
            this->n_PlogReaction++;
        }
        else
        {
            FatalErrorInFunction<< "unknown reaction type:"
                << reactionTypeName << exit(FatalError);
        }
    }

    {
        this->Itbr[0] = 0;
        this->Itbr[1] = this->n_NonEquilibriumThirdBodyReaction;
        this->Itbr[2] = this->Itbr[1] + this->n_ThirdBodyReaction;
        this->Itbr[3] = this->Itbr[2] + this->n_Fall_Off_Reaction;   
        this->Itbr[4] = this->Itbr[3] + this->n_ChemicallyActivated_Reaction;       
        this->Itbr[5] = this->Itbr[4] + this->n_NonEquilibriumThirdBodyReaction; 
    }

    {
        Ikf[0] = 0;
        Ikf[1] = this->n_Arrhenius;
        Ikf[2] = Ikf[1] + this->n_NonEquilibriumReversibleArrhenius;
        Ikf[3] = Ikf[2] + this->n_NonEquilibriumThirdBodyReaction;   
        Ikf[4] = Ikf[3] + this->n_ThirdBodyReaction;       
        Ikf[5] = Ikf[4] + this->n_Fall_Off_Reaction; 
        Ikf[6] = Ikf[5] + this->n_ChemicallyActivated_Reaction; 
        Ikf[7] = Ikf[6] + this->n_PlogReaction;
        Ikf[8] = Ikf[7] + this->n_Fall_Off_Reaction;   
        Ikf[9] = Ikf[8] + this->n_ChemicallyActivated_Reaction;   
        Ikf[10] = Ikf[9] + this->n_NonEquilibriumReversibleArrhenius;        
        Ikf[11] = Ikf[10] + this->n_NonEquilibriumThirdBodyReaction;
        Ikf[12] = Ikf[11] + this->n_PlogReaction;

    }
    this->offset_kinf = - Ikf[4] + Ikf[7];

    //this->n_Reactions                        = n_Reactions;
    this->nSpecies = speciesTable.size();

    this->A.resize(Ikf[12]);
    this->beta.resize(Ikf[12]);
    this->Ta.resize(Ikf[12]);
    this->lhsSpeciesIndex.resize(this->n_Reactions);
    this->rhsSpeciesIndex.resize(this->n_Reactions);
    this->lhsStoichCoeff.resize(this->n_Reactions);    
    this->rhsStoichCoeff.resize(this->n_Reactions);
    this->lhsReactionOrder.resize(this->n_Reactions);
    this->rhsReactionOrder.resize(this->n_Reactions);
    std::vector<std::vector<double>> ThirdBodyFactor(this->Itbr[5]);
    //ThirdBodyFactor.resize(Itbr[5]);
    this->alpha_.resize(0);
    this->alpha_.reserve(nTroe);
    this->Ts_.resize(0);
    this->Ts_.reserve(nTroe);    
    this->Tss_.resize(0);
    this->Tss_.reserve(nTroe);    
    this->Tsss_.resize(0);
    this->Tsss_.reserve(nTroe);
    this->a_.resize(0);
    this->b_.resize(0);
    this->c_.resize(0);
    this->d_.resize(0);
    this->e_.resize(0);
    this->a_.reserve(nSRI);
    this->b_.reserve(nSRI);
    this->c_.reserve(nSRI);
    this->d_.reserve(nSRI);
    this->e_.reserve(nSRI);
    this->HCoeffs.resize(this->nSpecies);
    this->LCoeffs.resize(this->nSpecies);
    this->Tlow.resize(this->nSpecies);
    this->Thigh.resize(this->nSpecies);
    this->Tcommon.resize(this->nSpecies);
    this->TcommonMin=0;
    this->TcommonMax=1e10;
    this->PtrCoeffs.resize(this->nSpecies);

    //this->W.resize(this->nSpecies);

    unsigned int alignSpecies = ((this->nSpecies+3)/4)*4;
    if (posix_memalign(reinterpret_cast<void**>(&this->W), 32, alignSpecies * sizeof(double)))
    {
        throw std::bad_alloc();
    }
    std::memset(this->W, 0, alignSpecies * sizeof(double));

    //this->invW.resize(this->nSpecies);
    if (posix_memalign(reinterpret_cast<void**>(&this->invW), 32, alignSpecies * sizeof(double)))
    {
        throw std::bad_alloc();
    }
    std::memset(this->invW, 0, alignSpecies * sizeof(double));

    if (posix_memalign(reinterpret_cast<void**>(&this->negGstdByRT), 32, alignSpecies * sizeof(double)))
    {
        throw std::bad_alloc();
    }
    std::memset(this->negGstdByRT, 0, alignSpecies * sizeof(double));

    if (posix_memalign(reinterpret_cast<void**>(&this->Hf), 32, alignSpecies * sizeof(double)))
    {
        throw std::bad_alloc();
    }
    std::memset(this->Hf, 0, alignSpecies * sizeof(double));

    this->isIrreversible.resize(this->n_Reactions,0); 
    this->isGlobal.resize(this->n_Reactions,0);

    scalar TcommonMax_ = 0;
    scalar TcommonMin_ = 1e10;
    for(int i = 0; i < speciesTable.size();i++)
    {
        const dictionary specieDict(physicalDict.subDict(speciesTable[i]));
        const dictionary thermodynamicsDict(specieDict.subDict("thermodynamics"));
        this->Tcommon[i] = thermodynamicsDict.lookup<scalar>("Tcommon");
        this->Tlow[i] = thermodynamicsDict.lookup<scalar>("Tlow");
        this->Thigh[i] = thermodynamicsDict.lookup<scalar>("Thigh");
        FixedList<scalar,7> temp1(thermodynamicsDict.lookup("highCpCoeffs"));

        const dictionary species(specieDict.subDict("specie"));
        this->W[i] = species.lookup<scalar>("molWeight");
        this->invW[i] = 1.0/this->W[i];
        for(unsigned int j = 0; j < 7; j ++)
        {
            this->HCoeffs[i][j] = temp1[j];
        }
        FixedList<scalar,7> temp2(thermodynamicsDict.lookup("lowCpCoeffs")); 
        for(unsigned int j = 0; j < 7; j ++)
        {
            this->LCoeffs[i][j] = temp2[j];
        }               
        TcommonMax_ = (this->Tcommon[i]>TcommonMax_)?this->Tcommon[i]:TcommonMax_;
        TcommonMin_ = (this->Tcommon[i]<TcommonMin_)?this->Tcommon[i]:TcommonMin_;
    }
    this->TcommonMin = TcommonMin_;
    this->TcommonMax = TcommonMax_;    

    // Find temperature independent Arrhenius reaction
    int iArrhenius = 0;
    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");

        bool isInteger = this->checkInteger(reactDict);

        
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
                {this->isIrreversible[iArrhenius]=1;}
                this->reactionType_.push_back(reactionTypeName);
                this->reactionName_.push_back(key);
                this->A[iArrhenius] = reactDict.lookup<scalar>("A");
                this->beta[iArrhenius] = reactDict.lookup<scalar>("beta");
                this->Ta[iArrhenius] = reactDict.lookup<scalar>("Ta");

                if(isInteger==true)
                {
                    this->readReactionInfo
                    (
                        this->lhsSpeciesIndex[iArrhenius],
                        this->rhsSpeciesIndex[iArrhenius],
                        reactDict,
                        speciesTable
                    );
                }
                else
                {
                    this->readReactionInfo
                    (
                        this->lhsSpeciesIndex[iArrhenius],
                        this->lhsStoichCoeff[iArrhenius],
                        this->lhsReactionOrder[iArrhenius],
                        this->rhsSpeciesIndex[iArrhenius],
                        this->rhsStoichCoeff[iArrhenius],
                        this->rhsReactionOrder[iArrhenius],
                        reactDict,
                        speciesTable
                    );
                    this->isGlobal[iArrhenius] = 1;
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
        bool isInteger = this->checkInteger(reactDict);

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
                {this->isIrreversible[iArrhenius]=1;}
                this->reactionType_.push_back(reactionTypeName);
                this->reactionName_.push_back(key);
                this->A[iArrhenius] = reactDict.lookup<scalar>("A");
                this->beta[iArrhenius] = reactDict.lookup<scalar>("beta");
                this->Ta[iArrhenius] = reactDict.lookup<scalar>("Ta");

                if(isInteger==true)
                {
                    this->readReactionInfo
                    (
                        this->lhsSpeciesIndex[iArrhenius],
                        this->rhsSpeciesIndex[iArrhenius],
                        reactDict,
                        speciesTable
                    );
                }
                else
                {
                    this->readReactionInfo
                    (
                        this->lhsSpeciesIndex[iArrhenius],
                        this->lhsStoichCoeff[iArrhenius],
                        this->lhsReactionOrder[iArrhenius],
                        this->rhsSpeciesIndex[iArrhenius],
                        this->rhsStoichCoeff[iArrhenius],
                        this->rhsReactionOrder[iArrhenius],
                        reactDict,
                        speciesTable
                    );
                    this->isGlobal[iArrhenius] = 1;                    
                }
                iArrhenius++;
            }
        }
    }


    auto j = this->Ikf[9];
    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = this->checkInteger(reactDict);

        if(reactionTypeName=="nonEquilibriumReversibleArrhenius")
        {
            this->isIrreversible[iArrhenius]=2;
            const dictionary& forwardDict = reactDict.subDict("forward");
            const dictionary& reverseDict = reactDict.subDict("reverse");

            this->reactionType_.push_back(reactionTypeName);            
            this->reactionName_.push_back(key);
            this->A[iArrhenius] = forwardDict.lookup<scalar>("A");       
            this->beta[iArrhenius] = forwardDict.lookup<scalar>("beta");
            this->Ta[iArrhenius] = forwardDict.lookup<scalar>("Ta");

            this->A[j] = reverseDict.lookup<scalar>("A");       
            this->beta[j] = reverseDict.lookup<scalar>("beta");
            this->Ta[j] = reverseDict.lookup<scalar>("Ta");
            if(isInteger==true)
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable
                );
            }
            else
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->lhsStoichCoeff[iArrhenius],
                    this->lhsReactionOrder[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    this->rhsStoichCoeff[iArrhenius],
                    this->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable
                );                
                    this->isGlobal[iArrhenius] = 1;                
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
        bool isInteger = this->checkInteger(reactDict);

        if(reactionTypeName=="nonEquilibriumReversibleThirdBodyArrhenius")
        {
            this->isIrreversible[iArrhenius]=2;
            const dictionary& forwardDict = reactDict.subDict("forward");
            const dictionary& reverseDict = reactDict.subDict("reverse");

            this->reactionType_.push_back(reactionTypeName); 
            this->reactionName_.push_back(key);

            this->A[iArrhenius] = forwardDict.lookup<scalar>("A");     
            this->beta[iArrhenius] = forwardDict.lookup<scalar>("beta");
            this->Ta[iArrhenius] = forwardDict.lookup<scalar>("Ta");

            this->A[j] = reverseDict.lookup<scalar>("A");        
            this->beta[j] = reverseDict.lookup<scalar>("beta");
            this->Ta[j] = reverseDict.lookup<scalar>("Ta");     

            List<Tuple2<word, scalar>> forwardCoeffs(forwardDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(forwardCoeffs.size());
            forAll(forwardCoeffs, n)
            {
                const int l = speciesTable[(forwardCoeffs[n].first())];
                const scalar ThirdBodyFactor_n = forwardCoeffs[n].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_n;
            }
            
            List<Tuple2<word, scalar>> reverseCoeffs(reverseDict.lookup("coeffs"));

            auto begin = k + this->Itbr[4];
            

            ThirdBodyFactor[begin].resize(reverseCoeffs.size());
            forAll(reverseCoeffs, n)
            {
                const int l = speciesTable[(reverseCoeffs[n].first())];
                const scalar ThirdBodyFactor_n = reverseCoeffs[n].second();
                ThirdBodyFactor[begin][l] = ThirdBodyFactor_n;
            }
            if(isInteger==true)
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable
                );
            }
            else
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->lhsStoichCoeff[iArrhenius],
                    this->lhsReactionOrder[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    this->rhsStoichCoeff[iArrhenius],
                    this->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable
                ); 
                    this->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;j++;k++;
        }
    }

    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = this->checkInteger(reactDict);

        if
        (
            reactionTypeName=="reversibleThirdBodyArrhenius"||
            reactionTypeName=="irreversibleThirdBodyArrhenius"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {this->isIrreversible[iArrhenius]=1;}
            this->reactionType_.push_back(reactionTypeName); 
            this->reactionName_.push_back(key);

            this->A[iArrhenius] = reactDict.lookup<scalar>("A");        
            this->beta[iArrhenius] = reactDict.lookup<scalar>("beta");
            this->Ta[iArrhenius] = reactDict.lookup<scalar>("Ta");

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
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable
                );
            }
            else
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->lhsStoichCoeff[iArrhenius],
                    this->lhsReactionOrder[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    this->rhsStoichCoeff[iArrhenius],
                    this->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable
                );
                    this->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;k++;
        }
    }

    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = this->checkInteger(reactDict);

        if
        (
            reactionTypeName=="reversibleArrheniusLindemannFallOff"||
            reactionTypeName=="irreversibleArrheniusLindemannFallOff"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {this->isIrreversible[iArrhenius]=1;}
            this->reactionType_.push_back(reactionTypeName);
            this->reactionName_.push_back(key);
            const dictionary& k0Dict = reactDict.subDict("k0");
            const dictionary& kInfDict = reactDict.subDict("kInf");

            const dictionary& thirdBodyEfficienciesDict =
            reactDict.subDict("thirdBodyEfficiencies");
    
            this->A[iArrhenius] = k0Dict.lookup<scalar>("A");        
            this->beta[iArrhenius] = k0Dict.lookup<scalar>("beta");
            this->Ta[iArrhenius] = k0Dict.lookup<scalar>("Ta");
    
            auto begin = iArrhenius - Ikf[4] + Ikf[7];
            this->A[begin] = (kInfDict.lookup<scalar>("A")) ;         
            this->beta[begin] = (kInfDict.lookup<scalar>("beta")) ;
            this->Ta[begin] = (kInfDict.lookup<scalar>("Ta")) ;            
    
            List<Tuple2<word, scalar>> coeffs(thirdBodyEfficienciesDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
    
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }
            this->Lindemann.push_back(iArrhenius);
            if(isInteger==true)
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable
                );
            }
            else
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->lhsStoichCoeff[iArrhenius],
                    this->lhsReactionOrder[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    this->rhsStoichCoeff[iArrhenius],
                    this->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable
                );
                    this->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;k++;
       }
    }
    
    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = this->checkInteger(reactDict);
        if
        (
            reactionTypeName=="reversibleArrheniusTroeFallOff"||
            reactionTypeName=="irreversibleArrheniusTroeFallOff"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {this->isIrreversible[iArrhenius]=1;}
            this->reactionType_.push_back(reactionTypeName);
            this->reactionName_.push_back(key);
            const dictionary& k0Dict = reactDict.subDict("k0");
            const dictionary& kInfDict = reactDict.subDict("kInf");
            const dictionary& FDict = reactDict.subDict("F");
            const dictionary& thirdBodyEfficienciesDict =
            reactDict.subDict("thirdBodyEfficiencies");
    
            this->A[iArrhenius] = k0Dict.lookup<scalar>("A");        
            this->beta[iArrhenius] = k0Dict.lookup<scalar>("beta");
            this->Ta[iArrhenius] = k0Dict.lookup<scalar>("Ta");
    
            auto begin = iArrhenius - Ikf[4] + Ikf[7];
            this->A[begin] = (kInfDict.lookup<scalar>("A")) ;         
            this->beta[begin] = (kInfDict.lookup<scalar>("beta")) ;
            this->Ta[begin] = (kInfDict.lookup<scalar>("Ta")) ;            
    
            List<Tuple2<word, scalar>> coeffs(thirdBodyEfficienciesDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
    
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }

            this->Troe.push_back(iArrhenius);
            this->alpha_.push_back(FDict.lookup<scalar>("alpha"));    
            this->Ts_.push_back(FDict.lookup<scalar>("Ts"));    
            this->Tss_.push_back(FDict.lookup<scalar>("Tss"));    
            this->Tsss_.push_back(FDict.lookup<scalar>("Tsss"));
            if(isInteger==true)
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable
                );
            }
            else
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->lhsStoichCoeff[iArrhenius],
                    this->lhsReactionOrder[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    this->rhsStoichCoeff[iArrhenius],
                    this->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable
                );
                this->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;k++;
       }
    }

    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = this->checkInteger(reactDict); 
        if
        (
            reactionTypeName=="reversibleArrheniusSRIFallOff"||
            reactionTypeName=="irreversibleArrheniusSRIFallOff"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {this->isIrreversible[iArrhenius]=1;}
            this->reactionType_.push_back(reactionTypeName);
            this->reactionName_.push_back(key);
            const dictionary& k0Dict = reactDict.subDict("k0");
            const dictionary& kInfDict = reactDict.subDict("kInf");
            const dictionary& FDict = reactDict.subDict("F");
            const dictionary& thirdBodyEfficienciesDict =
            reactDict.subDict("thirdBodyEfficiencies");
    
            this->A[iArrhenius] = k0Dict.lookup<scalar>("A");        
            this->beta[iArrhenius] = k0Dict.lookup<scalar>("beta");
            this->Ta[iArrhenius] = k0Dict.lookup<scalar>("Ta");
    
            auto begin = iArrhenius - Ikf[4] + Ikf[7];
            this->A[begin] = (kInfDict.lookup<scalar>("A")) ;         
            this->beta[begin] = (kInfDict.lookup<scalar>("beta")) ;
            this->Ta[begin] = (kInfDict.lookup<scalar>("Ta")) ;            
    
            List<Tuple2<word, scalar>> coeffs(thirdBodyEfficienciesDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
    
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }
            
            this->SRI.push_back(iArrhenius);
            this->a_.push_back(FDict.lookup<scalar>("a"));    
            this->b_.push_back(FDict.lookup<scalar>("b"));    
            this->c_.push_back(FDict.lookup<scalar>("c"));    
            this->d_.push_back(FDict.lookup<scalar>("d"));  
            this->e_.push_back(FDict.lookup<scalar>("e"));  
            
            if(isInteger==true)
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable
                );
            }
            else
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->lhsStoichCoeff[iArrhenius],
                    this->lhsReactionOrder[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    this->rhsStoichCoeff[iArrhenius],
                    this->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable
                );
                this->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;k++;
       }
    }

    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = this->checkInteger(reactDict);

        if(
            reactionTypeName=="reversibleArrheniusLindemannChemicallyActivated"||
            reactionTypeName=="irreversibleArrheniusLindemannChemicallyActivated"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {this->isIrreversible[iArrhenius]=1;}

            this->reactionType_.push_back(reactionTypeName); 
            this->reactionName_.push_back(key);   

            const dictionary& k0Dict = reactions.subDict("k0");
            const dictionary& kInfDict = reactions.subDict("kInf");

            const dictionary& thirdBodyEfficienciesDict = 
            reactions.subDict("thirdBodyEfficiencies");

            this->A[iArrhenius] = k0Dict.lookup<scalar>("A");
            this->beta[iArrhenius] = k0Dict.lookup<scalar>("beta");
            this->Ta[iArrhenius] = k0Dict.lookup<scalar>("Ta");

            auto begin = iArrhenius - this->Ikf[5] + this->Ikf[8];
            this->A[begin] = kInfDict.lookup<scalar>("A");
            this->beta[begin] = kInfDict.lookup<scalar>("beta");
            this->Ta[begin] = kInfDict.lookup<scalar>("Ta");

            List<Tuple2<word, scalar>> coeffs(thirdBodyEfficienciesDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }
            this->Lindemann.push_back(iArrhenius);
            if(isInteger==true)
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable
                );
            }
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->lhsStoichCoeff[iArrhenius],
                    this->lhsReactionOrder[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    this->rhsStoichCoeff[iArrhenius],
                    this->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable
                );
                this->isGlobal[iArrhenius] = 1;                
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
        bool isInteger = this->checkInteger(reactDict);

        if(
            reactionTypeName=="reversibleArrheniusTroeChemicallyActivated"||
            reactionTypeName=="irreversibleArrheniusTroeChemicallyActivated"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {this->isIrreversible[iArrhenius]=1;}

            this->reactionType_.push_back(reactionTypeName); 
            this->reactionName_.push_back(key);   

            const dictionary& k0Dict = reactions.subDict("k0");
            const dictionary& kInfDict = reactions.subDict("kInf");
            const dictionary& FDict = reactions.subDict("F");
            const dictionary& thirdBodyEfficienciesDict = 
            reactions.subDict("thirdBodyEfficiencies");

            this->A[iArrhenius] = k0Dict.lookup<scalar>("A");
            this->beta[iArrhenius] = k0Dict.lookup<scalar>("beta");
            this->Ta[iArrhenius] = k0Dict.lookup<scalar>("Ta");

            auto begin = iArrhenius - Ikf[5] + Ikf[8];
            this->A[begin] = kInfDict.lookup<scalar>("A");
            this->beta[begin] = kInfDict.lookup<scalar>("beta");
            this->Ta[begin] = kInfDict.lookup<scalar>("Ta");

            List<Tuple2<word, scalar>> coeffs(thirdBodyEfficienciesDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }

            this->Troe.push_back(iArrhenius);
            this->alpha_.push_back(FDict.lookup<scalar>("alpha"));    
            this->Ts_.push_back(FDict.lookup<scalar>("Ts"));    
            this->Tss_.push_back(FDict.lookup<scalar>("Tss"));    
            this->Tsss_.push_back(FDict.lookup<scalar>("Tsss"));                
            if(isInteger==true)
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable
                );
            }
            else
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->lhsStoichCoeff[iArrhenius],
                    this->lhsReactionOrder[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    this->rhsStoichCoeff[iArrhenius],
                    this->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable
                );
                this->isGlobal[iArrhenius] = 1;                
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
        bool isInteger = this->checkInteger(reactDict);

        if(
            reactionTypeName=="reversibleArrheniusSRIChemicallyActivated" ||
            reactionTypeName=="irreversibleArrheniusSRIChemicallyActivated" 
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {this->isIrreversible[iArrhenius]=1;}

            this->reactionType_.push_back(reactionTypeName); 
            this->reactionName_.push_back(key);   

            const dictionary& k0Dict = reactions.subDict("k0");
            const dictionary& kInfDict = reactions.subDict("kInf");
            const dictionary& FDict = reactions.subDict("F");
            const dictionary& thirdBodyEfficienciesDict = 
            reactions.subDict("thirdBodyEfficiencies");

            this->A[iArrhenius] = k0Dict.lookup<scalar>("A");
            this->beta[iArrhenius] = k0Dict.lookup<scalar>("beta");
            this->Ta[iArrhenius] = k0Dict.lookup<scalar>("Ta");

            auto begin = iArrhenius - Ikf[5] + Ikf[8];
            this->A[begin] = kInfDict.lookup<scalar>("A");
            this->beta[begin] = kInfDict.lookup<scalar>("beta");
            this->Ta[begin] = kInfDict.lookup<scalar>("Ta");

            List<Tuple2<word, scalar>> coeffs(thirdBodyEfficienciesDict.lookup("coeffs"));
            ThirdBodyFactor[k].resize(coeffs.size());
            forAll(coeffs, m)
            {
                const int l = speciesTable[(coeffs[m].first())];
                const scalar ThirdBodyFactor_m = coeffs[m].second();
                ThirdBodyFactor[k][l] = ThirdBodyFactor_m;
            }
            this->SRI.push_back(iArrhenius);
            this->a_.push_back(FDict.lookup<scalar>("a"));    
            this->b_.push_back(FDict.lookup<scalar>("b"));    
            this->c_.push_back(FDict.lookup<scalar>("c"));    
            this->d_.push_back(FDict.lookup<scalar>("d"));  
            this->e_.push_back(FDict.lookup<scalar>("e"));
            if(isInteger==true)
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable
                );
            }
            else
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->lhsStoichCoeff[iArrhenius],
                    this->lhsReactionOrder[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    this->rhsStoichCoeff[iArrhenius],
                    this->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable
                );
                this->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;
            k++;
        }
    }

   {
        unsigned int remain = 4 - this->nSpecies%4;

        this->AlignSpecies = this->nSpecies+remain;

        if (
            posix_memalign
            (
                reinterpret_cast<void**>(&this->ThirdBodyFactor1D), 
                32, 
                ThirdBodyFactor.size()*this->AlignSpecies*sizeof(double)
            )
            )
        {
            throw std::bad_alloc();
        }
        std::memset(this->ThirdBodyFactor1D, 0, ThirdBodyFactor.size()*this->AlignSpecies*sizeof(double));


        unsigned int count = 0;
        for(unsigned int i = 0; i < ThirdBodyFactor.size();i++)
        {
            for(unsigned int J = 0; J < ThirdBodyFactor[i].size();J++)
            {
                ThirdBodyFactor1D[count] = ThirdBodyFactor[i][J];
                count++;
            }
            for(unsigned int J = 0; J < remain;J++)
            {
                ThirdBodyFactor1D[count] = 0;
                count++;                
            }
        }
   }



    APlog.resize(this->n_PlogReaction);
    logAPlog.resize(this->n_PlogReaction);
    betaPlog.resize(this->n_PlogReaction);
    TaPlog.resize(this->n_PlogReaction);
    Prange.resize(this->n_PlogReaction);
    rDeltaP_.resize(this->n_PlogReaction);
    logPi.resize(this->n_PlogReaction);
    Pindex.resize(this->n_PlogReaction);

    unsigned int a = 0;
    forAllConstIter(dictionary, reactions, iter)
    {
        const word& key = iter().keyword();
        const dictionary& reactDict = reactions.subDict(key);
        const word reactionTypeName = reactDict.lookup("type");
        bool isInteger = this->checkInteger(reactDict);
        if
        (
            reactionTypeName=="reversibleArrheniusPLOG"||
            reactionTypeName=="irreversibleArrheniusPLOG"
        )
        {
            if(reactionTypeName.find("irreversible",0)!=std::string::npos)
            {
                this->isIrreversible[iArrhenius]=1;
            }
            this->reactionType_.push_back(reactionTypeName); 
            this->reactionName_.push_back(key);   

            List<List<double>> PlogData(reactDict.lookup("ArrheniusData"));
            unsigned int pSize = PlogData.size();

            APlog[a].resize(pSize);
            logAPlog[a].resize(pSize); 
            betaPlog[a].resize(pSize);
            TaPlog[a].resize(pSize);
            Prange[a].resize(pSize);
            rDeltaP_[a].resize(pSize-1);
            logPi[a].resize(pSize);

            for(unsigned int i = 0; i < pSize; i ++)
            {
                Prange[a][i] = PlogData[i][0];  
                APlog[a][i] = PlogData[i][1];
                logAPlog[a][i] = std::log(APlog[a][i]);
                betaPlog[a][i] = PlogData[i][2];
                TaPlog[a][i] = PlogData[i][3];
                logPi[a][i] = std::log(Prange[a][i]);
            }

            for(unsigned int i = 0; i < pSize-1; i ++)
            {
                rDeltaP_[a][i] = 1.0/(logPi[a][i+1]-logPi[a][i]);
            }
            if(isInteger==true)
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    reactDict,
                    speciesTable
                );
            }
            else
            {
                this->readReactionInfo
                (
                    this->lhsSpeciesIndex[iArrhenius],
                    this->lhsStoichCoeff[iArrhenius],
                    this->lhsReactionOrder[iArrhenius],
                    this->rhsSpeciesIndex[iArrhenius],
                    this->rhsStoichCoeff[iArrhenius],
                    this->rhsReactionOrder[iArrhenius],
                    reactDict,
                    speciesTable
                );
                this->isGlobal[iArrhenius] = 1;                
            }
            iArrhenius++;a++;
        }
    }

    List<int> sumVki(this->n_Reactions);

     for(unsigned int i = 0;i<this->n_Reactions;i++)
    {
        sumVki[i] = 0;
        for(unsigned int jj = 0; jj<this->rhsSpeciesIndex[i].size();jj++)
        {
            sumVki[i] = sumVki[i] + 1;
        }  
        for(unsigned int jj = 0; jj<this->lhsSpeciesIndex[i].size();jj++)
        {
            sumVki[i] = sumVki[i] - 1;
        }  
    } 

    this->Pow_pByRT_SumVki_I.insert({sumVki[0],0.0});
    for(int i = 1; i < sumVki.size();i++)
    {
        auto it = this->Pow_pByRT_SumVki_I.find(sumVki[i]);
        if(it==this->Pow_pByRT_SumVki_I.end())
        {this->Pow_pByRT_SumVki_I.insert({sumVki[i],0.0});}
    }

    this->Kf_.resize(this->Ikf[12]);
    this->dKfdT_.resize(this->Ikf[12]);
    this->dKfdC_.resize(this->Itbr[5]);
    this->tmp_M.resize(this->Itbr[5]);

    {
        tmp_ExpSize = (this->nSpecies + static_cast<unsigned int>(this->Troe.size())*3 + static_cast<unsigned int>(this->SRI.size())*2);
        size_t bytes = (this->nSpecies + this->Troe.size()*3 + this->SRI.size()*2)  * sizeof(double);
        if (posix_memalign(reinterpret_cast<void**>(&this->tmp_Exp), 32, bytes))
        {
            throw std::bad_alloc();
        }
        std::memset(this->tmp_Exp, 0, bytes);
    }



    this->invTs_.resize(this->Ts_.size());
    this->invTsss_.resize(this->Tsss_.size());
    for(unsigned int i = 0; i < this->Ts_.size();i++)
    {
        this->invTs_[i] = 1.0/this->Ts_[i];
    }
    for(unsigned int i = 0; i < this->Tsss_.size();i++)
    {
        this->invTsss_[i] = 1.0/this->Tsss_[i];
    }

    this->invc_.resize(this->c_.size());
    for(unsigned int i = 0; i < this->c_.size();i++)
    {
        this->invc_[i] = 1.0/this->c_[i];
    }    

    this->n_Temperature_Independent_Reaction =0;
    if(this->n_Arrhenius>0)
    {
        for(unsigned int ii = 0; ii < this->n_Arrhenius;ii++)
        {
            if(this->beta[ii]==0&&this->Ta[ii]==0)
            {this->n_Temperature_Independent_Reaction++;}
        }
    }

    for(unsigned int i0 = 0; i0 < this->n_Temperature_Independent_Reaction;i0++)
    {
        this->Kf_[i0]=this->A[i0];
        this->dKfdT_[i0] = 0;
    }

    this->n_ = this->nSpecies+1;
    
    unsigned int ArrSize = this->nSpecies+(4-this->nSpecies%4);
    size_t bytes = 8 * ArrSize * sizeof(double);
    if (posix_memalign(reinterpret_cast<void**>(&this->buffer), 32, bytes))
    {
        throw std::bad_alloc();
    }
    std::memset(this->buffer, 0, bytes);


    TlowMin=1e10;
    ThighMax=1;
    for(unsigned int i1 = 0; i1 < this->nSpecies;i1++)
    {
        if (TlowMin>Tlow[i1])
        {TlowMin=Tlow[i1];}
        if(ThighMax<Thigh[i1])
        {ThighMax=Thigh[i1];}
    }


    {
        unsigned int lhsAll=0;
        for(size_t i = 0; i < lhsSpeciesIndex.size();i++)
        {
            for(size_t J = 0; J < lhsSpeciesIndex[i].size();J++)
            {
                lhsAll++;
            }
        }
        lhsSpeciesIndex1D.resize(lhsAll);    
        lhsOffset.resize(lhsSpeciesIndex.size()+1);
        lhsAll=0;
        for(size_t i = 0; i < lhsSpeciesIndex.size();i++)
        {
            lhsOffset[i+1] = lhsOffset[i] + static_cast<unsigned int>(lhsSpeciesIndex[i].size());
            for(size_t J = 0; J < lhsSpeciesIndex[i].size();J++)
            {
                lhsSpeciesIndex1D[lhsAll] = lhsSpeciesIndex[i][J];
                lhsAll++;
            }
        }         
        lhsOffset[lhsSpeciesIndex.size()] = static_cast<unsigned int>(lhsSpeciesIndex1D.size());

        unsigned int rhsAll=0;
        for(size_t i = 0; i < rhsSpeciesIndex.size();i++)
        {
            for(size_t J = 0; J < rhsSpeciesIndex[i].size();J++)
            {
                rhsAll++;
            }
        }
        rhsSpeciesIndex1D.resize(rhsAll);    
        rhsOffset.resize(rhsSpeciesIndex.size()+1);
        rhsAll=0;
        for(size_t i = 0; i < rhsSpeciesIndex.size();i++)
        {
            rhsOffset[i+1] = rhsOffset[i] + static_cast<unsigned int>(rhsSpeciesIndex[i].size());
            for(size_t J = 0; J < rhsSpeciesIndex[i].size();J++)
            {
                rhsSpeciesIndex1D[rhsAll] = rhsSpeciesIndex[i][J];
                rhsAll++;
            }
        }       
        rhsOffset[rhsSpeciesIndex.size()] = static_cast<unsigned int>(rhsSpeciesIndex1D.size());
    }


    for(unsigned int i = 0; i < this->nSpecies; i++)
    {
        const double Tstd = 298.15;
        
        if(this->Tcommon[i]<Tstd)
        {
            auto& Coeff = this->HCoeffs[i];
            this->Hf[i] = (((((Coeff[4]*Tstd*0.2+Coeff[3]*0.25)*Tstd+Coeff[2]*(1.0/3.0))*Tstd+Coeff[1]*0.5)*Tstd+Coeff[0])*Tstd +Coeff[5])*this->Ru*this->invW[i];
        }
        else
        {
            auto& Coeff = this->LCoeffs[i];
            this->Hf[i] = (((((Coeff[4]*Tstd*0.2+Coeff[3]*0.25)*Tstd+Coeff[2]*(1.0/3.0))*Tstd+Coeff[1]*0.5)*Tstd+Coeff[0])*Tstd +Coeff[5])*this->Ru*this->invW[i];
        }
        
    }


}



OptReaction::~OptReaction
(
)
{
    free(this->buffer);
    free(this->W);
    free(this->invW);
    free(this->tmp_Exp);
    free(this->negGstdByRT);
    free(this->Hf);
    free(this->ThirdBodyFactor1D);
}
