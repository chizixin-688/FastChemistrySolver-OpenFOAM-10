/*---------------------------------------------------------------------------*\
  Description
      Header file for the gas-phase reaction object (class OptReaction).
      Declares the data members and the member functions used to
        - compute the molar reaction rates and the molar-concentration based
          Jacobian matrix,
        - evaluate the pressure dependent rate constants (Lindemann / Troe /
          SRI) and their partial derivatives,
        - evaluate the Plog rate constants,
        - compute the chemistry time scale (Tc).

      Member layout:
        - constructors / destructor / operators
        - public data: reaction information, reaction counters, index tables,
          Arrhenius and pressure-dependent parameters, temporary buffers,
          thermodynamic variables, function-pointer dispatch tables
        - member functions: chemical-equilibrium helpers, pressure-dependent
          kernels (incl. the AVX-512 variants), Jacobian / reaction-rate /
          time-scale functions

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. Standard C++ library headers
//---------------------------------
#include <vector>
#include <string>
#include <math.h>
#include <unordered_map>
#include <array>

//---------------------------------
// 2. FastChemistry headers
//---------------------------------
#include "Macro.H"
#include "vec_math_avx2.H"

//---------------------------------
// 3. SIMD / AVX2 headers
//---------------------------------
#include <immintrin.h>

//=============================================================================//

#ifndef OptReaction_H
#define OptReaction_H

namespace FastChemistry
{

/*---------------------------------------------------------------------------*\
                        Class Reaction Declaration
\*---------------------------------------------------------------------------*/
class OptReaction
{

public:

    // Constructors

        // default constructor
        OptReaction();

        // disable copy construction
        OptReaction(const OptReaction& reactions) = delete;

        // disable move constructor
        OptReaction(const OptReaction&& reaction) = delete;

    // Destructor

        ~OptReaction();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const OptReaction&) = delete;


    // Public data

        // Save the basic information of reactions and species.

            // Contain the species name, e.g. "CO", "H2O".
            mutable std::vector<std::string> speciesTable_;

            // Contain the reaction type name, e.g. "irreversibleArrhenius", "reversibleArrhenius".
            mutable std::vector<std::string> reactionType_;

            // Contain the reaction name, e.g. "un-named-reaction-1", "un-named-reaction-2".
            mutable std::vector<std::string> reactionName_;

            // Contain the chemical reaction equation, e.g. "O + H2 = H + OH","O + H2O2 = OH + HO2".
            mutable std::vector<std::string> reactionEquation_;

        // Data to store the numbers of reactions and species.

            // The number of all reactions.
            mutable unsigned int n_Reactions = 0;

            // The number of Arrhenius reversible reactions.
            mutable unsigned int n_Arrhenius = 0;

            // The number of reversible reactions with reverse Arrhenius rate constant.
            mutable unsigned int n_NonEquilibriumReversibleArrhenius = 0;

            // The number of third body reactions with reverse Arrhenius rate constant.
            mutable unsigned int n_NonEquilibriumThirdBodyReaction = 0;

            // The number of Lindemann reaction(including Falloff and ChemicallyActivated).
            mutable unsigned int n_Lindemann = 0;

            // The number of Lindemann reaction(Falloff).
            mutable unsigned int n_LindemannFO = 0;

            // The number of Lindemann reaction(ChemicallyActivated).
            mutable unsigned int n_LindemannCA = 0;

            // The number of Troe reaction(including Falloff and ChemicallyActivated).
            mutable unsigned int n_Troe = 0;

            // The number of Troe reaction(Falloff).
            mutable unsigned int n_TroeFO = 0;

            // The number of Troe reaction(ChemicallyActivated).
            mutable unsigned int n_TroeCA = 0;

            // The number of SRI reaction(including Falloff and ChemicallyActivated).
            mutable unsigned int n_SRI = 0;

            // The number of SRI reaction(Falloff).
            mutable unsigned int n_SRIFO = 0;

            // The number of SRI reaction(ChemicallyActivated).
            mutable unsigned int n_SRICA = 0;

            // The number of third body reactions.
            mutable unsigned int n_ThirdBodyReaction = 0;

            // The number of fall-off reactions.
            mutable unsigned int n_Fall_Off_Reaction = 0;

            // The number of chemically activated reactions.
            mutable unsigned int n_ChemicallyActivated_Reaction = 0;

            // The number of Plog reaction;
            mutable unsigned int n_PlogReaction = 0;

            // The number of temperature independent reactions in the type of
            // Arrhenius reversible reactions, it doesn't include
            // the third body reaction, pressure reaction or other type reaction
            mutable unsigned int n_Temperature_Independent_Reaction = 0;

            // Whether this reaction is reversible
            // 0: reversible using chemical equilibrium constant
            // 1: irreversible
            // 2: reversible using Arrhenius form
            mutable std::vector<char> isIrreversible;

            // Whether this reaction is not elementary reaction
            // (stoichiometric number is not an integer or reaction order is not an integer)
            mutable std::vector<char> isGlobal;

            // offset of kinf and k0 for pressure reaction
            mutable unsigned int offset_kinf = 0;

            // The size of ODE systems, n_ = number of species + 1 (Temperature)
            mutable unsigned int n_ = 0;

            // The size of ODE systems with blank byte padding for memory alignment. e.g.
            // n_ = 3, alignN = 4
            // n_ = 4, alignN = 4
            // n_ = 5, alignN = 8
            mutable unsigned int alignN = 0;

            // The number of chemical species, including the default species.
            mutable unsigned int nSpecies = 0;

            // The number of species with blank byte padding
            // nSpecies = 3, AlignSpecies = 4;
            // nSpecies = 4, AlignSpecies = 4;
            // nSpecies = 5, AlignSpecies = 8;
            mutable unsigned int AlignSpecies = 0;

            // The number of temporary array, this array is used to store data
            // for computing its exponential value
            mutable unsigned int tmp_ExpSize = 0;

            // Store the index of third body reaction in the reaction list
            // Itbr[1] = number of NonEquilibriumThirdBodyArrhenius
            // Itbr[2] = number of NonEquilibriumThirdBodyArrhenius + ThirdBodyArrhenius
            // ...
            mutable unsigned int Itbr[6];

            // Store the index of reaction in the reaction list, it is used
            // to know where the rate constant of different type of reaction is stored in the Kf_ array.
            // The order is as follows:
            // Ikf[1] = number of Arrhenius
            // Ikf[2] = number of NonEquilibriumReversibleArrhenius + Arrhenius
            // ...
            mutable unsigned int Ikf[13];

        // Data to save the index of falloff and chemicallyactivated reactions.

            // Store the index of the Lindemann form.
            mutable std::vector<unsigned int> LindemannFO;

            // Store the index of the Troe form.
            mutable std::vector<unsigned int> TroeFO;

            // Store the index of the SRI form.
            mutable std::vector<unsigned int> SRIFO;

            // Store the index of the Lindemann form.
            mutable std::vector<unsigned int> LindemannCA;

            // Store the index of the Troe form.
            mutable std::vector<unsigned int> TroeCA;

            // Store the index of the SRI form.
            mutable std::vector<unsigned int> SRICA;

        // Data to save the information of elementary reaction with integer stoichiometric number.
            // Save the index of reaction and species index of reaction
            // example:
            // reversibleR22allindex[0] = 5.    The 5-th reaction in reaction list is the 1-st reversible reaction
            // reversibleR22allindex[1] = 0.    its index of 1-st reactant is 0
            // reversibleR22allindex[2] = 9.    its index of 2-nd reactant is 9
            // reversibleR22allindex[3] = 1.    its index of 1-st product is 1
            // reversibleR22allindex[4] = 2.    its index of 2-nd product is 2
            mutable std::vector<unsigned int> RR11AllIndex;
            mutable std::vector<unsigned int> IR11AllIndex;
            mutable std::vector<unsigned int> NER11AllIndex;

            mutable std::vector<unsigned int> RR12AllIndex;
            mutable std::vector<unsigned int> RR12DupAllIndex;
            mutable std::vector<unsigned int> IR12AllIndex;
            mutable std::vector<unsigned int> IR12DupAllIndex;
            mutable std::vector<unsigned int> NER12AllIndex;
            mutable std::vector<unsigned int> NER12DupAllIndex;

            mutable std::vector<unsigned int> RR13AllIndex;
            mutable std::vector<unsigned int> RR13DupAllIndex;
            mutable std::vector<unsigned int> IR13AllIndex;
            mutable std::vector<unsigned int> IR13DupAllIndex;
            mutable std::vector<unsigned int> NER13AllIndex;
            mutable std::vector<unsigned int> NER13DupAllIndex;

            mutable std::vector<unsigned int> RR21AllIndex;
            mutable std::vector<unsigned int> RR21DupAllIndex;
            mutable std::vector<unsigned int> IR21AllIndex;
            mutable std::vector<unsigned int> IR21DupAllIndex;
            mutable std::vector<unsigned int> NER21AllIndex;
            mutable std::vector<unsigned int> NER21DupAllIndex;

            mutable std::vector<unsigned int> RR22AllIndex;
            mutable std::vector<unsigned int> RR22DupAllIndex;
            mutable std::vector<unsigned int> IR22AllIndex;
            mutable std::vector<unsigned int> IR22DupAllIndex;
            mutable std::vector<unsigned int> NER22AllIndex;
            mutable std::vector<unsigned int> NER22DupAllIndex;

            mutable std::vector<unsigned int> RR23AllIndex;
            mutable std::vector<unsigned int> RR23DupAllIndex;
            mutable std::vector<unsigned int> IR23AllIndex;
            mutable std::vector<unsigned int> IR23DupAllIndex;
            mutable std::vector<unsigned int> NER23AllIndex;
            mutable std::vector<unsigned int> NER23DupAllIndex;

            mutable std::vector<unsigned int> RR32AllIndex;
            mutable std::vector<unsigned int> RR32DupAllIndex;
            mutable std::vector<unsigned int> IR32AllIndex;
            mutable std::vector<unsigned int> IR32DupAllIndex;
            mutable std::vector<unsigned int> NER32AllIndex;
            mutable std::vector<unsigned int> NER32DupAllIndex;

            mutable std::vector<unsigned int> RR31AllIndex;
            mutable std::vector<unsigned int> RR31DupAllIndex;
            mutable std::vector<unsigned int> IR31AllIndex;
            mutable std::vector<unsigned int> IR31DupAllIndex;
            mutable std::vector<unsigned int> NER31AllIndex;
            mutable std::vector<unsigned int> NER31DupAllIndex;

            mutable std::vector<unsigned int> RR33AllIndex;
            mutable std::vector<unsigned int> RR33DupAllIndex;
            mutable std::vector<unsigned int> IR33AllIndex;
            mutable std::vector<unsigned int> IR33DupAllIndex;
            mutable std::vector<unsigned int> NER33AllIndex;
            mutable std::vector<unsigned int> NER33DupAllIndex;

            mutable unsigned int RR11Size;
            mutable unsigned int IR11Size;
            mutable unsigned int NER11Size;

            mutable unsigned int RR12Size;
            mutable unsigned int RR12DupSize;
            mutable unsigned int IR12Size;
            mutable unsigned int IR12DupSize;
            mutable unsigned int NER12Size;
            mutable unsigned int NER12DupSize;

            mutable unsigned int RR13Size;
            mutable unsigned int RR13DupSize;
            mutable unsigned int IR13Size;
            mutable unsigned int IR13DupSize;
            mutable unsigned int NER13Size;
            mutable unsigned int NER13DupSize;

            mutable unsigned int RR21Size;
            mutable unsigned int RR21DupSize;
            mutable unsigned int IR21Size;
            mutable unsigned int IR21DupSize;
            mutable unsigned int NER21Size;
            mutable unsigned int NER21DupSize;

            mutable unsigned int RR22Size;
            mutable unsigned int RR22DupSize;
            mutable unsigned int IR22Size;
            mutable unsigned int IR22DupSize;
            mutable unsigned int NER22Size;
            mutable unsigned int NER22DupSize;

            mutable unsigned int RR23Size;
            mutable unsigned int RR23DupSize;
            mutable unsigned int IR23Size;
            mutable unsigned int IR23DupSize;
            mutable unsigned int NER23Size;
            mutable unsigned int NER23DupSize;

            mutable unsigned int RR32Size;
            mutable unsigned int RR32DupSize;
            mutable unsigned int IR32Size;
            mutable unsigned int IR32DupSize;
            mutable unsigned int NER32Size;
            mutable unsigned int NER32DupSize;

            mutable unsigned int RR31Size;
            mutable unsigned int RR31DupSize;
            mutable unsigned int IR31Size;
            mutable unsigned int IR31DupSize;
            mutable unsigned int NER31Size;
            mutable unsigned int NER31DupSize;

            mutable unsigned int RR33Size;
            mutable unsigned int RR33DupSize;
            mutable unsigned int IR33Size;
            mutable unsigned int IR33DupSize;
            mutable unsigned int NER33Size;
            mutable unsigned int NER33DupSize;

            // Data to save the index of reactants for all reactions.
            mutable std::vector<std::vector<unsigned int>> lhsSpeciesIndex;
            mutable std::vector<std::vector<unsigned int>> rhsSpeciesIndex;

            // The index of global integer reaction
            // e.g. 2A+22B = 222C
            mutable std::vector<unsigned int> reactionGIindex;

            // The index of global non-integer reaction
            // e.g. 2A+2.2B = 2.22C
            mutable std::vector<unsigned int> reactionGNIindex;

            // 1D Data to save the index of reactants for all reactions.
            mutable std::vector<unsigned int> lhsSpeciesIndex1D;
            mutable std::vector<unsigned int> lhsOffset;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D;
            mutable std::vector<unsigned int> rhsOffset;


            // Data to save the stoichiometric number of reactants for reaction with non-integer
            // stoichiometric number
            mutable std::vector<std::vector<double>> lhsStoichCoeff;
            mutable std::vector<std::vector<double>> rhsStoichCoeff;

            // Data to save the reaction order of reaction with non-integer stoichiometric number.
            mutable std::vector<std::vector<double>> lhsReactionOrder;
            mutable std::vector<std::vector<double>> rhsReactionOrder;

        // Store the Arrhenius parameters for all reactions. For pressure dependent reaction,
        // these contains the parameters of low pressure limits and high pressure limits.

            // Arrhenius pre-exp factor.
            mutable std::vector<double>  A;

            // Arrhenius temperature corrector.
            mutable std::vector<double>  beta;

            // Arrhenius activation temperature.
            mutable std::vector<double>  Ta;

        // The data to save the temporary results

            // Data to save the Arrhenius forward rate constants for all reactions.
            mutable std::vector<double> Kf_;

            // Data to save the derivatives of forward rate constants w.r.t
            // temperature for all reactions.
            mutable std::vector<double> dKfdT_;

            // Data to save the additional Exp computation, including
            // species and pressure reaction
            mutable double* tmp_Exp;

            // Data to save the derivatives of forward rate constants w.r.t
            // species concentration for third body reactions and pressure
            // dependent reactions.
            mutable std::vector<double> dKfdC_;

            // Data to save the third body efficiency.
            mutable std::vector<double> tmp_M;

        // Store the third body efficiency for third body reactions,
        // falloff and chemicallyactivated reactions

            // Each row stores the three body efficiency factor of a reaction
            mutable double* ThirdBodyFactor1D;

        // Store the parameters of Troe form for falloff and chemicallyactivated reaction.

            mutable std::vector<double> alpha_;

            mutable std::vector<double> Ts_;

            mutable std::vector<double> invTs_;

            mutable std::vector<double> Tss_;

            mutable std::vector<double> Tsss_;

            mutable std::vector<double> invTsss_;

        // Store the parameters of SRI form for falloff and chemicallyactivated reaction.

            mutable std::vector<double> a_;

            mutable std::vector<double> b_;

            mutable std::vector<double> c_;

            mutable std::vector<double> invc_;

            mutable std::vector<double> d_;

            mutable std::vector<double> e_;

        // General pressure dependence using logarithmic interpolation

            mutable std::vector<std::vector<double>> APlog;

            mutable std::vector<std::vector<double>> logAPlog;

            mutable std::vector<std::vector<double>> betaPlog;

            mutable std::vector<std::vector<double>> TaPlog;

            mutable std::vector<std::vector<double>> Prange;

            mutable std::vector<std::vector<double>> rDeltaP_;

            mutable std::vector<std::vector<double>> logPi;

            mutable std::vector<unsigned int> Pindex;

            // Index of reaction that requires linear interpolation
            mutable std::vector<unsigned int> ActivePlogReactionIndex;

            // Number of reaction that requires linear interpolation
            mutable unsigned int ActivePlogReactionNumber = 0;

        // Thermodynamic parameter

            // The standard Gibbs energy divided by RT
            mutable double* negGstdByRT;

            // The standard Gibbs energy divided by RT
            mutable double* invNegGstdByRT;

            // Universal gas constant. [J/(kmol·K)]
            const static constexpr double Ru = 8314.470066505449722171761095523834228515625;

            // Inverse of universal gas constant. [J/(kmol·K)]
            const static constexpr double invRu = 1.0/Ru;

            // Standard pressure. [Pa]
            // consistent with OpenFOAM implementation
            const static constexpr double Pstd = 100000;
            //const static constexpr double Pstd = 101325;

        // Chemical equilibrium calculation

            // Table to save the temporary result for chemical equilibrium calculation.
            // For global reaction(stoichiometric number is an integer)
            mutable std::unordered_map<int,double> Pow_pByRT_SumVki_I;

            // Table to save the temporary result for chemical equilibrium calculation.
            // 5 element version of (Pstd/(Ru*T))^sumVki;
            mutable std::array<double,5> Pow_pByRT_SumVki;

        // Thermophysical variable

            // log(T).
            mutable double logT;

            // T
            mutable double T;

            // 1/T.
            mutable double invT;

            // log(P)
            mutable double logP;

            // log(Ru*T)
            mutable double logRuT;

        // Array
            //mutable double* buffer;


        // Pointer for reaction function
        using RFptrType = void(OptReaction::*)
        (
            const double* __restrict__ C,
            double*  __restrict__ dNdtByV,
            const double*  __restrict__ ExpNegGbyRT
        ) const noexcept;

        std::vector<RFptrType>  RFptr;

        // Pointer for jacobian function
        using JFptrType = void(OptReaction::*)
        (
            const double* __restrict__ C,
            double* __restrict__ dNdtByV,
            double* __restrict__ ddNdtByVdcTp,
            const double* __restrict__ ExpNegGbyRT,
            const double* __restrict__ dBdT
        ) const noexcept;

        std::vector<JFptrType>  JFptr;

    // Member function

        // Chemical equilibrium function

            // Update the variable (Pstd/(Ru*T))^sumVki for global reactions.
            inline void update_Pow_pByRT_SumVki(const double T) const;

            // Update the variable (Pstd/(Ru*T))^sumVki for element reactions.
            inline void update_Pow_pByRT_SumVki2(const double T) const;

        // Third body reaction function

            // Return the derivatives of third body efficiency w.r.t species concentration.
            inline double dMdC(const unsigned int ri, const unsigned int cj) const noexcept;

        // Pressure dependent reaction function

            void evalLindemannPartialDerivative() const noexcept;

            void evalLindemannRateConstant() const noexcept;

            void evalSRIPartialDerivative() const noexcept;

            void evalSRIRateConstant() const noexcept;

            void evalTroeRateConstant() const noexcept;

            void evalTroePartialDerivative() const noexcept;

            void findPlogPressureRange(double p) const noexcept;

            void evalPlogRateConstant() const noexcept;

            void evalPlogPartialDerivative() const noexcept;

        // Jacobian matrix function

            // Compute the Jacobian matrix and reaction rate
            void ddNdtByVdcTp
            (
                const double p,
                const double T,
                double* __restrict__ YTp,
                double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ dBdT,
                double* __restrict__ ddNdtByVdcT
            ) const noexcept;


            #include "JcFuncDeclaration.H"

        // Reaction rate function

            // Chemistry time scale(OpenFOAM-10)
            void Tc
            (
                int celli,
                double p,
                double T,
                double* C,
                double& sumW,
                double& sumWRateByCTot
            ) const noexcept;

            // Compute the reaction rate
            void dNdtByV
            (
                double p,
                double T,
                double* __restrict__ c,
                double* __restrict__ dNdtByV
            ) const noexcept;
            
            #include "ReactionFuncDeclaration.H"
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}// End namespace FastChemistry

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#include "ChemicalEquilibriumFunctionI.H"
#include "thirdBodyReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
