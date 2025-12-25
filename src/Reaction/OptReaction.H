/*---------------------------------------------------------------------------*\
  Description
      Header file for gas phase reaction objects. Contains class declarations 
      and function prototypes related to gas-phase chemistry reactions.

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

#ifndef OptReaction_H
#define OptReaction_H

//=============================================================================//

//---------------------------------
// 1. Standard C++ library headers
//---------------------------------
#include <vector>
#include <string>
#include <math.h>
#include <unordered_map>

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


namespace FastChemistry
{

/*---------------------------------------------------------------------------*\
                        Class Reaction Declaration
\*---------------------------------------------------------------------------*/
class OptReaction
{

public:

    // Constructors

        // construct from using OpenFOAM
        //OptReaction(const dictionary& chemistryDict,const dictionary& physicalDict);

        // default constructor
        OptReaction();

        // disable copy construction
        OptReaction(const OptReaction& reactions) = delete;

    // Destructor

        ~OptReaction();

    // Member function

        // Read information from OpenFOAM
        /*void readInfo
        (
            const dictionary& chemistryDict,
            const dictionary& physicalDict  
        );*/

        // Read reaction information from OpenFOAM, read reaction with integer number
        // e.g. A+B=C
        /*void readReactionInfo
        (
            std::vector<unsigned int>& lhsIndex,
            std::vector<unsigned int>& rhsIndex,
            const dictionary& nthreaction,
            const hashedWordList& speciesTable
        );*/

        // Read reaction information from OpenFOAM, read reaction with non-integer number
        // e.g. A+1.5B=C, A+B^1.5=B
        /*void readReactionInfo
        (
            std::vector<unsigned int>& lhsIndex,
            std::vector<double>& lhsStoichCoeff,
            std::vector<double>& lhsReactionOrder,
            std::vector<unsigned int>& rhsIndex,
            std::vector<double>& rhsStoichCoeff,
            std::vector<double>& rhsRectionOrder,
            const dictionary& nthreaction,
            const hashedWordList& speciesTable
        );*/

        // Find reaction with two reactant and two products
        /*bool findTwoTwoReaction
        (
            const dictionary& nthreaction,
            const hashedWordList& speciesTable
        );*/

        // Check whether the stoichiometric number is an integer, this class only treat 
        // integer.
        /*bool checkInteger
        (
            const dictionary& nthreaction
        );*/


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
            mutable std::vector<std::string> reactionTable_;  

        // Data to store the numbers of reactions and species.

            // The number of all reactions.
            mutable unsigned int n_Reactions=0;

            // The number of Arrhenius reversible reactions.
            mutable unsigned int n_Arrhenius=0;

            // The number of reversible reactions with reverse Arrhenius rate constant.
            mutable unsigned int n_NonEquilibriumReversibleArrhenius=0;       

            // The number of third body reactions with reverse Arrhenius rate constant.
            mutable unsigned int n_NonEquilibriumThirdBodyReaction=0;

            // The number of Lindemann reaction(including Falloff and ChemicallyActivated).
            mutable unsigned int n_Lindemann=0;

            // The number of Lindemann reaction(Falloff).
            mutable unsigned int n_LindemannFO=0;

            // The number of Lindemann reaction(ChemicallyActivated).
            mutable unsigned int n_LindemannCA=0;

            // The number of Troe reaction(including Falloff and ChemicallyActivated).
            mutable unsigned int n_Troe=0;

            // The number of Troe reaction(Falloff).
            mutable unsigned int n_TroeFO=0;

            // The number of Troe reaction(ChemicallyActivated).
            mutable unsigned int n_TroeCA=0;

            // The number of SRI reaction(including Falloff and ChemicallyActivated).
            mutable unsigned int n_SRI=0;

            // The number of SRI reaction(Falloff).
            mutable unsigned int n_SRIFO=0;

            // The number of SRI reaction(ChemicallyActivated).
            mutable unsigned int n_SRICA=0;

            // The number of third body reactions.
            mutable unsigned int n_ThirdBodyReaction=0;     

            // The number of fall-off reactions.
            mutable unsigned int n_Fall_Off_Reaction=0;

            // The number of chemically activated reactions.
            mutable unsigned int n_ChemicallyActivated_Reaction=0;

            // The number of Plog reaction;
            mutable unsigned int n_PlogReaction=0;

            // The number of temperature independent reactions in the type of
            // Arrhenius reversible reactions, it doesn't include 
            // the third body reaction, pressure reaction or other type reaction
            mutable unsigned int n_Temperature_Independent_Reaction=0;

            // Whether this reaction is reversible
            mutable std::vector<char> isIrreversible;

            // Whether this reaction is not elementary reaction
            // (stoichiometric number is not an integer or reaction order is not an integer)
            mutable std::vector<char> isGlobal;

            // offset of k inf
            mutable unsigned int offset_kinf=0;

            // The size of ODE systems, n_ = number of species + 1 (Temperature)
            mutable unsigned int n_ = 0;

            // The size of ODE systems with blank byte padding for memory alignment. e.g.
            // n_ = 3, alignN = 4
            // n_ = 4, alignN = 4
            // n_ = 5, alignN = 8
            mutable unsigned int alignN = 0;

            // The number of chemical species, including the default species.
            mutable unsigned int nSpecies=0;

            // The number of species with blank byte padding
            // nSpecies = 3, AlignSpecies = 4;
            // nSpecies = 4, AlignSpecies = 4;
            // nSpecies = 5, AlignSpecies = 8;
            mutable unsigned int AlignSpecies = 0;
            
            // The number of temporary array, this array is used to store data 
            // for computating its exponential value
            mutable unsigned int tmp_ExpSize = 0;

            // Store the number of different types of third body reaction 
            // Itbr[1] = number of NonEquilibriumThirdBodyArrhenius
            // Itbr[2] = number of NonEquilibriumThirdBodyArrhenius + ThirdBodyArrhenius
            // ...
            mutable unsigned int Itbr[6];

            // Store the number of different types of reaction 
            // Itbr[1] = number of Arrhenius
            // Itbr[2] = number of NonEquilibriumReversibleArrhenius + Arrhenius
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


            // Data to save the index of reactants for all reactions.
            mutable std::vector<std::vector<unsigned int>>     lhsSpeciesIndex;
            mutable std::vector<std::vector<unsigned int>>     rhsSpeciesIndex;

            //mutable std::vector<unsigned int> reaction11index;
            mutable std::vector<unsigned int> reversibleReaction11index;
            mutable std::vector<unsigned int> irreversibleReaction11index;
            mutable std::vector<unsigned int> nonEquilibriumReaction11index;

            //mutable std::vector<unsigned int> reaction12index;
            mutable std::vector<unsigned int> reversibleReaction12index;
            mutable std::vector<unsigned int> irreversibleReaction12index;
            mutable std::vector<unsigned int> nonEquilibriumReaction12index;

            //mutable std::vector<unsigned int> reaction13index;
            mutable std::vector<unsigned int> reversibleReaction13index;
            mutable std::vector<unsigned int> irreversibleReaction13index;
            mutable std::vector<unsigned int> nonEquilibriumReaction13index;


            //mutable std::vector<unsigned int> reaction21index;
            mutable std::vector<unsigned int> reversibleReaction21index;
            mutable std::vector<unsigned int> irreversibleReaction21index;
            mutable std::vector<unsigned int> nonEquilibriumReaction21index;


            //mutable std::vector<unsigned int> reaction22index;
            mutable std::vector<unsigned int> reversibleReaction22index;
            mutable std::vector<unsigned int> irreversibleReaction22index;
            mutable std::vector<unsigned int> nonEquilibriumReaction22index;

            //mutable std::vector<unsigned int> reaction23index;
            mutable std::vector<unsigned int> reversibleReaction23index;
            mutable std::vector<unsigned int> irreversibleReaction23index;
            mutable std::vector<unsigned int> nonEquilibriumReaction23index;

            //mutable std::vector<unsigned int> reaction31index;
            mutable std::vector<unsigned int> reversibleReaction31index;
            mutable std::vector<unsigned int> irreversibleReaction31index;
            mutable std::vector<unsigned int> nonEquilibriumReaction31index;

            //mutable std::vector<unsigned int> reaction32index;
            mutable std::vector<unsigned int> reversibleReaction32index;
            mutable std::vector<unsigned int> irreversibleReaction32index;
            mutable std::vector<unsigned int> nonEquilibriumReaction32index;


            //mutable std::vector<unsigned int> reaction33index;
            mutable std::vector<unsigned int> reversibleReaction33index;
            mutable std::vector<unsigned int> irreversibleReaction33index;
            mutable std::vector<unsigned int> nonEquilibriumReaction33index;


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


            // 1D Data to save the index of reactants for specific reactions.
            //mutable std::vector<unsigned int> lhsSpeciesIndex1D11;
            //mutable std::vector<unsigned int> rhsSpeciesIndex1D11;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D11RR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D11RR;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D11IR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D11IR;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D11NER;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D11NER;

            //mutable std::vector<unsigned int> lhsSpeciesIndex1D12;
            //mutable std::vector<unsigned int> rhsSpeciesIndex1D12;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D12RR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D12RR;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D12IR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D12IR;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D12NER;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D12NER;

            //mutable std::vector<unsigned int> lhsSpeciesIndex1D13;
            //mutable std::vector<unsigned int> rhsSpeciesIndex1D13;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D13RR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D13RR;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D13IR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D13IR;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D13NER;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D13NER;

            //mutable std::vector<unsigned int> lhsSpeciesIndex1D21;
            //mutable std::vector<unsigned int> rhsSpeciesIndex1D21;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D21RR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D21RR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D21IR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D21IR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D21NER;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D21NER;


            //mutable std::vector<unsigned int> lhsSpeciesIndex1D22;
            //mutable std::vector<unsigned int> rhsSpeciesIndex1D22;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D22RR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D22RR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D22IR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D22IR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D22NER;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D22NER;


            //mutable std::vector<unsigned int> lhsSpeciesIndex1D23;
            //mutable std::vector<unsigned int> rhsSpeciesIndex1D23;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D23RR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D23RR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D23IR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D23IR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D23NER;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D23NER;

            //mutable std::vector<unsigned int> lhsSpeciesIndex1D31;
            //mutable std::vector<unsigned int> rhsSpeciesIndex1D31;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D31RR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D31RR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D31IR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D31IR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D31NER;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D31NER;

            //mutable std::vector<unsigned int> lhsSpeciesIndex1D32;
            //mutable std::vector<unsigned int> rhsSpeciesIndex1D32;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D32RR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D32RR;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D32IR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D32IR;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D32NER;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D32NER;

            //mutable std::vector<unsigned int> lhsSpeciesIndex1D33;
            //mutable std::vector<unsigned int> rhsSpeciesIndex1D33;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D33RR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D33RR;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D33IR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D33IR;
            mutable std::vector<unsigned int> lhsSpeciesIndex1D33NER;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D33NER;

            // Data to save the stoichiometric number of reactants for reaction with non-integer
            // stoichiometric number. For integer stoichiometric number, value = 1.
            mutable std::vector<std::vector<double>>  lhsStoichCoeff;
            mutable std::vector<std::vector<double>>  rhsStoichCoeff;

            // Data to save the reaction order of reaction with non-integer stoichiometric number.
            // For integer stoichiometric number, value = 1.
            mutable std::vector<std::vector<double>>  lhsReactionOrder;
            mutable std::vector<std::vector<double>>  rhsReactionOrder;

        // Store the Arrhenius parameters for all reactions. For pressure depedent reaction,
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

        // Store the third dody efficiency for third body reactions, 
        // falloff and chemicallyactivated reactions

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

        // Pressure depedent reaction function

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

        // Math function

            // s = v0 + v1 + v2 + v3
            inline double hsum4(__m256d v) const
            {
                __m128d lo = _mm256_castpd256_pd128(v);
                __m128d hi = _mm256_extractf128_pd(v,1);
                lo = _mm_add_pd(lo,hi);               // [v0+v2, v1+v3]
                __m128d sh = _mm_unpackhi_pd(lo,lo);  // [v1+v3, v1+v3]
                __m128d s = _mm_add_sd(lo, sh);       // [sum, …]
                return _mm_cvtsd_f64(s);
            }
            inline __m256d hsum4x4(__m256d sum0, __m256d sum1,__m256d sum2, __m256d sum3) const
            {
                __m256d t0 = _mm256_hadd_pd(sum0, sum1);  // [a0+a1, b0+b1, a2+a3, b2+b3]
                __m256d t1 = _mm256_hadd_pd(sum2, sum3);  // [c0+c1, d0+d1, c2+c3, d2+d3]
                __m256d t2 = _mm256_permute4x64_pd (t0, _MM_SHUFFLE(3,1,2,0)); // [a0+a1, a2+a3, b0+b1,  b2+b3]
                __m256d t3 = _mm256_permute4x64_pd(t1, _MM_SHUFFLE(3,1,2,0)); // [c0+c1, c2+c3, d0+d1, d2+d3]
                __m256d t4 = _mm256_hadd_pd(t2, t3); // [a0+a1+a2+a3, c0+c1+c2+c3, b0+b1+b2+b3, d0+d1+d2+d3]
                __m256d t5 = _mm256_permute4x64_pd(t4, _MM_SHUFFLE(3,1,2,0));
                return t5;
            }

            inline void transpose4x4_pd(__m256d& v0,__m256d& v1,__m256d& v2,__m256d& v3)const
            {
                __m256d t0 = _mm256_unpacklo_pd(v0, v1); // a0 b0 a1 b1
                __m256d t1 = _mm256_unpackhi_pd(v0, v1); // a2 b2 a3 b3
                __m256d t2 = _mm256_unpacklo_pd(v2, v3); // c0 d0 c1 d1
                __m256d t3 = _mm256_unpackhi_pd(v2, v3); // c2 d2 c3 d3

                v0 = _mm256_permute2f128_pd(t0, t2, 0x20); // [a0 b0 c0 d0]
                v1 = _mm256_permute2f128_pd(t1, t3, 0x20); // [a1 b1 c1 d1]
                v2 = _mm256_permute2f128_pd(t0, t2, 0x31); // [a2 b2 c2 d2]
                v3 = _mm256_permute2f128_pd(t1, t3, 0x31); // [a3 b3 c3 d3]
            }
        
            inline double get_elem0(__m256d vec) const
            {
                return _mm256_cvtsd_f64(vec); 
            }

            inline double get_elem1(__m256d vec) const
            {
                __m128d low = _mm256_castpd256_pd128(vec); 
                return _mm_cvtsd_f64(_mm_unpackhi_pd(low, low));
            }

            inline double get_elem2(__m256d vec) const
            {
                __m128d high = _mm256_extractf128_pd(vec, 1); 
                return _mm_cvtsd_f64(high);
            }

            inline double get_elem3(__m256d vec) const
            {
                __m128d high = _mm256_extractf128_pd(vec, 1); 
                return _mm_cvtsd_f64(_mm_unpackhi_pd(high, high));
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}// End namespace FastChemistry

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#include "ChemicalEquilibriumFunctionI.H"
#include "thirdBodyReactionI.H"

/*#ifdef NoRepository
    #include "EvalTroePartialDerivative.C"
    #include "EvalTroeRateConstant.C"
    #include "EvalSRIPartialDerivative.C"
    #include "EvalSRIRateConstant.C"
    #include "EvalLindemannPartialDerivative.C"
    #include "EvalLindemannRateConstant.C"
    #include "FindPlogPressureRange.C"
    #include "EvalPlogPartialDerivative.C"
    #include "EvalPlogRateConstant.C"
    #include "OptReaction.C"
    #include "Jacobian.C"
    #include "ReactionRate.C"
    #include "Tc.C"
    #include "updateJacobian11.C"
    #include "updateJacobian12.C"
    #include "updateJacobian13.C"
    #include "updateJacobian21.C"
    #include "updateJacobian22.C"
    #include "updateJacobian23.C"
    #include "updateJacobian31.C"
    #include "updateJacobian32.C"
    #include "updateJacobian33.C"
    #include "updateJacobianGlobalIntegerReaction.C"
    #include "updateJacobianGlobalNonIntegerReaction.C"
    #include "update11Reaction.C"
    #include "update12Reaction.C"
    #include "update13Reaction.C"
    #include "update21Reaction.C"
    #include "update22Reaction.C"
    #include "update23Reaction.C"
    #include "update31Reaction.C"
    #include "update32Reaction.C"
    #include "update33Reaction.C"
    #include "updateGlobalIntegerReaction.C"
    #include "updateGlobalNonIntegerReaction.C"
#endif*/




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
