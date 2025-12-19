/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    OptReaction

Description
    A class to optimized the computation of reaction rate

SourceFiles
    Reaction.C

\*---------------------------------------------------------------------------*/

#ifndef OptReaction_H
#define OptReaction_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <immintrin.h>
#include <math.h>
#include <unordered_map>
#include <iomanip>
#include <array>
#include "PtrList.H"
#include "scalarField.H"
#include "fvCFD.H"
#include "HashTable.H"
#include "Macro.H"
#include "vec_math_avx2.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class Reaction Declaration
\*---------------------------------------------------------------------------*/


class OptReaction
{

public:

    // Constructors

        // construct from using OpenFOAM
        OptReaction(const dictionary& chemistryDict,const dictionary& physicalDict);

        // default constructor
        OptReaction();

        // disable copy construction
        OptReaction(const OptReaction& reactions) = delete;

    // Destructor

        ~OptReaction();

    // Member function

        // Read information from OpenFOAM
        void readInfo
        (
            const dictionary& chemistryDict,
            const dictionary& physicalDict  
        );

        // Read reaction information from OpenFOAM, read reaction with integer number
        // e.g. A+B=C
        void readReactionInfo
        (
            std::vector<unsigned int>& lhsIndex,
            std::vector<unsigned int>& rhsIndex,
            const dictionary& nthreaction,
            const hashedWordList& speciesTable
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
            const hashedWordList& speciesTable
        );

        // Find reaction with two reactant and two products
        /*bool findTwoTwoReaction
        (
            const dictionary& nthreaction,
            const hashedWordList& speciesTable
        );*/

        // Check whether the stoichiometric number is an integer, this class only treat 
        // integer.
        bool checkInteger
        (
            const dictionary& nthreaction
        );


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

            mutable std::vector<unsigned int> reaction11index;
            mutable std::vector<unsigned int> reaction12index;
            mutable std::vector<unsigned int> reaction13index;



            mutable std::vector<unsigned int> reaction21index;
            mutable std::vector<unsigned int> reversibleReaction21index;
            mutable std::vector<unsigned int> irreversibleReaction21index;
            mutable std::vector<unsigned int> nonEquilibriumReaction21index;


            mutable std::vector<unsigned int> reaction22index;

            mutable std::vector<unsigned int> reversibleReaction22index;
            mutable std::vector<unsigned int> irreversibleReaction22index;
            mutable std::vector<unsigned int> nonEquilibriumReaction22index;

            mutable std::vector<unsigned int> reaction23index;

            mutable std::vector<unsigned int> reversibleReaction23index;
            mutable std::vector<unsigned int> irreversibleReaction23index;
            mutable std::vector<unsigned int> nonEquilibriumReaction23index;

            mutable std::vector<unsigned int> reaction31index;
            mutable std::vector<unsigned int> reaction32index;
            mutable std::vector<unsigned int> reaction33index;

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
            mutable std::vector<unsigned int> lhsSpeciesIndex1D11;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D11;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D12;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D12;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D13;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D13;


            mutable std::vector<unsigned int> lhsSpeciesIndex1D21;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D21;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D21RR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D21RR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D21IR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D21IR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D21NER;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D21NER;


            mutable std::vector<unsigned int> lhsSpeciesIndex1D22;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D22;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D22RR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D22RR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D22IR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D22IR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D22NER;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D22NER;


            mutable std::vector<unsigned int> lhsSpeciesIndex1D23;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D23;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D23RR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D23RR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D23IR;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D23IR;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D23NER;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D23NER;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D31;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D31;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D32;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D32;

            mutable std::vector<unsigned int> lhsSpeciesIndex1D33;
            mutable std::vector<unsigned int> rhsSpeciesIndex1D33;

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

        // Data to store the parameters of Janf polynomial Fits.

            // Store the parameters of polynomial Fits at high temperature range
            // for all species.
            mutable std::vector<std::array<double,7>> HCoeffs;

            // Store the parameters of polynomial Fits at low temperature range
            // for all species.            
            mutable std::vector<std::array<double,7>> LCoeffs;

            // Store the low temperature limits of all species.
            mutable std::vector<double> Tlow;

            // Store the high temperature limits of all species.
            mutable std::vector<double> Thigh;

            // Store the common temperature limits of all species.           
            mutable std::vector<double> Tcommon;

            // Store the minimum value of Tcommon of all species.
            mutable double TcommonMin;

            // Store the maximum value of Tcommon of all species.
            mutable double TcommonMax;

            // Store the minimum value of Tlow of all species.
            mutable double TlowMin;

            // Store the maximum value of Tlow of all species.
            mutable double ThighMax;

            // Pointer to point the coefficients of all species.
            mutable std::vector<std::array<double,7>*> PtrCoeffs;

        // Thermodynamic parameter

            // Moleweight of species. [kg/kmol]
            mutable double* W;

            // The reciprocal of moleweight of species.
            mutable double* invW;

            // The standard Gibbs energy divided by RT
            mutable double* negGstdByRT;

            // The standard enthalpy of formation
            mutable double* Hf;

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



            // mixture density
            mutable double rhoM;

            mutable double logP;
      
        // Array
            mutable double* buffer; 

    // Member function

        // Thermodynamic function

            // Pointing the janaf coefficient using the given temperature.
            //inline void setPtrCoeffs (const double T) const noexcept;

            // Compute the thermodynamic variable for Jacobian matrix
                // YTp:             Array for species mass fraction, temperature and pressure.
                // c_:              Species concentration.                  [kmol/m^3]
                // dNdtByV:         Instantaneous reaction rate.            [kmol/m^3/s]
                // ExpNegGstdByRT:  exp(-Gstd/(Ru*T)).                      [dimLess]
                // Gstd:            Standard Gibbs free energy.             [J/kmol]
                // dBdT:            ddT(-Gstd)/Ru, Partial derivative Gstd 
                //                  w.r.t temperature divided by Ru.        [dimLess]
                // dCpdT:           Partial derivative of Specific heat 
                //                  capacity w.r.t temperature.             [dimLess]
                // Cp:              Specific heat capacity at constant 
                //                  pressure.                               [J/kg/k]
                // Ha:              Absolute enthalpy.                      [J/kg]
                // rhoMByRhoi:      Mixture density multiplied by specific  
                //                  volume (divided by density of i th species)
                //                                                          [dimLess].
                // WiByrhoM:        Molecular weight of i-th species divided
                //                  by mixture density                      [m^3/kmol].
            /*inline void JacobianThermo
            (
                double p,
                double T,
                double* __restrict__ YTp,
                double* __restrict__ c_,
                double* __restrict__ ExpNegGstdByRT,
                double* __restrict__ dBdT,
                double* __restrict__ dCpdT,
                double* __restrict__ Cp,
                double* __restrict__ Ha,
                double* __restrict__ rhoMByRhoi,
                double* __restrict__ WiByrhoM
            ) const noexcept;*/

            // Compute the thermodynamic variable for reaction rate.
            /*inline void CpHaExpNegGstdByRT
            (
                const double T,
                double* __restrict__ Cp,
                double* __restrict__ Ha,
                double* __restrict__ ExpNegGstdByRT
            )const noexcept;*/

            // Compute the exp(-Gstd/(Ru*T)).           
            //inline void ExpNegGstdByRT (const double T, double* ExpNegGstdByRT) const;

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

            void evalSRIPartialDerivative(double T) const noexcept;

            void evalSRIRateConstant(double T) const noexcept;

            // Compute the F using SRI form.
            /*inline double SRI_F
            (
                const double T,
                const double Pr,
                const unsigned i
            ) const;*/

            // Compute the F, the derivative of F w.r.t temperature
            // and the derivative of F w.r.t Pr.
            // Pr: reduced pressure.
            /*inline void SRI_F_dFdT_dFdPr
            (
                const double T,
                const double Pr,
                const unsigned int i,
                double& __restrict__ F,
                double& __restrict__ dFdT,
                double& __restrict__ dFdPr
            ) const;*/

            // Vectorized
            /*inline void SRI_F_dFdT_dFdPr
            (
                const double T,
                __m256d& __restrict__ Pr,
                __m256d& __restrict__ F,
                __m256d& __restrict__ dFdT,
                __m256d& __restrict__ dFdPr,
                const unsigned int i
            ) const;*/

            void evalTroeRateConstant() const noexcept;

            void evalTroePartialDerivative() const noexcept;

            //void Troe_F_1() const noexcept;

            //void Troe_F_2() const noexcept;         
            
            //void Troe_F_3() const noexcept;   

            // Compute the F, the derivative of F w.r.t temperature
            // and the derivative of F w.r.t Pr.
            // Pr: reduced pressure.

            //void Troe_Jac_1() const noexcept;

            //void Troe_Jac_2() const noexcept;

            //void Troe_Jac_3() const noexcept;


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


            // for one-one reactions. e.g. A<=>B.
            void updateJacobian11
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            /*void inline JF11
            (
                unsigned int i,
                double Kf,
                double dKfdT,
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;*/

            // for one-two reactions. e.g. A<=>B+B, A<=>B+C.
            void updateJacobian12
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;

            /*void inline JF12
            (
                unsigned int i,
                double Kf,
                double dKfdT,
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;*/

            // for one-three reactions. e.g. A<=>B+B+C, A<=>B+C+D.
            /*void inline JF13
            (
                unsigned int i,
                double Kf,
                double dKfdT,
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;*/
            void updateJacobian13
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;
            // for two-one reactions. e.g. A+A<=>B, A+B<=>C.
            /*void inline JF21
            (
                unsigned int i,
                double Kf,
                double dKfdT,
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;  */       
            void updateJacobian21
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;
            // for two-two reactions. e.g. A+A<=>B+C, A+B<=>C+D.
            /*void inline JF22
            (
                unsigned int i,
                double Kf,
                double dKfdT,
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;   */
            void updateJacobian22
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;
            void updateJacobian22ReversibleReaction
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;
            void updateJacobian22IrreversibleReaction
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;
            void updateJacobian22NonEquilibriumReaction
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;
            // for two-three reactions. e.g. A+A<=>B+C+D, A+B<=>C+D+D.
            /*void inline JF23
            (
                unsigned int i,
                double Kf,
                double dKfdT,
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;*/ 
            void updateJacobian23
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;
            // for three-one reactions. e.g. A+A+A<=>B, A+B+C<=>D.
            /*void inline JF31
            (
                unsigned int i,
                double Kf,
                double dKfdT,
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;*/  
            void updateJacobian31
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;
            // for three-two reactions. e.g. A+A+A<=>B+C, A+B+C<=>D+E.
            /*void inline JF32
            (
                unsigned int i,
                double Kf,
                double dKfdT,
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept; */ 
            void updateJacobian32
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;
            // for three-three reactions. e.g. A+A+A<=>B+C+D, A+B+C<=>D+E+F.
            /*void inline JF33
            (
                unsigned int i,
                double Kf,
                double dKfdT,
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;*/  
            void updateJacobian33
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;
            // for general reactions with integer number. e.g. A+A+A+A<=>B+C+D, A+B+C<=>D+E+F.
            /*void inline JFGI
            (
                unsigned int i,
                double Kf,
                double dKfdT,
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;*/  
            void updateJacobianGlobalIntegerReaction
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;
            // for general reactions with non-integer number. e.g. 3.5A<=>B+C+D.
            /*void inline JFGNI
            (
                unsigned int i,
                double Kf,
                double dKfdT,
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;*/  
            void updateJacobianGlobalNonIntegerReaction
            (
                const double* __restrict__ C,
                double* __restrict__ dNdtByV,
                double* __restrict__ ddNdtByVdcTp,
                const double* __restrict__ ExpNegGbyRT,
                const double* __restrict__ dBdT
            )const noexcept;
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
            /*void dNdtByV
            (
                double p,
                double T,
                double* __restrict__ c,
                double* __restrict__ dNdtByV,
                double* __restrict__ Cp,
                double* __restrict__ Ha
            ) const noexcept;*/

            // Compute the reaction rate
            void dNdtByV
            (
                double p,
                double T,
                double* __restrict__ c,
                double* __restrict__ dNdtByV
            ) const noexcept;

            // for one-one reactions. e.g. A<=>B.
            /*void inline RF11
            (
                unsigned int i,
                double Kf,
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            // for one-two reactions. e.g. A<=>B+B, A<=>B+C.
            void inline RF12
            (
                unsigned int i,
                double Kf,
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
         

            // for one-three reactions. e.g. A<=>B+B+C, A<=>B+C+D.
            void inline RF13
            (
                unsigned int i,
                double Kf,
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            // for two-one reactions. e.g. A+A<=>B, A+B<=>C.
            void inline RF21
            (
                unsigned int i,
                double Kf,
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            // for two-two reactions. e.g. A+A<=>B+C, A+B<=>C+D.
            void inline RF22
            (
                unsigned int i,
                double Kf,
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            // for two-three reactions. e.g. A+A<=>B+C+D, A+B<=>C+D+D.
            void inline RF23
            (
                unsigned int i,
                double Kf,
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            // for three-one reactions. e.g. A+A+A<=>B, A+B+C<=>D.
            void inline RF31
            (
                unsigned int i,
                double Kf,
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            // for three-two reactions. e.g. A+A+A<=>B+C, A+B+C<=>D+E.
            void inline RF32
            (
                unsigned int i,
                double Kf,
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            // for three-three reactions. e.g. A+A+A<=>B+C+D, A+B+C<=>D+E+F.
            void inline RF33
            (
                unsigned int i,
                double Kf,
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;*/


             // for reactions with reactants' or products' number larger than 4 
             // e.g. A+A+A+B<=>B+C+D, A+B+C<=>D+D+E+F.
            /*void inline RFGI
            (
                unsigned int i,
                double Kf,
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;*/

             // for reactions with non-integer number
             // e.g. 1.5A <=> 1B+0.5C.
            /*void inline RFGNI
            (
                unsigned int i,
                double Kf,
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;*/

            void updateGlobalNonIntegerReaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void updateGlobalIntegerReaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update11Reaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update12Reaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update13Reaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update21Reaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update21ReversibleReaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update21IrreversibleReaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update21NonEquilibriumReaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update22Reaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update22ReversibleReaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update22IrreversibleReaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update22NonEquilibriumReaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update23Reaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update23ReversibleReaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
            void update23IrreversibleReaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
            void update23NonEquilibriumReaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update31Reaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update32Reaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;

            void update33Reaction
            (
                const double* __restrict__ C,
                double*  __restrict__ dNdtByV,
                const double*  __restrict__ ExpNegGbyRT
            )const noexcept;
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

 // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "ThermoFunctionI.H"
#include "ChemicalEquilibriumFunctionI.H"
#include "thirdBodyReactionI.H"
#include "pressureReactionFunctionI.H"

#include "ReactionFunctionI.H"

#ifdef NoRepository

    #include "EvalTroePartialDerivative.C"
    #include "EvalTroeRateConstant.C"
    #include "EvalSRIPartialDerivative.C"
    #include "EvalSRIRateConstant.C"
    #include "EvalLindemannPartialDerivative.C"
    #include "EvalLindemannRateConstant.C"
    #include "JacobianFunction.C"
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
#endif




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
