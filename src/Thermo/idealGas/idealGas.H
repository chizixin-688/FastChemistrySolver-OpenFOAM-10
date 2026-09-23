/*---------------------------------------------------------------------------*\
  Description
      Header file for the ideal gas object (class idealGas).

      The class stores the NASA / Janaf 7-coefficient polynomial fits of all
      species (high / low temperature ranges) and provides the thermo
      properties required by the chemical solver:
        - Cp, Ha, standard Gibbs energy (-Gstd/(Ru*T)) and their temperature
          derivatives,
        - molar concentrations and mixture density / mole weight / volume,
        - the Jacobian (Y-T and N-T formulations) and derivative (dPhidt)
          thermo quantities.

      All kernels are implemented with AVX-2 (256-bit, 4 species per group)
      SIMD code and read the coefficients directly from the SoA (column-major)
      storage LCoeffsSoA / HCoeffsSoA. Temperature regime selection:
        - T <= TcommonMin  -> all species use the low-temperature fit
        - T >= TcommonMax  -> all species use the high-temperature fit
        - otherwise        -> per-species selection via vblendvpd (mixed)

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. Standard C++ library headers
//---------------------------------
#include <vector>
#include <array>

//---------------------------------
// 2. FastChemistry headers
//---------------------------------
#include "Macro.H"

//---------------------------------
// 3. SIMD / AVX2 headers
//---------------------------------
#include <immintrin.h>

//=============================================================================//

#ifndef idealGas_H
#define idealGas_H

namespace FastChemistry
{


// Use janaf and idealGas assumption
// This module mainly provides the function required by chemical reaction
class idealGas
{
public:

    // Constructors

        // Disable default constructor
        idealGas() = delete;

        // Disable copy construction
        idealGas(const idealGas& gas) = delete;

        // Move constructor
        idealGas(const idealGas&& gas) = delete;

        // Construct using number of species
        idealGas(const int n);

    // Destructor
        ~idealGas();

public:

    // Public Data

        //  Number of species
            mutable int nSpecies;

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
            mutable double TlowMax;

            // Store the minimum value of Thigh of all species.
            mutable double ThighMin;

            // Store the maximum value of Thigh of all species.
            mutable double ThighMax;

            // SoA (structure-of-arrays) coefficient storage.
            // Column-major layout: LCoeffsSoA[j][i] = LCoeffs[i][j].
            // Every segment holds alignSpecies (a multiple of 4) doubles and is
            // 64-byte aligned, so 256-bit aligned loads are always safe.
            mutable double* bufferSoA;              // 64-byte aligned, 15*alignSpecies doubles
            mutable double* LCoeffsSoA[7];          // low-temperature coefficients, column-major
            mutable double* HCoeffsSoA[7];          // high-temperature coefficients, column-major
            mutable double* TcommonSoA;             // per-species crossover temperature

        // Thermodynamic parameter

            // Data buffer
            mutable double* buffer;

            // Moleweight of species [kg/kmol]
            mutable double* W;

            // The reciprocal of moleweight of species [kmol/kg]
            mutable double* invW;

            // The standard enthalpy of formation [J/kg]
            mutable double* Hf;

            // The reciprocal of moleweight of species multiplied with Ru [kmol/kg]
            mutable double* invWRu;

            // Universal gas constant. [J/(kmol·K)]
            const static constexpr double Ru = 8314.470066505449722171761095523834228515625; 

            // Inverse of universal gas constant. [(kmol·K)/J]
            const static constexpr double invRu = 1.0/Ru;

            // Standard pressure [Pa]
            // consistent with OpenFOAM implementation
            const static constexpr double Pstd = 100000;

            // Standard temperature [K]
            const static constexpr double Tstd = 298.15;    

        // Thermophysical variable will be changed during computation

            // log(T)
            mutable double logT;

            // 1/T
            mutable double invT;

            // T*T
            mutable double sqrT;

            // log(P)
            mutable double logP;

            // mixture density [kg/m^3]
            mutable double rhoM;

            // specific volume of mixture  [m^3/kg]
            // = 1.0/rhoM
            mutable double vM;

            // mixture moleweight [kg/kmol]
            mutable double Mw;

            // mixture total volume [m^3]
            mutable double V;

            // 1.0/V [1/m^3]
            mutable double invV;

            mutable double invp;

            mutable double invNtot;

            mutable double Ntot;

            mutable double invRuT;
public:

    // Member Functions

        // Read Function

            //void readThermodynamicsInformation(const dictionary& physicalDict);

        // Member functions for chemistry

            // Compute the thermodynamic parameters required by chemical Jacobian matrix
            // d(dPhidt)dPhi. Phi=[Y0 Y1 Y2 ... YNs T]
            // Phi:             Array for all species mass fraction and temperature
            // c:               Species concentration.                  [kmol/m^3]
            // negGstdByRT:     (-Gstd/(Ru*T)).                         [dimLess]
            // Gstd:            Standard Gibbs free energy.             [J/kmol]
            // dBdT:            ddT(-Gstd)/Ru, Partial derivative Gstd 
            //                  w.r.t temperature divided by Ru.        [dimLess]
            // dCpdT:           Partial derivative of Specific heat 
            //                  capacity w.r.t temperature.             [dimLess]
            // Cp:              Specific heat capacity at constant 
            //                  pressure.                               [J/kg/k]
            // Ha:              Absolute enthalpy.                      [J/kg]
            // rhoMByRhoi:      Mixture density multiplied by specific  
            //                  volume (divided by density of i-th species)
            //                                                          [dimLess].
            // WiByrhoM:        Molecular weight of i-th species divided
            //                  by mixture density                      [m^3/kmol].
            void JacobianThermoYT
            (
                double p,
                double T,
                double* __restrict__ Phi,
                double* __restrict__ c,
                double* __restrict__ negGstdByRT,
                double* __restrict__ negGstdByRT2,
                double* __restrict__ dBdT,
                double* __restrict__ dCpdT,
                double* __restrict__ Cp,
                double* __restrict__ Ha,
                double* __restrict__ rhoMByRhoi,
                double* __restrict__ WiByrhoM
            ) const;

            void JacobianThermoNT
            (
                double p,
                double T,
                double* __restrict__ Phi,
                double* __restrict__ c,
                double* __restrict__ negGstdByRT,
                double* __restrict__ negGstdByRT2,
                double* __restrict__ dBdT,
                double* __restrict__ dcpdT,
                double* __restrict__ cp,
                double* __restrict__ ha
            ) const;

            // Compute the thermodynamic parameters required by chemical derivative
            // (dPhidt). Phi=[Y0 Y1 Y2 ... YNs T]
            // Phi:             Array for all species mass fraction and temperature
            // c:               Species concentration.                  [kmol/m^3]
            // negGstdByRT:     (-Gstd/(Ru*T)).                         [dimLess]
            // Gstd:            Standard Gibbs free energy.             [J/kmol]
            // Cp:              Specific heat capacity at constant 
            //                  pressure.                               [J/kg/k]
            // Ha:              Absolute enthalpy.                      [J/kg]
            // negGstdByRT:     (-Gstd/(Ru*T)).                         [dimLess]
            void DerivativeThermoYT
            (
                double T,
                const double p,
                double* __restrict__ Phi,
                double* __restrict__ c,
                double* __restrict__ Cp,
                double* __restrict__ Ha,
                double* __restrict__ negGstdByRT,
                double* __restrict__ negGstdByRT2
            ) const;


            // Compute the thermodynamic parameters required by chemical derivative
            // (dPhidt). Phi=[N0 N1 N2 ... NNs T]
            // Phi:             Array for all species mass fraction and temperature
            // c:               Species concentration.                  [kmol/m^3]
            // negGstdByRT:     (-Gstd/(Ru*T)).                         [dimLess]
            // Gstd:            Standard Gibbs free energy.             [J/kmol]
            // Cp:              Specific heat capacity at constant 
            //                  pressure.                               [J/kg/k]
            // Ha:              Absolute enthalpy.                      [J/kg]
            // negGstdByRT:     (-Gstd/(Ru*T)).                         [dimLess]
            void DerivativeThermoNT
            (
                double T,
                const double p,
                double* __restrict__ Phi,
                double* __restrict__ c,
                double* __restrict__ Cp,
                double* __restrict__ Ha,
                double* __restrict__ negGstdByRT,
                double* __restrict__ negGstdByRT2
            ) const;


            // Compute the thermodynamic variable required by chemistry time scale
            // Tc
            // negGstdByRT:     (-Gstd/(Ru*T)).                         [dimLess]
            void negGstdByRT
            (
                double T,
                double* __restrict__ negGstdByRT,
                double* __restrict__ negGstdByRT2
            )const;

    // Thermo inline function

            // Volumetric coefficient of thermal expansion [1/T]
            inline double alphav()
            {
                return this->invT;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif