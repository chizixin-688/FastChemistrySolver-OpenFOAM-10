


#ifndef idealGas_H
#define idealGas_H

#include "fvCFD.H"
#include "HashTable.H"
#include <vector>
#include <array>
#include <immintrin.h>
#include "Macro.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




namespace FastThermo
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
        idealGas(const idealGas&& gas);

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
            mutable double ThighMax;

            // Pointer to point the coefficients of all species.
            mutable std::vector<std::array<double,7>*> PtrCoeffs;

        // Thermodynamic parameter

            // Data buffer
            mutable double* buffer;

            // Moleweight of species [kg/kmol]
            mutable double* W;

            // The reciprocal of moleweight of species [kmol/kg]
            mutable double* invW;

            // The standard enthalpy of formation [J/kg]
            mutable double* Hf;

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

public:

    // Member Functions

        // Read Function

            //void readThermodynamicsInformation(const dictionary& physicalDict);

        // Member functions for chemistry

            // Pointing the janaf coefficient using the given temperature.
            inline void setPtrCoeffs (const double T) const;

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


            // Compute the thermodynamic variable required by chemistry derivatives
            // dPhidt. Phi=[Y0 Y1 Y2 ... YNs T]
            // Cp:              Specific heat capacity at constant 
            //                  pressure.                               [J/kg/k]
            // Ha:              Absolute enthalpy.                      [J/kg]
            // negGstdByRT:     (-Gstd/(Ru*T)).                         [dimLess]
            void CpHaNegGstdByRT
            (
                double T,
                double* __restrict__ Cp,
                double* __restrict__ Ha,
                double* __restrict__ negGstdByRT
            )const;

            // Compute the thermodynamic variable required by chemistry derivatives
            // dPhidt. Phi=[Y0 Y1 Y2 ... YNs T]
            // rhoM:            Mixture density                         [kg/m^3]
            // c:               Species concentration                   [kmol/m^3]
            void rhoMc
            (
                double T,
                const double p,
                double* __restrict__ Phi,
                double* __restrict__ c
            )const;

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

    // SIMD inline function

            // s = v0 + v1 + v2 + v3
            inline double hsum4(__m256d v) const;
            inline __m256d hsum4x4(__m256d sum0, __m256d sum1,__m256d sum2, __m256d sum3) const;
            inline void transpose4x4_pd(__m256d& v0,__m256d& v1,__m256d& v2,__m256d& v3)const;
        
            inline double get_elem0(__m256d vec) const;

            inline double get_elem1(__m256d vec) const;

            inline double get_elem2(__m256d vec) const;

            inline double get_elem3(__m256d vec) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif