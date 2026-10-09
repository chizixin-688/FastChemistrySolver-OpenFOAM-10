/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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
    Foam::FastChemistryModel

Description
    Chemistry model with an AVX-vectorised, cell-wise ODE solver for stiff
    combustion chemistry.  Two Jacobian options are provided: a fast
    diagonal approximation and the exact Jacobian (see FCjacobianType_).
    The mass-fraction (Y) and molar-number (N) Jacobian kernels are declared
    in Jacobian/JacobianFuncDeclaration.H and implemented in the
    Jacobian.C translation units included below.

SourceFiles
    FastChemistryModelI.H
    FastChemistryModel.C

Author
    Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. Standard C++ library headers
//---------------------------------
#include <chrono>

//---------------------------------
// 2. OpenFOAM library headers
//---------------------------------
#include "basicChemistryModel.H"

//---------------------------------
// 3. FastChemistry headers
//---------------------------------
#include "LUsolver.H"
#include "dataBlock.H"
#include "simpleDataBlock.H"
#include "idealGas.H"
#include "OptReaction.H"
#include "createIdealGasFromFoamDict.H"
#include "createGasPhaseReactionFromFoamDict.H"
#include "vec_math_avx2.H"
//=============================================================================//

#ifndef FastChemistryModel_H
#define FastChemistryModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class FastChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class UnusedThermo>
class FastChemistryModel
:
    public basicChemistryModel
{
    // Public data
public:
        //- Pointer for ideal gas
        std::unique_ptr<FastChemistry::idealGas> gas;

        bool massFraction;
    // Private data
private:
        //- Reference to the field of specie mass fractions
        const PtrList<volScalarField>& Yvf_;

        //- Number of species
        label nSpecie_;

        //- Type of the Jacobian to be calculated in Fast Chemistry Model
        // 0:Fast 1:Exact
        int FCjacobianType_;

        //- Gas phase reactions

        mutable std::unique_ptr<FastChemistry::OptReaction> GasReaction;

        //- List of reaction rate per specie [kg/m^3/s]
        PtrList<volScalarField> RR_;

        //- Kernel function pointers for the Jacobian blocks below; the actual
        //  implementation is selected when the model is constructed according
        //  to FCjacobianType_ (0 = fast approximation, 1 = exact).
        using ddYdtdYPtrType = void(FastChemistryModel<UnusedThermo>::*)
        (
            const double* __restrict__ ,
            const double* __restrict__ ,
            const double* __restrict__ ,
            const double* __restrict__ ,
            const double* __restrict__ ,
            const double* __restrict__ ,
            double* __restrict__ ,
            double* __restrict__ ,
            double rhoM
        ) const noexcept;
        ddYdtdYPtrType  ddYdtdYPtr = nullptr;

        using ddYdtdTPtrType = void(FastChemistryModel<UnusedThermo>::*)
        (
            const double* __restrict__ ddNdtByVdcT,
            const double* __restrict__ WiByrhoM,
            const double* __restrict__ c,
            double* __restrict__ dPhidt,
            double* __restrict__ Jac,
            double alphavM
        ) const noexcept;
        ddYdtdTPtrType  ddYdtdTPtr = nullptr;

        using ddTdtdYTPtrType = void(FastChemistryModel<UnusedThermo>::*)
        (
            const double* __restrict__ Cp,
            const double* __restrict__ dCpdT,
            const double* __restrict__ Ha,
            double* __restrict__ dPhidt,
            double* __restrict__ Jac
        ) const noexcept;
        ddTdtdYTPtrType  ddTdtdYTPtr = nullptr;

public:
        //- Data buffer
        mutable double*     buffer;
        mutable double*     YTpWork[12];
        mutable double*     YTpYTpWork[3];

        //- Size of the ODE system
        label n_;
        mutable unsigned int alignN=0;

private:

        //- Minimum reaction temperature
        const scalar Treact;

        //****************************MPI Communication*****************************//

        //- CPU time larger tave*DLBthreshold is busy process, and will send data to other process.
        scalar DLBthreshold = 1.0;

        //- Update the load-balancing trigger each iteration
        label MaxIter = 1;
        label CurrentIter = 0;

        //- CPU time for data transfer
        List<List<int64_t>> cpuLoadTransferTable;

        //- CPU time to solve chemistry for each cell
        std::vector<std::pair<int64_t,label>> CPUtimeField;

        //- CPU time spent solving chemistry per process
        List<int64_t> chemistryIntegrationTime;

        //- Size of send buffer
        List<std::streamsize> sendBufferSize_;

        //- Size of receive buffer
        List<std::streamsize> recvBufferSize_;

        //- Send buffer, containing the binary data block
        List<DynamicList<char>> sendBuffer_;

        //- Receive buffer, containing the binary data block
        List<DynamicList<char>> recvBuffer_;

        //- The receive position of receive buffer
        List<label> recvBufPos_;

        //- first time this process is executed
        bool firstTime = true;

        //- Skip this cell or solve?
        List<bool> skip;

        //- Determine whether this process is busy
        List<bool> IamBusyProcess;

        //- Whether to switch on the load balancing
        bool Balance;

        void exchange
        (
            const UList<DynamicList<char>>& sendBufs,
            const List<std::streamsize>& recvSizes,
            List<DynamicList<char>>& recvBufs,
            const int tag,
            const label comm,
            const bool block
        );

        //****************************MPI Communication*****************************//

    // Private Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<volScalarField>& RR();

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  Variable number of species added
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT){return 0;}

public:

    //- Runtime type information
    TypeName("FastChemistryModel");

    // Constructors

        //- Construct from thermo
        FastChemistryModel(const fluidReactionThermo& thermo);

        //- Disallow default bitwise copy construction
        FastChemistryModel(const FastChemistryModel&) = delete;

    //- Destructor
    virtual ~FastChemistryModel();

    // Member Functions

        //- The number of reactions
        virtual inline label nReaction() const;

        //- Calculates the reaction rates
        virtual void calculate();

        // Return the reaction object
        inline std::unique_ptr<FastChemistry::OptReaction>&  getReaction() const;

        //- The number of species
        inline virtual label nSpecie() const;

        //- Return const access to the chemical source terms for specie, i
        inline const volScalarField::Internal& RR(const label i) const;

        //- Return non const access to chemical source terms [kg/m^3/s]
        virtual volScalarField::Internal& RR(const label i);

        //- Return the reaction rate of species i in reactionI
        virtual tmp<volScalarField::Internal> calculateRR
        (
            const label reactionI,
            const label speciei
        ) const;

        //- Solve the reaction system for the given time step
        //  and return the characteristic time
        virtual scalar solve(const scalar deltaT);
        virtual scalar solve(const scalarField& deltaT);

        //- Return the chemical time scale
        virtual tmp<volScalarField> tc() const;

        //- Return the heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const;

    // ODE interface used by the cell-wise stiff-integration solvers

        virtual void derivatives
        (
            const scalar t,
            const label li,
            const double p,
            double* __restrict__ YT,
            double* __restrict__ dYTdt,
            double* __restrict__ Cp,
            double* __restrict__ Ha
        ) const;

        virtual void jacobian
        (
            const scalar t,
            const label li,
            const double p,
            double* __restrict__ Phi,
            double* __restrict__ dPhidt,
            double* __restrict__ Jac
        ) const;

        //      Compute the change rate of mass fraction and temperature during
        //  constant pressure adiabatic combustion
        // t                :Current time                                       [s]
        // li               :Current cell index                                 [dimLess]
        // c[i]             :Molar concentration                                [kmol/m^3]
        // Cp[i]            :Specific heat capacity                             [J/kg/K]
        // Ha[i]            :Specific absolute enthalpy                         [J/kg]
        // dPhidt[i]        :Mass fraction change rate                          [1/s]
        // Cp[nSpecie()]    :Mixture specific heat capacity                     [J/kg/K]
        // Phi[nSpecie()]   :Mixture temperature                                [K]
        // dPhidt[nSpecie()]:Mixture temperature change rate                    [K/s]
        virtual void getdYTdt
        (
            const scalar t,
            const label li,
            const double p,
            double* __restrict__ YT,
            double* __restrict__ dYTdt,
            double* __restrict__ c,
            double* __restrict__ Cp,
            double* __restrict__ Ha
        ) const;

        //      Compute the change rate of molar number and temperature during
        //  constant pressure adiabatic combustion
        // t                :Current time                                       [s]
        // li               :Current cell index                                 [dimLess]
        // c[i]             :Molar concentration                                [kmol/m^3]
        // cp[i]            :Specific molar heat capacity                       [J/kmol/K]
        // ha[i]            :Specific molar absolute enthalpy                   [J/kmol]
        // dPhidt[i]        :Molar number change rate                           [1/s]
        // cp[nSpecie()]    :Mixture specific molar heat capacity               [J/kmol/K]
        // Phi[nSpecie()]   :Mixture temperature                                [K]
        // dPhidt[nSpecie()]:Mixture temperature change rate                    [K/s]
        virtual void getdNTdt
        (
            const scalar t,
            const label li,
            const double p,
            double* __restrict__ NT,
            double* __restrict__ dNTdt,
            double* __restrict__ c,
            double* __restrict__ Cp,
            double* __restrict__ Ha
        ) const;

        // Compute the Jacobian matrix of change rate of mass fraction and
        // temperature during constant pressure adiabatic combustion
        //      ODE variable:
        // t                :Current time                                       [s]
        // li               :Current cell index                                 [dimLess]
        //      i-th species variables (size:nSpecie):
        // c[i]             :Molar concentration                                [kmol/m^3]
        // Cp[i]            :Specific heat capacity                             [J/kg/K]
        // Ha[i]            :Specific absolute enthalpy                         [J/kg]
        // dPhidt[i]        :Mass fraction change rate                          [1/s]
        // dBdT:            ddT(-Gstd)/Ru, Partial derivative Gstd
        //                  w.r.t temperature divided by Ru.                    [dimLess]
        // dCpdT:           Partial derivative of Specific heat
        //                  capacity w.r.t temperature.                         [J/kg]
        // rhoMByRhoi:      Mixture density multiplied by specific
        //                  volume (divided by density of i-th species)         [dimLess].
        // WiByrhoM:        Molecular weight of i-th species divided
        //                  by mixture density                                  [m^3/kmol].
        //      Mixture variable:
        // Cp[nSpecie()]    :Mixture specific heat capacity                     [J/kg/K]
        // Phi[nSpecie()]   :Mixture temperature                                [K]
        // dPhidt[nSpecie()]:Mixture temperature change rate                    [K/s]
        void getddYTdtdYT
        (
            const scalar t,
            const label li,
            const double p,
            double* __restrict__ c,
            double* __restrict__ YT,
            double* __restrict__ dYTdt,
            double* __restrict__ Cp,
            double* __restrict__ Ha,
            double* __restrict__ dBdT,
            double* __restrict__ dCpdT,
            double* __restrict__ rhoMByRhoi,
            double* __restrict__ WiByrhoM,
            double* __restrict__ ddNdtByVdcT,
            double* __restrict__ dcdY,
            double* __restrict__ ddYTdtdYT
        ) const;

        void getddNTdtdNT
        (
            const scalar t,
            const label li,
            const double p,
            double* __restrict__ c,
            double* __restrict__ NT,
            double* __restrict__ dNTdt,
            double* __restrict__ cp,
            double* __restrict__ ha,
            double* __restrict__ dBdT,
            double* __restrict__ dcpdT,
            double* __restrict__ ddNdtByVdcT,
            double* __restrict__ ddNTdtdNT
        ) const;

        void getddNTdtdNTwithMergeI
        (
            const scalar t,
            const label li,
            const double p,
            const double scale,
            double* __restrict__ c,
            double* __restrict__ NT,
            double* __restrict__ dNTdt,
            double* __restrict__ cp,
            double* __restrict__ ha,
            double* __restrict__ dBdT,
            double* __restrict__ dcpdT,
            double* __restrict__ ddNdtByVdcT,
            double* __restrict__ ddNTdtdNT
        ) const;

        //- Pure-virtual driver overloads (implemented by the concrete model)
        virtual void solve
        (
            scalar& p,
            scalar& T,
            scalarField& Y,
            const label li,
            scalar& deltaT,
            scalar& subDeltaT
        ) const = 0;

        virtual void solve
        (
            const label li,
            double T,
            double& __restrict__ deltaT,
            double& __restrict__ subDeltaT
        ) const = 0;

        // Jacobian functions
        #include "JacobianFuncDeclaration.H"

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const FastChemistryModel&) = delete;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FastChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FastChemistryModel.C"
    #include "FastddYdtdY_Vec44_0.C"
    #include "FastddYdtdY_Vec44_1.C"
    #include "FastddYdtdY_Vec44_2.C"
    #include "FastddYdtdY_Vec44_3.C"
    #include "ddYdtdY_Vec844_0.C"
    #include "ddYdtdY_Vec844_1.C"
    #include "ddYdtdY_Vec844_2.C"
    #include "ddYdtdY_Vec844_3.C"
    #include "ddYdtdY_Vec844_4.C"
    #include "ddYdtdY_Vec844_5.C"
    #include "ddYdtdY_Vec844_6.C"
    #include "ddYdtdY_Vec844_7.C"
    #include "ddYdtdT_Vec88_0.C"
    #include "ddYdtdT_Vec88_1.C"
    #include "ddYdtdT_Vec88_2.C"
    #include "ddYdtdT_Vec88_3.C"
    #include "ddYdtdT_Vec88_4.C"
    #include "ddYdtdT_Vec88_5.C"
    #include "ddYdtdT_Vec88_6.C"
    #include "ddYdtdT_Vec88_7.C"
    #include "ddTdtdYT_Vec88_0.C"
    #include "ddTdtdYT_Vec88_1.C"
    #include "ddTdtdYT_Vec88_2.C"
    #include "ddTdtdYT_Vec88_3.C"
    #include "ddTdtdYT_Vec88_4.C"
    #include "ddTdtdYT_Vec88_5.C"
    #include "ddTdtdYT_Vec88_6.C"
    #include "ddTdtdYT_Vec88_7.C"
    #include "ddNdtdNT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

