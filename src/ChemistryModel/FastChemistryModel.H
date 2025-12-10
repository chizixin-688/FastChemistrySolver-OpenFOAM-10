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

SourceFiles
    FastChemistryModelI.H
    FastChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef FastChemistryModel_H
#define FastChemistryModel_H


#include "basicChemistryModel.H"
#include "dataBlock.H"
#include "simpleDataBlock.H"
#include "OptReaction.H"

#include "idealGas.H"
#include "createIdealGasFromFoamDict.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class FastChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class FastChemistryModel
:
    public basicChemistryModel
{
    // Private data

        //- Pointer for ideal gas
        std::unique_ptr<FastThermo::idealGas> gas;

        //- Reference to the field of specie mass fractions
        const PtrList<volScalarField>& Yvf_;

        //- Number of species
        label nSpecie_;

        //- Type of the Jacobian to be calculated
        const jacobianType jacobianType_;

        //- Gas phase reactions
        mutable OptReaction reaction;

        //- List of reaction rate per specie [kg/m^3/s]
        PtrList<volScalarField> RR_;

        using ddYdtdYPtrType = void(FastChemistryModel<ThermoType>::*)
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

        using ddYdtdTPtrType = void(FastChemistryModel<ThermoType>::*)
        (
            const double* __restrict__ ddNdtByVdcT,
            const double* __restrict__ WiByrhoM,
            const double* __restrict__ c,
            double* __restrict__ dPhidt,
            double* __restrict__ Jac,
            double alphavM
        ) const noexcept;
        ddYdtdTPtrType  ddYdtdTPtr = nullptr;


        using ddTdtdYTPtrType = void(FastChemistryModel<ThermoType>::*)
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

            //- Updte load balancing per Iter
            label MaxIter = 1;
            label CurrentIter = 0;

            //- CPU time for data transfer
            List<List<int64_t>> cpuLoadTransferTable;

            //- CPU time to solve chemistry for each cell
            std::vector<std::pair<int64_t,label>> CPUtimeField;

            //- CPU time so solve chemistry for each process
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

            //- firt time to be executed
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
        inline OptReaction&  getReaction() const;

        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- The number of species
            inline virtual label nSpecie() const;

            //- Return const access to the chemical source terms for specie, i
            inline const volScalarField::Internal& RR
            (
                const label i
            ) const;

            //- Return non const access to chemical source terms [kg/m^3/s]
            virtual volScalarField::Internal& RR
            (
                const label i
            );

            //- Return reaction rate of the speciei in reactionI
            virtual tmp<volScalarField::Internal> calculateRR
            (
                const label reactionI,
                const label speciei
            ) const;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT) ;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT) ;

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return the heat release rate [kg/m/s^3]
            virtual tmp<volScalarField> Qdot() const;

        // ODE functions (overriding abstract functions in ODE.H)
            virtual void derivatives
            (
                const scalar t,
                const label li,
                const double p,
                double* __restrict__ Phi,
                double* __restrict__ dPhidt,
                double* __restrict__ Cp,
                double* __restrict__ Ha
            ) const ;

            virtual void jacobian
            (
                const scalar t,
                const label li,
                const double p,
                double* __restrict__ Phi,
                double* __restrict__ dPhidt,
                double* __restrict__ Jac
            ) const ;

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
            /*void FastddYdtdY_Vec44_0
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ rhoMByRhoi,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ dPhidt,
                const double* __restrict__ invW,
                const double* __restrict__ Phi,
                double* __restrict__ dcdY,
                double* __restrict__ Jac,
                double rhoM
            ) const noexcept;

            void FastddYdtdY_Vec44_1
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ rhoMByRhoi,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ dPhidt,
                const double* __restrict__ invW,
                const double* __restrict__ Phi,
                double* __restrict__ dcdY,
                double* __restrict__ Jac,
                double rhoM
            ) const noexcept;

            void FastddYdtdY_Vec44_2
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ rhoMByRhoi,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ dPhidt,
                const double* __restrict__ invW,
                const double* __restrict__ Phi,
                double* __restrict__ dcdY,
                double* __restrict__ Jac,
                double rhoM
            ) const noexcept;

            void FastddYdtdY_Vec44_3
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ rhoMByRhoi,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ dPhidt,
                const double* __restrict__ invW,
                const double* __restrict__ Phi,
                double* __restrict__ dcdY,
                double* __restrict__ Jac,
                double rhoM
            ) const noexcept;

            void ddYdtdY_Vec844_0
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ rhoMByRhoi,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ dPhidt,
                const double* __restrict__ invW,
                const double* __restrict__ Phi,
                double* __restrict__ dcdY,
                double* __restrict__ Jac,
                double rhoM
            ) const noexcept;

            void ddYdtdY_Vec844_1
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ rhoMByRhoi,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ dPhidt,
                const double* __restrict__ invW,
                const double* __restrict__ Phi,
                double* __restrict__ dcdY,
                double* __restrict__ Jac,
                double rhoM
            ) const noexcept;

            void ddYdtdY_Vec844_2
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ rhoMByRhoi,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ dPhidt,
                const double* __restrict__ invW,
                const double* __restrict__ Phi,
                double* __restrict__ dcdY,
                double* __restrict__ Jac,
                double rhoM
            ) const noexcept;

            void ddYdtdY_Vec844_3
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ rhoMByRhoi,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ dPhidt,
                const double* __restrict__ invW,
                const double* __restrict__ Phi,
                double* __restrict__ dcdY,
                double* __restrict__ Jac,
                double rhoM
            ) const noexcept;

            void ddYdtdY_Vec844_4
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ rhoMByRhoi,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ dPhidt,
                const double* __restrict__ invW,
                const double* __restrict__ Phi,
                double* __restrict__ dcdY,
                double* __restrict__ Jac,
                double rhoM
            ) const noexcept;

            void ddYdtdY_Vec844_5
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ rhoMByRhoi,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ dPhidt,
                const double* __restrict__ invW,
                const double* __restrict__ Phi,
                double* __restrict__ dcdY,
                double* __restrict__ Jac,
                double rhoM
            ) const noexcept;

            void ddYdtdY_Vec844_6
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ rhoMByRhoi,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ dPhidt,
                const double* __restrict__ invW,
                const double* __restrict__ Phi,
                double* __restrict__ dcdY,
                double* __restrict__ Jac,
                double rhoM
            ) const noexcept;

            void ddYdtdY_Vec844_7
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ rhoMByRhoi,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ dPhidt,
                const double* __restrict__ invW,
                const double* __restrict__ Phi,
                double* __restrict__ dcdY,
                double* __restrict__ Jac,
                double rhoM
            ) const noexcept;

            void ddYdtdT_Vec88_0
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ c,
                double* __restrict__ dPhidt,
                double* __restrict__ Jac,
                double alphavM
            ) const noexcept;

            void ddYdtdT_Vec88_1
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ c,
                double* __restrict__ dPhidt,
                double* __restrict__ Jac,
                double alphavM
            ) const noexcept;

            void ddYdtdT_Vec88_2
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ c,
                double* __restrict__ dPhidt,
                double* __restrict__ Jac,
                double alphavM
            ) const noexcept;

            void ddYdtdT_Vec88_3
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ c,
                double* __restrict__ dPhidt,
                double* __restrict__ Jac,
                double alphavM
            ) const noexcept;

            void ddYdtdT_Vec88_4
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ c,
                double* __restrict__ dPhidt,
                double* __restrict__ Jac,
                double alphavM
            ) const noexcept;

            void ddYdtdT_Vec88_5
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ c,
                double* __restrict__ dPhidt,
                double* __restrict__ Jac,
                double alphavM
            ) const noexcept;

            void ddYdtdT_Vec88_6
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ c,
                double* __restrict__ dPhidt,
                double* __restrict__ Jac,
                double alphavM
            ) const noexcept;

            void ddYdtdT_Vec88_7
            (
                const double* __restrict__ ddNdtByVdcT,
                const double* __restrict__ WiByrhoM,
                const double* __restrict__ c,
                double* __restrict__ dPhidt,
                double* __restrict__ Jac,
                double alphavM
            ) const noexcept;

            void ddTdtdYT_Vec88_0
            (
                const double* __restrict__ Cp,
                const double* __restrict__ dCpdT,
                const double* __restrict__ Ha,      
                double* __restrict__ dPhidt,  
                double* __restrict__ Jac
            ) const noexcept;

            void ddTdtdYT_Vec88_1
            (
                const double* __restrict__ Cp,
                const double* __restrict__ dCpdT,
                const double* __restrict__ Ha,      
                double* __restrict__ dPhidt,  
                double* __restrict__ Jac
            ) const noexcept;

            void ddTdtdYT_Vec88_2
            (
                const double* __restrict__ Cp,
                const double* __restrict__ dCpdT,
                const double* __restrict__ Ha,      
                double* __restrict__ dPhidt,  
                double* __restrict__ Jac
            ) const noexcept;

            void ddTdtdYT_Vec88_3
            (
                const double* __restrict__ Cp,
                const double* __restrict__ dCpdT,
                const double* __restrict__ Ha,      
                double* __restrict__ dPhidt,  
                double* __restrict__ Jac
            ) const noexcept;

            void ddTdtdYT_Vec88_4
            (
                const double* __restrict__ Cp,
                const double* __restrict__ dCpdT,
                const double* __restrict__ Ha,      
                double* __restrict__ dPhidt,  
                double* __restrict__ Jac
            ) const noexcept;

            void ddTdtdYT_Vec88_5
            (
                const double* __restrict__ Cp,
                const double* __restrict__ dCpdT,
                const double* __restrict__ Ha,      
                double* __restrict__ dPhidt,  
                double* __restrict__ Jac
            ) const noexcept;

            void ddTdtdYT_Vec88_6
            (
                const double* __restrict__ Cp,
                const double* __restrict__ dCpdT,
                const double* __restrict__ Ha,      
                double* __restrict__ dPhidt,  
                double* __restrict__ Jac
            ) const noexcept;

            void ddTdtdYT_Vec88_7
            (
                const double* __restrict__ Cp,
                const double* __restrict__ dCpdT,
                const double* __restrict__ Ha,      
                double* __restrict__ dPhidt,  
                double* __restrict__ Jac
            ) const noexcept;*/

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const FastChemistryModel&) = delete;


    // SIMD inline function

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
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
