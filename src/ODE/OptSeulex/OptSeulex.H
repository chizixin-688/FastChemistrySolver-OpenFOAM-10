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
    Foam::OptSeulex

Description
    An extrapolation-algorithm, based on the linearly implicit Euler method
    with step size control and order selection.

    Reference:
    \verbatim
        Hairer, E., Nørsett, S. P., & Wanner, G. (1996).
        Solving Ordinary Differential Equations II:
        Stiff and Differential-Algebraic Problems, second edition",
        Springer-Verlag, Berlin.
    \endverbatim
    
SourceFiles
    OptSeulex.C

\*---------------------------------------------------------------------------*/

#ifndef OptSeulex_H
#define OptSeulex_H

#include <iostream>
#include <vector>
#include <algorithm>

#include "chemistrySolver.H"
#include "simpleMatrix.H"


#include "LUsolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class OptSeulex Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class OptSeulex
:
    public chemistrySolver<ChemistryModel>
{
    // Private Data

    class stepState
    {
        public:

        const bool forward;
        scalar dxTry;
        scalar dxDid;
        bool first;
        bool last;
        bool reject;
        bool prevReject;

        stepState(const scalar dx)
        :
            forward(dx > 0 ? true : false),
            dxTry(dx),
            dxDid(0),
            first(true),
            last(false),
            reject(false),
            prevReject(false)
        {}
    };


        //- Coefficients dictionary
        dictionary coeffsDict_;

    // ODE solver

        //- Absolute convergence tolerance per step
        const scalar absTol_=1e-8;

        //- Relative convergence tolerance per step
        const scalar relTol_=1e-4;

        scalar logTol;

        //- The maximum number of sub-steps allowed for the integration step
        const label maxSteps_;


    // Seulex

        // Static constants

        static const label kMaxx_ = 12;
        static const label iMaxx_ = kMaxx_ + 1;
        constexpr static const scalar stepFactor1_ = 0.8;
        constexpr static const scalar stepFactor2_ = 0.93;
        constexpr static const scalar stepFactor3_ = 0.1;
        constexpr static const scalar stepFactor4_ = 4;
        constexpr static const scalar stepFactor5_ = 0.5;
        constexpr static const scalar kFactor1_ = 0.7;
        constexpr static const scalar kFactor2_ = 0.9;

        // Evaluated constants

            scalar jacRedo_;
            labelField nSeq_;
            scalarField cpu_;
            scalarField invCpu_;
            scalarSquareMatrix coeff_;

        // Temporary storage
        // held to avoid dynamic memory allocation between calls
        // and to transfer internal values between functions

            mutable scalar theta_;
            mutable label kTarg_;
            mutable scalarRectangularMatrix table_;

            mutable labelList pivotIndices_;

            // Fields space for "solve" function


            mutable scalarField dxOpt_, temp_;
            //mutable scalarField y0_, ySequence_, scale_;
            
            mutable double* y0_;
            mutable double* ySequence_;
            mutable double* scale_;

            // Fields used in "seul" function
            //mutable scalarField dy_, yTemp_, dydx_;

        mutable FastChemistry::LUsolver LU;




public:

    //- Runtime type information
    TypeName("OptSeulex");


    // Constructors

        //- Construct from thermo
        OptSeulex(const fluidReactionThermo& thermo);


    //- Destructor
    virtual ~OptSeulex();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalar& __restrict__ p,
            scalar& __restrict__ T,
            scalarField& c,
            const label li,
            scalar& __restrict__ deltaT,
            scalar& __restrict__ subDeltaT
        ) const;


        virtual void solve
        (
            const label li,
            double p,
            double& __restrict__ deltaT,
            double& __restrict__ subDeltaT
        ) const;
    // ODE solver

        void ODESolve
        (
            const scalar xEnd,
            const label li,
            scalar& dxTry,
            const double p,
            double* __restrict__ y0,
            double* __restrict__ yTemp_,
            double* __restrict__ Cp,
            double* __restrict__ Ha,
            double* __restrict__ dy_,
            double* __restrict__ dydx_,
            double* __restrict__ k9,
            double* __restrict__ k11,
            double* __restrict__ Jy,
            double* __restrict__ a
        )const;

        void SeulexSolve
        (
            scalar& x,
            const label li,
            stepState& step,
            const double p,
            double* __restrict__ y,
            double* __restrict__ yTemp__,
            double* __restrict__ Cp,
            double* __restrict__ Ha,
            double* __restrict__ dy_,
            double* __restrict__ dydx_,
            double* __restrict__ k9,
            double* __restrict__ k11,            
            double* __restrict__ Jy,
            double* __restrict__ a            
        )const;

    // Private Member Functions

        //- Computes the j-th line of the extrapolation table
        bool seul
        (
            const scalar x0,
            double* __restrict__ y0,
            const label li,
            const scalar dxTot,
            const label k,
            const double p,
            double* __restrict__ y,
            double* __restrict__ yTemp_,
            double* __restrict__ Cp,
            double* __restrict__ Ha,
            double* __restrict__ dy_,
            double* __restrict__ dydx_,
            double* __restrict__ k11,             
            double* __restrict__ Jy,
            double* __restrict__ a                
        ) const;

        //- Polynomial extrapolation
        void extrapolate
        (
            const label k,
            scalarRectangularMatrix& table,
            double* y
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "OptSeulex.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
