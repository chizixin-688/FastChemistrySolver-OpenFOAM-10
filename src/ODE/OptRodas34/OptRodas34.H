/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::rodas23

Description
    L-stable, stiffly-accurate embedded Rosenbrock ODE solver of order (2)3.

    References:
    \verbatim
        Sandu, A., Verwer, J. G., Blom, J. G., Spee, E. J., Carmichael, G. R.,
        & Potra, F. A. (1997).
        Benchmarking stiff ODE solvers for atmospheric chemistry problems II:
        Rosenbrock solvers.
        Atmospheric environment, 31(20), 3459-3472.
    \endverbatim

SourceFiles
    rodas34.C

\*---------------------------------------------------------------------------*/

#ifndef OptRodas34_H
#define OptRodas34_H

#include "chemistrySolver.H"
#include "simpleMatrix.H"
#include <iostream>
#include <vector>
#include <algorithm>

#include "LUsolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class OptRodas34 Declaration
\*---------------------------------------------------------------------------*/
template<class ChemistryModel>
class OptRodas34
:
    public chemistrySolver<ChemistryModel>
{
    // Private Data
    class stepState
    {
        public:


        scalar dxTry;
        bool last;

        stepState(const scalar dx)
        :
            dxTry(dx),
            last(false)
        {}
    };

        //- Coefficients dictionary
        dictionary coeffsDict_;

    // ODE solver

        //- Absolute convergence tolerance per step
        const scalar absTol_;

        //- Relative convergence tolerance per step
        const scalar relTol_;

        //- The maximum number of sub-steps allowed for the integration step
        const label maxSteps_;
        
    // Adaptive Solver

        //- Step-size adjustment controls
        constexpr static const scalar safeScale_ =     0.9;
        constexpr static const scalar alphaInc_  =     0.2;
        constexpr static const scalar alphaDec_ =     0.25;
        constexpr static const scalar minScale_ =      0.2;
        constexpr static const scalar maxScale_ =      10;

        //constexpr static const scalar ratio = std::pow(maxScale_/safeScale_, -1.0/alphaInc_);
        constexpr static const scalar ratio = 5.904900000E-06;
        constexpr static const scalar safeMaxScale = safeScale_*maxScale_;

    // Rodas34

        constexpr static const scalar c2 = 0.386;
        constexpr static const scalar c3 = 0.21;
        constexpr static const scalar c4 = 0.63;
        constexpr static const scalar d1 =  0.25;
        constexpr static const scalar d2 = -0.1043;
        constexpr static const scalar d3 =  0.1035;
        constexpr static const scalar d4 = -0.3620000000000023e-01;
        constexpr static const scalar a21 =  0.1544e1;
        constexpr static const scalar a31 =  0.9466785280815826;
        constexpr static const scalar a32 =  0.2557011698983284;
        constexpr static const scalar a41 =  0.3314825187068521e1;
        constexpr static const scalar a42 =  0.2896124015972201e1;
        constexpr static const scalar a43 =  0.9986419139977817;
        constexpr static const scalar a51 =  0.1221224509226641e1;
        constexpr static const scalar a52 =  0.6019134481288629e1;
        constexpr static const scalar a53 =  0.1253708332932087e2;
        constexpr static const scalar a54 = -0.6878860361058950;
        constexpr static const scalar c21 = -0.56688e1;
        constexpr static const scalar c31 = -0.2430093356833875e1;
        constexpr static const scalar c32 = -0.2063599157091915;
        constexpr static const scalar c41 = -0.1073529058151375;
        constexpr static const scalar c42 = -0.9594562251023355e1;
        constexpr static const scalar c43 = -0.2047028614809616e2;
        constexpr static const scalar c51 =  0.7496443313967647e1;
        constexpr static const scalar c52 = -0.1024680431464352e2;
        constexpr static const scalar c53 = -0.3399990352819905e2;
        constexpr static const scalar c54 =  0.1170890893206160e2;
        constexpr static const scalar c61 =  0.8083246795921522e1;
        constexpr static const scalar c62 = -0.7981132988064893e1;
        constexpr static const scalar c63 = -0.3152159432874371e2;
        constexpr static const scalar c64 =  0.1631930543123136e2;
        constexpr static const scalar c65 = -0.6058818238834054e1;
        constexpr static const scalar gamma = 0.25;
        constexpr static const scalar Invgamma= 1.0/gamma;



    // Rodas34

        mutable FastChemistry::LUsolver LU;

public:

    //- Runtime type information
    TypeName("OptRodas34");


    // Constructors

        //- Construct from thermo
        OptRodas34(const fluidReactionThermo& thermo);


    //- Destructor
    virtual ~OptRodas34();


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
            const double p,
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
            double* __restrict__ Phi0,
            double* __restrict__ PhiTemp,
            double* __restrict__ k1,
            double* __restrict__ k2,
            double* __restrict__ k3,
            double* __restrict__ k4,
            double* __restrict__ k5,
            double* __restrict__ dy,
            double* __restrict__ err,
            double* __restrict__ dydx,
            double* __restrict__ dfdx,
            double* __restrict__ a_
        )const;


    // adaptive solver

        void adaptiveSolve
        (
            scalar& __restrict__ x,
            const label li,
            scalar& __restrict__ dxTry,
            const double p,
            double* __restrict__ Phi0,
            double* __restrict__ PhiTemp,
            double* __restrict__ k1,
            double* __restrict__ k2,
            double* __restrict__ k3,
            double* __restrict__ k4,
            double* __restrict__ k5,
            double* __restrict__ dy,
            double* __restrict__ err,
            double* __restrict__ dydx,
            double* __restrict__ dfdx,
            double* __restrict__ a_
        )const;

        scalar Rodas34Solve
        (
            const scalar x0,
            const label li,
            const scalar dx,
            const scalar invdx,
            const double p,
            double* __restrict__ Phi0,
            double* __restrict__ PhiTemp,
            double* __restrict__ k1,
            double* __restrict__ k2,
            double* __restrict__ k3,
            double* __restrict__ k4,
            double* __restrict__ k5,
            double* __restrict__ dy,
            double* __restrict__ err,
            double* __restrict__ dydx,
            double* __restrict__ dfdx,
            double* __restrict__ a_
        )const;


        //- Return the normalised scalar error
        inline scalar normaliseError
        (
            const double* __restrict__ Phi0,
            const double* __restrict__ Phi,
            const double* __restrict__ err
        )const;         
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "OptRodas34.C"
#endif

#endif

// ************************************************************************* //
