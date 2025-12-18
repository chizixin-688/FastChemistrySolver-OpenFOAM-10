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
    Foam::OptRosenbrock34

Description
    L-stable embedded Rosenbrock ODE solver of order (3)4.

    \verbatim
        Hairer, E., Nørsett, S. P., & Wanner, G. (1996).
        Solving Ordinary Differential Equations II:
        Stiff and Differential-Algebraic Problems, second edition",
        Springer-Verlag, Berlin.
    \endverbatim

    The default constants are from:
    \verbatim
        Shampine, L. F. (1982).
        Implementation of Rosenbrock Methods.
        ACM Transactions on Mathematical Software, vol. 8, pp. 93–113.
    \endverbatim
    with which the scheme is more accurate than with the L-Stable coefficients
    for small step-size but less stable for large step-size.

    The L-Stable scheme constants are provided commented-out in Rosenbrock34.C
    
SourceFiles
    OptRosenbrock34.C

\*---------------------------------------------------------------------------*/

#ifndef OptRosenbrock34_H
#define OptRosenbrock34_H

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
                        Class OptRosenbrock34 Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class OptRosenbrock34
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

    // Rosenbrock

        constexpr static const scalar a21= 2;
        constexpr static const scalar a31= 48.0/25.0;
        constexpr static const scalar a32= 6.0/25.0;

        constexpr static const scalar c21= -8;
        constexpr static const scalar c31= 372.0/25.0;
        constexpr static const scalar c32= 12.0/5.0;

        constexpr static const scalar c41= -112.0/125.0;
        constexpr static const scalar c42= -54.0/125.0;
        constexpr static const scalar c43= -2.0/5.0;

        constexpr static const scalar b1= 19.0/9.0;
        constexpr static const scalar b2= 1.0/2.0;
        constexpr static const scalar b3= 25.0/108.0;
        constexpr static const scalar b4= 125.0/108.0;

        constexpr static const scalar e1= 34.0/108.0;
        constexpr static const scalar e2= 7.0/36.0;
        constexpr static const scalar e3= 0;
        constexpr static const scalar e4= 125.0/108.0;

        constexpr static const scalar gamma= 1.0/2.0;
        constexpr static const scalar Invgamma= 1.0/gamma;
        constexpr static const scalar c2= 1;
        constexpr static const scalar c3= 3.0/5.0;

        constexpr static const scalar d1= 1.0/2.0;
        constexpr static const scalar d2= -3.0/2.0;
        constexpr static const scalar d3= 605.0/250.0;
        constexpr static const scalar d4= 29.0/250.0;
        
    // Rosenbrock

        mutable FastChemistry::LUsolver LU;


public:

    //- Runtime type information
    TypeName("OptRosenbrock34");


    // Constructors

        //- Construct from thermo
        OptRosenbrock34(const fluidReactionThermo& thermo);

    //- Destructor
    virtual ~OptRosenbrock34();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalar& __restrict__ p,
            scalar& __restrict__ T,
            scalarField& Y,
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
            double* __restrict__ err,
            double* __restrict__ dydx,
            double* __restrict__ dfdx,
            double* __restrict__ a_
        )const;


        scalar Rosenbrock34Solve
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
    #include "OptRosenbrock34.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
