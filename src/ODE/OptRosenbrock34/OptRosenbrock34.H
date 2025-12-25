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

/*---------------------------------------------------------------------------*\
  Description
      Header file for Rosenbrock34 algorithm

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

#ifndef OptRosenbrock34_H
#define OptRosenbrock34_H

//=============================================================================//

//---------------------------------
// 1. Standard C++ library headers
//---------------------------------
#include <iostream>
#include <vector>
#include <algorithm>

//---------------------------------
// 2. OpenFOAM library headers
//---------------------------------
#include "chemistrySolver.H"

//---------------------------------
// 3. FastChemistry headers
//---------------------------------
#include "LUsolver.H"

//=============================================================================//

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


        double dxTry;
        bool last;
        stepState(const double dx)
        :
            dxTry(dx),
            last(false)
        {}
    };

        //- Coefficients dictionary
        dictionary coeffsDict_;

    // ODE solver
    
        //- Absolute convergence tolerance per step
        const double absTol_;

        //- Relative convergence tolerance per step
        const double relTol_;

        //- The maximum number of sub-steps allowed for the integration step
        const label maxSteps_;

    // Adaptive Solver

        //- Step-size adjustment controls
        constexpr static const double safeScale_ =     0.9;
        constexpr static const double alphaInc_  =     0.2;
        constexpr static const double alphaDec_ =     0.25;
        constexpr static const double minScale_ =      0.2;
        constexpr static const double maxScale_ =      10;

        //constexpr static const scalar ratio = std::pow(maxScale_/safeScale_, -1.0/alphaInc_);
        constexpr static const double ratio = 5.904900000E-06;

        constexpr static const double safeMaxScale = safeScale_*maxScale_;

    // Rosenbrock

        constexpr static const double a21= 2;
        constexpr static const double a31= 48.0/25.0;
        constexpr static const double a32= 6.0/25.0;

        constexpr static const double c21= -8;
        constexpr static const double c31= 372.0/25.0;
        constexpr static const double c32= 12.0/5.0;

        constexpr static const double c41= -112.0/125.0;
        constexpr static const double c42= -54.0/125.0;
        constexpr static const double c43= -2.0/5.0;

        constexpr static const double b1= 19.0/9.0;
        constexpr static const double b2= 1.0/2.0;
        constexpr static const double b3= 25.0/108.0;
        constexpr static const double b4= 125.0/108.0;

        constexpr static const double e1= 34.0/108.0;
        constexpr static const double e2= 7.0/36.0;
        constexpr static const double e3= 0;
        constexpr static const double e4= 125.0/108.0;

        constexpr static const double gamma= 1.0/2.0;
        constexpr static const double Invgamma= 1.0/gamma;
        constexpr static const double c2= 1;
        constexpr static const double c3= 3.0/5.0;

        constexpr static const double d1= 1.0/2.0;
        constexpr static const double d2= -3.0/2.0;
        constexpr static const double d3= 605.0/250.0;
        constexpr static const double d4= 29.0/250.0;
        
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
            double& __restrict__ p,
            double& __restrict__ T,
            scalarField& Y,
            const label li,
            double& __restrict__ deltaT,
            double& __restrict__ subDeltaT
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
            const double xEnd,
            const label li,
            double& dxTry,
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
            double& __restrict__ x,
            const label li,
            double& __restrict__ dxTry,
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


        double Rosenbrock34Solve
        (
            const double x0,
            const label li,
            const double dx,
            const double invdx,
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
        inline double normaliseError
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
