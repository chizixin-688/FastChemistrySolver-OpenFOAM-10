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

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
  Description
      Inplementation for Rosenbrock34 algorithm

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. Standard C++ library headers
//---------------------------------
#include <iomanip>

//---------------------------------
// 2. OpenFOAM library headers
//---------------------------------
#include "SubField.H"
#include "addToRunTimeSelectionTable.H"

//---------------------------------
// 3. FastChemistry headers
//---------------------------------
#include "OptRosenbrock34.H"

//---------------------------------
// 4. SIMD / AVX2 headers
//---------------------------------
#include <immintrin.h>  

//=============================================================================//


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::OptRosenbrock34<ChemistryModel>::OptRosenbrock34
(
    const fluidReactionThermo& thermo
)
:
    chemistrySolver<ChemistryModel>(thermo),
    coeffsDict_(this->subDict("OptRosenbrock34Coeffs")),
    absTol_(coeffsDict_.lookup<scalar>("absTol")),
    relTol_(coeffsDict_.lookup<scalar>("relTol")),
    maxSteps_(coeffsDict_.lookupOrDefault("maxSteps",10000)),
    LU(this->YTpYTpWork[1],this->n_)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::OptRosenbrock34<ChemistryModel>::~OptRosenbrock34()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::OptRosenbrock34<ChemistryModel>::solve
(
    scalar& __restrict__ p,
    scalar& __restrict__ T,
    scalarField& Phi,
    const label li,
    scalar& __restrict__ deltaT,
    scalar& __restrict__ subDeltaT
) const
{}

template<class ChemistryModel>
void Foam::OptRosenbrock34<ChemistryModel>::solve
(
    const label li,
    double p,
    double& __restrict__ deltaT,
    double& __restrict__ subDeltaT
) const
{
    double* __restrict__ Phi0     = this->YTpWork[1];
    double* __restrict__ PhiTemp  = this->YTpWork[2];
    double* __restrict__ k1     = this->YTpWork[3];
    double* __restrict__ k2     = this->YTpWork[4];
    double* __restrict__ k3     = this->YTpWork[5];
    double* __restrict__ k4     = this->YTpWork[6];
    double* __restrict__ err    = this->YTpWork[7];
    double* __restrict__ dydx   = this->YTpWork[8];
    double* __restrict__ dfdx   = this->YTpWork[9];
    double* __restrict__ Jac     = this->YTpYTpWork[1];

    this->ODESolve
    (
        deltaT,         // end( CFD flow time step)
        li,             // celli
        subDeltaT,      // sub time step
        p,
        Phi0,           // initial vector
        PhiTemp,        // predicted vector
        k1,
        k2,
        k3,
        k4,
        err,
        dydx,
        dfdx,
        Jac    
    );
}

template<class ChemistryModel>
void Foam::OptRosenbrock34<ChemistryModel>::ODESolve
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
    double* __restrict__ Jac
) const
{

    stepState step(dxTry);
    scalar x = 0;

    for (label nStep=0; nStep<maxSteps_; nStep++)
    {
        // Store previous iteration dxTry
        double dxTry0 = step.dxTry;

        // Check if this is a truncated step and set dxTry to integrate to xEnd
        if ((x + step.dxTry - xEnd)*(x + step.dxTry) > 0)
        {
            step.last = true;
            step.dxTry = xEnd - x;
        }

        // Integrate as far as possible up to step.dxTry
        //solve(x, y, li, step);


        adaptiveSolve
        (
            x,              //start
            li,             //celli
            step.dxTry,     //sub time step
            p,
            Phi0,
            PhiTemp,
            k1,
            k2,
            k3,
            k4,
            err,
            dydx,
            dfdx,
            Jac
        );

        // Check if reached xEnd
        if ((x - xEnd)*(xEnd) >= 0)
        {
            if (nStep > 0 && step.last)
            {
                step.dxTry = dxTry0;
            }

            dxTry = step.dxTry;
            return;
        }
    }
    FatalErrorInFunction
        << "Integration steps greater than maximum " << maxSteps_ << nl
        << exit(FatalError);

}

template<class ChemistryModel>
void Foam::OptRosenbrock34<ChemistryModel>::adaptiveSolve
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
    double* __restrict__ Jac
) const
{

    double dx = dxTry;
    double Err = 0.0;
    double invdx = 1.0/dx;

    // Loop over solver and adjust step-size as necessary
    // to achieve desired error
        Err = Rosenbrock34Solve
        (
            x,
            li,
            dx,
            invdx,
            p,
            Phi0,
            PhiTemp,
            k1,
            k2,
            k3,
            k4,
            err,
            dydx,
            dfdx,
            Jac
        );

    while(Err > 1)
    {

        scalar scale = max(safeScale_*pow(Err, -alphaDec_), minScale_);
        dx *= scale;
        invdx = 1.0/dx;
        Err = Rosenbrock34Solve
        (
            x,
            li,
            dx,
            invdx,
            p,
            Phi0,
            PhiTemp,
            k1,
            k2,
            k3,
            k4,
            err,
            dydx,
            dfdx,
            Jac
        );
    }

    // Update the state
    x += dx;
    for(label i = 0; i < this->n_;i++)
    {
        Phi0[i] = PhiTemp[i];
    }

    dxTry = (Err>ratio)?
            min(max(safeScale_*pow(Err, -alphaInc_), minScale_), maxScale_)*dx
            :safeMaxScale*dx;
}

template<class ChemistryModel>
Foam::scalar Foam::OptRosenbrock34<ChemistryModel>::Rosenbrock34Solve
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
    double* __restrict__ Jac
) const
{
    
    this->jacobian(x0, li, p, Phi0,  dfdx, Jac);
    
    {
        const unsigned int NN = this->alignN*this->n_;
        unsigned int  remain = NN%16;
        for(unsigned int  i = 0 ; i<NN-remain;i=i+16)
        {
            __m256d Av0 = load256d(&Jac[i+0]);
            __m256d Av1 = load256d(&Jac[i+4]);
            __m256d Av2 = load256d(&Jac[i+8]);
            __m256d Av3 = load256d(&Jac[i+12]);

            store256d(&Jac[i+0],-Av0);
            store256d(&Jac[i+4],-Av1);
            store256d(&Jac[i+8],-Av2);
            store256d(&Jac[i+12],-Av3);
        }
        for(unsigned int  i = NN-remain; i < NN;i++)
        {
            Jac[i] = -Jac[i];
        }
    }
    for (  label i=0; i<this->n_; i++)
    {
        Jac[i*this->alignN+i] += invdx*Invgamma;
    }

    LU.Block4LUDecompose();
    
    //const unsigned int remain = n_%4;

    // compute k1
    {
        const double dxd1 = dx*d1+1;
        __m256d dxd1v = _mm256_set1_pd(dxd1);
        for(unsigned int i = 0 ;i < this->alignN;i=i+4)
        {
            __m256d dfdxv = load256d(&dfdx[i+0]);
            __m256d k1v = _mm256_mul_pd(dxd1v,dfdxv);
            store256d(&k1[i+0],k1v);
        }
    }

    LU.xSolve(k1);
    
    {
        __m256d a21v = _mm256_set1_pd(a21);
        for(unsigned int i = 0 ;i < this->alignN;i=i+4)
        {
            __m256d y0v = load256d(&Phi0[i]);
            __m256d k1v = load256d(&k1[i]);
            __m256d yTempv = _mm256_fmadd_pd(a21v,k1v,y0v);
            store256d(&PhiTemp[i],yTempv);
        }
    }



    this->derivatives(x0 + c2*dx, li, p, PhiTemp, dydx, k3, k4);

    {
        double dxd2 = dx*d2;
        double c21invdx = c21*invdx;
        __m256d dxd2v = _mm256_set1_pd(dxd2);
        __m256d c21invdxv = _mm256_set1_pd(c21invdx);

        for(unsigned int i = 0 ;i < this->alignN;i=i+4)
        {
            k2[i] = dydx[i] + dxd2*dfdx[i] + c21invdx*k1[i];
            __m256d dydxv = load256d(&dydx[i]);
            __m256d dfdxv = load256d(&dfdx[i]);
            __m256d k1v = load256d(&k1[i]);
            __m256d k2v = _mm256_fmadd_pd(dxd2v,dfdxv,dydxv);
            k2v = _mm256_fmadd_pd(c21invdxv,k1v,k2v);
            store256d(&k2[i],k2v);
        } 
    }

    LU.xSolve(k2);

    {
        __m256d a31v = _mm256_set1_pd(a31);
        __m256d a32v = _mm256_set1_pd(a32);
        for(unsigned int i = 0 ;i < this->alignN;i=i+4)
        {
            __m256d k1v = load256d(&k1[i]);
            __m256d k2v = load256d(&k2[i]);
            __m256d y0v = load256d(&Phi0[i]);
            __m256d yTempv = _mm256_fmadd_pd(a31v,k1v,y0v);
            yTempv = _mm256_fmadd_pd(a32v,k2v,yTempv);
            store256d(&PhiTemp[i],yTempv);
        }   
    }

    this->derivatives(x0 + c2*dx, li, p, PhiTemp, dydx, k3, k4);
    
    {
        double dxd3 = dx*d3;
        double c31invdx = c31*invdx;
        double c32invdx = c32*invdx;
        __m256d dxd3v = _mm256_set1_pd(dxd3);
        __m256d c31invdxv = _mm256_set1_pd(c31invdx);
        __m256d c32invdxv = _mm256_set1_pd(c32invdx);
        for(unsigned int i = 0 ;i < this->alignN;i=i+4)
        {
            __m256d dydxv = load256d(&dydx[i]);
            __m256d dfdxv = load256d(&dfdx[i]);
            __m256d k1v = load256d(&k1[i]);
            __m256d k2v = load256d(&k2[i]);
            __m256d k3v = _mm256_fmadd_pd(dxd3v,dfdxv,dydxv);
            k3v = _mm256_fmadd_pd(c31invdxv,k1v,k3v);
            k3v = _mm256_fmadd_pd(c32invdxv,k2v,k3v);
            store256d(&k3[i],k3v);
        }          
    }



    LU.xSolve(k3);
    


    {
        double dxd4 = dx*d4;
        double c41invdx = c41*invdx;
        double c42invdx = c42*invdx;
        double c43invdx = c43*invdx;
        __m256d dxd4v = _mm256_set1_pd(dxd4);
        __m256d c41invdxv = _mm256_set1_pd(c41invdx);
        __m256d c42invdxv = _mm256_set1_pd(c42invdx);
        __m256d c43invdxv = _mm256_set1_pd(c43invdx);
        for(unsigned int i = 0 ;i < this->alignN;i=i+4)
        {
            __m256d dydxv = load256d(&dydx[i]);
            __m256d dfdxv = load256d(&dfdx[i]);
            __m256d k1v = load256d(&k1[i]);
            __m256d k2v = load256d(&k2[i]);
            __m256d k3v = load256d(&k3[i]);
            __m256d k4v = _mm256_fmadd_pd(dxd4v,dfdxv,dydxv);

            k4v = _mm256_fmadd_pd(c41invdxv,k1v,k4v);
            k4v = _mm256_fmadd_pd(c42invdxv,k2v,k4v);
            k4v = _mm256_fmadd_pd(c43invdxv,k3v,k4v);
            store256d(&k4[i],k4v);
        }      
    }


    LU.xSolve(k4);
    

    {
        __m256d b1v = _mm256_set1_pd(b1);
        __m256d b2v = _mm256_set1_pd(b2);
        __m256d b3v = _mm256_set1_pd(b3);
        __m256d b4v = _mm256_set1_pd(b4);
        __m256d e1v = _mm256_set1_pd(e1);
        __m256d e2v = _mm256_set1_pd(e2);
        __m256d e4v = _mm256_set1_pd(e4);        
        for(unsigned int i = 0 ;i < this->alignN;i=i+4)
        {
            __m256d y0v = load256d(&Phi0[i]);
            __m256d k1v = load256d(&k1[i]);
            __m256d k2v = load256d(&k2[i]);
            __m256d k3v = load256d(&k3[i]);
            __m256d k4v = load256d(&k4[i]);

            __m256d yTempv = _mm256_fmadd_pd(b1v,k1v,y0v);
            yTempv = _mm256_fmadd_pd(b2v,k2v,yTempv);
            yTempv = _mm256_fmadd_pd(b3v,k3v,yTempv);
            yTempv = _mm256_fmadd_pd(b4v,k4v,yTempv);

            __m256d errv = _mm256_mul_pd(e1v,k1v);
            errv = _mm256_fmadd_pd(e2v,k2v,errv);
            errv = _mm256_fmadd_pd(e4v,k4v,errv);

            store256d(&PhiTemp[i],yTempv);
            store256d(&err[i],errv);
        }         
    }

    double maxErr = 0;
    {
        __m256d maxErrv = _mm256_setzero_pd();
        __m256d absTolv = _mm256_set1_pd(absTol_);
        __m256d relTolv = _mm256_set1_pd(relTol_);
        for(unsigned int i = 0 ;i < this->alignN;i=i+4)
        {
            __m256d y0v = load256d(&Phi0[i]);
            __m256d yTempv = load256d(&PhiTemp[i]);
            __m256d errv = load256d(&err[i]);
            __m256d signMask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000));
            __m256d abserrv = _mm256_andnot_pd(signMask, errv);
            __m256d absy0v = _mm256_andnot_pd(signMask, y0v);
            __m256d absyTempv  = _mm256_andnot_pd(signMask, yTempv);

            __m256d maxY = _mm256_max_pd(absyTempv,absy0v);
            __m256d tolv = _mm256_fmadd_pd(relTolv,maxY,absTolv);

            __m256d errByTolv = _mm256_div_pd(abserrv,tolv);
            maxErrv = _mm256_max_pd(maxErrv,errByTolv);
        }

        __m128d low128  = _mm256_castpd256_pd128(maxErrv);           
        __m128d high128 = _mm256_extractf128_pd(maxErrv, 1);          
        __m128d max128  = _mm_max_pd(low128, high128);            

        __m128d shuffled = _mm_shuffle_pd(max128, max128, 0x1);   
        __m128d finalMax = _mm_max_sd(max128, shuffled);          

        maxErr = _mm_cvtsd_f64(finalMax);    
    }


    return maxErr;
}

template<class ChemistryModel>
inline Foam::scalar Foam::OptRosenbrock34<ChemistryModel>::normaliseError
(
    const double* __restrict__ Phi0,
    const double* __restrict__ Phi,
    const double* __restrict__ err
) const
{
    // Calculate the maximum error
    scalar maxErr = 0.0;
    for(unsigned int i = 0 ;i < this->alignN;i++)
    {
        scalar tol = absTol_ + relTol_*std::max(std::fabs(Phi0[i]), std::fabs(Phi[i]));
        maxErr = std::max(maxErr, std::fabs(err[i])/tol);
    }
    return maxErr;
}
// ************************************************************************* //
