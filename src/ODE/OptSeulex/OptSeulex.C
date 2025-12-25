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
      Note that this implementation is not consistent with OpenFOAM-10 version. 
      But same with Imren's version

      Reference

      [1] A. Imren, D. C. Haworth. On the merits of extrapolation-based stiff 
      ODE solvers for combustion CFD. Combustion and Flame 174 (2016) 1–15.


  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

//=============================================================================//

//---------------------------------
// 1. OpenFOAM library headers
//---------------------------------
#include "SubField.H"
#include "addToRunTimeSelectionTable.H"

//---------------------------------
// 2. FastChemistry headers
//---------------------------------
#include "OptSeulex.H"

//=============================================================================//

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::OptSeulex<ChemistryModel>::OptSeulex
(
    const fluidReactionThermo& thermo
)
:
    chemistrySolver<ChemistryModel>(thermo),
    coeffsDict_(this->subDict("OptSeulexCoeffs")),
    absTol_(coeffsDict_.lookup<scalar>("absTol")),
    relTol_(coeffsDict_.lookup<scalar>("relTol")),
    maxSteps_(coeffsDict_.lookupOrDefault("maxSteps",10000)),
    jacRedo_(min(1e-4, relTol_)),
    nSeq_(iMaxx_),
    cpu_(iMaxx_),
    invCpu_(iMaxx_),
    coeff_(iMaxx_, iMaxx_),
    theta_(2*jacRedo_),
    table_(kMaxx_,this->n_),
    pivotIndices_(this->n_),
    dxOpt_(iMaxx_),
    temp_(iMaxx_),
    y0_(nullptr),
    ySequence_(nullptr),
    scale_(nullptr),
    LU(this->YTpYTpWork[2],this->n_)
{

    if (posix_memalign(reinterpret_cast<void**>(&this->y0_), 32, this->alignN*sizeof(double)))
    {
        throw std::bad_alloc();
    }
    std::memset(this->y0_, 0, this->alignN);

    if (posix_memalign(reinterpret_cast<void**>(&this->ySequence_), 32, this->alignN*sizeof(double)))
    {
        throw std::bad_alloc();
    }
    std::memset(this->ySequence_, 0, this->alignN);

    if (posix_memalign(reinterpret_cast<void**>(&this->scale_), 32, this->alignN*sizeof(double)))
    {
        throw std::bad_alloc();
    }
    std::memset(this->scale_, 0, this->alignN);

    // The CPU time factors for the major parts of the algorithm
    const scalar cpuFunc = 1, cpuJac = 5, cpuLU = 1, cpuSolve = 1;

    nSeq_[0] = 2;
    nSeq_[1] = 3;

    for (int i=2; i<iMaxx_; i++)
    {
        nSeq_[i] = 2*nSeq_[i-2];
    }
    cpu_[0] = cpuJac + cpuLU + nSeq_[0]*(cpuFunc + cpuSolve);
    invCpu_[0] = 1.0/cpu_[0];

    for (int k=0; k<kMaxx_; k++)
    {
        cpu_[k+1] = cpu_[k] + (nSeq_[k+1]-1)*(cpuFunc + cpuSolve) + cpuLU;
        invCpu_[k+1] = 1.0/cpu_[k+1];
    }

    // Set the extrapolation coefficients array
    for (int k=0; k<iMaxx_; k++)
    {
        for (int l=0; l<k; l++)
        {
            scalar r = scalar(nSeq_[k])/nSeq_[l];
            coeff_(k, l) = 1/(r - 1);
        }
    }
    this->logTol = -log10(relTol_ + absTol_)*0.6 + 0.5;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::OptSeulex<ChemistryModel>::~OptSeulex()
{
    free(this->y0_);
    free(this->ySequence_);
    free(this->scale_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::OptSeulex<ChemistryModel>::solve
(
    scalar& __restrict__ p,
    scalar& __restrict__ T,
    scalarField& y,
    const label li,
    scalar& __restrict__ deltaT,
    scalar& __restrict__ subDeltaT
) const
{}

template<class ChemistryModel>
void Foam::OptSeulex<ChemistryModel>::solve
(
    const label li,
    double p,
    double& __restrict__ deltaT,
    double& __restrict__ subDeltaT
) const
{

    double* __restrict__ Phi0         = this->YTpWork[1];
    double* __restrict__ PhiTemp    = this->YTpWork[2];
    double* __restrict__ Cp         = this->YTpWork[3];
    double* __restrict__ Ha         = this->YTpWork[4];
    double* __restrict__ dy_       = this->YTpWork[6];
    double* __restrict__ dydx_     = this->YTpWork[7];
    double* __restrict__ k9         = this->YTpWork[9];
    double* __restrict__ k11         = this->YTpWork[11];
    double* __restrict__ Jac = this->YTpYTpWork[1];
    double* __restrict__ a = this->YTpYTpWork[2];  
     
    this->ODESolve
    (
        deltaT,
        li,
        subDeltaT,
        p,
        Phi0,
        PhiTemp,
        Cp,
        Ha,
        dy_,
        dydx_,
        k9,
        k11,
        Jac,
        a        
    );
}

template<class ChemistryModel>
void Foam::OptSeulex<ChemistryModel>::ODESolve
(
    const scalar xEnd,
    const label li,
    scalar& dxTry,
    const double p,
    double* __restrict__ Phi0,
    double* __restrict__ PhiTemp_,
    double* __restrict__ Cp,
    double* __restrict__ Ha,
    double* __restrict__ dy_,
    double* __restrict__ dydx_,
    double* __restrict__ k9,
    double* __restrict__ k11,
    double* __restrict__ Jy,
    double* __restrict__ a
) const
{
    stepState step(dxTry);
    scalar x = 0;

    for (label nStep=0; nStep<maxSteps_; nStep++)
    {
        // Store previous iteration dxTry
        scalar dxTry0 = step.dxTry;

        step.reject = false;

        // Check if this is a truncated step and set dxTry to integrate to xEnd
        if ((x + step.dxTry - xEnd)*(x + step.dxTry) > 0)
        {
            step.last = true;
            step.dxTry = xEnd - x;
        }

        // Integrate as far as possible up to step.dxTry
        SeulexSolve
        (
            x, 
            li, 
            step,
            p,
            Phi0,
            PhiTemp_,
            Cp,
            Ha,
            dy_,
            dydx_,
            k9,
            k11,
            Jy,
            a
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

        step.first = false;

        // If the step.dxTry was reject set step.prevReject
        if (step.reject)
        {
            step.prevReject = true;
        }
    }

    FatalErrorInFunction
        << "Integration steps greater than maximum " << maxSteps_ << nl
        << ", xEnd = " << xEnd
        << ", x = " << x << ", dxDid = " << step.dxDid << nl
        << "    y = " << Phi0
        << exit(FatalError);

}

template<class ChemistryModel>
void Foam::OptSeulex<ChemistryModel>::SeulexSolve
(
    scalar& x,
    const label li,
    stepState& step,
    const double p,
    double* __restrict__ Phi,
    double* __restrict__ PhiTemp_,
    double* __restrict__ Cp,
    double* __restrict__ Ha,
    double* __restrict__ dy_,
    double* __restrict__ dydx_,
    double* __restrict__ k9,
    double* __restrict__ k11,    
    double* __restrict__ Jy,
    double* __restrict__ a  
) const
{

    temp_[0] = great;
    scalar dx = step.dxTry;
    for(label i = 0; i < this->n_;i++)
    {
        y0_[i] = Phi[i];
    }
    dxOpt_[0] = mag(0.1*dx);

    if (step.first || step.prevReject)
    {
        theta_ = 2*jacRedo_;
    }

    if (step.first)
    {
        //scalar logTol = -log10(relTol_ + absTol_)*0.6 + 0.5;
        kTarg_ = max(1, min(kMaxx_ - 1, int(this->logTol)));
    }


    //forAll(scale_, i)
    //{
    //    scale_[i] = absTol_ + relTol_*mag(y[i]);
    //}
    for(label i = 0; i<this->n_;i++)
    {
        k11[i] = absTol_ + relTol_*mag(Phi[i]);
    }

    bool jacUpdated = false;

    if (theta_ > jacRedo_)
    {
        this->jacobian(x, li, p, Phi,  k9, Jy);
        jacUpdated = true;
    }

    int k;
    scalar dxNew = mag(dx);
    bool firstk = true;

    while (firstk || step.reject)
    {
        dx = step.forward ? dxNew : -dxNew;
        firstk = false;
        step.reject = false;

        scalar errOld = 0;

        for (k=0; k<=kTarg_+1; k++)
        {
            bool success = seul
            (
                x, 
                y0_, 
                li, 
                dx, 
                k,
                p,
                ySequence_, 
                PhiTemp_,
                Cp,
                Ha,
                dy_,
                dydx_,
                k11,
                Jy,
                a    
            );

            if (!success)
            {
                step.reject = true;
                dxNew = mag(dx)*stepFactor5_;
                break;
            }

            if (k == 0)
            {
                for(label i = 0; i < this->n_;i++)
                {
                    Phi[i] = ySequence_[i];
                }

            }
            else
            {
                
                for(label i = 0; i < this->n_;i++)                
                {
                    table_[k-1][i] = ySequence_[i];
                }
            }

            if (k != 0)
            {
                extrapolate(k, table_, Phi);
                scalar err = 0;
                //forAll(scale_, i) 
                //{
                //    scale_[i] = absTol_ + relTol_*mag(y0_[i]);
                //    err += sqr((y[i] - table_(0, i))/scale_[i]);
                //}
                for(label i = 0; i<this->n_;i++)
                {
                    k11[i] = absTol_ + relTol_*mag(y0_[i]);
                    err += sqr((Phi[i] - table_(0, i))/k11[i]);                    
                }
                err = sqrt(err/this->n_);
                if (err > 1/small || (k > 1 && err >= errOld))
                {
                    step.reject = true;
                    dxNew = mag(dx)*stepFactor5_;
                    break;
                }
                errOld = min(4*err, 1);
                scalar expo = 1.0/(k + 1);
                scalar facmin = pow(stepFactor3_, expo);
                scalar fac;
                if (err == 0)
                {
                    fac = 1.0/facmin;
                }
                else
                {
                    fac = stepFactor2_/pow(err/stepFactor1_, expo);
                    fac = max(facmin/stepFactor4_, min(1/facmin, fac));
                }
                dxOpt_[k] = mag(dx*fac);
                temp_[k] = cpu_[k]/dxOpt_[k];

                if ((step.first || step.last) && err <= 1.0)
                {
                    break;
                }
                if
                (
                    k == kTarg_ - 1
                 && !step.prevReject
                 && !step.first && !step.last
                )
                {
                    if (err <= 1.0)
                    {
                        break;
                    }
                    else if (err > nSeq_[kTarg_]*nSeq_[kTarg_ + 1]*4.0)
                    {
                        step.reject = true;
                        kTarg_ = k;
                        if (kTarg_>1 && temp_[k-1] < kFactor1_*temp_[k])
                        {
                            kTarg_--;
                        }
                        dxNew = dxOpt_[kTarg_];
                        break;
                    }
                }
                if (k == kTarg_)
                {
                    if (err <= 1.0)
                    {
                        break;
                    }
                    else if (err > nSeq_[k + 1]*2.0)
                    {
                        step.reject = true;
                        if (kTarg_>1 && temp_[k-1] < kFactor1_*temp_[k])
                        {
                            kTarg_--;
                        }
                        dxNew = dxOpt_[kTarg_];
                        break;
                    }
                }
                if (k == kTarg_+1)
                {
                    if (err > 1.0)
                    {
                        step.reject = true;
                        if
                        (
                            kTarg_ > 1
                         && temp_[kTarg_-1] < kFactor1_*temp_[kTarg_]
                        )
                        {
                            kTarg_--;
                        }
                        dxNew = dxOpt_[kTarg_];
                    }
                    break;
                }
            }
        }
        if (step.reject)
        {
            step.prevReject = true;
            if (!jacUpdated)
            {
                theta_ = 2.0*jacRedo_;
                if (theta_ > jacRedo_ && !jacUpdated)
                {
                    this->jacobian(x, li, p, Phi,  k9, Jy);
                    jacUpdated = true;
                }
            }
        }
    }

    jacUpdated = false;

    step.dxDid = dx;
    x += dx;

    label kopt;
    if (k == 1)
    {
        kopt = 2;
    }
    else if (k <= kTarg_)
    {
        kopt=k;
        if (temp_[k-1] < kFactor1_*temp_[k])
        {
            kopt = k - 1;
        }
        else if (temp_[k] < kFactor2_*temp_[k - 1])
        {
            kopt = min(k + 1, kMaxx_ - 1);
        }
    } 
    else
    {
        kopt = k - 1;
        if (k > 2 && temp_[k-2] < kFactor1_*temp_[k - 1])
        {
            kopt = k - 2;
        }
        if (temp_[k] < kFactor2_*temp_[kopt])
        {
            kopt = min(k, kMaxx_ - 1);
        }
    }

    if (step.prevReject)
    {
        kTarg_ = min(kopt, k);
        dxNew = min(mag(dx), dxOpt_[kTarg_]);
        step.prevReject = false;
    }
    else
    {
        if (kopt <= k)
        {
            dxNew = dxOpt_[kopt];
        }
        else
        {
            if (k < kTarg_ && temp_[k] < kFactor2_*temp_[k - 1])
            {

                dxNew = dxOpt_[k]*cpu_[kopt + 1]/cpu_[k];
            }
            else
            {
                dxNew = dxOpt_[k]*cpu_[kopt]/cpu_[k];
            }
        }
        kTarg_ = kopt;
    }

    step.dxTry = step.forward ? dxNew : -dxNew;

}

template<class ChemistryModel>
bool Foam::OptSeulex<ChemistryModel>::seul
(
    const scalar x0,
    double* __restrict__ y0,
    const label li,
    const scalar dxTot,
    const label k,
    const double p,
    double* __restrict__ ySequence,
    double* __restrict__ yTemp_,
    double* __restrict__ Cp,
    double* __restrict__ Ha,
    double* __restrict__ dy_,
    double* __restrict__ dydx_,
    double* __restrict__ k11,     
    double* __restrict__ Jac,
    double* __restrict__ a        
) const
{

    unsigned int nSteps = nSeq_[k];
    scalar dx = dxTot/nSteps;
    scalar invdx = nSteps/dxTot;

    {
        const unsigned int  NN = this->alignN*this->n_;
        unsigned int remain = NN%16;
        for(unsigned int i = 0 ; i<NN-remain;i=i+16)
        {
            __m256d Av0 = _mm256_loadu_pd(&Jac[i+0]);
            __m256d Av1 = _mm256_loadu_pd(&Jac[i+4]);
            __m256d Av2 = _mm256_loadu_pd(&Jac[i+8]);
            __m256d Av3 = _mm256_loadu_pd(&Jac[i+12]);

            _mm256_storeu_pd(&a[i+0],-Av0);
            _mm256_storeu_pd(&a[i+4],-Av1);
            _mm256_storeu_pd(&a[i+8],-Av2);
            _mm256_storeu_pd(&a[i+12],-Av3);
        }
        for(unsigned int i = NN-remain; i < NN;i++)
        {
            a[i] = -Jac[i];
        }
    }
    for (label i=0; i<this->n_; i++)
    {
        a[i*this->alignN+i] += invdx;
    }

    {
        LU.Block4LUDecompose();
    }



    scalar xnew = x0 + dx;

    this->derivatives(xnew, li, p, y0, dy_, Cp, Ha);

    {
        LU.xSolve(dy_);
    }

    for(label i = 0; i < this->n_;i++) 
    {
        yTemp_[i] = y0[i];
    }    
    for (unsigned int nn=1; nn<nSteps; nn++)
    {
        for(label i = 0; i < this->n_;i++) 
        {
            yTemp_[i] += dy_[i];
        }


        xnew += dx;

        if (nn == 1 && k<=1)
        {
            scalar dy1 = 0;
            for (label i=0; i<this->n_; i++) 
            {
                dy1 += sqr(dy_[i]/k11[i]);
            }
            dy1 = sqrt(dy1);

            //odes_.derivatives(x0 + dx, yTemp_, li, dydx_);
            this->derivatives(x0 + dx, li, p, yTemp_, dydx_, Cp, Ha);

            for (label i=0; i<this->n_; i++) 
            {
                dy_[i] = dydx_[i] - dy_[i]*invdx;
            }

            
            LU.xSolve(dy_);
            
            

            // This form from the original paper is unreliable
            // step size underflow for some cases
            // const scalar denom = max(1, dy1);

            // This form is reliable but limits how large the step size can be
            //const scalar denom = min(1, dy1 + small);

            scalar dy2 = 0.0;

	        for (label i=0; i<this->n_; i++)
	        {
	            if (mag(dy_[i]) > sqrt(vGreat) )
                { 
                    dy2 = vGreat; 
                    break;
                }
		        dy2 += sqr(dy_[i]/k11[i]);
	        }

	    	dy2 = sqrt(dy2);
		    theta_ = dy2/max(1.0, dy1 + small);
		    if (theta_ > 1.0 ) 
	    	{
		        return false;
	    	}
        }


        this->derivatives(xnew, li, p, yTemp_, dy_, Cp, Ha);

        {
            LU.xSolve(dy_);
        }
    }

    for (label i=0; i<this->n_; i++) 
    {
        ySequence[i] = yTemp_[i] + dy_[i];
    }

    return true;

}



template<class ChemistryModel>
void Foam::OptSeulex<ChemistryModel>::extrapolate
(
    const label k,
    scalarRectangularMatrix& table,
    double* y
) const
{
    for (int j=k-1; j>0; j--)
    {
        for (label i=0; i<this->n_; i++)
        {
            table[j-1][i] =
                table(j, i) + coeff_(k, j)*(table(j, i) - table[j-1][i]);
        }
    }

    for (label i=0; i<this->n_; i++) 
    {
        y[i] = table(0, i) + coeff_(k, 0)*(table(0, i) - y[i]);
    }
}


// ************************************************************************* //
