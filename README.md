# Chemistry Model 

A fast chemistry solver for OpenFOAM-10 chemistry computation

# Key Advantages

- **Mathematically consistent with OpenFOAM-10’s original formulation** —  

- **No third-party dependencies** —  

- **Significant performance improvement** —  

# Limitations

This optimized combustion solver is still under development.  
Some functions and reaction types may not work as expected:  
1. Doesn't support tabulation method or reduction method to  
 further accelerating computation.  
2. Only support ideal gas  
3. Support the following reaction type:
    irreversibleArrhenius  
    irreversibleArrheniusLindemannFallOff  
    irreversibleArrheniusSRIFallOff  
    irreversibleArrheniusTroeFallOff  
    irreversibleThirdBodyArrhenius  
    nonEquilibriumReversibleArrhenius  
    nonEquilibriumReversibleThirdBodyArrhenius  
    reversibleArrhenius  
    reversibleArrheniusLindemannFallOff  
    reversibleArrheniusSRIFallOff  
    reversibleArrheniusTroeFallOff  
    reversibleThirdBodyArrhenius  
4. Unsupport the following reaction type:
    irreversibleArrheniusLindemannChemicallyActivated  
    irreversibleArrheniusSRIChemicallyActivated  
    irreversibleArrheniusTroeChemicallyActivated  
    reversibleArrheniusLindemannChemicallyActivated  
    reversibleArrheniusSRIChemicallyActivated  
    reversibleArrheniusTroeChemicallyActivated  
    irreversibleJanev  
    irreversibleLandauTeller  
    irreversibleLangmuirHinshelwood  
    irreversiblePowerSeries  
    nonEquilibriumReversibleLandauTeller  
    reversibleJanev  
    reversibleLandauTeller  
    reversibleLangmuirHinshelwood  
    reversiblePowerSeries  
    irreversibleFluxLimitedLangmuirHinshelwood  
    irreversibleSurfaceArrhenius  

## Numerical considerations
This code includes the analytical derivatives for Falloff and 
ChemicallyActivated reactions, which are not considered by default during 
compilation in OpenFOAM-10. For reference, see the 
ChemicallyActivatedReactionRateI.H and FallOffReactionRateI.H files in OpenFOAM

# Prerequisites

1. cmake >= 3.10.2  
2. gcc >= 5.4.0  
3. OpenFOAM-10  

# Compilation

**Option 1: using wmake**

    If your glibc version is ≥ 2.22, you can build with the OpenFOAM wmake 
    command, since glibc-2.22 provides vectorized mathematical functions. 
    The installation steps are:

     1. Enter the `src` folder.

     2. In `Macro.H`, set `USE_LOCALFILE_  false`. 

     3. Source your OpenFOAM environment.  

     4. Run the `wmake`

**Option 2: using Cmake and make**

     1. Enter the `src` folder.

     2. Source your OpenFOAM environment.

     3. Create a `build` directory. 
        `mkdir build`

     4. Enter the `build` folder. 
        `cd build`

     5. Run CMake.
        `cmake ..`

     6. Build with Make.
        `make`  

     7. Copy the dynamic library `libFastChemistryModel.so` to the folder 
     `$FOAM_USER_LIBBIN`

     *  Note: When building on a compute cluster, CMake may pick the system’s 
        GCC by default, which can be too old (e.g., GCC 4.8.5 lacks C++14 
        support) and cause compilation failures. In that case, load a newer GCC 
        and run:
    `cmake -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++) ..`

**Vectorized math function**

    If your glibc version is **< 2.22**, the libmvec.so lacks and the vector 
    math function exp, pow or log cannot be used. In this situation, you can 
    set the macro `USE_LOCALFILE_ true` in Macro.H. The module will use 
    assembly code of vectorized math function of glibc 2.22 in asmCode folder


# Usage

**controlDict**

    libs
    (
        "libFastChemistryModel.so"
    );


**chemistryProperties**

    chemistryType
    {
        solver          OptRodas34;
        //solver          OptRosenbrock34;
        //solver          OptSeulex;
        method          FastChemistryModel;
    }


    // The option is consistent with OpenFOAM-10. The Jacobian matrix of OF-10 
    // is based on the mass fraction (Jy). The matrix dCdY is used to convert 
    // the Jc to Jy. If fast option is used, the non-diagonal element dCdY is 
    // zero, and the converting will be fast. However, this may influence the 
    // converging rate, especially the torelance is tight. The convergence rate 
    // can be evaluated through 0D ignition problems
    jacobian    exact;
    //jacobian    fast;

    // switch on the chemistry
    chemistry       on;
    //chemistry       off;

    // Solving chemistry for cells with temperature larger than Treact.
    Treact          0;
    //Treact          300;

    // Initial sub step of ODE system.
    initialChemicalTimeStep 1e-07;

    // Turn on the load balancing when parallel computing is performed.
    balance         on;

    // After a specified number of iterations, the DLB algorithm re evaluates 
    // the load. Default value is 1.
    Iter            1;

    // Only chemistry integration time large than Tave*DLBthreshold is 
    // identified as a high load process Tave is the average time of chemical 
    // reaction integration for all processes.
    DLBthreshold    1.0;


    OptRodas34Coeffs
    {
        absTol          1e-8;
        relTol          1e-4;
    }
    OptSeulexCoeffs
    {
        absTol          1e-8;
        relTol          1e-4;
    }

# PLOG reaction

    This chemistry solver supports Plog reaction, for more details about Plog 
    Reaction in OpenFOAM, see

    https://github.com/ZmengXu/PLOGArrheniusReactions
    https://github.com/yuchenzh/plogOF10


    You need to convert the CHEMKIN format to OF format manually or using other 
    utility. The following variable needs to be modified.
    1. pressure p
    2. prefactor A
    3. Ea->Ta
    The conversion of A should be careful. For reaction with one reactant. A is
    unchanged. For reaction with two reactant. A=A*1e-3. For reaction with 
    three reactant. A=A*1e-6

**CHEMKIN**

    HNO=H+NO        .18259e+21  -3.008  47880.0 
    PLOG /  0.100   .20121e+20  -3.021  47792.0   /
    PLOG /  1.000   .18259e+21  -3.008  47880.0   /
    PLOG /  10.00   .12762e+22  -2.959  48100.0   /
    PLOG /  100.0   .56445e+22  -2.855  48459.0   /
    PLOG /  1000.   .97111e+22  -2.642  48940.0   / 

    H+HNO2=NO+H2O   .3380E+10   1.070   5565.0 
    PLOG /  0.010   .3390E+10   1.070   5568.0    /
    PLOG /  0.100   .3390E+10   1.070   5567.0    /
    PLOG /  0.316   .3390E+10   1.070   5567.0    /
    PLOG /  1.000   .3380E+10   1.070   5565.0    /
    PLOG /  3.160   .3380E+10   1.070   5560.0    /
    PLOG /  10.00   .3400E+10   1.070   5546.0    /
    PLOG /  31.60   .4320E+10   1.040   5591.0    /
    PLOG /  100.0   .1270E+11   0.910   5968.0    /

**OpenFOAM**

    un-named-reaction-176
    {
        type            reversibleArrheniusPLOG;
        reaction        "HNO = H + NO";
        A               1.8259e+20;
        beta            -3.008;
        Ta              24092.79652; 
        ArrheniusData
        (
            // PLOG /p A beta Ta/
            // convert atm to Pa
            // convert Ea to Ta, Ta = Ea/R (8.314 J⋅K−1⋅mol−1), cal to J 
            // (4.18400 J/cal)
            (1E4  2.0121e+19   -3.021   24051.2     )   
            (1E5  1.8259e+20   -3.008   24095.4     )   
            (1E6  1.2762e+21   -2.959   24206.2     )   
            (1E7  5.6445e+21   -2.855   24386.9     )   
            (1E8  9.7111e+21   -2.642   24628.9     )   
        );
    }

    un-named-reaction-198
    {
        type            reversibleArrheniusPLOG;
        reaction        "H + HNO2 = NO + H2O";
        A               3380000;
        beta            1.07;
        Ta              2800.259244;
        ArrheniusData
        (
            // PLOG /p A beta Ta/
            // convert atm to Pa
            // convert Ea to Ta, Ta = Ea/R (8.314 J⋅K−1⋅mol−1), cal to J 
            // (4.18400 J/cal)
            (1.00E3 3390000  1.07 2802.082271 )
            (1.00E4 3390000  1.07 2801.579023 )
            (3.16E4 3390000  1.07 2801.579023 )
            (1.00E5 3380000  1.07 2800.572528 )
            (3.16E5 3380000  1.07 2798.056291 )
            (1.00E6 3400000  1.07 2791.010825 )
            (3.16E6 4320000  1.04 2813.656964 )
            (1.00E7 12700000 0.91 3003.381285 )
        );
    }

## License

    This OpenFOAM library is under the GNU General Public License.

## Contact

    Maintainer: Zixin Chi  
    Email: chizixin@buaa.edu.cn  
    Issues: Please use the [GitHub Issues] page for bug reports and questions.

  If you use this library in your research, please cite the following paper:
  [1] Z. Chi, X. Hui, K. Ji, B. Wang. A high-performance chemistry solver for 
  reactive flow simulations in OpenFOAM. Physics of Fluids 37 (2025) 125154.


