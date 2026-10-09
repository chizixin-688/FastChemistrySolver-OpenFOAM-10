/*---------------------------------------------------------------------------*\
  Description
      Some compile options and parameter settings for OpenFOAM-10

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

#ifndef MACRO_H
#define MACRO_H

//=============================================================================//

//---------------------------------
// 1. SIMD / AVX2 headers
//---------------------------------
#include <immintrin.h>

//=============================================================================//



//=============================================================================//

// Select vectorized math function source:
//   false -> system glibc implementation
//   true  -> project-provided glibc 2.22 implementation 
//
// IMPORTANT:
//   For systems with glibc < 2.22, set USE_LOCALFILE_ to true.
//
// NOTE:
//   1.If the system glibc version is lower than 2.22, the system library does
//     NOT provide the required vectorized math symbols. In this case, you MUST
//     set USE_LOCALFILE_ to true to avoid linker/runtime errors.
//   2.When USE_LOCALFILE_ is set to true, this project uses hand-written
//     assembly files for vectorized math functions. In this mode, OpenFOAM's
//     wmake build system is NOT supported.
//     Please use CMake + make to build the project instead.
#define USE_LOCALFILE_ false    
//#define USE_LOCALFILE_ true 

//=============================================================================//



//=============================================================================//

namespace FastChemistry
{
    
    constexpr const char* thermoDictName = "physicalProperties";
    //constexpr const char* thermoDictName = "thermoPhysicalProperties";

    constexpr const char* chemistryDictName = "chemistryProperties";
    //constexpr const char* chemistryDictName = "chemistryProperties";

    constexpr const double KcLimiter = 1.4901171103413047e-8;

    constexpr const double invKcLimiter = (1.0/KcLimiter);

    constexpr const double TroeLimiter = 2.220446049250313e-16;

    constexpr const double SRILimiter = 2.220446049250313e-16;

    constexpr const double LuLimiter = 2.220446049250313e-16;

    constexpr const double ConcentrationLimiter = 2.220446049250313e-16;
}


#endif

