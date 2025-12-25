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










//#define SIMD_LOAD_MODE 0    	//use the safe unaligned load/store instruction
#define SIMD_LOAD_MODE 1    	//use the aligned load/store instruction
//#define SIMD_LOAD_MODE 2    	//use the aligned load/store instruction with security check, for Debug purpose

inline __m256d load256d(const double* p)
{

#if SIMD_LOAD_MODE == 0
    return _mm256_loadu_pd(p);

#elif SIMD_LOAD_MODE == 1
    return _mm256_load_pd(p);

#elif SIMD_LOAD_MODE == 2
    uintptr_t addr = reinterpret_cast<uintptr_t>(p);
    if (addr & 31)
    {
        fprintf(stderr,
            "[ERROR] load_pd called with UNALIGNED pointer %p\n", p);
        abort();
    }
    return _mm256_load_pd(p);

#else
#error "Invalid SIMD_LOAD_MODE (must be 0,1,2)."
#endif
}

inline void store256d(double* p, __m256d v)
{


#if SIMD_LOAD_MODE == 0
    _mm256_storeu_pd(p, v);

#elif SIMD_LOAD_MODE == 1
    _mm256_store_pd(p, v);

#elif SIMD_LOAD_MODE == 2
    uintptr_t addr = reinterpret_cast<uintptr_t>(p);
    if (addr & 31)
    {
        fprintf(stderr,
            "[ERROR] store_pd called with UNALIGNED pointer %p\n", p);
        abort();
    }
    _mm256_store_pd(p, v);
#else
    #error "Invalid SIMD_STORE_MODE (must be 0,1,2)."
#endif
}


inline __m128d load128d(const double* p)
{

#if SIMD_LOAD_MODE == 0
    return _mm_loadu_pd(p);

#elif SIMD_LOAD_MODE == 1
    return _mm_load_pd(p);

#elif SIMD_LOAD_MODE == 2
    uintptr_t addr = reinterpret_cast<uintptr_t>(p);
    if (addr & 15)
    {
        fprintf(stderr,
            "[ERROR] load_pd called with UNALIGNED pointer %p\n", p);
        abort();
    }
    return _mm_load_pd(p);
#else
    #error "Invalid SIMD_LOAD_MODE (must be 0,1,2)."
#endif
}

inline void store128d(double* p, __m128d v)
{


#if SIMD_LOAD_MODE == 0
    _mm_storeu_pd(p, v);

#elif SIMD_LOAD_MODE == 1
    _mm_store_pd(p, v);

#elif SIMD_LOAD_MODE == 2
    uintptr_t addr = reinterpret_cast<uintptr_t>(p);
    if (addr & 15)
    {
        fprintf(stderr,
            "[ERROR] store_pd called with UNALIGNED pointer %p\n", p);
        abort();
    }
    _mm_store_pd(p, v);
#else
#error "Invalid SIMD_STORE_MODE (must be 0,1,2)."
#endif
}


#endif

