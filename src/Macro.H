#ifndef MACRO_H
#define MACRO_H

#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>





#define USE_LOCALFILE_ false    //use the vectorized math function in system library
//#define USE_LOCALFILE_ true   //use the vectorized math function in asmCode folder


//#define SIMD_LOAD_MODE 0    	//use the safe unaligned load/store instruction
#define SIMD_LOAD_MODE 1    	//use the aligned load/store instruction
//#define SIMD_LOAD_MODE 2    	//use the aligned load/store instruction with security check, for Debug purpose


#define KcLimiter (1.4901171103413047e-8)
//#define KcLimiter 0.0
#define invKcLimiter (1.0/KcLimiter)


#define TroeLimiter 2.220446049250313e-16

#define SRILimiter 2.220446049250313e-16

#define LuLimiter 2.220446049250313e-16






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

