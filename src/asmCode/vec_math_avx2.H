#ifndef VEC_MATH_AVX2_H
#define VEC_MATH_AVX2_H

#include <immintrin.h>

#include "Macro.H"
namespace FastChemistry
{
    #if USE_LOCALFILE_

        extern "C" 
        {
            __m256d _ZGVdN4v_gnuexp_avx2(__m256d);
            __m256d _ZGVdN4vv_gnupow_avx2(__m256d, __m256d);
            __m256d _ZGVdN4v_gnulog_avx2(__m256d);
        }

        static inline __m256d vec256_expd(__m256d x) 
        {
            return _ZGVdN4v_gnuexp_avx2(x);
        }
        static inline __m256d vec256_powd(__m256d x, __m256d y) 
        {
            return _ZGVdN4vv_gnupow_avx2(x, y);
        }
        static inline __m256d vec256_logd(__m256d x) 
        {
            return _ZGVdN4v_gnulog_avx2(x);
        }

    #else

        extern __m256d vec256_expd(__m256d) asm("_ZGVdN4v_exp");
        extern __m256d vec256_powd(__m256d, __m256d) asm("_ZGVdN4vv_pow");
        extern __m256d vec256_logd(__m256d) asm("_ZGVdN4v_log");

    #endif
}
#endif
