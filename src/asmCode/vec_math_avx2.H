/*---------------------------------------------------------------------------*\
  Description
      Vector math wrapper functions for double-precision SIMD operations
      based on the AVX2 instruction set (256-bit and 128-bit vectors).

      The wrappers provide:
        - aligned / unaligned load and store of full or partial vectors
        - broadcast and algebraic operations (add / sub / mul / div / FMA)
        - comparison and selection (cmp_*, select)
        - horizontal reductions (hsum*), 4x4 transpose
        - element extraction (get*), abs / max
        - vectorized exp / pow / log (vec256_*)

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/

#ifndef VEC_MATH_AVX2_H
#define VEC_MATH_AVX2_H

//=============================================================================//

//---------------------------------
// 1. SIMD / AVX2 headers (built-in instruction set)
//---------------------------------
#include <immintrin.h>

//---------------------------------
// 2. FastChemistry headers
//---------------------------------
#include "Macro.H"

//=============================================================================//

namespace FastChemistry
{
    // Vectorized transcendental functions (exp / pow / log).
    // When USE_LOCALFILE_ is true (e.g. glibc < 2.22), the GNU vectorized
    // symbols are declared in an extern "C" block and wrapped; otherwise the
    // compiler built-in vectorized symbols are linked through asm labels.
    #if USE_LOCALFILE_

        extern "C"
        {
            __m256d _ZGVdN4v_gnuexp_avx2(__m256d);
            __m256d _ZGVdN4vv_gnupow_avx2(__m256d, __m256d);
            __m256d _ZGVdN4v_gnulog_avx2(__m256d);
        }

        // Vectorized exponential: returns exp(x) element-wise.
        static inline __m256d vec256_expd(__m256d x)
        {
            return _ZGVdN4v_gnuexp_avx2(x);
        }

        // Vectorized power: returns x^y element-wise.
        static inline __m256d vec256_powd(__m256d x, __m256d y)
        {
            return _ZGVdN4vv_gnupow_avx2(x, y);
        }

        // Vectorized natural logarithm: returns log(x) element-wise.
        static inline __m256d vec256_logd(__m256d x)
        {
            return _ZGVdN4v_gnulog_avx2(x);
        }

    #else

        extern __m256d vec256_expd(__m256d) asm("_ZGVdN4v_exp");
        extern __m256d vec256_powd(__m256d, __m256d) asm("_ZGVdN4vv_pow");
        extern __m256d vec256_logd(__m256d) asm("_ZGVdN4v_log");

    #endif

    // ---------- load: 256-bit vectors ----------

    // Load 4 aligned doubles into a 256-bit vector.
    static inline __m256d load256d(const double* p)
    {
        return _mm256_load_pd(p);
    }

    // Load 4 unaligned doubles into a 256-bit vector.
    static inline __m256d loadu256d(const double* p)
    {
        return _mm256_loadu_pd(p);
    }

    // Load the first 3 doubles, the 4th element is set to 0.
    static inline __m256d load256d012(const double* ptr)
    {
        __m256i mask = _mm256_set_epi64x(0, -1, -1, -1);
        return _mm256_maskload_pd(ptr, mask);
    }

    // Load the first 3 doubles (unaligned safe), the 4th element is set to 0.
    static inline __m256d loadu256d012(const double* ptr)
    {
        __m256i mask = _mm256_set_epi64x(0, -1, -1, -1);
        return _mm256_maskload_pd(ptr, mask);
    }

    // Load the first 2 doubles, the remaining elements are set to 0.
    static inline __m256d load256d01(const double* ptr)
    {
        __m128d lo = _mm_load_pd(ptr);
        return _mm256_insertf128_pd(_mm256_castpd128_pd256(lo), _mm_setzero_pd(), 1);
    }

    // Load the first 2 doubles (unaligned safe), the remaining elements are set to 0.
    static inline __m256d loadu256d01(const double* ptr)
    {
        __m128d lo = _mm_loadu_pd(ptr);
        return _mm256_insertf128_pd(_mm256_castpd128_pd256(lo), _mm_setzero_pd(), 1);
    }

    // Load the first double, the remaining elements are set to 0.
    static inline __m256d load256d0(const double* ptr)
    {
        __m128d lo = _mm_load_sd(ptr);
        return _mm256_insertf128_pd(_mm256_castpd128_pd256(lo), _mm_setzero_pd(), 1);
    }

    // ---------- store: 256-bit vectors ----------

    // Store 4 doubles to an aligned address.
    static inline void store256d(double* p, __m256d v)
    {
        _mm256_store_pd(p, v);
    }

    // Store 4 doubles to an unaligned address.
    static inline void storeu256d(double* p, __m256d v)
    {
        _mm256_storeu_pd(p, v);
    }

    // Store the first 3 doubles.
    static inline void store256d012(double* ptr, __m256d v)
    {
        __m256i mask = _mm256_set_epi64x(0, -1, -1, -1);
        _mm256_maskstore_pd(ptr, mask, v);
    }

    // Store the first 3 doubles (unaligned safe).
    static inline void storeu256d012(double* ptr, __m256d v)
    {
        __m256i mask = _mm256_set_epi64x(0, -1, -1, -1);
        _mm256_maskstore_pd(ptr, mask, v);
    }

    // Store the first 2 doubles.
    static inline void store256d01(double* ptr, __m256d v)
    {
        __m128d lo = _mm256_castpd256_pd128(v);
        _mm_store_pd(ptr, lo);
    }

    // Store the first 2 doubles (unaligned safe).
    static inline void storeu256d01(double* ptr, __m256d v)
    {
        __m128d lo = _mm256_castpd256_pd128(v);
        _mm_storeu_pd(ptr, lo);
    }

    // Store the first double.
    static inline void store256d0(double* ptr, __m256d v)
    {
        __m128d lo = _mm256_castpd256_pd128(v);
        _mm_store_sd(ptr, lo);
    }

    // Store the first double (unaligned safe).
    static inline void storeu256d0(double* ptr, __m256d v)
    {
        __m128d lo = _mm256_castpd256_pd128(v);
        _mm_store_sd(ptr, lo);
    }

    // ---------- load / store: 128-bit vectors ----------

    // Load 2 aligned doubles into a 128-bit vector.
    static inline __m128d load128d(const double* p)
    {
        return _mm_load_pd(p);
    }

    // Load 2 unaligned doubles into a 128-bit vector.
    static inline __m128d loadu128d(const double* p)
    {
        return _mm_loadu_pd(p);
    }

    // Load the first double, the second element is set to 0.
    static inline __m128d load128d0(const double* ptr)
    {
        return _mm_load_sd(ptr);
    }

    // Store 2 doubles to an aligned address.
    static inline void store128d(double* p, __m128d v)
    {
        _mm_store_pd(p, v);
    }

    // Store 2 doubles to an unaligned address.
    static inline void storeu128d(double* p, __m128d v)
    {
        _mm_storeu_pd(p, v);
    }

    // Store the first double.
    static inline void store128d0(double* p, __m128d v)
    {
        _mm_store_sd(p, v);
    }

    // ---------- broadcast ----------

    // Broadcast a scalar to all 2 lanes of a 128-bit vector.
    static inline __m128d set1_128d(double v)
    {
        return _mm_set1_pd(v);
    }

    // Broadcast a scalar to all 4 lanes of a 256-bit vector.
    static inline __m256d set1_256d(double v)
    {
        return _mm256_set1_pd(v);
    }

    // ---------- algebraic operations: 128-bit ----------

    // a + b, where the second operand can be a vector or a memory location.
    static inline __m128d add128d(__m128d a, __m128d b)
    {
        return _mm_add_pd(a, b);
    }
    static inline __m128d add128d(__m128d a, const double* b)
    {
        return _mm_add_pd(a, load128d(b));
    }

    // a - b, where both operands can be a vector or a memory location.
    static inline __m128d sub128d(__m128d a, __m128d b)
    {
        return _mm_sub_pd(a, b);
    }
    static inline __m128d sub128d(__m128d a, const double* b)
    {
        return _mm_sub_pd(a, load128d(b));
    }
    static inline __m128d sub128d(const double* a, __m128d b)
    {
        return _mm_sub_pd(load128d(a), b);
    }

    // a * b, where both operands can be a vector or a memory location.
    static inline __m128d mul128d(__m128d a, __m128d b)
    {
        return _mm_mul_pd(a, b);
    }
    static inline __m128d mul128d(__m128d a, const double* b)
    {
        return _mm_mul_pd(a, load128d(b));
    }
    static inline __m128d mul128d(const double* a, __m128d b)
    {
        return _mm_mul_pd(load128d(a), b);
    }

    // a / b, where both operands can be a vector or a memory location.
    static inline __m128d div(__m128d a, __m128d b)
    {
        return _mm_div_pd(a, b);
    }
    static inline __m128d div(__m128d a, const double* b)
    {
        return _mm_div_pd(a, load128d(b));
    }
    static inline __m128d div(const double* a, __m128d b)
    {
        return _mm_div_pd(load128d(a), b);
    }

    // Fused multiply-add: a * b + c, where c can be a vector or a memory location.
    static inline __m128d fmadd128d(__m128d a, __m128d b, __m128d c)
    {
        return _mm_fmadd_pd(a, b, c);
    }
    static inline __m128d fmadd128d(__m128d a, __m128d b, const double* c)
    {
        return _mm_fmadd_pd(a, b, load128d(c));
    }

    // Fused multiply-subtract: a * b - c, where c can be a vector or a memory location.
    static inline __m128d fmsub128d(__m128d a, __m128d b, __m128d c)
    {
        return _mm_fmsub_pd(a, b, c);
    }
    static inline __m128d fmsub128d(__m128d a, __m128d b, const double* c)
    {
        return _mm_fmsub_pd(a, b, load128d(c));
    }

    // Fused negative multiply-add: -a * b + c.
    static inline __m128d fnmadd128d(__m128d a, __m128d b, __m128d c)
    {
        return _mm_fnmadd_pd(a, b, c);
    }

    // Fused negative multiply-subtract: -a * b - c, where c can be a vector or a memory location.
    static inline __m128d fnmsub128d(__m128d a, __m128d b, __m128d c)
    {
        return _mm_fnmsub_pd(a, b, c);
    }
    static inline __m128d fnmsub128d(__m128d a, __m128d b, const double* c)
    {
        return _mm_fnmsub_pd(a, b, load128d(c));
    }

    // ---------- algebraic operations: 256-bit ----------

    // a + b, where both operands can be a vector or a memory location.
    static inline __m256d add256d(__m256d a, __m256d b)
    {
        return _mm256_add_pd(a, b);
    }
    static inline __m256d add256d(__m256d a, const double* b)
    {
        return _mm256_add_pd(a, load256d(b));
    }
    static inline __m256d add256d(const double* a, __m256d b)
    {
        return _mm256_add_pd(load256d(a), b);
    }

    // a - b, where both operands can be a vector or a memory location.
    static inline __m256d sub256d(__m256d a, __m256d b)
    {
        return _mm256_sub_pd(a, b);
    }
    static inline __m256d sub256d(__m256d a, const double* b)
    {
        return _mm256_sub_pd(a, load256d(b));
    }
    static inline __m256d sub256d(const double* a, __m256d b)
    {
        return _mm256_sub_pd(load256d(a), b);
    }

    // a * b, where both operands can be a vector or a memory location.
    static inline __m256d mul256d(__m256d a, __m256d b)
    {
        return _mm256_mul_pd(a, b);
    }
    static inline __m256d mul256d(__m256d a, const double* b)
    {
        return _mm256_mul_pd(a, load256d(b));
    }
    static inline __m256d mul256d(const double* a, __m256d b)
    {
        return _mm256_mul_pd(load256d(a), b);
    }

    // a / b, where both operands can be a vector or a memory location.
    static inline __m256d div256d(__m256d a, __m256d b)
    {
        return _mm256_div_pd(a, b);
    }
    static inline __m256d div256d(__m256d a, const double* b)
    {
        return _mm256_div_pd(a, load256d(b));
    }
    static inline __m256d div256d(const double* a, __m256d b)
    {
        return _mm256_div_pd(load256d(a), b);
    }

    // Fused multiply-add: a * b + c, where c can be a vector or a memory location.
    static inline __m256d fmadd256d(__m256d a, __m256d b, __m256d c)
    {
        return _mm256_fmadd_pd(a, b, c);
    }


    // Fused multiply-subtract: a * b - c, where c can be a vector or a memory location.
    static inline __m256d fmsub256d(__m256d a, __m256d b, __m256d c)
    {
        return _mm256_fmsub_pd(a, b, c);
    }


    // Fused negative multiply-add: -a * b + c, where c can be a vector or a memory location.
    static inline __m256d fnmadd256d(__m256d a, __m256d b, __m256d c)
    {
        return _mm256_fnmadd_pd(a, b, c);
    }

    // Fused negative multiply-subtract: -a * b - c, where c can be a vector or a memory location.
    static inline __m256d fnmsub256d(__m256d a, __m256d b, __m256d c)
    {
        return _mm256_fnmsub_pd(a, b, c);
    }

    // ---------- comparison and selection ----------

    // Select lanes from true_val where mask is set, otherwise from false_val.
    static inline __m256d select(__m256d mask, __m256d true_val, __m256d false_val)
    {
        return _mm256_blendv_pd(false_val, true_val, mask);
    }

    // Element-wise comparisons; each lane of the result is 0xFFFFFFFFFFFFFFFF
    // (true) or 0 (false). cmp_eq: a == b.
    static inline __m256d cmp_eq(__m256d a, __m256d b)
    {
        return _mm256_cmp_pd(a, b, _CMP_EQ_OQ);
    }
    // a != b.
    static inline __m256d cmp_neq(__m256d a, __m256d b)
    {
        return _mm256_cmp_pd(a, b, _CMP_NEQ_OQ);
    }
    // a > b.
    static inline __m256d cmp_gt(__m256d a, __m256d b)
    {
        return _mm256_cmp_pd(a, b, _CMP_GT_OQ);
    }
    // a >= b.
    static inline __m256d cmp_ge(__m256d a, __m256d b)
    {
        return _mm256_cmp_pd(a, b, _CMP_GE_OQ);
    }
    // a < b.
    static inline __m256d cmp_lt(__m256d a, __m256d b)
    {
        return _mm256_cmp_pd(a, b, _CMP_LT_OQ);
    }
    // a <= b.
    static inline __m256d cmp_le(__m256d a, __m256d b)
    {
        return _mm256_cmp_pd(a, b, _CMP_LE_OQ);
    }

    // ---------- horizontal reduction ----------

    // Sum of all 4 lanes of a 256-bit vector.
    static inline double hsum4(__m256d v)
    {
        __m128d lo = _mm256_castpd256_pd128(v);
        __m128d hi = _mm256_extractf128_pd(v, 1);
        lo = _mm_add_pd(lo, hi);
        __m128d sh = _mm_unpackhi_pd(lo, lo);
        __m128d s = _mm_add_sd(lo, sh);
        return _mm_cvtsd_f64(s);
    }

    // Pairwise sums of two vectors: returns [sum(a), sum(b)] as a 128-bit vector.
    static inline __m128d hsum4x2(__m256d a, __m256d b)
    {
        __m256d t0 = _mm256_hadd_pd(a, b);                              // [a0+a1, b0+b1, a2+a3, b2+b3]
        t0 = _mm256_permute4x64_pd(t0, _MM_SHUFFLE(3,1,2,0));           // [a0+a1, a2+a3, b0+b1, b2+b3]
        t0 = _mm256_hadd_pd(t0, t0);                                    // [a_sum, a_sum, b_sum, b_sum]
        __m128d low  = _mm256_castpd256_pd128(t0);                      // [a_sum, a_sum]
        __m128d high = _mm256_extractf128_pd(t0, 1);                    // [b_sum, b_sum]
        return _mm_unpacklo_pd(low, high);                              // [a_sum, b_sum]
    }

    // Sum of 4 vectors lane-by-lane: returns the element-wise total of a, b, c, d.
    static inline __m256d hsum4x4(__m256d sum0, __m256d sum1, __m256d sum2, __m256d sum3)
    {
        __m256d t0 = _mm256_hadd_pd(sum0, sum1);                        // [a0+a1, b0+b1, a2+a3, b2+b3]
        __m256d t1 = _mm256_hadd_pd(sum2, sum3);                        // [c0+c1, d0+d1, c2+c3, d2+d3]
        __m256d t2 = _mm256_permute4x64_pd(t0, _MM_SHUFFLE(3,1,2,0));   // [a0+a1, a2+a3, b0+b1, b2+b3]
        __m256d t3 = _mm256_permute4x64_pd(t1, _MM_SHUFFLE(3,1,2,0));   // [c0+c1, c2+c3, d0+d1, d2+d3]
        __m256d t4 = _mm256_hadd_pd(t2, t3);                            // [a0+a1+a2+a3, c0+c1+c2+c3, b0+b1+b2+b3, d0+d1+d2+d3]
        __m256d t5 = _mm256_permute4x64_pd(t4, _MM_SHUFFLE(3,1,2,0));
        return t5;
    }

    // ----------  transpose ----------

    // In-place 4x4 transpose of four 256-bit vectors (each row is one vector).
    static inline void transpose4x4_pd(__m256d& v0, __m256d& v1, __m256d& v2, __m256d& v3)
    {
        __m256d t0 = _mm256_unpacklo_pd(v0, v1);                        // [a0 b0 a1 b1]
        __m256d t1 = _mm256_unpackhi_pd(v0, v1);                        // [a2 b2 a3 b3]
        __m256d t2 = _mm256_unpacklo_pd(v2, v3);                        // [c0 d0 c1 d1]
        __m256d t3 = _mm256_unpackhi_pd(v2, v3);                        // [c2 d2 c3 d3]

        v0 = _mm256_permute2f128_pd(t0, t2, 0x20);                      // [a0 b0 c0 d0]
        v1 = _mm256_permute2f128_pd(t1, t3, 0x20);                      // [a1 b1 c1 d1]
        v2 = _mm256_permute2f128_pd(t0, t2, 0x31);                      // [a2 b2 c2 d2]
        v3 = _mm256_permute2f128_pd(t1, t3, 0x31);                      // [a3 b3 c3 d3]
    }

    // In-place 2x2 transpose of two 128-bit vectors (each row is one vector).
    static inline void transpose2x2_pd(__m128d& v0, __m128d& v1)
    {
        __m128d t0 = _mm_shuffle_pd(v0, v1, 0);  // [a0, b0]
        __m128d t1 = _mm_shuffle_pd(v0, v1, 3);  // [a1, b1]
        
        v0 = t0;
        v1 = t1;
    }

    // In-place 2x4 transpose of two 256-bit vectors (each row is one vector).
    static inline void transpose2x4_pd(__m256d& v0, __m256d& v1, __m128d& w0, __m128d& w1, __m128d& w2, __m128d& w3)
    {
        // [a0 a1 a2 a3]
        // [b0 b1 b2 b3]
        __m128d v0_lo = _mm256_castpd256_pd128(v0);          // [a0, a1]
        __m128d v0_hi = _mm256_extractf128_pd(v0, 1);        // [a2, a3]
        __m128d v1_lo = _mm256_castpd256_pd128(v1);          // [b0, b1]
        __m128d v1_hi = _mm256_extractf128_pd(v1, 1);        // [b2, b3]

        w0 = _mm_unpacklo_pd(v0_lo, v1_lo);   // [a0, b0]
        w1 = _mm_unpackhi_pd(v0_lo, v1_lo);   // [a1, b1]
        w2 = _mm_unpacklo_pd(v0_hi, v1_hi);   // [a2, b2]
        w3 = _mm_unpackhi_pd(v0_hi, v1_hi);   // [a3, b3]
    }

    // In-place 2x4 transpose of four 128-bit vectors (each row is one vector).
    static inline void transpose4x2_pd(__m256d& v0, __m256d& v1, __m128d& w0, __m128d& w1, __m128d& w2, __m128d& w3)
    {
        __m128d a01 = _mm_unpacklo_pd(w0, w1);   // [a0, a1]
        __m128d b01 = _mm_unpackhi_pd(w0, w1);   // [b0, b1]
        __m128d a23 = _mm_unpacklo_pd(w2, w3);   // [a2, a3]
        __m128d b23 = _mm_unpackhi_pd(w2, w3);   // [b2, b3]

        v0 = _mm256_set_m128d(a23, a01);         
        v1 = _mm256_set_m128d(b23, b01);         
    }

    // ---------- element extraction ----------

    // Extract lane 0 of a 256-bit vector as a scalar.
    static inline double get0(__m256d vec)
    {
        return _mm256_cvtsd_f64(vec);
    }
    // Extract lane 0 of a 128-bit vector as a scalar.
    static inline double get0(__m128d v)
    {
        return _mm_cvtsd_f64(v);
    }
    // Extract lane 1 of a 256-bit vector as a scalar.
    static inline double get1(__m256d vec)
    {
        __m128d low = _mm256_castpd256_pd128(vec);
        return _mm_cvtsd_f64(_mm_unpackhi_pd(low, low));
    }
    // Extract lane 1 of a 128-bit vector as a scalar.
    static inline double get1(__m128d v)
    {
        return _mm_cvtsd_f64(_mm_unpackhi_pd(v, v));
    }
    // Extract lane 2 of a 256-bit vector as a scalar.
    static inline double get2(__m256d vec)
    {
        __m128d high = _mm256_extractf128_pd(vec, 1);
        return _mm_cvtsd_f64(high);
    }
    // Extract lane 3 of a 256-bit vector as a scalar.
    static inline double get3(__m256d vec)
    {
        __m128d high = _mm256_extractf128_pd(vec, 1);
        return _mm_cvtsd_f64(_mm_unpackhi_pd(high, high));
    }
    // Extract lane 34 of a 256-bit vector as a 128-bit vector.
    static inline __m128d get23(__m256d vec)
    {
        return _mm256_extractf128_pd(vec, 1);
    }
    // Extract lane 12 of a 256-bit vector as a 128-bit vector.
    static inline __m128d get01(__m256d vec)
    {
        return _mm256_castpd256_pd128(vec);
    }
    // ---------- misc ----------

    // Maximum absolute element of a 256-bit vector.
    static inline double maxAbsElement(__m256d v)
    {
        const __m256d signMask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7fffffffffffffffLL));
        __m256d absV = _mm256_and_pd(v, signMask);
        __m256d max1 = _mm256_max_pd(absV, _mm256_permute4x64_pd(absV, 0b01001110));
        __m256d max2 = _mm256_max_pd(max1, _mm256_permute_pd(max1, 0b0101));
        return _mm256_cvtsd_f64(max2);
    }

    // Element-wise absolute value of a 256-bit vector.
    static inline __m256d abs256d(__m256d v)
    {
        const __m256i mask = _mm256_set1_epi64x(0x7FFFFFFFFFFFFFFFLL);
        return _mm256_and_pd(v, _mm256_castsi256_pd(mask));
    }

    // Broadcast the maximum absolute element of a 256-bit vector to all lanes.
    static inline __m256d broadcastMax(__m256d v)
    {
        __m256d max1 = _mm256_max_pd(v, _mm256_permute4x64_pd(v, 0b01001110));
        __m256d max2 = _mm256_max_pd(max1, _mm256_permute_pd(max1, 0b0101));
        return max2;
    }
}
#endif
