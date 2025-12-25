/*---------------------------------------------------------------------------*\
  Description
      Inline function file for LU decomposition of matrix.

  Author
      Zixin Chi <chizixin@buaa.edu.cn>
\*---------------------------------------------------------------------------*/


inline double& FastChemistry::LUsolver::operator()
(
    const unsigned int i,
    const unsigned int j
)
{
    return v_[i*N + j];
}


inline double FastChemistry::LUsolver::operator() 
(
    const unsigned int i,
    const unsigned int j
)const
{
    return v_[i*N + j];
}

inline unsigned int& FastChemistry::LUsolver::pivotIndice
(
    const unsigned int i
)
{
    return this->pivotIndice_[i];
}

inline unsigned int FastChemistry::LUsolver::pivotIndice
(
    const unsigned int i
) const
{
    return this->pivotIndice_[i];
}

inline double FastChemistry::LUsolver::hsum4(__m256d v)
{
    __m128d lo = _mm256_castpd256_pd128(v);
    __m128d hi = _mm256_extractf128_pd(v,1);
    lo = _mm_add_pd(lo,hi);               
    __m128d sh = _mm_unpackhi_pd(lo,lo);  
    __m128d s = _mm_add_sd(lo, sh);       
    return _mm_cvtsd_f64(s);
}

inline void FastChemistry::LUsolver::transpose4x4_pd(__m256d& v0,__m256d& v1,__m256d& v2,__m256d& v3)
{
    __m256d t0 = _mm256_unpacklo_pd(v0, v1); // a0 b0 a1 b1
    __m256d t1 = _mm256_unpackhi_pd(v0, v1); // a2 b2 a3 b3
    __m256d t2 = _mm256_unpacklo_pd(v2, v3); // c0 d0 c1 d1
    __m256d t3 = _mm256_unpackhi_pd(v2, v3); // c2 d2 c3 d3

    v0 = _mm256_permute2f128_pd(t0, t2, 0x20); // [a0 b0 c0 d0]
    v1 = _mm256_permute2f128_pd(t1, t3, 0x20); // [a1 b1 c1 d1]
    v2 = _mm256_permute2f128_pd(t0, t2, 0x31); // [a2 b2 c2 d2]
    v3 = _mm256_permute2f128_pd(t1, t3, 0x31); // [a3 b3 c3 d3]
}

