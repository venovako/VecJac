#include "csplit.h"

fint csplit_(const fnat m[static restrict 1], const fnat n[static restrict 1], const float complex A[static restrict VSL_2], const fnat ldA[static restrict 1], float Ar[static restrict VSL], const fnat ldAr[static restrict 1], float Ai[static restrict VSL], const fnat ldAi[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -9;
  if (*m & VSL__2)
    return -1;
  if (IS_NOT_ALIGNED(A))
    return -3;
  if (*ldA < *m)
    return -4;
  if (*ldA & VSL__2)
    return -4;
  if (IS_NOT_ALIGNED(Ar))
    return -5;
  if (*ldAr < *m)
    return -6;
  if (*ldAr & VSL_1)
    return -6;
  if (IS_NOT_ALIGNED(Ai))
    return -7;
  if (*ldAi < *m)
    return -8;
  if (*ldAi & VSL_1)
    return -8;
#endif /* !NDEBUG */

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,Ar,ldAr,Ai,ldAi)
  for (fnat j = 0u; j < *n; ++j) {
    register const VI idx = _mm512_set_epi32(15, 13, 11, 9, 7, 5, 3, 1, 14, 12, 10, 8, 6, 4, 2, 0);
    const float complex *const Aj = A + j * (size_t)(*ldA);
    float *const Arj = Ar + j * (size_t)(*ldAr);
    float *const Aij = Ai + j * (size_t)(*ldAi);
    for (fnat i = 0u; i < *m; i += VSL_2) {
      register const VS ri = _mm512_permutexvar_ps(idx, _mm512_load_ps(Aj + i));
#ifdef __AVX512DQ__
      _mm256_store_ps((Arj + i), _mm512_extractf32x8_ps(ri, 0x00u));
      _mm256_store_ps((Aij + i), _mm512_extractf32x8_ps(ri, 0x01u));
#else /* !__AVX512DQ__ */
      _mm256_store_ps((Arj + i), _mm512_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(ri), 0x00u)));
      _mm256_store_ps((Aij + i), _mm512_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(ri), 0x01u)));
#endif /* ?__AVX512DQ__ */
    }
  }
  return 1;
#else /* !_OPENMP */
  register const VI idx = _mm512_set_epi32(15, 13, 11, 9, 7, 5, 3, 1, 14, 12, 10, 8, 6, 4, 2, 0);
  for (fnat j = 0u; j < *n; ++j) {
    const float complex *const Aj = A + j * (size_t)(*ldA);
    float *const Arj = Ar + j * (size_t)(*ldAr);
    float *const Aij = Ai + j * (size_t)(*ldAi);
    for (fnat i = 0u; i < *m; i += VSL_2) {
      register const VS ri = _mm512_permutexvar_ps(idx, _mm512_load_ps(Aj + i)); VSP(ri);
#ifdef __AVX512DQ__
      _mm256_store_ps((Arj + i), _mm512_extractf32x8_ps(ri, 0x00u));
      _mm256_store_ps((Aij + i), _mm512_extractf32x8_ps(ri, 0x01u));
#else /* !__AVX512DQ__ */
      _mm256_store_ps((Arj + i), _mm512_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(ri), 0x00u)));
      _mm256_store_ps((Aij + i), _mm512_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(ri), 0x01u)));
#endif /* ?__AVX512DQ__ */
    }
  }
  return 0;
#endif /* ?_OPENMP */
}
