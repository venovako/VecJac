#include "cmerge.h"

fint cmerge_(const fnat m[static restrict 1], const fnat n[static restrict 1], const float Ar[static restrict VSL], const fnat ldAr[static restrict 1], const float Ai[static restrict VSL], const fnat ldAi[static restrict 1], float complex A[static restrict VSL_2], const fnat ldA[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -9;
  if (*m & VSL__2)
    return -1;
  if (IS_NOT_ALIGNED(Ar))
    return -3;
  if (*ldAr < *m)
    return -4;
  if (*ldAr & VSL_1)
    return -4;
  if (IS_NOT_ALIGNED(Ai))
    return -5;
  if (*ldAi < *m)
    return -6;
  if (*ldAi & VSL_1)
    return -6;
  if (IS_NOT_ALIGNED(A))
    return -7;
  if (*ldA < *m)
    return -8;
  if (*ldA & VSL__2)
    return -8;
#endif /* !NDEBUG */

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,Ar,ldAr,Ai,ldAi)
  for (fnat j = 0u; j < *n; ++j) {
    register const VI idx = _mm512_set_epi32(15, 7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 8, 0);
    float complex *const Aj = A + j * (size_t)(*ldA);
    const float *const Arj = Ar + j * (size_t)(*ldAr);
    const float *const Aij = Ai + j * (size_t)(*ldAi);
    for (fnat i = 0u; i < *m; i += VSL_2) {
#ifdef __AVX512DQ__
      _mm512_store_ps((Aj + i), _mm512_permutexvar_ps(idx, _mm512_insertf32x8(_mm512_zextps256_ps512(_mm256_load_ps(Arj + i)), _mm256_load_ps(Aij + i), 0x01u)));
#else /* !__AVX512DQ__ */
      _mm512_store_ps((Aj + i), _mm512_permutexvar_ps(idx, _mm512_castpd_ps(_mm512_insertf64x4(_mm512_castps_pd(_mm512_zextps256_ps512(_mm256_load_ps(Arj + i))), _mm256_castps_pd(_mm256_load_ps(Aij + i)), 0x01u))));
#endif /* ?__AVX512DQ__ */
    }
  }
  return 1;
#else /* !_OPENMP */
  register const VI idx = _mm512_set_epi32(15, 7, 14, 6, 13, 5, 12, 4, 11, 3, 10, 2, 9, 1, 8, 0);
  for (fnat j = 0u; j < *n; ++j) {
    float complex *const Aj = A + j * (size_t)(*ldA);
    const float *const Arj = Ar + j * (size_t)(*ldAr);
    const float *const Aij = Ai + j * (size_t)(*ldAi);
    for (fnat i = 0u; i < *m; i += VSL_2) {
#ifdef __AVX512DQ__
      _mm512_store_ps((Aj + i), _mm512_permutexvar_ps(idx, _mm512_insertf32x8(_mm512_zextps256_ps512(_mm256_load_ps(Arj + i)), _mm256_load_ps(Aij + i), 0x01u)));
#else /* !__AVX512DQ__ */
      _mm512_store_ps((Aj + i), _mm512_permutexvar_ps(idx, _mm512_castpd_ps(_mm512_insertf64x4(_mm512_castps_pd(_mm512_zextps256_ps512(_mm256_load_ps(Arj + i))), _mm256_castps_pd(_mm256_load_ps(Aij + i)), 0x01u))));
#endif /* ?__AVX512DQ__ */
    }
  }
  return 0;
#endif /* ?_OPENMP */
}
