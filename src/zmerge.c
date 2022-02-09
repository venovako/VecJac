#include "zmerge.h"

fint zmerge_(const fnat m[static restrict 1], const fnat n[static restrict 1], const double Ar[static restrict VDL], const fnat ldAr[static restrict 1], const double Ai[static restrict VDL], const fnat ldAi[static restrict 1], double complex A[static restrict VDL_2], const fnat ldA[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -9;
  if (*m & VDL__2)
    return -1;
  if (IS_NOT_ALIGNED(Ar))
    return -3;
  if (*ldAr < *m)
    return -4;
  if (*ldAr & VDL_1)
    return -4;
  if (IS_NOT_ALIGNED(Ai))
    return -5;
  if (*ldAi < *m)
    return -6;
  if (*ldAi & VDL_1)
    return -6;
  if (IS_NOT_ALIGNED(A))
    return -7;
  if (*ldA < *m)
    return -8;
  if (*ldA & VDL__2)
    return -8;
#endif /* !NDEBUG */

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,Ar,ldAr,Ai,ldAi)
  for (fnat j = 0u; j < *n; ++j) {
    register const VI idx = _mm512_set_epi64(7, 3, 6, 2, 5, 1, 4, 0);
    double complex *const Aj = A + j * (size_t)(*ldA);
    const double *const Arj = Ar + j * (size_t)(*ldAr);
    const double *const Aij = Ai + j * (size_t)(*ldAi);
    for (fnat i = 0u; i < *m; i += VDL_2)
      _mm512_store_pd((Aj + i), _mm512_permutexvar_pd(idx, _mm512_insertf64x4(_mm512_zextpd256_pd512(_mm256_load_pd(Arj + i)), _mm256_load_pd(Aij + i), 0x01u)));
  }
  return 1;
#else /* !_OPENMP */
  register const VI idx = _mm512_set_epi64(7, 3, 6, 2, 5, 1, 4, 0);
  for (fnat j = 0u; j < *n; ++j) {
    double complex *const Aj = A + j * (size_t)(*ldA);
    const double *const Arj = Ar + j * (size_t)(*ldAr);
    const double *const Aij = Ai + j * (size_t)(*ldAi);
    for (fnat i = 0u; i < *m; i += VDL_2)
      _mm512_store_pd((Aj + i), _mm512_permutexvar_pd(idx, _mm512_insertf64x4(_mm512_zextpd256_pd512(_mm256_load_pd(Arj + i)), _mm256_load_pd(Aij + i), 0x01u)));
  }
  return 0;
#endif /* ?_OPENMP */
}
