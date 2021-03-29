#include "zmerge.h"

fint zmerge_(const fnat m[static restrict 1], const fnat n[static restrict 1], const double Ar[static restrict VDL], const fnat ldAr[static restrict 1], const double Ai[static restrict VDL], const fnat ldAi[static restrict 1], double complex A[static restrict VDL_2], const fnat ldA[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -2;
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
    register const VI ri = _mm512_set_epi64(7, 3, 6, 2, 5, 1, 4, 0);
    register const VI ii = _mm512_set_epi64(3, 7, 2, 6, 1, 5, 0, 4);
    double complex *const Aj = A + j * (*ldA);
    const double *const Arj = Ar + j * (*ldAr);
    const double *const Aij = Ai + j * (*ldAi);
    for (fnat i = 0u; i < *m; i += VDL_2)
      _mm512_store_pd((Aj + i), _mm512_mask_blend_pd(0xAAu, _mm512_permutexvar_pd(ri, _mm512_zextpd256_pd512(_mm256_load_pd(Arj + i))), _mm512_permutexvar_pd(ii, _mm512_zextpd256_pd512(_mm256_load_pd(Aij + i)))));
  }

  return omp_get_max_threads();
#else /* !_OPENMP */
  register const VI ri = _mm512_set_epi64(7, 3, 6, 2, 5, 1, 4, 0);
  register const VI ii = _mm512_set_epi64(3, 7, 2, 6, 1, 5, 0, 4);

  for (fnat j = 0u; j < *n; ++j) {
    double complex *const Aj = A + j * (*ldA);
    const double *const Arj = Ar + j * (*ldAr);
    const double *const Aij = Ai + j * (*ldAi);
    for (fnat i = 0u; i < *m; i += VDL_2)
      _mm512_store_pd((Aj + i), _mm512_mask_blend_pd(0xAAu, _mm512_permutexvar_pd(ri, _mm512_zextpd256_pd512(_mm256_load_pd(Arj + i))), _mm512_permutexvar_pd(ii, _mm512_zextpd256_pd512(_mm256_load_pd(Aij + i)))));
  }

  return 0;
#endif /* ?_OPENMP */
}
