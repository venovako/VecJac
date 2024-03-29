#include "aalloc.h"

fint salloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], float **const restrict A, fnat ldA[static restrict 1])
{
  if (A)
    *A = (float*)NULL;
  *ldA = 0u;
  if (!*m)
    return 0;
  if (!*n)
    return 0;

  const fnat k = *m & VSL_1;
  *ldA = (k ? (*m + (VSL - k)) : *m);

  if (A) {
    const size_t s = (*n) * ((*ldA) * sizeof(float));
    *A = (float*)aligned_alloc(VA, s);
    if (!*A)
      return -3;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,A,ldA)
#endif /* ?_OPENMP */
    for (fnat j = 0u; j < *n; ++j) {
      register const VS z = _mm512_setzero_ps();
      float *const Aj = *A + j * (size_t)(*ldA);
      for (fnat i = 0u; i < *ldA; i += VSL)
        _mm512_store_ps((Aj + i), z);
    }
  }

#ifdef _OPENMP
  return 1;
#else /* !_OPENMP */
  return 0;
#endif /* ?_OPENMP */
}

fint dalloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], double **const restrict A, fnat ldA[static restrict 1])
{
  if (A)
    *A = (double*)NULL;
  *ldA = 0u;
  if (!*m)
    return 0;
  if (!*n)
    return 0;

  const fnat k = *m & VDL_1;
  *ldA = (k ? (*m + (VDL - k)) : *m);

  if (A) {
    const size_t s = (*n) * ((*ldA) * sizeof(double));
    *A = (double*)aligned_alloc(VA, s);
    if (!*A)
      return -3;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,A,ldA)
#endif /* ?_OPENMP */
    for (fnat j = 0u; j < *n; ++j) {
      register const VD z = _mm512_setzero_pd();
      double *const Aj = *A + j * (size_t)(*ldA);
      for (fnat i = 0u; i < *ldA; i += VDL)
        _mm512_store_pd((Aj + i), z);
    }
  }

#ifdef _OPENMP
  return 1;
#else /* !_OPENMP */
  return 0;
#endif /* ?_OPENMP */
}

fint calloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], float complex **const restrict A, fnat ldA[static restrict 1], float **const restrict Ar, fnat ldAr[static restrict 1], float **const restrict Ai, fnat ldAi[static restrict 1])
{
  if (A)
    *A = (float complex*)NULL;
  *ldA = 0u;
  if (Ar)
    *Ar = (float*)NULL;
  *ldAr = 0u;
  if (Ai)
    *Ai = (float*)NULL;
  *ldAi = 0u;
  if (!*m)
    return 0;
  if (!*n)
    return 0;

  const fnat k = *m & VSL__2;
  *ldA = (k ? (*m + (VSL_2 - k)) : *m);

  if (A) {
    const size_t s = (*n) * ((*ldA) * sizeof(float complex));
    *A = (float complex*)aligned_alloc(VA, s);
    if (!*A)
      return -3;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,A,ldA)
#endif /* ?_OPENMP */
    for (fnat j = 0u; j < *n; ++j) {
      register const VS z = _mm512_setzero_ps();
      float complex *const Aj = *A + j * (size_t)(*ldA);
      for (fnat i = 0u; i < *ldA; i += VSL_2)
        _mm512_store_ps((Aj + i), z);
    }
  }

  if (salloc2_(m, n, Ar, ldAr) < 0)
    return -5;

  if (salloc2_(m, n, Ai, ldAi) < 0)
    return -7;

#ifdef _OPENMP
  return 1;
#else /* !_OPENMP */
  return 0;
#endif /* ?_OPENMP */
}

fint zalloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], double complex **const restrict A, fnat ldA[static restrict 1], double **const restrict Ar, fnat ldAr[static restrict 1], double **const restrict Ai, fnat ldAi[static restrict 1])
{
  if (A)
    *A = (double complex*)NULL;
  *ldA = 0u;
  if (Ar)
    *Ar = (double*)NULL;
  *ldAr = 0u;
  if (Ai)
    *Ai = (double*)NULL;
  *ldAi = 0u;
  if (!*m)
    return 0;
  if (!*n)
    return 0;

  const fnat k = *m & VDL__2;
  *ldA = (k ? (*m + (VDL_2 - k)) : *m);

  if (A) {
    const size_t s = (*n) * ((*ldA) * sizeof(double complex));
    *A = (double complex*)aligned_alloc(VA, s);
    if (!*A)
      return -3;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,A,ldA)
#endif /* ?_OPENMP */
    for (fnat j = 0u; j < *n; ++j) {
      register const VD z = _mm512_setzero_pd();
      double complex *const Aj = *A + j * (size_t)(*ldA);
      for (fnat i = 0u; i < *ldA; i += VDL_2)
        _mm512_store_pd((Aj + i), z);
    }
  }

  if (dalloc2_(m, n, Ar, ldAr) < 0)
    return -5;

  if (dalloc2_(m, n, Ai, ldAi) < 0)
    return -7;

#ifdef _OPENMP
  return 1;
#else /* !_OPENMP */
  return 0;
#endif /* ?_OPENMP */
}

void czfree_(void **const A)
{
  if (A) {
    if (*A) {
      free(*A);
      *A = NULL;
    }
  }
}
