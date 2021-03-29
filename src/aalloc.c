#include "aalloc.h"

void salloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], float *A[static restrict 1], fnat ldA[static restrict 1])
{
  *A = (float*)NULL;
  *ldA = 0u;
  if (!*m)
    return;
  if (!*n)
    return;

  const fnat k = *m & VSL_1;
  *ldA = (k ? (*m + (VSL - k)) : *m);
  const size_t s = (*n) * ((*ldA) * sizeof(float));
  *A = (float*)aligned_alloc(VA, s);
  if (!*A)
    return;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,A,ldA)
  for (fnat j = 0u; j < *n; ++j) {
    register const VS z = _mm512_setzero_ps();
    float *const Aj = *A + j * (*ldA);
    for (fnat i = 0u; i < *ldA; i += VSL)
      _mm512_store_ps((Aj + i), z);
  }
#else /* !_OPENMP */
  (void)memset(*A, 0, s);
#endif /* ?_OPENMP */
}

void dalloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], double *A[static restrict 1], fnat ldA[static restrict 1])
{
  *A = (double*)NULL;
  *ldA = 0u;
  if (!*m)
    return;
  if (!*n)
    return;

  const fnat k = *m & VDL_1;
  *ldA = (k ? (*m + (VDL - k)) : *m);
  const size_t s = (*n) * ((*ldA) * sizeof(double));
  *A = (double*)aligned_alloc(VA, s);
  if (!*A)
    return;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,A,ldA)
  for (fnat j = 0u; j < *n; ++j) {
    register const VD z = _mm512_setzero_pd();
    double *const Aj = *A + j * (*ldA);
    for (fnat i = 0u; i < *ldA; i += VDL)
      _mm512_store_pd((Aj + i), z);
  }
#else /* !_OPENMP */
  (void)memset(*A, 0, s);
#endif /* ?_OPENMP */
}

void calloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], float complex *A[static restrict 1], fnat ldA[static restrict 1], float *Ar[static restrict 1], fnat ldAr[static restrict 1], float *Ai[static restrict 1], fnat ldAi[static restrict 1])
{
  *A = (float complex*)NULL;
  *ldA = 0u;
  *Ar = (float*)NULL;
  *ldAr = 0u;
  *Ai = (float*)NULL;
  *ldAi = 0u;
  if (!*m)
    return;
  if (!*n)
    return;

  const fnat k = *m & VSL__2;
  *ldA = (k ? (*m + (VSL_2 - k)) : *m);
  const size_t s = (*n) * ((*ldA) * sizeof(float complex));
  *A = (float complex*)aligned_alloc(VA, s);
  if (!*A)
    return;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,A,ldA)
  for (fnat j = 0u; j < *n; ++j) {
    register const VS z = _mm512_setzero_ps();
    float complex *const Aj = *A + j * (*ldA);
    for (fnat i = 0u; i < *ldA; i += VSL_2)
      _mm512_store_ps((Aj + i), z);
  }
#else /* !_OPENMP */
  (void)memset(*A, 0, s);
#endif /* ?_OPENMP */

  salloc2_(m, n, Ar, ldAr);
  if (!*Ar)
    return;
  salloc2_(m, n, Ai, ldAi);
}

void zalloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], double complex *A[static restrict 1], fnat ldA[static restrict 1], double *Ar[static restrict 1], fnat ldAr[static restrict 1], double *Ai[static restrict 1], fnat ldAi[static restrict 1])
{
  *A = (double complex*)NULL;
  *ldA = 0u;
  *Ar = (double*)NULL;
  *ldAr = 0u;
  *Ai = (double*)NULL;
  *ldAi = 0u;
  if (!*m)
    return;
  if (!*n)
    return;

  const fnat k = *m & VDL__2;
  *ldA = (k ? (*m + (VDL_2 - k)) : *m);
  const size_t s = (*n) * ((*ldA) * sizeof(double complex));
  *A = (double complex*)aligned_alloc(VA, s);
  if (!*A)
    return;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,A,ldA)
  for (fnat j = 0u; j < *n; ++j) {
    register const VD z = _mm512_setzero_pd();
    double complex *const Aj = *A + j * (*ldA);
    for (fnat i = 0u; i < *ldA; i += VDL_2)
      _mm512_store_pd((Aj + i), z);
  }
#else /* !_OPENMP */
  (void)memset(*A, 0, s);
#endif /* ?_OPENMP */

  dalloc2_(m, n, Ar, ldAr);
  if (!*Ar)
    return;
  dalloc2_(m, n, Ai, ldAi);
}
