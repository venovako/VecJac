#include "dnormx.h"

double dnormx_(const fnat m[static restrict 1], const fnat n[static restrict 1], const double A[static restrict VDL], const fnat ldA[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -2.0;
  if (*m & VDL_1)
    return -1.0;
  if (IS_NOT_ALIGNED(A))
    return -3.0;
  if (*ldA < *m)
    return -4.0;
#endif /* !NDEBUG */

#ifdef _OPENMP
  double y = -HUGE_VAL;
#ifdef NDEBUG
#pragma omp parallel for default(none) shared(m,n,A,ldA) reduction(max:y)
#else /* !NDEBUG */
#pragma omp parallel for default(none) shared(m,n,A,ldA,stderr) reduction(max:y)
#endif /* ?NDEBUG */
  for (fnat j = 0u; j < *n; ++j) {
    register const VD inf = _mm512_set1_pd(HUGE_VAL);
    register VD x = _mm512_set1_pd(-HUGE_VAL);
    const double *const Aj = A + j * (*ldA);
    for (fnat i = 0u; i < *m; i += VDL) {
      x = _mm512_max_pd(_mm512_min_pd(_mm512_abs_pd(_mm512_load_pd(Aj + i)), inf), x);
#ifndef NDEBUG
      char s[28] = { '\0' };
      (void)sprintf(s, "x(%10u,%10u)%3u", (unsigned)i, (unsigned)j, (unsigned)omp_get_thread_num());
#pragma omp critical
      (void)VDprintf(stderr, s, x);
#endif /* !NDEBUG */
    }
    y = fmax(_mm512_reduce_max_pd(x), y);
  }

  return y;
#else /* !_OPENMP */  
  register const VD inf = _mm512_set1_pd(HUGE_VAL);
  register VD x = _mm512_set1_pd(-HUGE_VAL);

  for (fnat j = 0u; j < *n; ++j) {
    const double *const Aj = A + j * (*ldA);
    for (fnat i = 0u; i < *m; i += VDL) {
      x = _mm512_max_pd(_mm512_min_pd(_mm512_abs_pd(_mm512_load_pd(Aj + i)), inf), x); VP(x);
    }
  }

  return _mm512_reduce_max_pd(x);
#endif /* ?_OPENMP */
}
