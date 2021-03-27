#include "znormx.h"

double znormx_(const fnat m[static restrict 1], const fnat n[static restrict 1], const double Ar[static restrict VDL], const double Ai[static restrict VDL], const fnat ldA[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -2.0;
  if (*m & VDL_1)
    return -1.0;
  if (IS_NOT_ALIGNED(Ar))
    return -3.0;
  if (IS_NOT_ALIGNED(Ai))
    return -4.0;
  if (*ldA < *m)
    return -5.0;
#endif /* !NDEBUG */

#ifdef _OPENMP
  double y = -HUGE_VAL;

#ifdef NDEBUG
#pragma omp parallel for default(none) shared(m,n,Ar,ldA) reduction(max:y)
#else /* !NDEBUG */
#pragma omp parallel for default(none) shared(m,n,Ar,ldA,stderr) reduction(max:y)
#endif /* ?NDEBUG */
  for (fnat j = 0u; j < *n; ++j) {
    register const VD inf = _mm512_set1_pd(HUGE_VAL);
    register VD x = _mm512_set1_pd(-HUGE_VAL);
    const double *const Aj = Ar + j * (*ldA);
    for (fnat i = 0u; i < *m; i += VDL) {
      x = _mm512_max_pd(_mm512_min_pd(_mm512_abs_pd(_mm512_load_pd(Aj + i)), inf), x);
#ifndef NDEBUG
      char s[28] = { '\0' };
      (void)sprintf(s, "r(%10u,%10u)%3u", (unsigned)i, (unsigned)j, (unsigned)omp_get_thread_num());
#pragma omp critical
      (void)VDprintf(stderr, s, x);
#endif /* !NDEBUG */
    }
    y = fmax(_mm512_reduce_max_pd(x), y);
  }

#ifdef NDEBUG
#pragma omp parallel for default(none) shared(m,n,Ai,ldA) reduction(max:y)
#else /* !NDEBUG */
#pragma omp parallel for default(none) shared(m,n,Ai,ldA,stderr) reduction(max:y)
#endif /* ?NDEBUG */
  for (fnat j = 0u; j < *n; ++j) {
    register const VD inf = _mm512_set1_pd(HUGE_VAL);
    register VD x = _mm512_set1_pd(-HUGE_VAL);
    const double *const Aj = Ai + j * (*ldA);
    for (fnat i = 0u; i < *m; i += VDL) {
      x = _mm512_max_pd(_mm512_min_pd(_mm512_abs_pd(_mm512_load_pd(Aj + i)), inf), x);
#ifndef NDEBUG
      char s[28] = { '\0' };
      (void)sprintf(s, "i(%10u,%10u)%3u", (unsigned)i, (unsigned)j, (unsigned)omp_get_thread_num());
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
    const double *const Aj = Ar + j * (*ldA);
    for (fnat i = 0u; i < *m; i += VDL) {
      x = _mm512_max_pd(_mm512_min_pd(_mm512_abs_pd(_mm512_load_pd(Aj + i)), inf), x); VDP(x);
    }
  }
  for (fnat j = 0u; j < *n; ++j) {
    const double *const Aj = Ai + j * (*ldA);
    for (fnat i = 0u; i < *m; i += VDL) {
      x = _mm512_max_pd(_mm512_min_pd(_mm512_abs_pd(_mm512_load_pd(Aj + i)), inf), x); VDP(x);
    }
  }

  return _mm512_reduce_max_pd(x);
#endif /* ?_OPENMP */
}
