#include "dnormx.h"

#include "dznrmx.h"

double dnormx_(const fnat m[static restrict 1], const fnat n[static restrict 1], const double A[static restrict VDL], const fnat ldA[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -5.0;
  if (*m & VDL_1)
    return -1.0;
  if (IS_NOT_ALIGNED(A))
    return -3.0;
  if (*ldA < *m)
    return -4.0;
  if (*ldA & VDL_1)
    return -4.5;
#endif /* !NDEBUG */

#ifdef _OPENMP
  double y = 0.0;

#pragma omp parallel for default(none) shared(m,n,A,ldA) reduction(max:y)
  DZNRMX_LOOP(A,ldA);

  return y;
#else /* !_OPENMP */
  register const VD _zero = _mm512_set1_pd(-0.0);
  register const VD inf = _mm512_set1_pd(HUGE_VAL);
  register VD x = _mm512_setzero_pd();

  DZNRMX_LOOP(A,ldA);

  return _mm512_reduce_max_pd(x);
#endif /* ?_OPENMP */
}
