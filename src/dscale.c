#include "dscale.h"

#include "dzscal.h"

int dscale_(const fnat m[static restrict 1], const fnat n[static restrict 1], double A[static restrict VDL], const fnat ldA[static restrict 1], const double e[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -2;
  if (*m & VDL_1)
    return -1;
  if (IS_NOT_ALIGNED(A))
    return -3;
  if (*ldA < *m)
    return -4;
  if (*ldA & VDL_1)
    return -4;
#endif /* !NDEBUG */

#ifdef _OPENMP
  int t = 0;

#pragma omp parallel for default(none) shared(m,n,A,ldA,e) reduction(max:t)
  DZSCAL_LOOP(A,ldA);

  return (t + 1);
#else /* !_OPENMP */
  register const VD s = _mm512_set1_pd(*e);

  DZSCAL_LOOP(A,ldA);

  return 0;
#endif /* ?_OPENMP */
}
