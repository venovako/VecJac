#include "sscale.h"

#include "scscal.h"

fint sscale_(const fnat m[static restrict 1], const fnat n[static restrict 1], float A[static restrict VSL], const fnat ldA[static restrict 1], const fint e[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -6;
  if (*m & VDL_1)
    return -1;
  if (IS_NOT_ALIGNED(A))
    return -3;
  if (*ldA < *m)
    return -4;
  if (*ldA & VDL_1)
    return -4;
#endif /* !NDEBUG */

  if (!*e)
    return 0;
  const float e_ = (float)*e;

#ifdef _OPENMP
  fint t = 0;

#pragma omp parallel for default(none) shared(m,n,A,ldA,e_) reduction(max:t)
  SCSCAL_LOOP(A,ldA);

  return (t + 1);
#else /* !_OPENMP */
  register const VS s = _mm512_set1_ps(e_);

  SCSCAL_LOOP(A,ldA);

  return 0;
#endif /* ?_OPENMP */
}
