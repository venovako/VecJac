#include "dscale.h"

#include "dnormx.h"
#include "dzscal.h"

fint dscale_(const fnat m[static restrict 1], const fnat n[static restrict 1], double A[static restrict VDL], const fnat ldA[static restrict 1], const fint e[static restrict 1])
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
  const double e_ = (double)*e;

#ifdef _OPENMP
  fint t = 0;

#pragma omp parallel for default(none) shared(m,n,A,ldA,e_) reduction(max:t)
  DZSCAL_LOOP(A,ldA);

  return (t + 1);
#else /* !_OPENMP */
  register const VD s = _mm512_set1_pd(e_);

  DZSCAL_LOOP(A,ldA);

  return 0;
#endif /* ?_OPENMP */
}

fint dlscal_(const fnat m[static restrict 1], const fnat n[static restrict 1], double A[static restrict VDL], const fnat ldA[static restrict 1], const fnat l[static restrict 1], double M[static restrict 1], fint e[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -8;
  if (*m & VDL_1)
    return -1;
  if (IS_NOT_ALIGNED(A))
    return -3;
  if (*ldA < *m)
    return -4;
  if (*ldA & VDL_1)
    return -4;
#endif /* !NDEBUG */

  *M = dnormx_(m, n, A, ldA);
  *e = (isfinite(*M) ? s(*M, *l) : (fint)0);
  return dscale_(m, n, A, ldA, e);
}
