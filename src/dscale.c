#include "dscale.h"

#include "dnormx.h"
#include "dzscal.h"

int dscale_(const fnat m[static restrict 1], const fnat n[static restrict 1], double A[static restrict VDL], const fnat ldA[static restrict 1], const fint e[static restrict 1])
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

  if (!*e)
    return 0;
  const double e_ = (double)*e;

#ifdef _OPENMP
  int t = 0;

#pragma omp parallel for default(none) shared(m,n,A,ldA,e_) reduction(max:t)
  DZSCAL_LOOP(A,ldA);

  return (t + 1);
#else /* !_OPENMP */
  register const VD s = _mm512_set1_pd(e_);

  DZSCAL_LOOP(A,ldA);

  return 0;
#endif /* ?_OPENMP */
}

#ifdef EDOMEGA
#error EDOMEGA already defined
#else /* !EDOMEGA */
#define EDOMEGA 1024
#endif /* ?EDOMEGA */

static inline fint s(const double M, const fint l)
{
  int e
#ifndef NDEBUG
    = 0
#endif /* !NDEBUG */
    ;
  (void)frexp(M, &e);
  return (EDOMEGA - l - e);
}

int dlscal_(const fnat m[static restrict 1], const fnat n[static restrict 1], double A[static restrict VDL], const fnat ldA[static restrict 1], const fnat l[static restrict 1], double M[static restrict 1], fint e[static restrict 1])
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

  *M = dnormx_(m, n, A, ldA);
  *e = s(*M, *l);
  return dscale_(m, n, A, ldA, e);
}
