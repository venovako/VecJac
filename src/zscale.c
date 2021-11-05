#include "zscale.h"

#include "znormx.h"
#include "dzscal.h"

fint zscale_(const fnat m[static restrict 1], const fnat n[static restrict 1], double Ar[static restrict VDL], const fnat ldAr[static restrict 1], double Ai[static restrict VDL], const fnat ldAi[static restrict 1], const fint e[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -8;
  if (*m & VDL_1)
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
#endif /* !NDEBUG */

  if (!*e)
    return 0;
  const double e_ = (double)*e;

#ifdef _OPENMP
  fint t = 0;

#pragma omp parallel for default(none) shared(m,n,Ar,ldAr,e_) reduction(max:t)
  DZSCAL_LOOP(Ar,ldAr);
#pragma omp parallel for default(none) shared(m,n,Ai,ldAi,e_) reduction(max:t)
  DZSCAL_LOOP(Ai,ldAi);

  return (t + 1);
#else /* !_OPENMP */
  register const VD s = _mm512_set1_pd(e_);

  DZSCAL_LOOP(Ar,ldAr);
  DZSCAL_LOOP(Ai,ldAi);

  return 0;
#endif /* ?_OPENMP */
}

fint zlscal_(const fnat m[static restrict 1], const fnat n[static restrict 1], double Ar[static restrict VDL], const fnat ldAr[static restrict 1], double Ai[static restrict VDL], const fnat ldAi[static restrict 1], const fnat l[static restrict 1], double M[static restrict 1], fint e[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -10;
  if (*m & VDL_1)
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
#endif /* !NDEBUG */

  *M = znormx_(m, n, Ar, ldAr, Ai, ldAi);
  *e = (isfinite(*M) ? s(*M, *l) : (fint)0);
  return zscale_(m, n, Ar, ldAr, Ai, ldAi, e);
}
