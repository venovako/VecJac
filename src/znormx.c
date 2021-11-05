#include "znormx.h"

#include "dznrmx.h"

double znormx_(const fnat m[static restrict 1], const fnat n[static restrict 1], const double Ar[static restrict VDL], const fnat ldAr[static restrict 1], const double Ai[static restrict VDL], const fnat ldAi[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -7.0;
  if (*m & VDL_1)
    return -1.0;
  if (IS_NOT_ALIGNED(Ar))
    return -3.0;
  if (*ldAr < *m)
    return -4.0;
  if (*ldAr & VDL_1)
    return -4.5;
  if (IS_NOT_ALIGNED(Ai))
    return -5.0;
  if (*ldAi < *m)
    return -6.0;
  if (*ldAi & VDL_1)
    return -6.5;
#endif /* !NDEBUG */

#ifdef _OPENMP
  double y = -HUGE_VAL;

#pragma omp parallel for default(none) shared(m,n,Ar,ldAr) reduction(max:y)
  DZNRMX_LOOP(Ar,ldAr);
#pragma omp parallel for default(none) shared(m,n,Ai,ldAi) reduction(max:y)
  DZNRMX_LOOP(Ai,ldAi);

  return y;
#else /* !_OPENMP */  
  register const VD inf = _mm512_set1_pd(HUGE_VAL);
  register VD x = _mm512_set1_pd(-HUGE_VAL);

  DZNRMX_LOOP(Ar,ldAr);
  DZNRMX_LOOP(Ai,ldAi);

  return _mm512_reduce_max_pd(x);
#endif /* ?_OPENMP */
}
