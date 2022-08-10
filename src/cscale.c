#include "cscale.h"

#include "scscal.h"

fint cscale_(const fnat m[static restrict 1], const fnat n[static restrict 1], float Ar[static restrict VSL], const fnat ldAr[static restrict 1], float Ai[static restrict VSL], const fnat ldAi[static restrict 1], const fint e[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -8;
  if (*m & VSL_1)
    return -1;
  if (IS_NOT_ALIGNED(Ar))
    return -3;
  if (*ldAr < *m)
    return -4;
  if (*ldAr & VSL_1)
    return -4;
  if (IS_NOT_ALIGNED(Ai))
    return -5;
  if (*ldAi < *m)
    return -6;
  if (*ldAi & VSL_1)
    return -6;
#endif /* !NDEBUG */

  if (!*e)
    return 0;
  const float e_ = (float)*e;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,Ar,ldAr,e_)
  SCSCAL_LOOP(Ar,ldAr);
#pragma omp parallel for default(none) shared(m,n,Ai,ldAi,e_)
  SCSCAL_LOOP(Ai,ldAi);
  return 1;
#else /* !_OPENMP */
  register const VS s = _mm512_set1_ps(e_);
  SCSCAL_LOOP(Ar,ldAr);
  SCSCAL_LOOP(Ai,ldAi);
  return 0;
#endif /* ?_OPENMP */
}
