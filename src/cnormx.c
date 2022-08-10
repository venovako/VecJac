#include "cnormx.h"

#include "scnrmx.h"

float cnormx_(const fnat m[static restrict 1], const fnat n[static restrict 1], const float Ar[static restrict VSL], const fnat ldAr[static restrict 1], const float Ai[static restrict VSL], const fnat ldAi[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -7.0f;
  if (*m & VSL_1)
    return -1.0f;
  if (IS_NOT_ALIGNED(Ar))
    return -3.0f;
  if (*ldAr < *m)
    return -4.0f;
  if (*ldAr & VSL_1)
    return -4.5f;
  if (IS_NOT_ALIGNED(Ai))
    return -5.0f;
  if (*ldAi < *m)
    return -6.0f;
  if (*ldAi & VSL_1)
    return -6.5f;
#endif /* !NDEBUG */

#ifdef _OPENMP
  float y = 0.0f;

#pragma omp parallel for default(none) shared(m,n,Ar,ldAr) reduction(max:y)
  SCNRMX_LOOP(Ar,ldAr);
#pragma omp parallel for default(none) shared(m,n,Ai,ldAi) reduction(max:y)
  SCNRMX_LOOP(Ai,ldAi);

  return y;
#else /* !_OPENMP */
  register const VS _zerof = _mm512_set1_ps(-0.0f);
  register const VS inff = _mm512_set1_ps(HUGE_VALF);
  register VS x = _mm512_setzero_ps();

  SCNRMX_LOOP(Ar,ldAr);
  SCNRMX_LOOP(Ai,ldAi);

  return _mm512_reduce_max_ps(x);
#endif /* ?_OPENMP */
}
