#include "snormx.h"

#include "scnrmx.h"

float snormx_(const fnat m[static restrict 1], const fnat n[static restrict 1], const float A[static restrict VSL], const fnat ldA[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -5.0f;
  if (*m & VSL_1)
    return -1.0f;
  if (IS_NOT_ALIGNED(A))
    return -3.0f;
  if (*ldA < *m)
    return -4.0f;
  if (*ldA & VSL_1)
    return -4.5f;
#endif /* !NDEBUG */

#ifdef _OPENMP
  float y = 0.0f;

#pragma omp parallel for default(none) shared(m,n,A,ldA) reduction(max:y)
  SCNRMX_LOOP(A,ldA);

  return y;
#else /* !_OPENMP */
  register const VS _zerof = _mm512_set1_ps(-0.0f);
  register const VS inff = _mm512_set1_ps(HUGE_VALF);
  register VS x = _mm512_setzero_ps();

  SCNRMX_LOOP(A,ldA);

  return _mm512_reduce_max_ps(x);
#endif /* ?_OPENMP */
}
