#include "sdpscl.h"

float sdpscl_(const fnat m[static restrict 1], const float x[static restrict VSL], const float y[static restrict VSL], const float e[static restrict 2], const float f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*m & VSL_1)
    return NAN;
  if (IS_NOT_ALIGNED(x))
    return NAN;
  if (IS_NOT_ALIGNED(y))
    return NAN;
  if (!(e[0u] < HUGE_VALF))
    return NAN;
  if (!(e[1u] < HUGE_VALF))
    return NAN;
  if (!(f[0u] >= 1.0f) || !(f[0u] < 2.0f))
    return NAN;
  if (!(f[1u] >= 1.0f) || !(f[0u] < 2.0f))
    return NAN;
#endif /* !NDEBUG */

  if (!(e[0u] > -HUGE_VALF))
    return 0.0f;
  if (!(e[1u] > -HUGE_VALF))
    return 0.0f;

  register const VS xe = _mm512_set1_ps(-(e[0u]));
  register const VS ye = _mm512_set1_ps(-(e[1u]));
  register VS p = _mm512_setzero_ps();

  for (fnat i = 0u; i < *m; i += VSL) {
    p = _mm512_fmadd_ps(_mm512_scalef_ps(_mm512_load_ps(x + i), xe), _mm512_scalef_ps(_mm512_load_ps(y + i), ye), p); VSP(p);
  }

  return (_mm512_reduce_add_ps(p) / (f[0u] * f[1u]));
}
