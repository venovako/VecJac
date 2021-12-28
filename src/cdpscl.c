#include "cdpscl.h"

float complex cdpscl_(const fnat m[static restrict 1], const float xr[static restrict VSL], const float xi[static restrict VSL], const float yr[static restrict VSL], const float yi[static restrict VSL], const float e[static restrict 2], const float f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*m & VSL_1)
    return NAN;
  if (IS_NOT_ALIGNED(xr))
    return NAN;
  if (IS_NOT_ALIGNED(xi))
    return NAN;
  if (IS_NOT_ALIGNED(yr))
    return NAN;
  if (IS_NOT_ALIGNED(yi))
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

  register VS pr = _mm512_setzero_ps();
  for (fnat i = 0u; i < *m; i += VSL) {
    pr = _mm512_fmadd_ps(_mm512_scalef_ps(_mm512_load_ps(xr + i), xe), _mm512_scalef_ps(_mm512_load_ps(yr + i), ye), pr); VSP(pr);
  }
  for (fnat i = 0u; i < *m; i += VSL) {
    pr = _mm512_fmadd_ps(_mm512_scalef_ps(_mm512_load_ps(xi + i), xe), _mm512_scalef_ps(_mm512_load_ps(yi + i), ye), pr); VSP(pr);
  }

  register VS pi = _mm512_setzero_ps();
  for (fnat i = 0u; i < *m; i += VSL) {
    pi = _mm512_fmadd_ps(_mm512_scalef_ps(_mm512_load_ps(xr + i), xe), _mm512_scalef_ps(_mm512_load_ps(yi + i), ye), pi); VSP(pi);
  }
  for (fnat i = 0u; i < *m; i += VSL) {
    pi = _mm512_fnmadd_ps(_mm512_scalef_ps(_mm512_load_ps(xi + i), xe), _mm512_scalef_ps(_mm512_load_ps(yr + i), ye), pi); VSP(pi);
  }

  const float ff = (f[0u] * f[1u]);
  const float r = (_mm512_reduce_add_ps(pr) / ff);
  const float i = (_mm512_reduce_add_ps(pi) / ff);
  const float complex z = CMPLXF(r, i);
  return z;
}
