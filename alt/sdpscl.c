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

  const float ex = e[0u];
  if (!(ex > -HUGE_VALF))
    return 0.0f;

  const float ey = e[1u];
  if (!(ey > -HUGE_VALF))
    return 0.0f;

  register const VS xe = _mm512_set1_ps(-ex);
  register const VS ye = _mm512_set1_ps(-ey);
  register VS p = _mm512_setzero_ps();
  register VS c = _mm512_setzero_ps();
  register VS d = _mm512_setzero_ps();
  register VS s = _mm512_setzero_ps();

  for (fnat i = 0u; i < *m; i += VSL) {
    register VS xi = _mm512_load_ps(x + i);
    register VS yi = _mm512_load_ps(y + i);
    xi = _mm512_scalef_ps(xi, xe); VSP(xi);
    yi = _mm512_scalef_ps(yi, ye); VSP(yi);
    c = _mm512_mul_round_ps(xi, yi, _MM_FROUND_TO_NEG_INF); VSP(c);
    d = _mm512_fmsub_ps(xi, yi, c); VSP(d);
    p = _mm512_add_ps(p, c); VSP(p);
    s = _mm512_add_ps(s, d); VSP(s);
  }

  const float fx = f[0u];
  const float fy = f[1u];
  return ((_mm512_reduce_add_ps(p) + _mm512_reduce_add_ps(s)) / (fx * fy));
}
