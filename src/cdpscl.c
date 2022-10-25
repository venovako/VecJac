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
  if (!(e[0u] <= FLT_MAX))
    return NAN;
  if (!(e[1u] <= FLT_MAX))
    return NAN;
  if (!(f[0u] >= 1.0f) || !(f[0u] < 2.0f))
    return NAN;
  if (!(f[1u] >= 1.0f) || !(f[0u] < 2.0f))
    return NAN;
#endif /* !NDEBUG */

  const float ex = e[0u];
  if (!(ex >= -FLT_MAX))
    return 0.0f;

  const float ey = e[1u];
  if (!(ey >= -FLT_MAX))
    return 0.0f;

  register const VS xe = _mm512_set1_ps(-ex);
  register const VS ye = _mm512_set1_ps(-ey);

  register VS pr = _mm512_setzero_ps();
  register VS pi = _mm512_setzero_ps();

  for (fnat i = 0u; i < *m; i += VSL) {
    register VS xri = _mm512_load_ps(xr + i);
    register VS yri = _mm512_load_ps(yr + i);
    register VS yii = _mm512_load_ps(yi + i);
    register VS xii = _mm512_load_ps(xi + i);
    xri = _mm512_scalef_ps(xri, xe); VSP(xri);
    yri = _mm512_scalef_ps(yri, ye); VSP(yri);
    yii = _mm512_scalef_ps(yii, ye); VSP(yii);
    xii = _mm512_scalef_ps(xii, xe); VSP(xii);
    pr = _mm512_fmadd_ps(xri, yri, pr); VSP(pr);
    pi = _mm512_fmadd_ps(xri, yii, pi); VSP(pi);
    pr = _mm512_fmadd_ps(xii, yii, pr); VSP(pr);
    pi = _mm512_fnmadd_ps(xii, yri, pi); VSP(pi);
  }

  const float fx = f[0u];
  const float fy = f[1u];
  alignas(16u) float complex z[2u];
  _mm_store_ps((float*)z, _mm_div_ps(_mm_set_ps(0.0f, 1.0f, _mm512_reduce_add_ps(pi), _mm512_reduce_add_ps(pr)), _mm_set1_ps(fx * fy)));
  return *z;
}
