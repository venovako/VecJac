#include "sgsscl.h"

#include "vecdef.h"

float sgsscl_(const fint m[static restrict 1], const float t[static restrict 1], float x[static restrict VSL], float y[static restrict VSL], const float e[static restrict 2], const float f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -7.0f;
  if (*m & VSL_1)
    return -1.0f;
  if (!(fabsf(*t) <= FLT_MAX))
    return -2.0f;
  if (IS_NOT_ALIGNED(x))
    return -3.0f;
  if (IS_NOT_ALIGNED(y))
    return -4.0f;
  if (!(e[0u] <= FLT_MAX))
    return -5.0f;
  if (!(e[1u] <= FLT_MAX))
    return -5.5f;
  if (!(f[0u] >= 1.0f) || !(f[0u] < 2.0f))
    return -6.0f;
  if (!(f[1u] >= 1.0f) || !(f[0u] < 2.0f))
    return -6.5f;
#endif /* !NDEBUG */

  if (!*m)
    return 0.0f;

  const float ex = e[0u];
  if (!(ex >= -FLT_MAX))
    return 0.0f;

  const float ey = e[1u];
  if (!(ey >= -FLT_MAX))
    return 0.0f;

  const float fx = f[0u];
  const float fy = f[1u];

  register const VS x_e = _mm512_set1_ps(-ex);
  register const VS y_e = _mm512_set1_ps(-ey);
  register const VS _zerof = _mm512_set1_ps(-0.0f);
  register VS M = _mm512_setzero_ps();

  if (*m < 0) { // transform x and permute
    register const VS tf = _mm512_set1_ps(-*t * (fx / fy));
    register const VS xe = _mm512_set1_ps(ex);
    const fnat _m = (fnat)-*m;

    for (fnat i = 0u; i < _m; i += VSL) {
      float *const x_i = x + i;
      float *const y_i = y + i;
      register VS xi = _mm512_load_ps(x_i);
      register VS yi = _mm512_load_ps(y_i);
      xi = _mm512_scalef_ps(xi, x_e); VSP(xi);
      _mm512_store_ps(x_i, yi);
      yi = _mm512_scalef_ps(yi, y_e); VSP(yi);
      xi = _mm512_scalef_ps(_mm512_fmadd_ps(tf, yi, xi), xe); VSP(xi);
      yi = VSABS(xi); VSP(yi);
      _mm512_store_ps(y_i, xi);
      M = _mm512_max_ps(M, yi); VSP(M);
    }
  }
  else { // transform y
    register const VS tf = _mm512_set1_ps(-*t * (fy / fx));
    register const VS ye = _mm512_set1_ps(ey);
    const fnat _m = (fnat)*m;

    for (fnat i = 0u; i < _m; i += VSL) {
      float *const y_i = y + i;
      register VS xi = _mm512_load_ps(x + i);
      register VS yi = _mm512_load_ps(y_i);
      xi = _mm512_scalef_ps(xi, x_e); VSP(xi);
      yi = _mm512_scalef_ps(yi, y_e); VSP(yi);
      yi = _mm512_scalef_ps(_mm512_fmadd_ps(tf, xi, yi), ye); VSP(yi);
      xi = VSABS(yi); VSP(xi);
      _mm512_store_ps(y_i, yi);
      M = _mm512_max_ps(M, xi); VSP(M);
    }
  }

  return _mm512_reduce_max_ps(M);
}
