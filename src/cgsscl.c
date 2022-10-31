#include "cgsscl.h"

#include "vecdef.h"

float cgsscl_(const fint m[static restrict 1], const float tr[static restrict 1], const float ti[static restrict 1], float xr[static restrict VSL], float xi[static restrict VSL], float yr[static restrict VSL], float yi[static restrict VSL], const float e[static restrict 2], const float f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*m & VSL_1)
    return NAN;
  if (!(fabsf(*tr) <= FLT_MAX))
    return NAN;
  if (!(fabsf(*ti) <= FLT_MAX))
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
    const float fx_y = (fx / fy);
    // -t * fx_y
    register const VS tfr = _mm512_set1_ps(-*tr * fx_y);
    register const VS tfi = _mm512_set1_ps(-*ti * fx_y);
    register const VS xe = _mm512_set1_ps(ex);
    const fnat _m = (fnat)-*m;

    for (fnat i = 0u; i < _m; i += VSL) {
      float *const xr_i = xr + i;
      float *const xi_i = xi + i;
      float *const yr_i = yr + i;
      float *const yi_i = yi + i;
      register VS yri = _mm512_load_ps(yr_i); VSP(yri);
      register VS yii = _mm512_load_ps(yi_i); VSP(yii);
      register const VS xri = _mm512_load_ps(xr_i); VSP(xri);
      register const VS xii = _mm512_load_ps(xi_i); VSP(xii);
      _mm512_store_ps(xr_i, yri);
      _mm512_store_ps(xi_i, yii);
      register VS yri_ = _mm512_scalef_ps(yri, y_e); VSP(yri_);
      register VS yii_ = _mm512_scalef_ps(yii, y_e); VSP(yii_);
      register VS xri_ = _mm512_scalef_ps(xri, x_e); VSP(xri_);
      register VS xii_ = _mm512_scalef_ps(xii, x_e); VSP(xii_);
      VCFMA(yri,yii,tfr,tfi,yri_,yii_,xri_,xii_);
      xri_ = _mm512_scalef_ps(yri, xe); VSP(xri_);
      xii_ = _mm512_scalef_ps(yii, xe); VSP(xii_);
      _mm512_store_ps(yr_i, xri_);
      _mm512_store_ps(yi_i, xii_);
      xri_ = VSABS(xri_);
      xii_ = VSABS(xii_);
      M = _mm512_max_ps(M, xri_);
      M = _mm512_max_ps(M, xii_);
    }
  }
  else { // transform y
    const float fy_x = (fy / fx);
    // -conjf(t) * fy_x
    register const VS tfr = _mm512_set1_ps(-*tr * fy_x);
    register const VS tfi = _mm512_set1_ps( *ti * fy_x);
    register const VS ye = _mm512_set1_ps(ey);
    const fnat _m = (fnat)*m;

    for (fnat i = 0u; i < _m; i += VSL) {
      float *const yr_i = yr + i;
      float *const yi_i = yi + i;
      register VS xri = _mm512_load_ps(xr + i); VSP(xri);
      register VS xii = _mm512_load_ps(xi + i); VSP(xii);
      register VS xri_ = _mm512_scalef_ps(xri, x_e); VSP(xri_);
      register VS xii_ = _mm512_scalef_ps(xii, x_e); VSP(xii_);
      register const VS yri = _mm512_load_ps(yr_i); VSP(yri);
      register const VS yii = _mm512_load_ps(yi_i); VSP(yii);
      register VS yri_ = _mm512_scalef_ps(yri, y_e); VSP(yri_);
      register VS yii_ = _mm512_scalef_ps(yii, y_e); VSP(yii_);
      VCFMA(xri,xii,tfr,tfi,xri_,xii_,yri_,yii_);
      yri_ = _mm512_scalef_ps(xri, ye); VSP(yri_);
      yii_ = _mm512_scalef_ps(xii, ye); VSP(yii_);
      _mm512_store_ps(yr_i, yri_);
      _mm512_store_ps(yi_i, yii_);
      yri_ = VSABS(yri_);
      yii_ = VSABS(yii_);
      M = _mm512_max_ps(M, yri_);
      M = _mm512_max_ps(M, yii_);
    }
  }

  return _mm512_reduce_max_ps(M);
}
