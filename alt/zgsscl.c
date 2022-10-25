#include "zgsscl.h"

#include "vecdef.h"

double zgsscl_(const fint m[static restrict 1], const double tr[static restrict 1], const double ti[static restrict 1], double xr[static restrict VDL], double xi[static restrict VDL], double yr[static restrict VDL], double yi[static restrict VDL], const double e[static restrict 2], const double f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*m & VDL_1)
    return NAN;
  if (!(fabs(*tr) <= DBL_MAX))
    return NAN;
  if (!(fabs(*ti) <= DBL_MAX))
    return NAN;
  if (IS_NOT_ALIGNED(xr))
    return NAN;
  if (IS_NOT_ALIGNED(xi))
    return NAN;
  if (IS_NOT_ALIGNED(yr))
    return NAN;
  if (IS_NOT_ALIGNED(yi))
    return NAN;
  if (!(e[0u] <= DBL_MAX))
    return NAN;
  if (!(e[1u] <= DBL_MAX))
    return NAN;
  if (!(f[0u] >= 1.0) || !(f[0u] < 2.0))
    return NAN;
  if (!(f[1u] >= 1.0) || !(f[0u] < 2.0))
    return NAN;
#endif /* !NDEBUG */

  if (!*m)
    return -0.0;

  const double ex = e[0u];
  if (!(ex >= -DBL_MAX))
    return -0.0;

  const double ey = e[1u];
  if (!(ey >= -DBL_MAX))
    return -0.0;

  const double fx = f[0u];
  const double fy = f[1u];

  register const VD x_e = _mm512_set1_pd(-ex);
  register const VD y_e = _mm512_set1_pd(-ey);
  register const VD _zero = _mm512_set1_pd(-0.0);
  register VD M = _mm512_setzero_pd();
  register MD ne = MDXOR(ne, ne);

  if (*m < 0) { // transform x and permute
    const double fx_y = (fx / fy);
    // -t * fx_y
    register const VD tfr = _mm512_set1_pd(-*tr * fx_y);
    register const VD tfi = _mm512_set1_pd(-*ti * fx_y);
    register const VD xe = _mm512_set1_pd(ex);
    const fnat _m = (fnat)-*m;

    for (fnat i = 0u; i < _m; i += VDL) {
      double *const xr_i = xr + i;
      double *const xi_i = xi + i;
      double *const yr_i = yr + i;
      double *const yi_i = yi + i;
      register VD yri = _mm512_load_pd(yr_i); VDP(yri);
      register VD yii = _mm512_load_pd(yi_i); VDP(yii);
      register const VD xri = _mm512_load_pd(xr_i); VDP(xri);
      register const VD xii = _mm512_load_pd(xi_i); VDP(xii);
      _mm512_store_pd(xr_i, yri);
      _mm512_store_pd(xi_i, yii);
      register VD yri_ = _mm512_scalef_pd(yri, y_e); VDP(yri_);
      register VD yii_ = _mm512_scalef_pd(yii, y_e); VDP(yii_);
      register VD xri_ = _mm512_scalef_pd(xri, x_e); VDP(xri_);
      register VD xii_ = _mm512_scalef_pd(xii, x_e); VDP(xii_);
      VDFMA(yri,yii,tfr,tfi,yri_,yii_,xri_,xii_);
      xri_ = _mm512_scalef_pd(yri, xe); VDP(xri_);
      xii_ = _mm512_scalef_pd(yii, xe); VDP(xii_);
      ne = MDOR(ne, _mm512_cmpneq_pd_mask(xri, xri_));
      ne = MDOR(ne, _mm512_cmpneq_pd_mask(xii, xii_));
      _mm512_store_pd(yr_i, xri_);
      _mm512_store_pd(yi_i, xii_);
      xri_ = VDABS(xri_);
      xii_ = VDABS(xii_);
      M = _mm512_max_pd(M, xri_);
      M = _mm512_max_pd(M, xii_);
    }
  }
  else { // transform y
    const double fy_x = (fy / fx);
    // -conj(t) * fy_x
    register const VD tfr = _mm512_set1_pd(-*tr * fy_x);
    register const VD tfi = _mm512_set1_pd( *ti * fy_x);
    register const VD ye = _mm512_set1_pd(ey);
    const fnat _m = (fnat)*m;

    for (fnat i = 0u; i < _m; i += VDL) {
      double *const yr_i = yr + i;
      double *const yi_i = yi + i;
      register VD xri = _mm512_load_pd(xr + i); VDP(xri);
      register VD xii = _mm512_load_pd(xi + i); VDP(xii);
      register VD xri_ = _mm512_scalef_pd(xri, x_e); VDP(xri_);
      register VD xii_ = _mm512_scalef_pd(xii, x_e); VDP(xii_);
      register const VD yri = _mm512_load_pd(yr_i); VDP(yri);
      register const VD yii = _mm512_load_pd(yi_i); VDP(yii);
      register VD yri_ = _mm512_scalef_pd(yri, y_e); VDP(yri_);
      register VD yii_ = _mm512_scalef_pd(yii, y_e); VDP(yii_);
      VDFMA(xri,xii,tfr,tfi,xri_,xii_,yri_,yii_);
      yri_ = _mm512_scalef_pd(xri, ye); VDP(yri_);
      yii_ = _mm512_scalef_pd(xii, ye); VDP(yii_);
      ne = MDOR(ne, _mm512_cmpneq_pd_mask(yri, yri_));
      ne = MDOR(ne, _mm512_cmpneq_pd_mask(yii, yii_));
      _mm512_store_pd(yr_i, yri_);
      _mm512_store_pd(yi_i, yii_);
      yri_ = VDABS(yri_);
      yii_ = VDABS(yii_);
      M = _mm512_max_pd(M, yri_);
      M = _mm512_max_pd(M, yii_);
    }
  }

  const double r = _mm512_reduce_max_pd(M);
  return (MD2U(ne) ? r : -r);
}
