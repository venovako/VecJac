#include "zjrotf.h"

#include "djrotf.h"
#include "vecdef.h"

fint zjrotf_(const fint n[static restrict 1], double xr[static restrict VDL], double xi[static restrict VDL], double yr[static restrict VDL], double yi[static restrict VDL], const double c[static restrict 1], const double cat[static restrict 1], const double sat[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -9;
  if (*n & VDL_1)
    return -1;
  if (IS_NOT_ALIGNED(xr))
    return -2;
  if (IS_NOT_ALIGNED(xi))
    return -3;
  if (IS_NOT_ALIGNED(yr))
    return -4;
  if (IS_NOT_ALIGNED(yi))
    return -5;
#endif /* !NDEBUG */

  if (!*n)
    return 0;

  if (*sat == 0.0) {
    // real rotation
    if (djrotf_(n, xr, yr, c, cat))
      return -10;
    if (djrotf_(n, xi, yi, c, cat))
      return -11;
    return 0;
  }

  register const VD cta = _mm512_set1_pd(*cat);
  register const VD sta = _mm512_set1_pd(*sat);
  register const VD cna = _mm512_set1_pd(-*cat);
  register const VD c_ = _mm512_set1_pd(*c);

  register VD
    _r
#ifndef NDEBUG
    = _mm512_setzero_pd()
#endif /* !NDEBUG */
    , _i
#ifndef NDEBUG
    = _mm512_setzero_pd()
#endif /* !NDEBUG */
    ;

  if (*n < 0) { // permute
    const fnat n_ = (fnat)-*n;
    if (*c == 1.0) {
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xri = xr + i;
        double *const xii = xi + i;
        double *const yri = yr + i;
        double *const yii = yi + i;

        register const VD x_r = _mm512_load_pd(xri);
        register const VD x_i = _mm512_load_pd(xii);
        register const VD y_r = _mm512_load_pd(yri);
        register const VD y_i = _mm512_load_pd(yii);

        VZFMA(_r,_i,y_r,y_i,cta,sta,x_r,x_i);
        _mm512_store_pd(yri, _r);
        _mm512_store_pd(yii, _i);

        VZFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        _mm512_store_pd(xri, _r);
        _mm512_store_pd(xii, _i);
      }
    }
    else { // cos < 1
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xri = xr + i;
        double *const xii = xi + i;
        double *const yri = yr + i;
        double *const yii = yi + i;

        register const VD x_r = _mm512_load_pd(xri);
        register const VD x_i = _mm512_load_pd(xii);
        register const VD y_r = _mm512_load_pd(yri);
        register const VD y_i = _mm512_load_pd(yii);

        VZFMA(_r,_i,y_r,y_i,cta,sta,x_r,x_i);
        register const VD yr_ = VDMULDIV(_r, c_);
        register const VD yi_ = VDMULDIV(_i, c_);
        _mm512_store_pd(yri, yr_);
        _mm512_store_pd(yii, yi_);

        VZFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        register const VD xr_ = VDMULDIV(_r, c_);
        register const VD xi_ = VDMULDIV(_i, c_);
        _mm512_store_pd(xri, xr_);
        _mm512_store_pd(xii, xi_);
      }
    }
  }
  else { // no permute
    const fnat n_ = (fnat)*n;
    if (*c == 1.0) {
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xri = xr + i;
        double *const xii = xi + i;
        double *const yri = yr + i;
        double *const yii = yi + i;

        register const VD x_r = _mm512_load_pd(xri);
        register const VD x_i = _mm512_load_pd(xii);
        register const VD y_r = _mm512_load_pd(yri);
        register const VD y_i = _mm512_load_pd(yii);

        VZFMA(_r,_i,y_r,y_i,cta,sta,x_r,x_i);
        _mm512_store_pd(xri, _r);
        _mm512_store_pd(xii, _i);

        VZFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        _mm512_store_pd(yri, _r);
        _mm512_store_pd(yii, _i);
      }
    }
    else { // cos < 1
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xri = xr + i;
        double *const xii = xi + i;
        double *const yri = yr + i;
        double *const yii = yi + i;

        register const VD x_r = _mm512_load_pd(xri);
        register const VD x_i = _mm512_load_pd(xii);
        register const VD y_r = _mm512_load_pd(yri);
        register const VD y_i = _mm512_load_pd(yii);

        VZFMA(_r,_i,y_r,y_i,cta,sta,x_r,x_i);
        register const VD xr_ = VDMULDIV(_r, c_);
        register const VD xi_ = VDMULDIV(_i, c_);
        _mm512_store_pd(xri, xr_);
        _mm512_store_pd(xii, xi_);

        VZFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        register const VD yr_ = VDMULDIV(_r, c_);
        register const VD yi_ = VDMULDIV(_i, c_);
        _mm512_store_pd(yri, yr_);
        _mm512_store_pd(yii, yi_);
      }
    }
  }

  return 0;
}
