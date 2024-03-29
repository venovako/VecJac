#include "zjrot.h"

#include "djrot.h"
#include "vecdef.h"

double zjrot_(const fint n[static restrict 1], double xr[static restrict VDL], double xi[static restrict VDL], double yr[static restrict VDL], double yi[static restrict VDL], const double c[static restrict 1], const double cat[static restrict 1], const double sat[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -9.0;
  if (*n & VDL_1)
    return -1.0;
  if (IS_NOT_ALIGNED(xr))
    return -2.0;
  if (IS_NOT_ALIGNED(xi))
    return -3.0;
  if (IS_NOT_ALIGNED(yr))
    return -4.0;
  if (IS_NOT_ALIGNED(yi))
    return -5.0;
#endif /* !NDEBUG */

  if (!*n)
    return 0.0;

  if (*sat == 0.0) {
    // real rotation
    const double r_r = djrot_(n, xr, yr, c, cat);
    if (!(r_r >= 0.0))
      return -10.0;
    const double r_i = djrot_(n, xi, yi, c, cat);
    if (!(r_i >= 0.0))
      return -11.0;
    return fmax(r_r, r_i);
  }

  register const VD cta = _mm512_set1_pd(*cat);
  register const VD sta = _mm512_set1_pd(*sat);
  register const VD cna = _mm512_set1_pd(-*cat);
  register const VD c_ = _mm512_set1_pd(*c);

  register const VD _zero = _mm512_set1_pd(-0.0);
  register VD mx = _mm512_setzero_pd()
    , _r
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

        _r = VDABS(_r);
        _i = VDABS(_i);
        mx = _mm512_max_pd(mx, _mm512_max_pd(_r, _i));

        VZFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        _mm512_store_pd(xri, _r);
        _mm512_store_pd(xii, _i);

        _r = VDABS(_r);
        _i = VDABS(_i);
        mx = _mm512_max_pd(mx, _mm512_max_pd(_r, _i));
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

        _r = VDABS(yr_);
        _i = VDABS(yi_);
        mx = _mm512_max_pd(mx, _mm512_max_pd(_r, _i));

        VZFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        register const VD xr_ = VDMULDIV(_r, c_);
        register const VD xi_ = VDMULDIV(_i, c_);
        _mm512_store_pd(xri, xr_);
        _mm512_store_pd(xii, xi_);

        _r = VDABS(xr_);
        _i = VDABS(xi_);
        mx = _mm512_max_pd(mx, _mm512_max_pd(_r, _i));
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

        _r = VDABS(_r);
        _i = VDABS(_i);
        mx = _mm512_max_pd(mx, _mm512_max_pd(_r, _i));

        VZFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        _mm512_store_pd(yri, _r);
        _mm512_store_pd(yii, _i);

        _r = VDABS(_r);
        _i = VDABS(_i);
        mx = _mm512_max_pd(mx, _mm512_max_pd(_r, _i));
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

        _r = VDABS(xr_);
        _i = VDABS(xi_);
        mx = _mm512_max_pd(mx, _mm512_max_pd(_r, _i));

        VZFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        register const VD yr_ = VDMULDIV(_r, c_);
        register const VD yi_ = VDMULDIV(_i, c_);
        _mm512_store_pd(yri, yr_);
        _mm512_store_pd(yii, yi_);

        _r = VDABS(yr_);
        _i = VDABS(yi_);
        mx = _mm512_max_pd(mx, _mm512_max_pd(_r, _i));
      }
    }
  }

  return _mm512_reduce_max_pd(mx);
}
