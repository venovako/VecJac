#include "cjrotf.h"

#include "sjrotf.h"
#include "vecdef.h"

float cjrotf_(const fint n[static restrict 1], float xr[static restrict VSL], float xi[static restrict VSL], float yr[static restrict VSL], float yi[static restrict VSL], const float c[static restrict 1], const float cat[static restrict 1], const float sat[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -9.0f;
  if (*n & VSL_1)
    return -1.0f;
  if (IS_NOT_ALIGNED(xr))
    return -2.0f;
  if (IS_NOT_ALIGNED(xi))
    return -3.0f;
  if (IS_NOT_ALIGNED(yr))
    return -4.0f;
  if (IS_NOT_ALIGNED(yi))
    return -5.0f;
#endif /* !NDEBUG */

  if (!*n)
    return 0.0f;

  if (*sat == 0.0f) {
    // real rotation
    if (sjrotf_(n, xr, yr, c, cat))
      return -10.0f;
    if (sjrotf_(n, xi, yi, c, cat))
      return -11.0f;
    return 0.0f;
  }

  register const VS cta = _mm512_set1_ps(*cat);
  register const VS sta = _mm512_set1_ps(*sat);
  register const VS cna = _mm512_set1_ps(-*cat);
  register const VS c_ = _mm512_set1_ps(*c);

  register VS
    _r
#ifndef NDEBUG
    = _mm512_setzero_ps()
#endif /* !NDEBUG */
    , _i
#ifndef NDEBUG
    = _mm512_setzero_ps()
#endif /* !NDEBUG */
    ;

  if (*n < 0) { // permute
    const fnat n_ = (fnat)-*n;
    if (*c == 1.0f) {
      for (fnat i = 0u; i < n_; i += VSL) {
        float *const xri = xr + i;
        float *const xii = xi + i;
        float *const yri = yr + i;
        float *const yii = yi + i;

        register const VS x_r = _mm512_load_ps(xri);
        register const VS x_i = _mm512_load_ps(xii);
        register const VS y_r = _mm512_load_ps(yri);
        register const VS y_i = _mm512_load_ps(yii);

        VCFMA(_r,_i,y_r,y_i,cta,sta,x_r,x_i);
        _mm512_store_ps(yri, _r);
        _mm512_store_ps(yii, _i);

        VCFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        _mm512_store_ps(xri, _r);
        _mm512_store_ps(xii, _i);
      }
    }
    else { // cos < 1
      for (fnat i = 0u; i < n_; i += VSL) {
        float *const xri = xr + i;
        float *const xii = xi + i;
        float *const yri = yr + i;
        float *const yii = yi + i;

        register const VS x_r = _mm512_load_ps(xri);
        register const VS x_i = _mm512_load_ps(xii);
        register const VS y_r = _mm512_load_ps(yri);
        register const VS y_i = _mm512_load_ps(yii);

        VCFMA(_r,_i,y_r,y_i,cta,sta,x_r,x_i);
        register const VS yr_ = VSMULDIV(_r, c_);
        register const VS yi_ = VSMULDIV(_i, c_);
        _mm512_store_ps(yri, yr_);
        _mm512_store_ps(yii, yi_);

        VCFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        register const VS xr_ = VSMULDIV(_r, c_);
        register const VS xi_ = VSMULDIV(_i, c_);
        _mm512_store_ps(xri, xr_);
        _mm512_store_ps(xii, xi_);
      }
    }
  }
  else { // no permute
    const fnat n_ = (fnat)*n;
    if (*c == 1.0f) {
      for (fnat i = 0u; i < n_; i += VSL) {
        float *const xri = xr + i;
        float *const xii = xi + i;
        float *const yri = yr + i;
        float *const yii = yi + i;

        register const VS x_r = _mm512_load_ps(xri);
        register const VS x_i = _mm512_load_ps(xii);
        register const VS y_r = _mm512_load_ps(yri);
        register const VS y_i = _mm512_load_ps(yii);

        VCFMA(_r,_i,y_r,y_i,cta,sta,x_r,x_i);
        _mm512_store_ps(xri, _r);
        _mm512_store_ps(xii, _i);

        VCFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        _mm512_store_ps(yri, _r);
        _mm512_store_ps(yii, _i);
      }
    }
    else { // cos < 1
      for (fnat i = 0u; i < n_; i += VSL) {
        float *const xri = xr + i;
        float *const xii = xi + i;
        float *const yri = yr + i;
        float *const yii = yi + i;

        register const VS x_r = _mm512_load_ps(xri);
        register const VS x_i = _mm512_load_ps(xii);
        register const VS y_r = _mm512_load_ps(yri);
        register const VS y_i = _mm512_load_ps(yii);

        VCFMA(_r,_i,y_r,y_i,cta,sta,x_r,x_i);
        register const VS xr_ = VSMULDIV(_r, c_);
        register const VS xi_ = VSMULDIV(_i, c_);
        _mm512_store_ps(xri, xr_);
        _mm512_store_ps(xii, xi_);

        VCFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        register const VS yr_ = VSMULDIV(_r, c_);
        register const VS yi_ = VSMULDIV(_i, c_);
        _mm512_store_ps(yri, yr_);
        _mm512_store_ps(yii, yi_);
      }
    }
  }

  return 0.0f;
}
