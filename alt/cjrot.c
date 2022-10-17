#include "cjrot.h"

#include "sjrot.h"
#include "vecdef.h"

float cjrot_(const fint n[static restrict 1], float xr[static restrict VSL], float xi[static restrict VSL], float yr[static restrict VSL], float yi[static restrict VSL], const float c[static restrict 1], const float cat[static restrict 1], const float sat[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*n & VSL_1)
    return NAN;
  if (IS_NOT_ALIGNED(xr))
    return NAN;
  if (IS_NOT_ALIGNED(xi))
    return NAN;
  if (IS_NOT_ALIGNED(yr))
    return NAN;
  if (IS_NOT_ALIGNED(yi))
    return NAN;
#endif /* !NDEBUG */

  if (!*n)
    return -0.0f;

  if (*sat == 0.0f) {
    // real rotation
    const float r_r = sjrot_(n, xr, yr, c, cat);
#ifndef NDEBUG
    if (!(r_r == r_r))
      return NAN;
#endif /* !NDEBUG */
    const float r_i = sjrot_(n, xi, yi, c, cat);
#ifndef NDEBUG
    if (!(r_i == r_i))
      return NAN;
#endif /* !NDEBUG */
    const float r = fmaxf(fabsf(r_r), fabsf(r_i));
    return (((copysignf(1.0f, r_r) == -1.0f) && (copysignf(1.0f, r_i) == -1.0f)) ? -r : r);
  }

  register const VS cta = _mm512_set1_ps(*cat);
  register const VS sta = _mm512_set1_ps(*sat);
  register const VS cna = _mm512_set1_ps(-*cat);
  register const VS c_ = _mm512_set1_ps(*c);

  register const VS _zerof = _mm512_set1_ps(-0.0f);
  register VS mx = _mm512_setzero_ps()
    , _r
#ifndef NDEBUG
    = _mm512_setzero_ps()
#endif /* !NDEBUG */
    , _i
#ifndef NDEBUG
    = _mm512_setzero_ps()
#endif /* !NDEBUG */
    ;
  register MS ne = _mm512_kxor(ne, ne);

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

        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(x_r, _r));
        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(x_i, _i));

        _r = VSABS(_r);
        _i = VSABS(_i);
        mx = _mm512_max_ps(mx, _mm512_max_ps(_r, _i));

        VCFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        _mm512_store_ps(xri, _r);
        _mm512_store_ps(xii, _i);

        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(y_r, _r));
        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(y_i, _i));

        _r = VSABS(_r);
        _i = VSABS(_i);
        mx = _mm512_max_ps(mx, _mm512_max_ps(_r, _i));
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

        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(x_r, yr_));
        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(x_i, yi_));

        _r = VSABS(yr_);
        _i = VSABS(yi_);
        mx = _mm512_max_ps(mx, _mm512_max_ps(_r, _i));

        VCFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        register const VS xr_ = VSMULDIV(_r, c_);
        register const VS xi_ = VSMULDIV(_i, c_);
        _mm512_store_ps(xri, xr_);
        _mm512_store_ps(xii, xi_);

        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(y_r, xr_));
        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(y_i, xi_));

        _r = VSABS(xr_);
        _i = VSABS(xi_);
        mx = _mm512_max_ps(mx, _mm512_max_ps(_r, _i));
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

        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(x_r, _r));
        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(x_i, _i));

        _r = VSABS(_r);
        _i = VSABS(_i);
        mx = _mm512_max_ps(mx, _mm512_max_ps(_r, _i));

        VCFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        _mm512_store_ps(yri, _r);
        _mm512_store_ps(yii, _i);

        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(y_r, _r));
        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(y_i, _i));

        _r = VSABS(_r);
        _i = VSABS(_i);
        mx = _mm512_max_ps(mx, _mm512_max_ps(_r, _i));
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

        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(x_r, xr_));
        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(x_i, xi_));

        _r = VSABS(xr_);
        _i = VSABS(xi_);
        mx = _mm512_max_ps(mx, _mm512_max_ps(_r, _i));

        VCFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        register const VS yr_ = VSMULDIV(_r, c_);
        register const VS yi_ = VSMULDIV(_i, c_);
        _mm512_store_ps(yri, yr_);
        _mm512_store_ps(yii, yi_);

        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(y_r, yr_));
        ne = _mm512_kor(ne, _mm512_cmpneq_ps_mask(y_i, yi_));

        _r = VSABS(yr_);
        _i = VSABS(yi_);
        mx = _mm512_max_ps(mx, _mm512_max_ps(_r, _i));
      }
    }
  }

  const float r = _mm512_reduce_max_ps(mx);
  return (MS2U(ne) ? r : -r);
}
