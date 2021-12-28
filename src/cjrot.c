#include "cjrot.h"

#include "sjrot.h"
#include "vecdef.h"

float cjrot_(const fint n[static restrict 1], float xr[static restrict VSL], float xi[static restrict VSL], float yr[static restrict VSL], float yi[static restrict VSL], const float t[static restrict 1], const float c[static restrict 1], const float ca[static restrict 1], const float sa[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -10.0f;
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
  if ((*t == 0.0f) && (*c != 1.0f))
    return -7.0f;
#endif /* !NDEBUG */

  if (!*n)
    return 0.0f;

  if ((*sa == 0.0f) && (fabsf(*ca) == 1.0f)) {
    // real rotation
    const float t_ = ((*ca == 1.0f) ? *t : -*t);
    const float r_r = sjrot_(n, xr, yr, &t_, c);
    if (!(r_r >= 0.0f))
      return -11.0f;
    const float r_i = sjrot_(n, xi, yi, &t_, c);
    if (!(r_i >= 0.0f))
      return -12.0f;
    return fmaxf(r_r, r_i);
  }

  const float cat = ((*ca) * (*t));
  const float sat = ((*sa) * (*t));

  register const VS cta = _mm512_set1_ps(cat);
  register const VS sta = _mm512_set1_ps(sat);
  register const VS cna = _mm512_set1_ps(-cat);
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

  if (*n < 0) { // permute
    const fnat n_ = (fnat)-*n;
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
      register const VS yr_ = _mm512_mul_ps(_r, c_);
      register const VS yi_ = _mm512_mul_ps(_i, c_);
      _mm512_store_ps(yri, yr_);
      _mm512_store_ps(yii, yi_);

      _r = VSABS(yr_);
      _i = VSABS(yi_);
      mx = _mm512_max_ps(mx, _mm512_max_ps(_r, _i));

      VCFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
      register const VS xr_ = _mm512_mul_ps(_r, c_);
      register const VS xi_ = _mm512_mul_ps(_i, c_);
      _mm512_store_ps(xri, xr_);
      _mm512_store_ps(xii, xi_);

      _r = VSABS(xr_);
      _i = VSABS(xi_);
      mx = _mm512_max_ps(mx, _mm512_max_ps(_r, _i));
    }
  }
  else { // no permute
    const fnat n_ = (fnat)*n;
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
      register const VS xr_ = _mm512_mul_ps(_r, c_);
      register const VS xi_ = _mm512_mul_ps(_i, c_);
      _mm512_store_ps(xri, xr_);
      _mm512_store_ps(xii, xi_);

      _r = VSABS(xr_);
      _i = VSABS(xi_);
      mx = _mm512_max_ps(mx, _mm512_max_ps(_r, _i));

      VCFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
      register const VS yr_ = _mm512_mul_ps(_r, c_);
      register const VS yi_ = _mm512_mul_ps(_i, c_);
      _mm512_store_ps(yri, yr_);
      _mm512_store_ps(yii, yi_);

      _r = VSABS(yr_);
      _i = VSABS(yi_);
      mx = _mm512_max_ps(mx, _mm512_max_ps(_r, _i));
    }
  }

  return _mm512_reduce_max_ps(mx);
}
