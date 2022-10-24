#include "djrotf.h"

#include "vecdef.h"

double djrotf_(const fint n[static restrict 1], double x[static restrict VDL], double y[static restrict VDL], const double c[static restrict 1], const double at[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*n & VDL_1)
    return NAN;
  if (IS_NOT_ALIGNED(x))
    return NAN;
  if (IS_NOT_ALIGNED(y))
    return NAN;
#endif /* !NDEBUG */

  if (!*n)
    return 0.0;

  register const VD _zero = _mm512_set1_pd(-0.0);
  register VD mx = _mm512_setzero_pd();

  if (*n < 0) { // permute
    const fnat n_ = (fnat)-*n;
    if (*at == 0.0) {
      if (*c == 1.0) {
        for (fnat i = 0u; i < n_; i += VDL) {
          double *const xi = x + i;
          double *const yi = y + i;
          register const VD x_ = _mm512_load_pd(xi);
          register const VD y_ = _mm512_load_pd(yi);
          _mm512_store_pd(yi, x_);
          _mm512_store_pd(xi, y_);
        }
        return 0.0;
      }
      // should never happen
      return NAN;
    }
    if (*c == 1.0) {
      register const VD t_ = _mm512_set1_pd(*at);
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xi = x + i;
        double *const yi = y + i;
        register const VD x_ = _mm512_load_pd(xi);
        register const VD y_ = _mm512_load_pd(yi);
        register const VD x_r = _mm512_fnmadd_pd(x_, t_, y_);
        register const VD y_r = _mm512_fmadd_pd(y_, t_, x_);
        mx = _mm512_max_pd(mx, _mm512_max_pd(VDABS(x_r), VDABS(y_r)));
        _mm512_store_pd(xi, x_r);
        _mm512_store_pd(yi, y_r);
      }
    }
    else { // cos < 1
      register const VD t_ = _mm512_set1_pd(*at);
      register const VD c_ = _mm512_set1_pd(*c);
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xi = x + i;
        double *const yi = y + i;
        register const VD x_ = _mm512_load_pd(xi);
        register const VD y_ = _mm512_load_pd(yi);
        register const VD x_r = VDMULDIV(_mm512_fnmadd_pd(x_, t_, y_), c_);
        register const VD y_r = VDMULDIV(_mm512_fmadd_pd(y_, t_, x_), c_);
        mx = _mm512_max_pd(mx, _mm512_max_pd(VDABS(x_r), VDABS(y_r)));
        _mm512_store_pd(xi, x_r);
        _mm512_store_pd(yi, y_r);
      }
    }
  }
  else { // no permute
    const fnat n_ = (fnat)*n;
    if (*at == 0.0)
      return ((*c == 1.0) ? 0.0 : (double)NAN);
    if (*c == 1.0) {
      register const VD t_ = _mm512_set1_pd(*at);
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xi = x + i;
        double *const yi = y + i;
        register const VD x_ = _mm512_load_pd(xi);
        register const VD y_ = _mm512_load_pd(yi);
        register const VD x_r = _mm512_fmadd_pd(y_, t_, x_);
        register const VD y_r = _mm512_fnmadd_pd(x_, t_, y_);
        mx = _mm512_max_pd(mx, _mm512_max_pd(VDABS(x_r), VDABS(y_r)));
        _mm512_store_pd(xi, x_r);
        _mm512_store_pd(yi, y_r);
      }
    }
    else { // cos < 1
      register const VD t_ = _mm512_set1_pd(*at);
      register const VD c_ = _mm512_set1_pd(*c);
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xi = x + i;
        double *const yi = y + i;
        register const VD x_ = _mm512_load_pd(xi);
        register const VD y_ = _mm512_load_pd(yi);
        register const VD x_r = VDMULDIV(_mm512_fmadd_pd(y_, t_, x_), c_);
        register const VD y_r = VDMULDIV(_mm512_fnmadd_pd(x_, t_, y_), c_);
        mx = _mm512_max_pd(mx, _mm512_max_pd(VDABS(x_r), VDABS(y_r)));
        _mm512_store_pd(xi, x_r);
        _mm512_store_pd(yi, y_r);
      }
    }
  }

  return _mm512_reduce_max_pd(mx);
}
