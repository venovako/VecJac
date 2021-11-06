#include "djrot.h"

double djrot_(const fint n[static restrict 1], double x[static restrict VDL], double y[static restrict VDL], const double t[static restrict 1], const double c[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -6.0;
  if (*n & VDL_1)
    return -1.0;
  if (IS_NOT_ALIGNED(x))
    return -2.0;
  if (IS_NOT_ALIGNED(y))
    return -3.0;
#endif /* !NDEBUG */

  if (!*n)
    return 0.0;

  register VD mx = _mm512_setzero_pd();
  if (*n < 0) { // permute
    const fnat n_ = (fnat)-*n;
    if (*t == 0.0) {
      if (*c == 1.0) {
        for (fnat i = 0u; i < n_; i += VDL) {
          double *const xi = x + i;
          double *const yi = y + i;
          register const VD x_ = _mm512_load_pd(xi);
          register const VD y_ = _mm512_load_pd(yi);
          mx = _mm512_max_pd(mx, _mm512_max_pd(_mm512_abs_pd(x_), _mm512_abs_pd(y_)));
          _mm512_store_pd(xi, y_);
          _mm512_store_pd(yi, x_);
        }
        return _mm512_reduce_max_pd(mx);
      }
      // should never happen
      return -5.0;
    }
    register const VD t_ = _mm512_set1_pd(*t);
    register const VD c_ = _mm512_set1_pd(*c);
    for (fnat i = 0u; i < n_; i += VDL) {
      double *const xi = x + i;
      double *const yi = y + i;
      register const VD x_ = _mm512_load_pd(xi);
      register const VD y_ = _mm512_load_pd(yi);
      register const VD x_r = _mm512_mul_pd(_mm512_fnmadd_pd(x_, t_, y_), c_);
      register const VD y_r = _mm512_mul_pd(_mm512_fmadd_pd(y_, t_, x_), c_);
      mx = _mm512_max_pd(mx, _mm512_max_pd(_mm512_abs_pd(x_r), _mm512_abs_pd(y_r)));
      _mm512_store_pd(xi, x_r);
      _mm512_store_pd(yi, y_r);
    }
  }
  else { // no permute
    const fnat n_ = (fnat)*n;
    if (*t == 0.0) {
      if (*c == 1.0) {
        for (fnat i = 0u; i < n_; i += VDL)
          mx = _mm512_max_pd(mx, _mm512_max_pd(_mm512_abs_pd(_mm512_load_pd(x + i)), _mm512_abs_pd(_mm512_load_pd(y + i))));
        return _mm512_reduce_max_pd(mx);
      }
      // should never happen
      return -5.0;
    }
    register const VD t_ = _mm512_set1_pd(*t);
    register const VD c_ = _mm512_set1_pd(*c);
    for (fnat i = 0u; i < n_; i += VDL) {
      double *const xi = x + i;
      double *const yi = y + i;
      register const VD x_ = _mm512_load_pd(xi);
      register const VD y_ = _mm512_load_pd(yi);
      register const VD x_r = _mm512_mul_pd(_mm512_fmadd_pd(y_, t_, x_), c_);
      register const VD y_r = _mm512_mul_pd(_mm512_fnmadd_pd(x_, t_, y_), c_);
      mx = _mm512_max_pd(mx, _mm512_max_pd(_mm512_abs_pd(x_r), _mm512_abs_pd(y_r)));
      _mm512_store_pd(xi, x_r);
      _mm512_store_pd(yi, y_r);
    }
  }

  return _mm512_reduce_max_pd(mx);
}
