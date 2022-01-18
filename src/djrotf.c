#include "djrotf.h"

fint djrotf_(const fint n[static restrict 1], double x[static restrict VDL], double y[static restrict VDL], const double c[static restrict 1], const double at[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -6;
  if (*n & VDL_1)
    return -1;
  if (IS_NOT_ALIGNED(x))
    return -2;
  if (IS_NOT_ALIGNED(y))
    return -3;
#endif /* !NDEBUG */

  if (!*n)
    return 0;

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
        return 0;
      }
      // should never happen
      return -4;
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
        register const VD x_r = _mm512_mul_pd(_mm512_fnmadd_pd(x_, t_, y_), c_);
        register const VD y_r = _mm512_mul_pd(_mm512_fmadd_pd(y_, t_, x_), c_);
        _mm512_store_pd(xi, x_r);
        _mm512_store_pd(yi, y_r);
      }
    }
  }
  else { // no permute
    const fnat n_ = (fnat)*n;
    if (*at == 0.0)
      return ((*c == 1.0) ? 0 : -4);
    if (*c == 1.0) {
      register const VD t_ = _mm512_set1_pd(*at);
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xi = x + i;
        double *const yi = y + i;
        register const VD x_ = _mm512_load_pd(xi);
        register const VD y_ = _mm512_load_pd(yi);
        register const VD x_r = _mm512_fmadd_pd(y_, t_, x_);
        register const VD y_r = _mm512_fnmadd_pd(x_, t_, y_);
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
        register const VD x_r = _mm512_mul_pd(_mm512_fmadd_pd(y_, t_, x_), c_);
        register const VD y_r = _mm512_mul_pd(_mm512_fnmadd_pd(x_, t_, y_), c_);
        _mm512_store_pd(xi, x_r);
        _mm512_store_pd(yi, y_r);
      }
    }
  }

  return 0;
}
