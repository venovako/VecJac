#include "sjrot.h"

#include "vecdef.h"

float sjrot_(const fint n[static restrict 1], float x[static restrict VSL], float y[static restrict VSL], const float c[static restrict 1], const float at[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -6.0f;
  if (*n & VSL_1)
    return -1.0f;
  if (IS_NOT_ALIGNED(x))
    return -2.0f;
  if (IS_NOT_ALIGNED(y))
    return -3.0f;
#endif /* !NDEBUG */

  if (!*n)
    return 0.0f;

  register const VS _zerof = _mm512_set1_ps(-0.0f);
  register VS mx = _mm512_setzero_ps();

  if (*n < 0) { // permute
    const fnat n_ = (fnat)-*n;
    if (*at == 0.0f) {
      if (*c == 1.0f) {
        for (fnat i = 0u; i < n_; i += VSL) {
          float *const xi = x + i;
          float *const yi = y + i;
          register const VS x_ = _mm512_load_ps(xi);
          register const VS y_ = _mm512_load_ps(yi);
          _mm512_store_ps(yi, x_);
          _mm512_store_ps(xi, y_);
        }
        return 0.0f;
      }
      // should never happen
      return -4.0f;
    }
    if (*c == 1.0f) {
      register const VS t_ = _mm512_set1_ps(*at);
      for (fnat i = 0u; i < n_; i += VSL) {
        float *const xi = x + i;
        float *const yi = y + i;
        register const VS x_ = _mm512_load_ps(xi);
        register const VS y_ = _mm512_load_ps(yi);
        register const VS x_r = _mm512_fnmadd_ps(x_, t_, y_);
        register const VS y_r = _mm512_fmadd_ps(y_, t_, x_);
        mx = _mm512_max_ps(mx, _mm512_max_ps(VSABS(x_r), VSABS(y_r)));
        _mm512_store_ps(xi, x_r);
        _mm512_store_ps(yi, y_r);
      }
    }
    else { // cos < 1
      register const VS t_ = _mm512_set1_ps(*at);
      register const VS c_ = _mm512_set1_ps(*c);
      for (fnat i = 0u; i < n_; i += VSL) {
        float *const xi = x + i;
        float *const yi = y + i;
        register const VS x_ = _mm512_load_ps(xi);
        register const VS y_ = _mm512_load_ps(yi);
        register const VS x_r = VSMULDIV(_mm512_fnmadd_ps(x_, t_, y_), c_);
        register const VS y_r = VSMULDIV(_mm512_fmadd_ps(y_, t_, x_), c_);
        mx = _mm512_max_ps(mx, _mm512_max_ps(VSABS(x_r), VSABS(y_r)));
        _mm512_store_ps(xi, x_r);
        _mm512_store_ps(yi, y_r);
      }
    }
  }
  else { // no permute
    const fnat n_ = (fnat)*n;
    if (*at == 0.0f)
      return ((*c == 1.0f) ? 0.0f : -4.0f);
    if (*c == 1.0f) {
      register const VS t_ = _mm512_set1_ps(*at);
      for (fnat i = 0u; i < n_; i += VSL) {
        float *const xi = x + i;
        float *const yi = y + i;
        register const VS x_ = _mm512_load_ps(xi);
        register const VS y_ = _mm512_load_ps(yi);
        register const VS x_r = _mm512_fmadd_ps(y_, t_, x_);
        register const VS y_r = _mm512_fnmadd_ps(x_, t_, y_);
        mx = _mm512_max_ps(mx, _mm512_max_ps(VSABS(x_r), VSABS(y_r)));
        _mm512_store_ps(xi, x_r);
        _mm512_store_ps(yi, y_r);
      }
    }
    else { // cos < 1
      register const VS t_ = _mm512_set1_ps(*at);
      register const VS c_ = _mm512_set1_ps(*c);
      for (fnat i = 0u; i < n_; i += VSL) {
        float *const xi = x + i;
        float *const yi = y + i;
        register const VS x_ = _mm512_load_ps(xi);
        register const VS y_ = _mm512_load_ps(yi);
        register const VS x_r = VSMULDIV(_mm512_fmadd_ps(y_, t_, x_), c_);
        register const VS y_r = VSMULDIV(_mm512_fnmadd_ps(x_, t_, y_), c_);
        mx = _mm512_max_ps(mx, _mm512_max_ps(VSABS(x_r), VSABS(y_r)));
        _mm512_store_ps(xi, x_r);
        _mm512_store_ps(yi, y_r);
      }
    }
  }

  return _mm512_reduce_max_ps(mx);
}
