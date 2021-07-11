#include "djrot.h"

int djrot_(const fint n[static restrict 1], double x[static restrict VDL], double y[static restrict VDL], const double t[static restrict 1], const double c[static restrict 1])
{
  if (IS_NOT_VFPENV)
    return -6;
  if (*n & VDL_1)
    return -1;
  if (IS_NOT_ALIGNED(x))
    return -2;
  if (IS_NOT_ALIGNED(y))
    return -3;

  if (!*n)
    return 0;

  const double tc[2] = { *t, *c };

  if (*n < 0) { // permute
    const fnat n_ = (fnat)-*n;
    if (*t == 0.0) {
      if (*c == 1.0) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n_,x,y)
#endif /* _OPENMP */
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
      return -5;
    }
#ifdef _OPENMP
#pragma omp parallel default(none) shared(n_,tc,x,y)
#endif /* _OPENMP */
    {
      register const VD t_ = _mm512_set1_pd(tc[0]);
      register const VD c_ = _mm512_set1_pd(tc[1]);
#ifdef _OPENMP
#pragma omp for
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xi = x + i;
        double *const yi = y + i;
        register const VD x_ = _mm512_load_pd(xi);
        register const VD y_ = _mm512_load_pd(yi);
        _mm512_store_pd(yi, _mm512_mul_pd(_mm512_fmadd_pd(y_, t_, x_), c_));
        _mm512_store_pd(xi, _mm512_mul_pd(_mm512_fnmadd_pd(x_, t_, y_), c_));
      }
    }
  }
  else { // no permute
    if (*t == 0.0) {
      if (*c == 1.0)
        return 0;
      // should never happen
      return -5;
    }
    const fnat n_ = (fnat)*n;
#ifdef _OPENMP
#pragma omp parallel default(none) shared(n_,tc,x,y)
#endif /* _OPENMP */
    {
      register const VD t_ = _mm512_set1_pd(tc[0]);
      register const VD c_ = _mm512_set1_pd(tc[1]);
#ifdef _OPENMP
#pragma omp for
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xi = x + i;
        double *const yi = y + i;
        register const VD x_ = _mm512_load_pd(xi);
        register const VD y_ = _mm512_load_pd(yi);
        _mm512_store_pd(xi, _mm512_mul_pd(_mm512_fmadd_pd(y_, t_, x_), c_));
        _mm512_store_pd(yi, _mm512_mul_pd(_mm512_fnmadd_pd(x_, t_, y_), c_));
      }
    }
  }

  return 0;
}
