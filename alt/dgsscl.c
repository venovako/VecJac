#include "dgsscl.h"

#include "vecdef.h"

double dgsscl_(const fint m[static restrict 1], const double t[static restrict 1], double x[static restrict VDL], double y[static restrict VDL], const double e[static restrict 2], const double f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*m & VDL_1)
    return NAN;
  if (!(fabs(*t) <= DBL_MAX))
    return NAN;
  if (IS_NOT_ALIGNED(x))
    return NAN;
  if (IS_NOT_ALIGNED(y))
    return NAN;
  if (!(e[0u] <= DBL_MAX))
    return NAN;
  if (!(e[1u] <= DBL_MAX))
    return NAN;
  if (!(f[0u] >= 1.0) || !(f[0u] < 2.0))
    return NAN;
  if (!(f[1u] >= 1.0) || !(f[0u] < 2.0))
    return NAN;
#endif /* !NDEBUG */

  if (!*m)
    return -0.0;

  const double ex = e[0u];
  if (!(ex >= -DBL_MAX))
    return -0.0;

  const double ey = e[1u];
  if (!(ey >= -DBL_MAX))
    return -0.0;

  const double fx = f[0u];
  const double fy = f[1u];

  register const VD x_e = _mm512_set1_pd(-ex);
  register const VD y_e = _mm512_set1_pd(-ey);
  register const VD _zero = _mm512_set1_pd(-0.0);
  register VD M = _mm512_setzero_pd();
  register MD ne = MDXOR(ne, ne);

  if (*m < 0) { // transform x and permute
    register const VD tf = _mm512_set1_pd(-*t * (fx / fy));
    register const VD xe = _mm512_set1_pd(ex);
    const fnat _m = (fnat)-*m;

    for (fnat i = 0u; i < _m; i += VDL) {
      double *const x_i = x + i;
      double *const y_i = y + i;
      register const VD xi = _mm512_load_pd(x_i); VDP(xi);
      register VD yi = _mm512_load_pd(y_i); VDP(yi);
      _mm512_store_pd(x_i, yi);
      yi = _mm512_scalef_pd(yi, y_e);
      yi = _mm512_scalef_pd(_mm512_fmadd_pd(tf, yi, _mm512_scalef_pd(xi, x_e)), xe);
      ne = MDOR(ne, _mm512_cmpneq_pd_mask(xi, yi));
      _mm512_store_pd(y_i, yi);
      yi = VDABS(yi);
      M = _mm512_max_pd(M, yi); VDP(M);
    }
  }
  else { // transform y
    register const VD tf = _mm512_set1_pd(-*t * (fy / fx));
    register const VD ye = _mm512_set1_pd(ey);
    const fnat _m = (fnat)*m;

    for (fnat i = 0u; i < _m; i += VDL) {
      double *const y_i = y + i;
      register VD xi = _mm512_load_pd(x + i); VDP(xi);
      register const VD yi = _mm512_load_pd(y_i); VDP(yi);
      xi = _mm512_scalef_pd(xi, x_e);
      xi = _mm512_scalef_pd(_mm512_fmadd_pd(tf, xi, _mm512_scalef_pd(yi, y_e)), ye);
      ne = MDOR(ne, _mm512_cmpneq_pd_mask(yi, xi));
      _mm512_store_pd(y_i, xi);
      xi = VDABS(xi);
      M = _mm512_max_pd(M, xi); VDP(M);
    }
  }

  const double r = _mm512_reduce_max_pd(M);
  return (MD2U(ne) ? r : -r);
}
