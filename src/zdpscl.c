#include "zdpscl.h"

double complex zdpscl_(const fnat m[static restrict 1], const double xr[static restrict VDL], const double xi[static restrict VDL], const double yr[static restrict VDL], const double yi[static restrict VDL], const double e[static restrict 2], const double f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*m & VDL_1)
    return NAN;
  if (IS_NOT_ALIGNED(xr))
    return NAN;
  if (IS_NOT_ALIGNED(xi))
    return NAN;
  if (IS_NOT_ALIGNED(yr))
    return NAN;
  if (IS_NOT_ALIGNED(yi))
    return NAN;
  if (!(e[0u] < HUGE_VAL))
    return NAN;
  if (!(e[1u] < HUGE_VAL))
    return NAN;
  if (!(f[0u] >= 1.0) || !(f[0u] < 2.0))
    return NAN;
  if (!(f[1u] >= 1.0) || !(f[0u] < 2.0))
    return NAN;
#endif /* !NDEBUG */

  const double ex = e[0u];
  if (!(ex > -HUGE_VAL))
    return 0.0;

  const double ey = e[1u];
  if (!(ey > -HUGE_VAL))
    return 0.0;

  register const VD xe = _mm512_set1_pd(-ex);
  register const VD ye = _mm512_set1_pd(-ey);

  register VD pr = _mm512_setzero_pd();
  register VD pi = _mm512_setzero_pd();

  for (fnat i = 0u; i < *m; i += VDL) {
    register VD xri = _mm512_load_pd(xr + i);
    register VD yri = _mm512_load_pd(yr + i);
    register VD yii = _mm512_load_pd(yi + i);
    register VD xii = _mm512_load_pd(xi + i);
    xri = _mm512_scalef_pd(xri, xe); VDP(xri);
    yri = _mm512_scalef_pd(yri, ye); VDP(yri);
    yii = _mm512_scalef_pd(yii, ye); VDP(yii);
    xii = _mm512_scalef_pd(xii, xe); VDP(xii);
    pr = _mm512_fmadd_pd(xri, yri, pr); VDP(pr);
    pi = _mm512_fmadd_pd(xri, yii, pi); VDP(pi);
    pr = _mm512_fmadd_pd(xii, yii, pr); VDP(pr);
    pi = _mm512_fnmadd_pd(xii, yri, pi); VDP(pi);
  }

  const double fx = f[0u];
  const double fy = f[1u];
  alignas(16u) double complex z;
  _mm_store_pd((double*)&z, _mm_div_pd(_mm_set_pd(_mm512_reduce_add_pd(pi), _mm512_reduce_add_pd(pr)), _mm_set1_pd(fx * fy)));
  return z;
}
