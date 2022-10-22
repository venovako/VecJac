#include "ddpscl.h"

double ddpscl_(const fnat m[static restrict 1], const double x[static restrict VDL], const double y[static restrict VDL], const double e[static restrict 2], const double f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*m & VDL_1)
    return NAN;
  if (IS_NOT_ALIGNED(x))
    return NAN;
  if (IS_NOT_ALIGNED(y))
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
  register VD p = _mm512_setzero_pd();

  for (fnat i = 0u; i < *m; i += VDL) {
    register VD xi = _mm512_load_pd(x + i);
    register VD yi = _mm512_load_pd(y + i);
    xi = _mm512_scalef_pd(xi, xe); VDP(xi);
    yi = _mm512_scalef_pd(yi, ye); VDP(yi);
    p = _mm512_fmadd_pd(xi, yi, p); VDP(p);
  }

  const double fx = f[0u];
  const double fy = f[1u];
  return (_mm512_reduce_add_pd(p) / (fx * fy));
}
