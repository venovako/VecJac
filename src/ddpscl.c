#include "ddpscl.h"

double ddpscl_(const fnat m[static restrict 1], const double x[static restrict VDL], const double y[static restrict VDL], const double e[static restrict 2], const double f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -5.0;
  if (*m & VDL_1)
    return -1.0;
  if (IS_NOT_ALIGNED(x))
    return -2.0;
  if (IS_NOT_ALIGNED(y))
    return -3.0;
  if (!(e[0u] > -HUGE_VAL))
    return -4.0;
  if (!(e[1u] > -HUGE_VAL))
    return -4.0;
#endif /* !NDEBUG */

  register const VD xe = _mm512_set1_pd(-(e[0u]));
  register const VD ye = _mm512_set1_pd(-(e[1u]));
  register VD p = _mm512_setzero_pd();

  for (fnat i = 0u; i < *m; i += VDL) {
    p = _mm512_fmadd_pd(_mm512_scalef_pd(_mm512_load_pd(x + i), xe), _mm512_scalef_pd(_mm512_load_pd(y + i), ye), p); VDP(p);
  }

  return (_mm512_reduce_add_pd(p) / (f[0u] * f[1u]));
}
