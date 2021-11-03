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

  if (!(e[0u] > -HUGE_VAL))
    return 0.0;
  if (!(e[1u] > -HUGE_VAL))
    return 0.0;

  register const VD xe = _mm512_set1_pd(-(e[0u]));
  register const VD ye = _mm512_set1_pd(-(e[1u]));

  register VD pr = _mm512_setzero_pd();
  for (fnat i = 0u; i < *m; i += VDL) {
    pr = _mm512_fmadd_pd(_mm512_scalef_pd(_mm512_load_pd(xr + i), xe), _mm512_scalef_pd(_mm512_load_pd(yr + i), ye), pr); VDP(pr);
  }
  for (fnat i = 0u; i < *m; i += VDL) {
    pr = _mm512_fmadd_pd(_mm512_scalef_pd(_mm512_load_pd(xi + i), xe), _mm512_scalef_pd(_mm512_load_pd(yi + i), ye), pr); VDP(pr);
  }

  register VD pi = _mm512_setzero_pd();
  for (fnat i = 0u; i < *m; i += VDL) {
    pi = _mm512_fmadd_pd(_mm512_scalef_pd(_mm512_load_pd(xr + i), xe), _mm512_scalef_pd(_mm512_load_pd(yi + i), ye), pi); VDP(pi);
  }
  for (fnat i = 0u; i < *m; i += VDL) {
    pi = _mm512_fnmadd_pd(_mm512_scalef_pd(_mm512_load_pd(xi + i), xe), _mm512_scalef_pd(_mm512_load_pd(yr + i), ye), pi); VDP(pi);
  }

  alignas(16u) double complex z
#ifndef NDEBUG
    = CMPLX(0.0, 0.0)
#endif /* !NDEBUG */
    ;
  _mm_store_pd((double*)&z, _mm_div_pd(_mm_set_pd(_mm512_reduce_add_pd(pi), _mm512_reduce_add_pd(pr)), _mm_set1_pd(f[0u] * f[1u])));
  return z;
}
