#include "dswp.h"

fint dswp_(const fnat n[static restrict 1], double x[static restrict VDL], double y[static restrict VDL])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -4;
  if (*n & VDL_1)
    return -1;
  if (IS_NOT_ALIGNED(x))
    return -2;
  if (IS_NOT_ALIGNED(y))
    return -3;
#endif /* !NDEBUG */

  for (fnat i = 0u; i < *n; i += VDL) {
    double *const x_i = x + i;
    double *const y_i = y + i;
    register const VD xi = _mm512_load_pd(x_i);
    register const VD yi = _mm512_load_pd(y_i);
    _mm512_store_pd(y_i, xi);
    _mm512_store_pd(x_i, yi);
  }

  return 0;
}
