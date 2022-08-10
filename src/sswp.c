#include "sswp.h"

fint sswp_(const fnat n[static restrict 1], float x[static restrict VSL], float y[static restrict VSL])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -4;
  if (*n & VSL_1)
    return -1;
  if (IS_NOT_ALIGNED(x))
    return -2;
  if (IS_NOT_ALIGNED(y))
    return -3;
#endif /* !NDEBUG */

  for (fnat i = 0u; i < *n; i += VSL) {
    float *const x_i = x + i;
    float *const y_i = y + i;
    register const VS xi = _mm512_load_ps(x_i);
    register const VS yi = _mm512_load_ps(y_i);
    _mm512_store_ps(y_i, xi);
    _mm512_store_ps(x_i, yi);
  }

  return 0;
}
