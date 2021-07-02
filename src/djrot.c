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
  if (*t == 0.0) {
    if (*c == 1.0)
      return 0;
    // should never happen
    return -5;
  }

  if (*n < 0) { // permute
    const fint n_ = -*n;
  }
  else { // no permute
    const fint n_ = *n;
  }

  return 0;
}
