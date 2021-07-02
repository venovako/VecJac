#include "zjrot.h"

#include "djrot.h"

int zjrot_(const fint n[static restrict 1], double xr[static restrict VDL], double xi[static restrict VDL], double yr[static restrict VDL], double yi[static restrict VDL], const double t[static restrict 1], const double c[static restrict 1], const double ca[static restrict 1], const double sa[static restrict 1])
{
  if (IS_NOT_VFPENV)
    return -10;
  if (*n & VDL_1)
    return -1;
  if (IS_NOT_ALIGNED(xr))
    return -2;
  if (IS_NOT_ALIGNED(xi))
    return -3;
  if (IS_NOT_ALIGNED(yr))
    return -4;
  if (IS_NOT_ALIGNED(yi))
    return -5;

  if (!*n)
    return 0;
  if (*t == 0.0) {
    if (*c == 1.0)
      return 0;
    // should never happen
    return -7;
  }

  if (*sa == 0.0) {
    // real rotation
    if (*ca == 1.0) {
      const double t_ = *t;
      return imax(djrot_(n, xr, yr, &t_, c), djrot_(n, xi, yi, &t_, c));
    }
    if (*ca == -1.0) {
      const double t_ = -*t;
      return imax(djrot_(n, xr, yr, &t_, c), djrot_(n, xi, yi, &t_, c));
    }
  }

  if (*n < 0) { // permute
    const fint n_ = -*n;
  }
  else { // no permute
    const fint n_ = *n;
  }

  return 0;
}
