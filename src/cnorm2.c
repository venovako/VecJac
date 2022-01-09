#include "cnorm2.h"

#include "scnrm2.h"

float cnorm2_(const fnat m[static restrict 1], const float zr[static restrict 1], const float zi[static restrict 1], float e0[static restrict 1], float f0[static restrict 1], float e1[static restrict 1], float f1[static restrict 1])
{
  static const fint incx = 1;
  const float dr = BLAS_S(nrm2)((const fint*)m, zr, &incx);
  const float di = BLAS_S(nrm2)((const fint*)m, zi, &incx);
  const float d = hypotf(dr, di);
  flt2ef(d, e0, f0);
  sqeff(e0, f0, e1, f1);
  return d;
}
