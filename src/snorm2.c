#include "snorm2.h"

#include "scnrm2.h"

float snorm2_(const fnat m[static restrict 1], const float x[static restrict 1], float e0[static restrict 1], float f0[static restrict 1], float e1[static restrict 1], float f1[static restrict 1])
{
  static const fint incx = 1;
  const float d = BLAS_S(nrm2)((const fint*)&m, x, &incx);
  flt2ef(d, e0, f0);
  sqeff(e0, f0, e1, f1);
  return d;
}
