#include "dnorm2.h"

#include "dznrm2.h"

double dnorm2_(const fnat m[static restrict 1], const double x[static restrict 1], double e0[static restrict 1], double f0[static restrict 1], double e1[static restrict 1], double f1[static restrict 1])
{
  static const fint incx = 1;
  const double d = BLAS_D(nrm2)((const fint*)&m, x, &incx);
  dbl2ef(d, e0, f0);
  dbl2ef((d * d), e1, f1);
  return d;
}
