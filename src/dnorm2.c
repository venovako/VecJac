#include "dnorm2.h"

#include "dznrm2.h"

#ifndef USE_MKL
extern double BLAS_D(nrm2)(const fint n[static 1], const double x[static 1], const fint incx[static 1]);
#endif /* !USE_MKL */

double dnorm2_(const fnat m[static restrict 1], const double x[static restrict 1], double e0[static restrict 1], double f0[static restrict 1], double e1[static restrict 1], double f1[static restrict 1])
{
  static const fint incx = 1;
  const double d = BLAS_D(nrm2)((const fint*)m, x, &incx);
  dbl2ef(d, e0, f0);
  sqef(e0, f0, e1, f1);
  return d;
}
