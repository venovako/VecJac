#include "znorm2.h"

#include "dznrm2.h"

#ifndef USE_MKL
extern double BLAS_D(nrm2)(const fint n[static 1], const double x[static 1], const fint incx[static 1]);
#endif /* !USE_MKL */

double znorm2_(const fnat m[static restrict 1], const double zr[static restrict 1], const double zi[static restrict 1], double e0[static restrict 1], double f0[static restrict 1], double e1[static restrict 1], double f1[static restrict 1])
{
  static const fint incx = 1;
  const double dr = BLAS_D(nrm2)((const fint*)m, zr, &incx);
  const double di = BLAS_D(nrm2)((const fint*)m, zi, &incx);
  const double d = cr_hypot(dr, di);
  dbl2ef(d, e0, f0);
  sqef(e0, f0, e1, f1);
  return d;
}
