#include "cnorm2.h"

#include "scnrm2.h"

#ifndef USE_MKL
extern float BLAS_S(nrm2)(const fint n[static 1], const float x[static 1], const fint incx[static 1]);
#endif /* !USE_MKL */

float cnorm2_(const fnat m[static restrict 1], const float cr[static restrict 1], const float ci[static restrict 1], float e0[static restrict 1], float f0[static restrict 1], float e1[static restrict 1], float f1[static restrict 1])
{
  static const fint incx = 1;
  const float fr = BLAS_S(nrm2)((const fint*)m, cr, &incx);
  const float fi = BLAS_S(nrm2)((const fint*)m, ci, &incx);
  const float d = hypotf(fr, fi);
  flt2ef(d, e0, f0);
  sqeff(e0, f0, e1, f1);
  return d;
}
