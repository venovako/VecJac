#ifndef WDP_H
#define WDP_H

#include "vec.h"

#ifndef USE_MKL
extern double BLAS_D(nrm2)(const fint n[static 1], const double x[static 1], const fint incx[static 1]);
extern double BLAS_D(dot)(const fint n[static 1], const double x[static 1], const fint incx[static 1], const double y[static 1], const fint incy[static 1]);
#endif /* !USE_MKL */

extern double wsq(const fnat n, const double x[static restrict 1]);
extern double dn2(const fnat n, const double x[static restrict 1]);
extern double ddp(const fnat n, const double x[static restrict 1]);
extern double xdp(const fnat n, const double x[static restrict 1]);

extern double dne(const fnat n, const double x[static restrict VDL]);
extern double dnf(const fnat n, const double x[static restrict VDL]);

extern double dre(const double c, const double e);

#endif /* !WDP_H */
