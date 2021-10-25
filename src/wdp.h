#ifndef WDP_H
#define WDP_H

#include "vec.h"

// call the reference BLAS
extern double dnb(const fnat n, const double x[static restrict 1]);

// call the MKL
extern double ddp(const fnat n, const double x[static restrict 1]);
extern double dn2(const fnat n, const double x[static restrict 1]);

// use long double
extern double dnf(const fnat n, const double x[static restrict 1]);

extern double dne(const fnat n, const double x[static restrict VDL]);
extern double dns(const fnat n, const double x[static restrict VDL]);

extern double dnc(const fnat n, const double x[static restrict VDL]);
extern double dnd(const fnat n, const double x[static restrict VDL]);

extern double wsq(const fnat n, const double x[static restrict 1]);
extern double dre(const double c, const double e);

#endif /* !WDP_H */
