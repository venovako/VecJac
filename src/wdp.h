#ifndef WDP_H
#define WDP_H

#include "vec.h"

extern double wsq(const fnat n, const double x[static restrict 1]);
extern double dn2(const fnat n, const double x[static restrict 1]);
extern double ddp(const fnat n, const double x[static restrict 1]);
extern double xdp(const fnat n, const double x[static restrict 1]);

extern double dne(const fnat n, const double x[static restrict VDL]);
extern double dnf(const fnat n, const double x[static restrict VDL]);

extern double dre(const double c, const double e);

#endif /* !WDP_H */
