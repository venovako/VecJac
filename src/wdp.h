#ifndef WDP_H
#define WDP_H

#include "vec.h"

extern wide wdsq(const fnat n, const double x[static restrict 1]);
extern wide wdn2(const fnat n, const double x[static restrict 1]);
extern wide wddp(const fnat n, const double x[static restrict 1]);

extern wide wdne(const fnat n, const double x[static restrict VDL]);
extern long double qdnre(const wide c, const wide e);

#endif /* !WDP_H */
