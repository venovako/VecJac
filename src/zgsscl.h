#ifndef ZGSSCL_H
#define ZGSSCL_H

#include "vec.h"

extern double zgsscl_(const fint m[static restrict 1], const double tr[static restrict 1], const double ti[static restrict 1], double xr[static restrict VDL], double xi[static restrict VDL], double yr[static restrict VDL], double yi[static restrict VDL], const double e[static restrict 2], const double f[static restrict 2]);

#endif /* !ZGSSCL_H */
