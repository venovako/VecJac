#ifndef DGSSCL_H
#define DGSSCL_H

#include "vec.h"

extern double dgsscl_(const fint m[static restrict 1], const double t[static restrict 1], double x[static restrict VDL], double y[static restrict VDL], const double e[static restrict 2], const double f[static restrict 2]);

#endif /* !DGSSCL_H */
