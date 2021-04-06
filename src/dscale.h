#ifndef DSCALE_H
#define DSCALE_H

#include "vec.h"

extern int dscale_(const fnat m[static restrict 1], const fnat n[static restrict 1], double A[static restrict VDL], const fnat ldA[static restrict 1], const fint e[static restrict 1]);
extern int dlscal_(const fnat m[static restrict 1], const fnat n[static restrict 1], double A[static restrict VDL], const fnat ldA[static restrict 1], const fnat l[static restrict 1], double M[static restrict 1], int e[static restrict 1]);

#endif /* !DSCALE_H */
