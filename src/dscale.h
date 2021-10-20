#ifndef DSCALE_H
#define DSCALE_H

#include "vec.h"

extern fint dscale_(const fnat m[static restrict 1], const fnat n[static restrict 1], double A[static restrict VDL], const fnat ldA[static restrict 1], const fint e[static restrict 1]);
extern fint dlscal_(const fnat m[static restrict 1], const fnat n[static restrict 1], double A[static restrict VDL], const fnat ldA[static restrict 1], const fnat l[static restrict 1], double M[static restrict 1], fint e[static restrict 1]);

#endif /* !DSCALE_H */
