#ifndef ZSCALE_H
#define ZSCALE_H

#include "vec.h"

extern int zscale_(const fnat m[static restrict 1], const fnat n[static restrict 1], double Ar[static restrict VDL], const fnat ldAr[static restrict 1], double Ai[static restrict VDL], const fnat ldAi[static restrict 1], const fint e[static restrict 1]);
extern int zlscal_(const fnat m[static restrict 1], const fnat n[static restrict 1], double Ar[static restrict VDL], const fnat ldAr[static restrict 1], double Ai[static restrict VDL], const fnat ldAi[static restrict 1], const fnat l[static restrict 1], double M[static restrict 1], fint e[static restrict 1]);

#endif /* !ZSCALE_H */
