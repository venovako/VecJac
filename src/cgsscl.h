#ifndef CGSSCL_H
#define CGSSCL_H

#include "vec.h"

extern float cgsscl_(const fint m[static restrict 1], const float tr[static restrict 1], const float ti[static restrict 1], float xr[static restrict VSL], float xi[static restrict VSL], float yr[static restrict VSL], float yi[static restrict VSL], const float e[static restrict 2], const float f[static restrict 2]);

#endif /* !CGSSCL_H */
