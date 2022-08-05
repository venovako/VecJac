#ifndef CSPLIT_H
#define CSPLIT_H

#include "vec.h"

extern fint csplit_(const fnat m[static restrict 1], const fnat n[static restrict 1], const float complex A[static restrict VSL_2], const fnat ldA[static restrict 1], float Ar[static restrict VSL], const fnat ldAr[static restrict 1], float Ai[static restrict VSL], const fnat ldAi[static restrict 1]);

#endif /* !CSPLIT_H */
