#ifndef CMERGE_H
#define CMERGE_H

#include "vec.h"

extern fint cmerge_(const fnat m[static restrict 1], const fnat n[static restrict 1], const float Ar[static restrict VSL], const fnat ldAr[static restrict 1], const float Ai[static restrict VSL], const fnat ldAi[static restrict 1], float complex A[static restrict VSL_2], const fnat ldA[static restrict 1]);

#endif /* !CMERGE_H */
