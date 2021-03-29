#ifndef AALLOC_H
#define AALLOC_H
#include "vec.h"
extern void salloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], float *A[static restrict 1], fnat ldA[static restrict 1]);
extern void dalloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], double *A[static restrict 1], fnat ldA[static restrict 1]);
extern void calloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], float complex *A[static restrict 1], fnat ldA[static restrict 1], float *Ar[static restrict 1], fnat ldAr[static restrict 1], float *Ai[static restrict 1], fnat ldAi[static restrict 1]);
extern void zalloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], double complex *A[static restrict 1], fnat ldA[static restrict 1], double *Ar[static restrict 1], fnat ldAr[static restrict 1], double *Ai[static restrict 1], fnat ldAi[static restrict 1]);
#endif /* !AALLOC_H */
