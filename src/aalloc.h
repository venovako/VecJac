#ifndef AALLOC_H
#define AALLOC_H
#include "vec.h"
extern fint salloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], float **const restrict A, fnat ldA[static restrict 1]);
extern fint dalloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], double **const restrict A, fnat ldA[static restrict 1]);
extern fint calloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], float complex **const restrict A, fnat ldA[static restrict 1], float **const restrict Ar, fnat ldAr[static restrict 1], float **const restrict Ai, fnat ldAi[static restrict 1]);
extern fint zalloc2_(const fnat m[static restrict 1], const fnat n[static restrict 1], double complex **const restrict A, fnat ldA[static restrict 1], double **const restrict Ar, fnat ldAr[static restrict 1], double **const restrict Ai, fnat ldAi[static restrict 1]);
#endif /* !AALLOC_H */