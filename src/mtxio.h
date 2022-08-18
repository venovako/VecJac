#ifndef MTXIO_H
#define MTXIO_H

#include "common.h"

extern int open_ro_(const char bn[static restrict 1], const char ex[static restrict 1]);
extern int open_wo_(const char bn[static restrict 1], const char ex[static restrict 1]);
extern int resizef_(const int fd[static restrict 1], const size_t sz[static restrict 1]);

extern int sread2_(const fnat m[static restrict 1], const fnat n[static restrict 1], float A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1]);
extern int dread2_(const fnat m[static restrict 1], const fnat n[static restrict 1], double A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1]);
extern int cread2_(const fnat m[static restrict 1], const fnat n[static restrict 1], float complex A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1]);
extern int zread2_(const fnat m[static restrict 1], const fnat n[static restrict 1], double complex A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1]);

extern int wwrite1_(const fnat m[static restrict 1], const wide w[static restrict 1], const int fd[static restrict 1]);

extern int swrite2_(const fnat m[static restrict 1], const fnat n[static restrict 1], const float A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1]);
extern int dwrite2_(const fnat m[static restrict 1], const fnat n[static restrict 1], const double A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1]);
extern int cwrite2_(const fnat m[static restrict 1], const fnat n[static restrict 1], const float complex A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1]);
extern int zwrite2_(const fnat m[static restrict 1], const fnat n[static restrict 1], const double complex A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1]);

#endif /* !MTXIO_H */
