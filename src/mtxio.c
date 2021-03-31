#include "mtxio.h"

#ifdef FNCAT
#error FNCAT already defined
#else /* !FNCAT */
#define FNCAT(nn,bn,ex)                                         \
  char *const nn = (char*)alloca(strlen(bn) + strlen(ex) + 2u); \
  strcat(strcat(strcpy(nn, bn), "."), ex)
#endif /* ?FNCAT */

int open_ro_(const char bn[static restrict 1], const char ex[static restrict 1])
{
  FNCAT(fn, bn, ex);
  return open(fn, (O_RDONLY | O_LARGEFILE));
}

int open_wo_(const char bn[static restrict 1], const char ex[static restrict 1])
{
  FNCAT(fn, bn, ex);
  return open(fn, (O_WRONLY | O_CREAT | O_LARGEFILE), (S_IRUSR | S_IWUSR));
}

int resizef_(const int fd[static restrict 1], const size_t sz[static restrict 1])
{
  return ((*fd >= 0) ? (ftruncate(*fd, *sz) || fsync(*fd)) : *fd);
}

int sread2_(const fnat m[static restrict 1], const fnat n[static restrict 1], float A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1])
{
  if (*ldA < *m)
    return -4;
  if (*fd < 0)
    return -5;
  if (!*m)
    return 0;
  if (!*n)
    return 0;

  const size_t sz = *m * sizeof(float);
  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd,sz) reduction(+:r)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    float *const Aj = A + j * (*ldA);
    const size_t of = j * (*m * sizeof(float));
    r += ((ssize_t)sz != pread(*fd, Aj, sz, of));
  }

  return r;
}

int dread2_(const fnat m[static restrict 1], const fnat n[static restrict 1], double A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1])
{
  if (*ldA < *m)
    return -4;
  if (*fd < 0)
    return -5;
  if (!*m)
    return 0;
  if (!*n)
    return 0;

  const size_t sz = *m * sizeof(double);
  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd,sz) reduction(+:r)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    double *const Aj = A + j * (*ldA);
    const size_t of = j * (*m * sizeof(double));
    r += ((ssize_t)sz != pread(*fd, Aj, sz, of));
  }

  return r;
}

int cread2_(const fnat m[static restrict 1], const fnat n[static restrict 1], float complex A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1])
{
  if (*ldA < *m)
    return -4;
  if (*fd < 0)
    return -5;
  if (!*m)
    return 0;
  if (!*n)
    return 0;

  const size_t sz = *m * sizeof(float complex);
  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd,sz) reduction(+:r)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    float complex *const Aj = A + j * (*ldA);
    const size_t of = j * (*m * sizeof(float complex));
    r += ((ssize_t)sz != pread(*fd, Aj, sz, of));
  }

  return r;
}

int zread2_(const fnat m[static restrict 1], const fnat n[static restrict 1], double complex A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1])
{
  if (*ldA < *m)
    return -4;
  if (*fd < 0)
    return -5;
  if (!*m)
    return 0;
  if (!*n)
    return 0;

  const size_t sz = *m * sizeof(double complex);
  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd,sz) reduction(+:r)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    double complex *const Aj = A + j * (*ldA);
    const size_t of = j * (*m * sizeof(double complex));
    r += ((ssize_t)sz != pread(*fd, Aj, sz, of));
  }

  return r;
}

int swrite2_(const fnat m[static restrict 1], const fnat n[static restrict 1], const float A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1])
{
  if (*ldA < *m)
    return -4;
  if (*fd < 0)
    return -5;
  if (!*m)
    return 0;
  if (!*n)
    return 0;

  const size_t sz = *m * sizeof(float);
  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd,sz) reduction(+:r)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    const float *const Aj = A + j * (*ldA);
    const size_t of = j * (*m * sizeof(float));
    r += ((ssize_t)sz != pwrite(*fd, Aj, sz, of));
  }

  return r;
}

int dwrite2_(const fnat m[static restrict 1], const fnat n[static restrict 1], const double A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1])
{
  if (*ldA < *m)
    return -4;
  if (*fd < 0)
    return -5;
  if (!*m)
    return 0;
  if (!*n)
    return 0;

  const size_t sz = *m * sizeof(double);
  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd,sz) reduction(+:r)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    const double *const Aj = A + j * (*ldA);
    const size_t of = j * (*m * sizeof(double));
    r += ((ssize_t)sz != pwrite(*fd, Aj, sz, of));
  }

  return r;
}

int cwrite2_(const fnat m[static restrict 1], const fnat n[static restrict 1], const float complex A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1])
{
  if (*ldA < *m)
    return -4;
  if (*fd < 0)
    return -5;
  if (!*m)
    return 0;
  if (!*n)
    return 0;

  const size_t sz = *m * sizeof(float complex);
  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd,sz) reduction(+:r)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    const float complex *const Aj = A + j * (*ldA);
    const size_t of = j * (*m * sizeof(float complex));
    r += ((ssize_t)sz != pwrite(*fd, Aj, sz, of));
  }

  return r;
}

int zwrite2_(const fnat m[static restrict 1], const fnat n[static restrict 1], const double complex A[static restrict 1], const fnat ldA[static restrict 1], const int fd[static restrict 1])
{
  if (*ldA < *m)
    return -4;
  if (*fd < 0)
    return -5;
  if (!*m)
    return 0;
  if (!*n)
    return 0;

  const size_t sz = *m * sizeof(double complex);
  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd,sz) reduction(+:r)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    const double complex *const Aj = A + j * (*ldA);
    const size_t of = j * (*m * sizeof(double complex));
    r += ((ssize_t)sz != pwrite(*fd, Aj, sz, of));
  }

  return r;
}
