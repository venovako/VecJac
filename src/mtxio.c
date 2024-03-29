#include "mtxio.h"

#ifdef FNCAT
#error FNCAT already defined
#else /* !FNCAT */
#define FNCAT(nn,bn,ex)                   \
  char nn[strlen(bn) + strlen(ex) + 2u];  \
  strcat(strcat(strcpy(nn, bn), "."), ex)
#endif /* ?FNCAT */

int open_ro_(const char bn[static restrict 1], const char ex[static restrict 1])
{
  FNCAT(fn, bn, ex);
  int oflag = O_RDONLY;
#ifdef _LARGEFILE64_SOURCE
  oflag |= O_LARGEFILE;
#endif /* _LARGEFILE64_SOURCE */
  return open(fn, oflag);
}

int open_wo_(const char bn[static restrict 1], const char ex[static restrict 1])
{
  FNCAT(fn, bn, ex);
  int oflag = (O_WRONLY | O_CREAT);
#ifdef _LARGEFILE64_SOURCE
  oflag |= O_LARGEFILE;
#endif /* _LARGEFILE64_SOURCE */
  return open(fn, oflag, (S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH));
}

int resizef_(const int fd[static restrict 1], const size_t sz[static restrict 1])
{
  return ((*fd >= 0) ? (ftruncate(*fd, *sz) || fsync(*fd)) : *fd);
}

#ifdef READ_LOOP
#error READ_LOOP already defined
#else /* !READ_LOOP */
#define READ_LOOP(T)                              \
  for (fnat j = 0u; j < *n; ++j) {                \
    T *const Aj = A + j * (size_t)(*ldA);         \
    const size_t sz = *m * sizeof(T);             \
    const size_t of = j * sz;                     \
    r += ((ssize_t)sz != pread(*fd, Aj, sz, of)); \
  }
#endif /* ?READ_LOOP */

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

  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd) reduction(+:r)
#endif /* _OPENMP */
  READ_LOOP(float);

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

  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd) reduction(+:r)
#endif /* _OPENMP */
  READ_LOOP(double);

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

  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd) reduction(+:r)
#endif /* _OPENMP */
  READ_LOOP(float complex);

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

  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd) reduction(+:r)
#endif /* _OPENMP */
  READ_LOOP(double complex);

  return r;
}

int wwrite1_(const fnat m[static restrict 1], const wide w[static restrict 1], const int fd[static restrict 1])
{
  const size_t sz = (*m * sizeof(wide));
  return (sz ? ((ssize_t)sz != write(*fd, w, sz)) : 0);
}

#ifdef WRITE_LOOP
#error WRITE_LOOP already defined
#else /* !WRITE_LOOP */
#define WRITE_LOOP(T)                              \
  for (fnat j = 0u; j < *n; ++j) {                 \
    const T *const Aj = A + j * (size_t)(*ldA);    \
    const size_t sz = *m * sizeof(T);              \
    const size_t of = j * sz;                      \
    r += ((ssize_t)sz != pwrite(*fd, Aj, sz, of)); \
  }
#endif /* ?WRITE_LOOP */

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

  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd) reduction(+:r)
#endif /* _OPENMP */
  WRITE_LOOP(float);

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

  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd) reduction(+:r)
#endif /* _OPENMP */
  WRITE_LOOP(double);

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

  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd) reduction(+:r)
#endif /* _OPENMP */
  WRITE_LOOP(float complex);

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

  int r = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,A,ldA,fd) reduction(+:r)
#endif /* _OPENMP */
  WRITE_LOOP(double complex);

  return r;
}
