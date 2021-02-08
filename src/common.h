#ifndef COMMON_H
#define COMMON_H

#ifdef __ICC
#include <mathimf.h>
#ifndef USE_EXTENDED
extern __float128 __fmaq(__float128, __float128, __float128);
extern __float128 __fmaxq(__float128, __float128);
extern __float128 __fminq(__float128, __float128);
extern __float128 __hypotq(__float128, __float128);
extern __float128 __scalbq(__float128, __float128);
extern __float128 __invsqrtq(__float128);
extern __float128 __sincosq(__float128, __float128*, __float128*);
#endif /* !USE_EXTENDED */
#else /* !__ICC */
#include <complex.h>
#include <math.h>
#ifndef USE_EXTENDED
#define USE_EXTENDED
#endif /* !USE_EXTENDED */
#endif /* ?__ICC */

#ifndef CMPLXF
#define CMPLXF(r,i) ((float)(r) + I * (float)(i))
#endif /* !CMPLXF */
#ifndef CMPLX
#define CMPLX(r,i) ((double)(r) + I * (double)(i))
#endif /* !CMPLX */
#ifndef CMPLXL
#define CMPLXL(r,i) ((long double)(r) + I * (long double)(i))
#endif /* !CMPLXL */

#ifdef fmaw
#error fmaw already defined
#endif /* fmaw */
#ifdef fmaxw
#error fmaxw already defined
#endif /* fmaxw */
#ifdef fminw
#error fminw already defined
#endif /* fminw */
#ifdef hypotw
#error hypotw already defined
#endif /* hypotw */
#ifdef scalbw
#error scalbw already defined
#endif /* scalbw */
#ifdef invsqrtw
#error invsqrtw already defined
#endif /* invsqrtw */
#ifdef sincosw
#error sincosw already defined
#endif /* sincosw */
#ifdef W_ONE
#error W_ONE already defined
#endif /* W_ONE */
#ifdef W_ZERO
#error W_ZERO already defined
#endif /* W_ZERO */
#ifdef W_MONE
#error W_MONE already defined
#endif /* W_MONE */
#ifdef W_PI
#error W_PI already defined
#endif /* W_PI */
#ifdef CMPLXW
#error CMPLXW already defined
#endif /* CMPLXW */

#ifdef USE_EXTENDED
typedef long double wide;
#define fmaw     fmal
#define fmaxw    fmaxl
#define fminw    fminl
#define hypotw   hypotl
#define scalbw   scalbl
#define invsqrtw invsqrtl
#define sincosw  sincosl
#define W_ONE     1.0L
#define W_ZERO    0.0L
#define W_MONE   -1.0L
#define W_PI      3.1415926535897932384626433832795L
#define CMPLXW(r,i) CMPLXL((r),(i))
#else /* USE_QUAD */
typedef __float128  wide;
#define fmaw     __fmaq
#define fmaxw    __fmaxq
#define fminw    __fminq
#define hypotw   __hypotq
#define scalbw   __scalbq
#define invsqrtw __invsqrtq
#define sincosw  __sincosq
#define W_ONE     1.0q
#define W_ZERO    0.0q
#define W_MONE   -1.0q
#define W_PI      3.1415926535897932384626433832795q
#define CMPLXW(r,i) ((wide)(r) + I * (wide)(i))
#endif /* ?USE_EXTENDED */

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fenv.h>
#include <float.h>
#include <limits.h>
#include <signal.h>
#include <stdalign.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#ifdef FINT_C
#error FINT_C already defined
#endif /* FINT_C */
#ifdef FNAT_C
#error FNAT_C already defined
#endif /* FNAT_C */

#ifdef MKL_ILP64
typedef  int64_t fint;
typedef uint64_t fnat;
#define FINT_C(x) INT64_C(x)
#define FNAT_C(x) UINT64_C(x)
#else /* LP64 */
typedef  int32_t fint;
typedef uint32_t fnat;
#define FINT_C(x) INT32_C(x)
#define FNAT_C(x) UINT32_C(x)
#endif /* ?MKL_ILP64 */

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

#ifdef BLAS_S
#error BLAS_S already defined
#endif /* BLAS_S */
#ifdef BLAS_D
#error BLAS_D already defined
#endif /* BLAS_D */
#ifdef BLAS_C
#error BLAS_C already defined
#endif /* BLAS_C */
#ifdef BLAS_Z
#error BLAS_Z already defined
#endif /* BLAS_Z */

#ifdef USE_MKL
#ifndef MKL_Complex8
#define MKL_Complex8 float _Complex
#endif /* !MKL_Complex8 */
#ifndef MKL_Complex16
#define MKL_Complex16 double _Complex
#endif /* !MKL_Complex16 */
#include <mkl.h>
#define BLAS_S(name) s##name
#define BLAS_D(name) d##name
#define BLAS_C(name) c##name
#define BLAS_Z(name) z##name
#else /* some other Fortran-compatible BLAS */
#define BLAS_S(name) s##name##_
#define BLAS_D(name) d##name##_
#define BLAS_C(name) c##name##_
#define BLAS_Z(name) z##name##_
#endif /* ?USE_MKL */

#ifdef LAPACK_S
#error LAPACK_S already defined
#else /* !LAPACK_S */
#define LAPACK_S(name) s##name##_
#endif /* ?LAPACK_S */
#ifdef LAPACK_D
#error LAPACK_D already defined
#else /* !LAPACK_D */
#define LAPACK_D(name) d##name##_
#endif /* ?LAPACK_D */
#ifdef LAPACK_C
#error LAPACK_C already defined
#else /* !LAPACK_C */
#define LAPACK_C(name) c##name##_
#endif /* ?LAPACK_C */
#ifdef LAPACK_Z
#error LAPACK_Z already defined
#else /* !LAPACK_Z */
#define LAPACK_Z(name) z##name##_
#endif /* ?LAPACK_Z */

extern size_t atoz(const char *const s);

extern char *stoa(char s[static 17], const float x);
extern char *dtoa(char s[static 26], const double x);
extern char *xtoa(char s[static 31], const long double x);

#endif /* !COMMON_H */
