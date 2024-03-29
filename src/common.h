#ifndef COMMON_H
#define COMMON_H

#if (defined(__ICC) || defined(__INTEL_COMPILER) || defined(__INTEL_CLANG_COMPILER) || defined(__INTEL_LLVM_COMPILER))
#include <mathimf.h>
#ifndef USE_EXTENDED
extern __float128 __copysignq(__float128, __float128);
extern __float128 __fabsq(__float128);
extern __float128 __fmaq(__float128, __float128, __float128);
extern __float128 __fmaxq(__float128, __float128);
extern __float128 __fminq(__float128, __float128);
extern __float128 __frexpq(__float128, int*);
extern __float128 __hypotq(__float128, __float128);
extern __float128 __scalbq(__float128, __float128);
extern __float128 __sqrtq(__float128);
#endif /* !USE_EXTENDED */
#else /* !__ICC */
#include <complex.h>
#include <math.h>
#ifndef USE_EXTENDED
#define USE_EXTENDED
#endif /* !USE_EXTENDED */
#endif /* ?__ICC */

#ifdef USE_CR_MATH
extern float cr_hypotf(float x, float y);
extern double cr_hypot(double x, double y);
#else /* !USE_CR_MATH */
#ifndef cr_hypotf
#define cr_hypotf hypotf
#endif /* !cr_hypotf */
#ifndef cr_hypot
#define cr_hypot hypot
#endif /* !cr_hypot */
#endif /* ?USE_CR_MATH */

#ifndef CMPLXF
#define CMPLXF(r,i) ((float)(r) + I * (float)(i))
#endif /* !CMPLXF */
#ifndef CMPLX
#define CMPLX(r,i) ((double)(r) + I * (double)(i))
#endif /* !CMPLX */
#ifndef CMPLXL
#define CMPLXL(r,i) ((long double)(r) + I * (long double)(i))
#endif /* !CMPLXL */

#ifdef copysignw
#error copysignw already defined
#endif /* copysignw */
#ifdef fabsw
#error fabsw already defined
#endif /* fabsw */
#ifdef fmaw
#error fmaw already defined
#endif /* fmaw */
#ifdef fmaxw
#error fmaxw already defined
#endif /* fmaxw */
#ifdef fminw
#error fminw already defined
#endif /* fminw */
#ifdef frexpw
#error frexpw already defined
#endif /* frexpw */
#ifdef hypotw
#error hypotw already defined
#endif /* hypotw */
#ifdef scalbw
#error scalbw already defined
#endif /* scalbw */
#ifdef sqrtw
#error sqrtw already defined
#endif /* sqrtw */

#ifdef W_ONE
#error W_ONE already defined
#endif /* W_ONE */
#ifdef W_ZERO
#error W_ZERO already defined
#endif /* W_ZERO */
#ifdef W_SQRT2
#error W_SQRT2 already defined
#endif /* W_SQRT2 */
#ifdef CMPLXW
#error CMPLXW already defined
#endif /* CMPLXW */

#ifdef USE_EXTENDED
typedef long double wide;
#define copysignw   copysignl
#define fabsw       fabsl
#define fmaw        fmal
#define fmaxw       fmaxl
#define fminw       fminl
#define frexpw      frexpl
#define hypotw      hypotl
#define scalbw      scalbl
#define sqrtw       sqrtl
#define W_ONE        1.0L
#define W_ZERO       0.0L
#define W_SQRT2      1.4142135623730950488016887242097L
#define CMPLXW(r,i) CMPLXL((r),(i))
#else /* USE_QUAD */
typedef __float128  wide;
#define copysignw   __copysignq
#define fabsw       __fabsq
#define fmaw        __fmaq
#define fmaxw       __fmaxq
#define fminw       __fminq
#define frexpw      __frexpq
#define hypotw      __hypotq
#define scalbw      __scalbq
#define sqrtw       __sqrtq
#define W_ONE        1.0q
#define W_ZERO       0.0q
#define W_SQRT2      1.4142135623730950488016887242097q
#define CMPLXW(r,i) ((wide)(r) + I * (wide)(i))
#endif /* ?USE_EXTENDED */

#include <alloca.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

#ifdef BLAS_I
#error BLAS_I already defined
#endif /* BLAS_I */
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
#ifdef MKL_INT
typedef  MKL_INT fint;
#else /* !MKL_INT */
#error MKL_INT not defined
#endif /* ?MKL_INT */
#ifdef MKL_UINT
typedef MKL_UINT fnat;
#else /* !MKL_UINT */
#error MKL_UINT not defined
#endif /* !MKL_UINT */
#define BLAS_I(name) i##name
#define BLAS_S(name) s##name
#define BLAS_D(name) d##name
#define BLAS_C(name) c##name
#define BLAS_Z(name) z##name
#else /* some other Fortran-compatible BLAS */
// take the same definitions as MKL
#ifdef MKL_ILP64
typedef          long long fint;
typedef unsigned long long fnat;
#else /* LP64 */
typedef      int fint;
typedef unsigned fnat;
#endif /* ?MKL_ILP64 */
#define BLAS_I(name) i##name##_
#define BLAS_S(name) s##name##_
#define BLAS_D(name) d##name##_
#define BLAS_C(name) c##name##_
#define BLAS_Z(name) z##name##_
#endif /* ?USE_MKL */

#ifdef LAPACK_I
#error LAPACK_I already defined
#else /* !LAPACK_I */
#define LAPACK_I(name) i##name##_
#endif /* ?LAPACK_I */
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

static inline fint imax(const fint a, const fint b)
{
  return ((a >= b) ? a : b);
}

extern size_t atoz(const char *const s);

extern char *stoa(char s[static restrict 17], const float x);
extern char *dtoa(char s[static restrict 26], const double x);
extern char *xtoa(char s[static restrict 31], const long double x);

extern int set_cbwr();

#endif /* !COMMON_H */
