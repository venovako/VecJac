#ifndef COMMON_H
#define COMMON_H

#if (defined(__ICC) || defined(__INTEL_COMPILER) || defined(__INTEL_CLANG_COMPILER) || defined(__INTEL_LLVM_COMPILER))
#include <mathimf.h>
#else /* !__ICC */
#ifdef __cplusplus
#include <cmath>
#else /* !__cplusplus */
#include <math.h>
#endif /* ?__cplusplus */
#endif /* ?__ICC */

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

#ifdef __cplusplus
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cfenv>
#include <cfloat>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#else /* !__cplusplus */
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fenv.h>
#include <float.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#endif /* ?__cplusplus */

#include <fcntl.h>
#ifdef _LARGEFILE64_SOURCE
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#else /* !_LARGEFILE64_SOURCE */
#include <io.h>
#include <malloc.h>
#include <windows.h>
#endif /* ?_LARGEFILE64_SOURCE */

#ifdef EXTERN
#error EXTERN already defined
#else /* !EXTERN */
#ifdef __cplusplus
#define EXTERN extern "C"
#else /* !__cplusplus */
#define EXTERN extern
#endif /* ?__cplusplus */
#endif /* ?EXTERN */

EXTERN size_t atoz(const char *const s);

EXTERN char *stoa(char *const s, const float x);
EXTERN char *dtoa(char *const s, const double x);
EXTERN char *xtoa(char *const s, const long double x);

#endif /* !COMMON_H */
