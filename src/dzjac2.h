#ifndef DZJAC2_H
#define DZJAC2_H

#include "vecdef.h"

#ifdef DBL_BIG_EXP
#error DBL_BIG_EXP already defined
#else /* !DBL_BIG_EXP */
// (double)(DBL_MAX_EXP - 4)
#define DBL_BIG_EXP 1020.0
#endif /* ?DBL_BIG_EXP */

#ifdef DBL_SQRT_HUGE
#error DBL_SQRT_HUGE already defined
#else /* !DBL_SQRT_HUGE */
#define DBL_SQRT_HUGE 1.34078079299425956E+154
#endif /* ?DBL_SQRT_HUGE */

#ifdef DZJAC2_PARAMS
#error DZJAC2_PARAMS already defined
#else /* !DZJAC2_PARAMS */
#define DZJAC2_PARAMS                                     \
  register const VD  zero = _mm512_setzero_pd();          \
  register const VD _zero = _mm512_set1_pd(-0.0);         \
  register const VD   one = _mm512_set1_pd(1.0);          \
  register const VD    sh = _mm512_set1_pd(DBL_SQRT_HUGE)
#endif /* ?DZJAC2_PARAMS */

#endif /* !DZJAC2_H */
