#ifndef DZJAC2_H
#define DZJAC2_H

#include "vecdef.h"

#ifdef SQRT_HUGE
#error SQRT_HUGE already defined
#else /* !SQRT_HUGE */
#define SQRT_HUGE 1.34078079299425956E+154
#endif /* ?SQRT_HUGE */

#ifdef DZJAC2_PARAMS
#error DZJAC2_PARAMS already defined
#else /* !DZJAC2_PARAMS */
#define DZJAC2_PARAMS                               \
  register const VD zero = _mm512_setzero_pd();     \
  register const VD m0 = _mm512_set1_pd(-0.0);      \
  register const VD one = _mm512_set1_pd(1.0);      \
  register const VD huge = _mm512_set1_pd(DBL_MAX); \
  register const VD sh = _mm512_set1_pd(SQRT_HUGE)
#endif /* ?DZJAC2_PARAMS */

extern int dzjac2_pp(const fnat n, const double *const restrict s, const double *const restrict c, const double l1[static restrict 1], const double l2[static restrict 1], const unsigned p[static restrict 1], double *const restrict L1, double *const restrict L2);

#endif /* !DZJAC2_H */
