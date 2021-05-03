#ifndef DZJAC2_H
#define DZJAC2_H

#ifdef SQRT_HUGE
#error SQRT_HUGE already defined
#else /* !SQRT_HUGE */
#define SQRT_HUGE 1.34078079299425956E+154
#endif /* ?SQRT_HUGE */

#ifdef BIG_EXP
#error BIG_EXP already defined
#else /* !BIG_EXP */
// (double)(DBL_MAX_EXP - 3)
#define BIG_EXP 1021.0
#endif /* ?BIG_EXP */

#ifdef DZJAC2_PARAMS
#error DZJAC2_PARAMS already defined
#else /* !DZJAC2_PARAMS */
#define DZJAC2_PARAMS                               \
  register const VD zero = _mm512_setzero_pd();     \
  register const VD m0 = _mm512_set1_pd(-0.0);      \
  register const VD one = _mm512_set1_pd(1.0);      \
  register const VD huge = _mm512_set1_pd(DBL_MAX); \
  register const VD sh = _mm512_set1_pd(SQRT_HUGE); \
  register const VD be = _mm512_set1_pd(BIG_EXP)
#endif /* ?DZJAC2_PARAMS */

#endif /* !DZJAC2_H */
