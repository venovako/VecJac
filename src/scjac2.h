#ifndef SCJAC2_H
#define SCJAC2_H

#include "vecdef.h"

#ifdef FLT_BIG_EXP
#error FLT_BIG_EXP already defined
#else /* !FLT_BIG_EXP */
// (float)(FLT_MAX_EXP - 4)
#define FLT_BIG_EXP 124.0f
#endif /* ?FLT_BIG_EXP */

#ifdef FLT_SQRT_HUGE
#error FLT_SQRT_HUGE already defined
#else /* !FLT_SQRT_HUGE */
#define FLT_SQRT_HUGE 1.844674297e+19f
#endif /* ?FLT_SQRT_HUGE */

#ifdef SCJAC2_PARAMS
#error SCJAC2_PARAMS already defined
#else /* !SCJAC2_PARAMS */
#define SCJAC2_PARAMS                                      \
  register const VS  zerof = _mm512_setzero_ps();          \
  register const VS _zerof = _mm512_set1_ps(-0.0f);        \
  register const VS   onef = _mm512_set1_ps(1.0f);         \
  register const VS    shf = _mm512_set1_ps(FLT_SQRT_HUGE)
#endif /* ?SCJAC2_PARAMS */

#endif /* !SCJAC2_H */
