#include "serial.h"

#ifdef __ICC
#include <mathimf.h>
#else /* !__ICC */
#ifdef __cplusplus
#include <cmath>
#else /* !__cplusplus */
#include <math.h>
#endif /* ?__cplusplus */
#endif /* ?__ICC */
#ifdef __cplusplus
#include <cfloat>
#include <climits>
#else /* !__cplusplus */
#include <float.h>
#include <limits.h>
#endif /* ?__cplusplus */

#ifdef FLT_BIG_EXP
#error FLT_BIG_EXP already defined
#else /* FLT_MAX_EXP - 3 */
#define FLT_BIG_EXP 125
#endif /* ?FLT_BIG_EXP */

#ifdef FLT_SQRT_HUGE
#error FLT_SQRT_HUGE already defined
#else /* sqrtf(FLT_MAX) */
#define FLT_SQRT_HUGE 1.844674297e+19f
#endif /* ?FLT_SQRT_HUGE */

#ifdef DBL_BIG_EXP
#error DBL_BIG_EXP already defined
#else /* DBL_MAX_EXP - 3 */
#define DBL_BIG_EXP 1021
#endif /* ?DBL_BIG_EXP */

#ifdef DBL_SQRT_HUGE
#error DBL_SQRT_HUGE already defined
#else /* sqrt(DBL_MAX) */
#define DBL_SQRT_HUGE 1.34078079299425956E+154
#endif /* ?DBL_SQRT_HUGE */

#ifdef LDBL_BIG_EXP
#error LDBL_BIG_EXP already defined
#else /* LDBL_MAX_EXP - 3 */
#define LDBL_BIG_EXP 16381
#endif /* ?LDBL_BIG_EXP */

#ifdef LDBL_SQRT_HUGE
#error LDBL_SQRT_HUGE already defined
#else /* sqrtl(LDBL_MAX) */
#define LDBL_SQRT_HUGE 1.090748135619415929404E+2466
#endif /* ?LDBL_SQRT_HUGE */

/* in the routines below, -es could be returned instead of es */

int csjac2(const float a11, const float a22, const float a21r, const float a21i, float *const cs, float *const snr, float *const sni, float *const l1, float *const l2)
{
#ifndef NDEBUG
  if (!isfinite(a11))
    return INT_MIN;
  if (!isfinite(a22))
    return INT_MIN;
  if (!isfinite(a21r))
    return INT_MIN;
  if (!isfinite(a21i))
    return INT_MIN;
#endif /* !NDEBUG */
  int e1, e2, er, ei;
  (void)frexpf(fmaxf(fabsf(a11), FLT_TRUE_MIN), &e1); e1 = FLT_BIG_EXP - e1;
  (void)frexpf(fmaxf(fabsf(a22), FLT_TRUE_MIN), &e2); e2 = FLT_BIG_EXP - e2;
  (void)frexpf(fmaxf(fabsf(a21r), FLT_TRUE_MIN), &er); er = FLT_BIG_EXP - er;
  (void)frexpf(fmaxf(fabsf(a21i), FLT_TRUE_MIN), &ei); ei = FLT_BIG_EXP - ei;
  e1 = ((e1 <= e2) ? e1 : e2);
  e2 = ((er <= ei) ? er : ei);
  int es = ((e1 <= e2) ? e1 : e2);
  const float
    a1 = scalbnf(a11, es),
    a2 = scalbnf(a22, es),
    ar = scalbnf(a21r, es),
    ai = scalbnf(a21i, es);
  float
    ar_ = fabsf(ar),
    ai_ = fabsf(ai);
  const float
    am = fminf(ar_, ai_),
    aM = fmaxf(ar_, ai_);
  float aa = fmaxf(am / aM, 0.0f);
  aa = sqrtf(fmaf(aa, aa, 1.0f)) * aM;
  ar_ = copysignf(fminf(ar_ / aa, 1.0f), ar);
  ai_ = ai / fmaxf(aa, FLT_TRUE_MIN);
  const float
    an = (aa * 2.0f),
    ad = (a1 - a2),
    t2 = copysignf(fminf(fmaxf(an / fabsf(ad), 0.0f), FLT_SQRT_HUGE), ad),
    t1 = (t2 / (1.0f + sqrtf(fmaf(t2, t2, 1.0f)))),
    s2 = fmaf(t1, t1, 1.0f),
    s1 = sqrtf(s2);
  *cs = 1.0f / s1;
  *snr = (ar_ * t1) / s1;
  *sni = (ai_ * t1) / s1;
  *l1 = fmaf(t1, fmaf(a2, t1,  an), a1) / s2;
  *l2 = fmaf(t1, fmaf(a1, t1, -an), a2) / s2;
  er = ((es << 1) | (*l1 < *l2));
  es = -es;
  *l1 = scalbnf(*l1, es);
  *l2 = scalbnf(*l2, es);
  return er;
}

int zsjac2(const double a11, const double a22, const double a21r, const double a21i, double *const cs, double *const snr, double *const sni, double *const l1, double *const l2)
{
#ifndef NDEBUG
  if (!isfinite(a11))
    return INT_MIN;
  if (!isfinite(a22))
    return INT_MIN;
  if (!isfinite(a21r))
    return INT_MIN;
  if (!isfinite(a21i))
    return INT_MIN;
#endif /* !NDEBUG */
  int e1, e2, er, ei;
  (void)frexp(fmax(fabs(a11), DBL_TRUE_MIN), &e1); e1 = DBL_BIG_EXP - e1;
  (void)frexp(fmax(fabs(a22), DBL_TRUE_MIN), &e2); e2 = DBL_BIG_EXP - e2;
  (void)frexp(fmax(fabs(a21r), DBL_TRUE_MIN), &er); er = DBL_BIG_EXP - er;
  (void)frexp(fmax(fabs(a21i), DBL_TRUE_MIN), &ei); ei = DBL_BIG_EXP - ei;
  e1 = ((e1 <= e2) ? e1 : e2);
  e2 = ((er <= ei) ? er : ei);
  int es = ((e1 <= e2) ? e1 : e2);
  const double
    a1 = scalbn(a11, es),
    a2 = scalbn(a22, es),
    ar = scalbn(a21r, es),
    ai = scalbn(a21i, es);
  double
    ar_ = fabs(ar),
    ai_ = fabs(ai);
  const double
    am = fmin(ar_, ai_),
    aM = fmax(ar_, ai_);
  double aa = fmax(am / aM, 0.0);
  aa = sqrt(fma(aa, aa, 1.0)) * aM;
  ar_ = copysign(fmin(ar_ / aa, 1.0), ar);
  ai_ = ai / fmax(aa, DBL_TRUE_MIN);
  const double
    an = (aa * 2.0),
    ad = (a1 - a2),
    t2 = copysign(fmin(fmax(an / fabs(ad), 0.0), DBL_SQRT_HUGE), ad),
    t1 = (t2 / (1.0 + sqrt(fma(t2, t2, 1.0)))),
    s2 = fma(t1, t1, 1.0),
    s1 = sqrt(s2);
  *cs = 1.0 / s1;
  *snr = (ar_ * t1) / s1;
  *sni = (ai_ * t1) / s1;
  *l1 = fma(t1, fma(a2, t1,  an), a1) / s2;
  *l2 = fma(t1, fma(a1, t1, -an), a2) / s2;
  er = ((es << 1) | (*l1 < *l2));
  es = -es;
  *l1 = scalbn(*l1, es);
  *l2 = scalbn(*l2, es);
  return er;
}

int wsjac2(const long double a11, const long double a22, const long double a21r, const long double a21i, long double *const cs, long double *const snr, long double *const sni, long double *const l1, long double *const l2)
{
#ifndef NDEBUG
  if (!isfinite(a11))
    return INT_MIN;
  if (!isfinite(a22))
    return INT_MIN;
  if (!isfinite(a21r))
    return INT_MIN;
  if (!isfinite(a21i))
    return INT_MIN;
#endif /* !NDEBUG */
  int e1, e2, er, ei;
  (void)frexpl(fmaxl(fabsl(a11), LDBL_TRUE_MIN), &e1); e1 = LDBL_BIG_EXP - e1;
  (void)frexpl(fmaxl(fabsl(a22), LDBL_TRUE_MIN), &e2); e2 = LDBL_BIG_EXP - e2;
  (void)frexpl(fmaxl(fabsl(a21r), LDBL_TRUE_MIN), &er); er = LDBL_BIG_EXP - er;
  (void)frexpl(fmaxl(fabsl(a21i), LDBL_TRUE_MIN), &ei); ei = LDBL_BIG_EXP - ei;
  e1 = ((e1 <= e2) ? e1 : e2);
  e2 = ((er <= ei) ? er : ei);
  int es = ((e1 <= e2) ? e1 : e2);
  const long double
    a1 = scalbnl(a11, es),
    a2 = scalbnl(a22, es),
    ar = scalbnl(a21r, es),
    ai = scalbnl(a21i, es);
  long double
    ar_ = fabsl(ar),
    ai_ = fabsl(ai);
  const long double
    am = fminl(ar_, ai_),
    aM = fmaxl(ar_, ai_);
  long double aa = fmaxl(am / aM, 0.0L);
  aa = sqrtl(fmal(aa, aa, 1.0L)) * aM;
  ar_ = copysignl(fminl(ar_ / aa, 1.0L), ar);
  ai_ = ai / fmaxl(aa, LDBL_TRUE_MIN);
  const long double
    an = (aa * 2.0L),
    ad = (a1 - a2),
    t2 = copysignl(fminl(fmaxl(an / fabsl(ad), 0.0L), LDBL_SQRT_HUGE), ad),
    t1 = (t2 / (1.0L + sqrtl(fmal(t2, t2, 1.0L)))),
    s2 = fmal(t1, t1, 1.0L),
    s1 = sqrtl(s2);
  *cs = 1.0L / s1;
  *snr = (ar_ * t1) / s1;
  *sni = (ai_ * t1) / s1;
  *l1 = fmal(t1, fmal(a2, t1,  an), a1) / s2;
  *l2 = fmal(t1, fmal(a1, t1, -an), a2) / s2;
  er = ((es << 1) | (*l1 < *l2));
  es = -es;
  *l1 = scalbnl(*l1, es);
  *l2 = scalbnl(*l2, es);
  return er;
}

int ssjac2(const float a11, const float a22, const float a21, float *const cs, float *const sn, float *const l1, float *const l2)
{
#ifndef NDEBUG
  if (!isfinite(a11))
    return INT_MIN;
  if (!isfinite(a22))
    return INT_MIN;
  if (!isfinite(a21))
    return INT_MIN;
#endif /* !NDEBUG */
  int e1, e2, er;
  (void)frexpf(fmaxf(fabsf(a11), FLT_TRUE_MIN), &e1); e1 = FLT_BIG_EXP - e1;
  (void)frexpf(fmaxf(fabsf(a22), FLT_TRUE_MIN), &e2); e2 = FLT_BIG_EXP - e2;
  (void)frexpf(fmaxf(fabsf(a21), FLT_TRUE_MIN), &er); er = FLT_BIG_EXP - er;
  e1 = ((e1 <= e2) ? e1 : e2);
  int es = ((e1 <= er) ? e1 : er);
  const float
    a1 = scalbnf(a11, es),
    a2 = scalbnf(a22, es),
    ar = scalbnf(a21, es),
    aa = fabsf(ar),
    as = copysignf(1.0f, ar),
    an = (aa * 2.0f),
    ad = (a1 - a2),
    t2 = copysignf(fminf(fmaxf(an / fabsf(ad), 0.0f), FLT_SQRT_HUGE), ad),
    t1 = (t2 / (1.0f + sqrtf(fmaf(t2, t2, 1.0f)))),
    s2 = fmaf(t1, t1, 1.0f),
    s1 = sqrtf(s2);
  *cs = 1.0f / s1;
  *sn = (t1 * as) / s1;
  *l1 = fmaf(t1, fmaf(a2, t1,  an), a1) / s2;
  *l2 = fmaf(t1, fmaf(a1, t1, -an), a2) / s2;
  er = ((es << 1) | (*l1 < *l2));
  es = -es;
  *l1 = scalbnf(*l1, es);
  *l2 = scalbnf(*l2, es);
  return er;
}

int dsjac2(const double a11, const double a22, const double a21, double *const cs, double *const sn, double *const l1, double *const l2)
{
#ifndef NDEBUG
  if (!isfinite(a11))
    return INT_MIN;
  if (!isfinite(a22))
    return INT_MIN;
  if (!isfinite(a21))
    return INT_MIN;
#endif /* !NDEBUG */
  int e1, e2, er;
  (void)frexp(fmax(fabs(a11), DBL_TRUE_MIN), &e1); e1 = DBL_BIG_EXP - e1;
  (void)frexp(fmax(fabs(a22), DBL_TRUE_MIN), &e2); e2 = DBL_BIG_EXP - e2;
  (void)frexp(fmax(fabs(a21), DBL_TRUE_MIN), &er); er = DBL_BIG_EXP - er;
  e1 = ((e1 <= e2) ? e1 : e2);
  int es = ((e1 <= er) ? e1 : er);
  const double
    a1 = scalbn(a11, es),
    a2 = scalbn(a22, es),
    ar = scalbn(a21, es),
    aa = fabs(ar),
    as = copysign(1.0, ar),
    an = (aa * 2.0),
    ad = (a1 - a2),
    t2 = copysign(fmin(fmax(an / fabs(ad), 0.0), DBL_SQRT_HUGE), ad),
    t1 = (t2 / (1.0 + sqrt(fma(t2, t2, 1.0)))),
    s2 = fma(t1, t1, 1.0),
    s1 = sqrt(s2);
  *cs = 1.0 / s1;
  *sn = (t1 * as) / s1;
  *l1 = fma(t1, fma(a2, t1,  an), a1) / s2;
  *l2 = fma(t1, fma(a1, t1, -an), a2) / s2;
  er = ((es << 1) | (*l1 < *l2));
  es = -es;
  *l1 = scalbn(*l1, es);
  *l2 = scalbn(*l2, es);
  return er;
}

int xsjac2(const long double a11, const long double a22, const long double a21, long double *const cs, long double *const sn, long double *const l1, long double *const l2)
{
#ifndef NDEBUG
  if (!isfinite(a11))
    return INT_MIN;
  if (!isfinite(a22))
    return INT_MIN;
  if (!isfinite(a21))
    return INT_MIN;
#endif /* !NDEBUG */
  int e1, e2, er;
  (void)frexpl(fmaxl(fabsl(a11), LDBL_TRUE_MIN), &e1); e1 = LDBL_BIG_EXP - e1;
  (void)frexpl(fmaxl(fabsl(a22), LDBL_TRUE_MIN), &e2); e2 = LDBL_BIG_EXP - e2;
  (void)frexpl(fmaxl(fabsl(a21), LDBL_TRUE_MIN), &er); er = LDBL_BIG_EXP - er;
  e1 = ((e1 <= e2) ? e1 : e2);
  int es = ((e1 <= er) ? e1 : er);
  const long double
    a1 = scalbnl(a11, es),
    a2 = scalbnl(a22, es),
    ar = scalbnl(a21, es),
    aa = fabsl(ar),
    as = copysignl(1.0L, ar),
    an = (aa * 2.0L),
    ad = (a1 - a2),
    t2 = copysignl(fminl(fmaxl(an / fabsl(ad), 0.0L), LDBL_SQRT_HUGE), ad),
    t1 = (t2 / (1.0L + sqrtl(fmal(t2, t2, 1.0L)))),
    s2 = fmal(t1, t1, 1.0L),
    s1 = sqrtl(s2);
  *cs = 1.0L / s1;
  *sn = (t1 * as) / s1;
  *l1 = fmal(t1, fmal(a2, t1,  an), a1) / s2;
  *l2 = fmal(t1, fmal(a1, t1, -an), a2) / s2;
  er = ((es << 1) | (*l1 < *l2));
  es = -es;
  *l1 = scalbnl(*l1, es);
  *l2 = scalbnl(*l2, es);
  return er;
}
