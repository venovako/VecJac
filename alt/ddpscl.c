#include "ddpscl.h"

// M\o ller's 2Sum
#ifdef TwoSum
#error TwoSum already defined
#else /* !TwoSum */
#define TwoSum(a,b,a_,b_,s,t) \
  s = _mm512_add_pd(a, b);    \
  a_ = _mm512_sub_pd(s, b);   \
  b_ = _mm512_sub_pd(s, a_);  \
  a_ = _mm512_sub_pd(a, a_);  \
  b_ = _mm512_sub_pd(b, b_);  \
  t = _mm512_add_pd(a_, b_)
#endif /* ?TwoSum */

double ddpscl_(const fnat m[static restrict 1], const double x[static restrict VDL], const double y[static restrict VDL], const double e[static restrict 2], const double f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*m & VDL_1)
    return NAN;
  if (IS_NOT_ALIGNED(x))
    return NAN;
  if (IS_NOT_ALIGNED(y))
    return NAN;
  if (!(e[0u] < HUGE_VAL))
    return NAN;
  if (!(e[1u] < HUGE_VAL))
    return NAN;
  if (!(f[0u] >= 1.0) || !(f[0u] < 2.0))
    return NAN;
  if (!(f[1u] >= 1.0) || !(f[0u] < 2.0))
    return NAN;
#endif /* !NDEBUG */

  const double ex = e[0u];
  if (!(ex > -HUGE_VAL))
    return 0.0;

  const double ey = e[1u];
  if (!(ey > -HUGE_VAL))
    return 0.0;

  register const VD xe = _mm512_set1_pd(-ex);
  register const VD ye = _mm512_set1_pd(-ey);
  register VD spd = _mm512_setzero_pd();
  register VD sdn = _mm512_setzero_pd();
  register VD tpd = _mm512_setzero_pd();
  register VD tdn = _mm512_setzero_pd();
  register VD pd = _mm512_setzero_pd();
  register VD dn = _mm512_setzero_pd();
  register VD a = _mm512_setzero_pd();
  register VD b = _mm512_setzero_pd();
  register VD a_ = _mm512_setzero_pd();
  register VD b_ = _mm512_setzero_pd();

  // Kahan + Graillat & al.
  for (fnat i = 0u; i < *m; i += VDL) {
    register VD xi = _mm512_load_pd(x + i);
    register VD yi = _mm512_load_pd(y + i);
    xi = _mm512_scalef_pd(xi, xe); VDP(xi);
    yi = _mm512_scalef_pd(yi, ye); VDP(yi);
    pd = _mm512_mul_round_pd(xi, yi, (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)); VDP(pd);
    dn = _mm512_fmsub_pd(xi, yi, pd); VDP(dn);
    a = spd;
    b = _mm512_add_pd(pd, tpd);
    TwoSum(a,b,a_,b_,spd,tpd);
    a = sdn;
    b = _mm512_add_pd(dn, tdn);
    TwoSum(a,b,a_,b_,sdn,tdn);
  }

  const double fx = f[0u];
  const double fy = f[1u];
  return ((_mm512_reduce_add_pd(spd) + (_mm512_reduce_add_pd(tpd) + (_mm512_reduce_add_pd(sdn) + _mm512_reduce_add_pd(tdn)))) / (fx * fy));
}
