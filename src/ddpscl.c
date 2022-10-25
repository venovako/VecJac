#include "ddpscl.h"

#ifdef USE_2SUM
#include "d2sum.h"
#endif /* USE_2SUM */

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

  register const VD xe = _mm512_set1_pd(-ex); // 1
  register const VD ye = _mm512_set1_pd(-ey); // 1
  register VD spd = _mm512_setzero_pd();
#ifdef USE_2SUM
  register VD sdn = _mm512_setzero_pd();
  register VD tpd = _mm512_setzero_pd();
  register VD tdn = _mm512_setzero_pd();
#endif /* USE_2SUM */

  // USE_2SUM: Kahan + Graillat & al.
  for (fnat i = 0u; i < *m; i += VDL) {
    register VD xi = _mm512_load_pd(x + i);
    register VD yi = _mm512_load_pd(y + i);
    xi = _mm512_scalef_pd(xi, xe); VDP(xi); // 1
    yi = _mm512_scalef_pd(yi, ye); VDP(yi); // 1
#ifdef USE_2SUM
    register const VD pd = _mm512_mul_round_pd(xi, yi, (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)); VDP(pd); // 1
    register const VD dn = _mm512_fmsub_pd(xi, yi, pd); VDP(dn); // 1
    register VD a = spd;
    register VD b = _mm512_add_pd(pd, tpd); // 1
    register VD a_, b_;
    TwoSum(a,b,a_,b_,spd,tpd); // 6
    a = sdn;
    b = _mm512_add_pd(dn, tdn); // 1
    TwoSum(a,b,a_,b_,sdn,tdn); // 6
#else /* !USE_2SUM */
    spd = _mm512_fmadd_pd(xi, yi, spd); VDP(spd); // 1
#endif /* ?USE_2SUM */
  }

  const double fx = f[0u];
  const double fy = f[1u];
#ifdef USE_2SUM
  const double nu = (_mm512_reduce_add_pd(spd) + (_mm512_reduce_add_pd(tpd) + (_mm512_reduce_add_pd(sdn) + _mm512_reduce_add_pd(tdn)))); // 31
#else /* !USE_2SUM */
  const double nu = _mm512_reduce_add_pd(spd); // 7
#endif /* ?USE_2SUM */
  const double de = (fx * fy); // 1
  return (nu / de); // 1
}
