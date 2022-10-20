#include "sdpscl.h"

// M\o ller's 2Sum
#ifdef TwoSum
#error TwoSum already defined
#else /* !TwoSum */
#define TwoSum(a,b,a_,b_,s,t) \
  s = _mm512_add_ps(a, b);    \
  a_ = _mm512_sub_ps(s, b);   \
  b_ = _mm512_sub_ps(s, a_);  \
  a_ = _mm512_sub_ps(a, a_);  \
  b_ = _mm512_sub_ps(b, b_);  \
  t = _mm512_add_ps(a_, b_)
#endif /* ?TwoSum */

float sdpscl_(const fnat m[static restrict 1], const float x[static restrict VSL], const float y[static restrict VSL], const float e[static restrict 2], const float f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*m & VSL_1)
    return NAN;
  if (IS_NOT_ALIGNED(x))
    return NAN;
  if (IS_NOT_ALIGNED(y))
    return NAN;
  if (!(e[0u] < HUGE_VALF))
    return NAN;
  if (!(e[1u] < HUGE_VALF))
    return NAN;
  if (!(f[0u] >= 1.0f) || !(f[0u] < 2.0f))
    return NAN;
  if (!(f[1u] >= 1.0f) || !(f[0u] < 2.0f))
    return NAN;
#endif /* !NDEBUG */

  const float ex = e[0u];
  if (!(ex > -HUGE_VALF))
    return 0.0f;

  const float ey = e[1u];
  if (!(ey > -HUGE_VALF))
    return 0.0f;

  register const VS xe = _mm512_set1_ps(-ex);
  register const VS ye = _mm512_set1_ps(-ey);
  register VS spd = _mm512_setzero_ps();
  register VS sdn = _mm512_setzero_ps();
  register VS tpd = _mm512_setzero_ps();
  register VS tdn = _mm512_setzero_ps();
  register VS pd = _mm512_setzero_ps();
  register VS dn = _mm512_setzero_ps();
  register VS a = _mm512_setzero_ps();
  register VS b = _mm512_setzero_ps();
  register VS a_ = _mm512_setzero_ps();
  register VS b_ = _mm512_setzero_ps();

  // Kahan + Graillat & al.
  for (fnat i = 0u; i < *m; i += VSL) {
    register VS xi = _mm512_load_ps(x + i);
    register VS yi = _mm512_load_ps(y + i);
    xi = _mm512_scalef_ps(xi, xe); VSP(xi);
    yi = _mm512_scalef_ps(yi, ye); VSP(yi);
    pd = _mm512_mul_round_ps(xi, yi, (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)); VSP(pd);
    dn = _mm512_fmsub_ps(xi, yi, pd); VSP(dn);
    a = spd;
    b = _mm512_add_ps(pd, tpd);
    TwoSum(a,b,a_,b_,spd,tpd);
    a = sdn;
    b = _mm512_add_ps(dn, tdn);
    TwoSum(a,b,a_,b_,sdn,tdn);
  }

  const float fx = f[0u];
  const float fy = f[1u];
  return ((_mm512_reduce_add_ps(spd) + (_mm512_reduce_add_ps(tpd) + (_mm512_reduce_add_ps(sdn) + _mm512_reduce_add_ps(tdn)))) / (fx * fy));
}
