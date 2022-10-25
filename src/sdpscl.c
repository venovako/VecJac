#include "sdpscl.h"

#ifdef USE_2SUM
#include "s2sum.h"
#endif /* USE_2SUM */

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
  if (!(e[0u] <= FLT_MAX))
    return NAN;
  if (!(e[1u] <= FLT_MAX))
    return NAN;
  if (!(f[0u] >= 1.0f) || !(f[0u] < 2.0f))
    return NAN;
  if (!(f[1u] >= 1.0f) || !(f[0u] < 2.0f))
    return NAN;
#endif /* !NDEBUG */

  const float ex = e[0u];
  if (!(ex >= -FLT_MAX))
    return 0.0f;

  const float ey = e[1u];
  if (!(ey >= -FLT_MAX))
    return 0.0f;

  register const VS xe = _mm512_set1_ps(-ex);
  register const VS ye = _mm512_set1_ps(-ey);
  register VS spd = _mm512_setzero_ps();
#ifdef USE_2SUM
  register VS sdn = _mm512_setzero_ps();
  register VS tpd = _mm512_setzero_ps();
  register VS tdn = _mm512_setzero_ps();
#endif /* USE_2SUM */

  // USE_2SUM: Kahan + Graillat & al.
  for (fnat i = 0u; i < *m; i += VSL) {
    register VS xi = _mm512_load_ps(x + i);
    register VS yi = _mm512_load_ps(y + i);
    xi = _mm512_scalef_ps(xi, xe); VSP(xi);
    yi = _mm512_scalef_ps(yi, ye); VSP(yi);
#ifdef USE_2SUM
    register const VS pd = _mm512_mul_round_ps(xi, yi, (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)); VSP(pd);
    register const VS dn = _mm512_fmsub_ps(xi, yi, pd); VSP(dn);
    register VS a = spd;
    register VS b = _mm512_add_ps(pd, tpd);
    register VS a_, b_;
    TwoSum(a,b,a_,b_,spd,tpd);
    a = sdn;
    b = _mm512_add_ps(dn, tdn);
    TwoSum(a,b,a_,b_,sdn,tdn);
#else /* !USE_2SUM */
    spd = _mm512_fmadd_ps(xi, yi, spd); VSP(spd);
#endif /* ?USE_2SUM */
  }

  const float fx = f[0u];
  const float fy = f[1u];
#ifdef USE_2SUM
  const float nu = (_mm512_reduce_add_ps(spd) + (_mm512_reduce_add_ps(tpd) + (_mm512_reduce_add_ps(sdn) + _mm512_reduce_add_ps(tdn))));
#else /* !USE_2SUM */
  const float nu = _mm512_reduce_add_ps(spd);
#endif /* ?USE_2SUM */
  const float de = (fx * fy);
  return (nu / de);
}
