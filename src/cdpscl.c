#include "cdpscl.h"

#ifdef USE_2SUM
#include "s2sum.h"
#include "vecdef.h"
#ifdef INMULSUM
#error INMULSUM already defined
#else /* !INMULSUM */
#define INMULSUM(x,y,s,t,s_,t_)       \
  xy = _mm512_mul_round_ps(x, y, rni);\
  a = s;                              \
  xy_ = _mm512_fmsub_ps(x, y, xy);    \
  b = _mm512_add_ps(xy, t);           \
  TwoSum(a,b,a_,b_,s,t);              \
  a = s_;                             \
  b = _mm512_add_ps(xy_, t_);         \
  TwoSum(a,b,a_,b_,s_,t_)
#endif /* ?INMULSUM */
#ifdef OUTSUM
#error OUTSUM already defined
#else /* !OUTSUM */
#define OUTSUM(x,y)    \
  a = x;               \
  b = y;               \
  TwoSum(a,b,a_,b_,x,y)
#endif /* ?OUTSUM */
#endif /* USE_2SUM */

float complex cdpscl_(const fnat m[static restrict 1], const float xr[static restrict VSL], const float xi[static restrict VSL], const float yr[static restrict VSL], const float yi[static restrict VSL], const float e[static restrict 2], const float f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*m & VSL_1)
    return NAN;
  if (IS_NOT_ALIGNED(xr))
    return NAN;
  if (IS_NOT_ALIGNED(xi))
    return NAN;
  if (IS_NOT_ALIGNED(yr))
    return NAN;
  if (IS_NOT_ALIGNED(yi))
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
  // USE_2SUM: Kahan + Graillat & al.
#ifdef USE_2SUM
  register const VS _zerof = _mm512_set1_ps(-0.0f);
  register VS sXrYr = _mm512_setzero_ps();
  register VS tXrYr = _mm512_setzero_ps();
  register VS sXiYi = _mm512_setzero_ps();
  register VS tXiYi = _mm512_setzero_ps();
  register VS sXrYi = _mm512_setzero_ps();
  register VS tXrYi = _mm512_setzero_ps();
  register VS sXiYr = _mm512_setzero_ps();
  register VS tXiYr = _mm512_setzero_ps();
  register VS sXrYr_ = _mm512_setzero_ps();
  register VS tXrYr_ = _mm512_setzero_ps();
  register VS sXiYi_ = _mm512_setzero_ps();
  register VS tXiYi_ = _mm512_setzero_ps();
  register VS sXrYi_ = _mm512_setzero_ps();
  register VS tXrYi_ = _mm512_setzero_ps();
  register VS sXiYr_ = _mm512_setzero_ps();
  register VS tXiYr_ = _mm512_setzero_ps();
  register VS a, b, a_, b_, xy, xy_;
  const int rni = (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
  for (fnat i = 0u; i < *m; i += VSL) {
    register VS xri = _mm512_load_ps(xr + i);
    register VS yri = _mm512_load_ps(yr + i);
    register VS yii = _mm512_load_ps(yi + i);
    register VS xii = _mm512_load_ps(xi + i);
    xri = _mm512_scalef_ps(xri, xe); VSP(xri);
    yri = _mm512_scalef_ps(yri, ye); VSP(yri);
    yii = _mm512_scalef_ps(yii, ye); VSP(yii);
    xii = _mm512_scalef_ps(xii, xe); VSP(xii);
    INMULSUM(xri,yri,sXrYr,tXrYr,sXrYr_,tXrYr_);
    INMULSUM(xii,yii,sXiYi,tXiYi,sXiYi_,tXiYi_);
    xii = VSNEG(xii); VSP(xii);
    INMULSUM(xri,yii,sXrYi,tXrYi,sXrYi_,tXrYi_);
    INMULSUM(xii,yri,sXiYr,tXiYr,sXiYr_,tXiYr_);
  }
#else /* !USE_2SUM */
  register VS pr = _mm512_setzero_ps();
  register VS pi = _mm512_setzero_ps();
  for (fnat i = 0u; i < *m; i += VSL) {
    register VS xri = _mm512_load_ps(xr + i);
    register VS yri = _mm512_load_ps(yr + i);
    register VS yii = _mm512_load_ps(yi + i);
    register VS xii = _mm512_load_ps(xi + i);
    xri = _mm512_scalef_ps(xri, xe); VSP(xri);
    yri = _mm512_scalef_ps(yri, ye); VSP(yri);
    yii = _mm512_scalef_ps(yii, ye); VSP(yii);
    xii = _mm512_scalef_ps(xii, xe); VSP(xii);
    pr = _mm512_fmadd_ps(xri, yri, pr); VSP(pr);
    pi = _mm512_fmadd_ps(xri, yii, pi); VSP(pi);
    pr = _mm512_fmadd_ps(xii, yii, pr); VSP(pr);
    pi = _mm512_fnmadd_ps(xii, yri, pi); VSP(pi);
  }
#endif /* ?USE_2SUM */

  const float fx = f[0u];
  const float fy = f[1u];
#ifdef USE_2SUM
  OUTSUM(tXrYr_,tXiYi_);
  OUTSUM(sXrYr_,sXiYi_);
  OUTSUM(tXrYr,tXiYi);
  OUTSUM(sXrYr,sXiYi);
  const float nur =
    _mm512_reduce_add_ps(sXrYr) + (_mm512_reduce_add_ps(sXiYi) + (_mm512_reduce_add_ps(tXrYr) + (_mm512_reduce_add_ps(tXiYi) + (_mm512_reduce_add_ps(sXrYr_) + (_mm512_reduce_add_ps(sXiYi_) + (_mm512_reduce_add_ps(tXrYr_) + _mm512_reduce_add_ps(tXiYi_)))))));
  OUTSUM(tXrYi_,tXiYr_);
  OUTSUM(sXrYi_,sXiYr_);
  OUTSUM(tXrYi,tXiYr);
  OUTSUM(sXrYi,sXiYr);
  const float nui =
    _mm512_reduce_add_ps(sXrYi) + (_mm512_reduce_add_ps(sXiYr) + (_mm512_reduce_add_ps(tXrYi) + (_mm512_reduce_add_ps(tXiYr) + (_mm512_reduce_add_ps(sXrYi_) + (_mm512_reduce_add_ps(sXiYr_) + (_mm512_reduce_add_ps(tXrYi_) + _mm512_reduce_add_ps(tXiYr_)))))));
#else /* !USE_2SUM */
  const float nur = _mm512_reduce_add_ps(pr);
  const float nui = _mm512_reduce_add_ps(pi);
#endif /* ?USE_2SUM */
  const float de = (fx * fy);
  alignas(16u) float complex z[2u];
  _mm_store_ps((float*)z, _mm_div_ps(_mm_set_ps(0.0f, 1.0f, nui, nur), _mm_set1_ps(de)));
  return *z;
}
