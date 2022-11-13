#include "zdpscl.h"

#ifdef USE_2SUM
#include "d2sum.h"
#include "vecdef.h"
#ifdef INMULSUM
#error INMULSUM already defined
#else /* !INMULSUM */
#define INMULSUM(x,y,s,t,s_,t_)       \
  xy = _mm512_mul_round_pd(x, y, rni);\
  a = s;                              \
  xy_ = _mm512_fmsub_pd(x, y, xy);    \
  b = _mm512_add_pd(xy, t);           \
  TwoSum(a,b,a_,b_,s,t);              \
  a = s_;                             \
  b = _mm512_add_pd(xy_, t_);         \
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

double complex zdpscl_(const fnat m[static restrict 1], const double xr[static restrict VDL], const double xi[static restrict VDL], const double yr[static restrict VDL], const double yi[static restrict VDL], const double e[static restrict 2], const double f[static restrict 2])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return NAN;
  if (*m & VDL_1)
    return NAN;
  if (IS_NOT_ALIGNED(xr))
    return NAN;
  if (IS_NOT_ALIGNED(xi))
    return NAN;
  if (IS_NOT_ALIGNED(yr))
    return NAN;
  if (IS_NOT_ALIGNED(yi))
    return NAN;
  if (!(e[0u] <= DBL_MAX))
    return NAN;
  if (!(e[1u] <= DBL_MAX))
    return NAN;
  if (!(f[0u] >= 1.0) || !(f[0u] < 2.0))
    return NAN;
  if (!(f[1u] >= 1.0) || !(f[0u] < 2.0))
    return NAN;
#endif /* !NDEBUG */

  const double ex = e[0u];
  if (!(ex >= -DBL_MAX))
    return 0.0;

  const double ey = e[1u];
  if (!(ey >= -DBL_MAX))
    return 0.0;

  register const VD xe = _mm512_set1_pd(-ex);
  register const VD ye = _mm512_set1_pd(-ey);

  // USE_2SUM: Kahan + Graillat & al.
#ifdef USE_2SUM
  register const VD _zero = _mm512_set1_pd(-0.0);
  register VD sXrYr = _mm512_setzero_pd();
  register VD tXrYr = _mm512_setzero_pd();
  register VD sXiYi = _mm512_setzero_pd();
  register VD tXiYi = _mm512_setzero_pd();
  register VD sXrYi = _mm512_setzero_pd();
  register VD tXrYi = _mm512_setzero_pd();
  register VD sXiYr = _mm512_setzero_pd();
  register VD tXiYr = _mm512_setzero_pd();
  register VD sXrYr_ = _mm512_setzero_pd();
  register VD tXrYr_ = _mm512_setzero_pd();
  register VD sXiYi_ = _mm512_setzero_pd();
  register VD tXiYi_ = _mm512_setzero_pd();
  register VD sXrYi_ = _mm512_setzero_pd();
  register VD tXrYi_ = _mm512_setzero_pd();
  register VD sXiYr_ = _mm512_setzero_pd();
  register VD tXiYr_ = _mm512_setzero_pd();
  register VD a, b, a_, b_, xy, xy_;
  const int rni = (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
  for (fnat i = 0u; i < *m; i += VDL) {
    register VD xri = _mm512_load_pd(xr + i);
    register VD yri = _mm512_load_pd(yr + i);
    register VD yii = _mm512_load_pd(yi + i);
    register VD xii = _mm512_load_pd(xi + i);
    xri = _mm512_scalef_pd(xri, xe); VDP(xri);
    yri = _mm512_scalef_pd(yri, ye); VDP(yri);
    yii = _mm512_scalef_pd(yii, ye); VDP(yii);
    xii = _mm512_scalef_pd(xii, xe); VDP(xii);
    INMULSUM(xri,yri,sXrYr,tXrYr,sXrYr_,tXrYr_);
    INMULSUM(xii,yii,sXiYi,tXiYi,sXiYi_,tXiYi_);
    xii = VDNEG(xii); VDP(xii);
    INMULSUM(xri,yii,sXrYi,tXrYi,sXrYi_,tXrYi_);
    INMULSUM(xii,yri,sXiYr,tXiYr,sXiYr_,tXiYr_);
  }
#else /* !USE_2SUM */
  register VD pr = _mm512_setzero_pd();
  register VD pi = _mm512_setzero_pd();
  for (fnat i = 0u; i < *m; i += VDL) {
    register VD xri = _mm512_load_pd(xr + i);
    register VD yri = _mm512_load_pd(yr + i);
    register VD yii = _mm512_load_pd(yi + i);
    register VD xii = _mm512_load_pd(xi + i);
    xri = _mm512_scalef_pd(xri, xe); VDP(xri);
    yri = _mm512_scalef_pd(yri, ye); VDP(yri);
    yii = _mm512_scalef_pd(yii, ye); VDP(yii);
    xii = _mm512_scalef_pd(xii, xe); VDP(xii);
    pr = _mm512_fmadd_pd(xri, yri, pr); VDP(pr);
    pi = _mm512_fmadd_pd(xri, yii, pi); VDP(pi);
    pr = _mm512_fmadd_pd(xii, yii, pr); VDP(pr);
    pi = _mm512_fnmadd_pd(xii, yri, pi); VDP(pi);
  }
#endif /* ?USE_2SUM */

  const double fx = f[0u];
  const double fy = f[1u];
#ifdef USE_2SUM
  OUTSUM(tXrYr_,tXiYi_);
  OUTSUM(sXrYr_,sXiYi_);
  OUTSUM(tXrYr,tXiYi);
  OUTSUM(sXrYr,sXiYi);
  const double nur =
    _mm512_reduce_add_pd(sXrYr) + (_mm512_reduce_add_pd(sXiYi) + (_mm512_reduce_add_pd(tXrYr) + (_mm512_reduce_add_pd(tXiYi) + (_mm512_reduce_add_pd(sXrYr_) + (_mm512_reduce_add_pd(sXiYi_) + (_mm512_reduce_add_pd(tXrYr_) + _mm512_reduce_add_pd(tXiYi_)))))));
  OUTSUM(tXrYi_,tXiYr_);
  OUTSUM(sXrYi_,sXiYr_);
  OUTSUM(tXrYi,tXiYr);
  OUTSUM(sXrYi,sXiYr);
  const double nui =
    _mm512_reduce_add_pd(sXrYi) + (_mm512_reduce_add_pd(sXiYr) + (_mm512_reduce_add_pd(tXrYi) + (_mm512_reduce_add_pd(tXiYr) + (_mm512_reduce_add_pd(sXrYi_) + (_mm512_reduce_add_pd(sXiYr_) + (_mm512_reduce_add_pd(tXrYi_) + _mm512_reduce_add_pd(tXiYr_)))))));
#else /* !USE_2SUM */
  const double nur = _mm512_reduce_add_pd(pr);
  const double nui = _mm512_reduce_add_pd(pi);
#endif /* ?USE_2SUM */
  const double de = (fx * fy);
  alignas(16u) double complex z;
  _mm_store_pd((double*)&z, _mm_div_pd(_mm_set_pd(nui, nur), _mm_set1_pd(de)));
  return z;
}
