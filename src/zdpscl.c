#include "zdpscl.h"

#ifdef USE_2SUM
#include "d2sum.h"
#include "vecdef.h"
#ifdef MUL2SUM
#error MUL2SUM already defined
#else /* !MUL2SUM */
#define MUL2SUM(x,y,s,t,s_,t_)                                                \
  xy = _mm512_mul_round_pd(x, y, (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));\
  a = s;                                                                      \
  xy_ = _mm512_fmsub_pd(x, y, xy);                                            \
  b = _mm512_add_pd(xy, t);                                                   \
  TwoSum(a,b,a_,b_,s,t);                                                      \
  a = s_;                                                                     \
  b = _mm512_add_pd(xy_, t_);                                                 \
  TwoSum(a,b,a_,b_,s_,t_)
#endif /* ?MUL2SUM */
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
  register VD sr = _mm512_setzero_pd();
  register VD tr = _mm512_setzero_pd();
  register VD sr_ = _mm512_setzero_pd();
  register VD tr_ = _mm512_setzero_pd();
  register VD si = _mm512_setzero_pd();
  register VD ti = _mm512_setzero_pd();
  register VD si_ = _mm512_setzero_pd();
  register VD ti_ = _mm512_setzero_pd();
  for (fnat i = 0u; i < *m; i += VDL) {
    register VD xri = _mm512_load_pd(xr + i);
    register VD yri = _mm512_load_pd(yr + i);
    register VD yii = _mm512_load_pd(yi + i);
    register VD xii = _mm512_load_pd(xi + i);
    xri = _mm512_scalef_pd(xri, xe); VDP(xri);
    yri = _mm512_scalef_pd(yri, ye); VDP(yri);
    yii = _mm512_scalef_pd(yii, ye); VDP(yii);
    xii = _mm512_scalef_pd(xii, xe); VDP(xii);
    register VD a, b, a_, b_, xy, xy_;
    MUL2SUM(xri,yri,sr,tr,sr_,tr_);
    MUL2SUM(xii,yii,sr,tr,sr_,tr_);
    xii = VDNEG(xii); VDP(xii);
    MUL2SUM(xri,yii,si,ti,si_,ti_);
    MUL2SUM(xii,yri,si,ti,si_,ti_);
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
  const double nur = (_mm512_reduce_add_pd(sr) + (_mm512_reduce_add_pd(tr) + (_mm512_reduce_add_pd(sr_) + _mm512_reduce_add_pd(tr_))));
  const double nui = (_mm512_reduce_add_pd(si) + (_mm512_reduce_add_pd(ti) + (_mm512_reduce_add_pd(si_) + _mm512_reduce_add_pd(ti_))));
#else /* !USE_2SUM */
  const double nur = _mm512_reduce_add_pd(pr);
  const double nui = _mm512_reduce_add_pd(pi);
#endif /* ?USE_2SUM */
  const double de = (fx * fy);
  alignas(16u) double complex z;
  _mm_store_pd((double*)&z, _mm_div_pd(_mm_set_pd(nui, nur), _mm_set1_pd(de)));
  return z;
}
