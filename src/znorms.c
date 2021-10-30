#include "znorms.h"

#include "dznrms.h"

double znorms_(const fnat m[static restrict 1], const double zr[static restrict VDL], const double zi[static restrict VDL], double e0[static restrict 1], double f0[static restrict 1], double e1[static restrict 1], double f1[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -5.0;
  if (*m & VDL_1)
    return -1.0;
  if (IS_NOT_ALIGNED(zr))
    return -2.0;
  if (IS_NOT_ALIGNED(zi))
    return -3.0;
#endif /* !NDEBUG */

#ifdef USE_SLEEF
  register const VD _zero = _mm512_set1_pd(-0.0);
  __float128 rq[VDL];

  Sleef_quadx8 rv = Sleef_splatq8_avx512f(Q_ZERO);
  for (fnat i = 0u; i < *m; i += VDL) {
    register VD xi = VDABS(_mm512_load_pd(zr + i)); VDP(xi);
    VDSORT(xi); VDP(xi);
    const Sleef_quadx8 qi = Sleef_cast_from_doubleq8_avx512f(xi);
    rv = Sleef_fmaq8_u05avx512f(qi, qi, rv);
  }

  Sleef_quadx8 iv = Sleef_splatq8_avx512f(Q_ZERO);
  for (fnat i = 0u; i < *m; i += VDL) {
    register VD xi = VDABS(_mm512_load_pd(zi + i)); VDP(xi);
    VDSORT(xi); VDP(xi);
    const Sleef_quadx8 qi = Sleef_cast_from_doubleq8_avx512f(xi);
    iv = Sleef_fmaq8_u05avx512f(qi, qi, iv);
  }

  rv = Sleef_addq8_u05avx512f(rv, iv);
  Sleef_storeq8_avx512f(rq, rv);

  *rq += rq[1u];
  *rq += rq[2u];
  *rq += rq[3u];
  *rq += rq[4u];
  *rq += rq[5u];
  *rq += rq[6u];
  *rq += rq[7u];

  pquad2ef(rq, e1, f1);
  *rq = __sqrtq(*rq);
  pquad2ef(rq, e0, f0);
  return (double)*rq;
#else /* !USE_SLEEF */
  *e1 = *e0 = -HUGE_VAL;
  *f1 = *f0 = 1.0;
  return -0.0;
#endif /* ?USE_SLEEF */
}
