#include "dnorms.h"

#include "dznrms.h"

double dnorms_(const fnat m[static restrict 1], const double x[static restrict VDL], double e0[static restrict 1], double f0[static restrict 1], double e1[static restrict 1], double f1[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -7.0;
  if (*m & VDL_1)
    return -1.0;
  if (IS_NOT_ALIGNED(x))
    return -2.0;
#endif /* !NDEBUG */

  register const VD _zero = _mm512_set1_pd(-0.0);
  __float128 rq[VDL];

  Sleef_quadx8 rv = Sleef_splatq8_avx512f(Q_ZERO);
  for (fnat i = 0u; i < *m; i += VDL) {
    register VD xi = VDABS(_mm512_load_pd(x + i)); VDP(xi);
    VDSORT(xi); VDP(xi);
    const Sleef_quadx8 qi = Sleef_cast_from_doubleq8_avx512f(xi);
    rv = Sleef_fmaq8_u05avx512f(qi, qi, rv);
  }

  Sleef_storeq8_avx512f(rq, rv);

  *rq += rq[1u];
  *rq += rq[2u];
  *rq += rq[3u];
  *rq += rq[4u];
  *rq += rq[5u];
  *rq += rq[6u];
  *rq += rq[7u];

  pquad2ef(rq, e1, f1);
#if (defined(__ICC) || defined(__INTEL_COMPILER) || defined(__INTEL_CLANG_COMPILER) || defined(__INTEL_LLVM_COMPILER))
  *rq = __sqrtq(*rq);
#else /* !Intel */
  *rq = sqrtq(*rq);
#endif /* ?Intel */
  pquad2ef(rq, e0, f0);
  return (double)*rq;
}
