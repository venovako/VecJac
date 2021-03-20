#include "proba.h"

static double dnorme(double ef[static restrict 2], const double x[static restrict 1], const size_t m)
{
  ASSERT_VFPENV;

  ef[0u] = -INFINITY;
  ef[1u] = 1.0;

  if (m & (size_t)VDL_1)
    return -0.0;
  if (!m)
    return 0.0;

  register VD re = _mm512_set1_pd(ef[0u]);
  register VD rf = _mm512_set1_pd(ef[1u]);
  register const VD _inf = re;

  for (size_t i = 0u; i < m; i += VDL) {
    register const VD xi = _mm512_load_pd(x + i);
#ifndef NDEBUG
    (void)VDprintf(stderr, "xi", xi);
#endif /* !NDEBUG */
    register const VD fi = _mm512_getmant_pd(xi, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero);
#ifndef NDEBUG
    (void)VDprintf(stderr, "fi", fi);
#endif /* !NDEBUG */
    register const VD ei = _mm512_getexp_pd(xi);
#ifndef NDEBUG
    (void)VDprintf(stderr, "ei", ei);
#endif /* !NDEBUG */

    register const VD reh = VDLSB(re);
#ifndef NDEBUG
    (void)VDprintf(stderr, "reh", reh);
#endif /* !NDEBUG */
    register const VD rep = _mm512_sub_pd(re, reh);
#ifndef NDEBUG
    (void)VDprintf(stderr, "rep", rep);
#endif /* !NDEBUG */
    register const VD rfp = _mm512_scalef_pd(rf, reh);
#ifndef NDEBUG
    (void)VDprintf(stderr, "rfp", rfp);
#endif /* !NDEBUG */

    register const VD eh = _mm512_scalef_pd(ei, _mm512_set1_pd(1.0));
#ifndef NDEBUG
    (void)VDprintf(stderr, "eh", eh);
#endif /* !NDEBUG */
    register const VD emax = _mm512_max_pd(eh, rep);
#ifndef NDEBUG
    (void)VDprintf(stderr, "emax", emax);
#endif /* !NDEBUG */
    register const VD ehp = _mm512_max_pd(_mm512_sub_pd(eh, emax), _inf);
#ifndef NDEBUG
    (void)VDprintf(stderr, "ehp", ehp);
#endif /* !NDEBUG */
    register const VD repp = _mm512_max_pd(_mm512_sub_pd(rep, emax), _inf);
#ifndef NDEBUG
    (void)VDprintf(stderr, "repp", repp);
#endif /* !NDEBUG */
    register const VD ehpp = _mm512_scalef_pd(ehp, _mm512_set1_pd(-1.0));
#ifndef NDEBUG
    (void)VDprintf(stderr, "ehpp", ehpp);
#endif /* !NDEBUG */

    register const VD fh = _mm512_scalef_pd(fi, ehpp);
#ifndef NDEBUG
    (void)VDprintf(stderr, "fh", fh);
#endif /* !NDEBUG */
    register const VD rfhp = _mm512_scalef_pd(rfp, repp);
#ifndef NDEBUG
    (void)VDprintf(stderr, "rfhp", rfhp);
#endif /* !NDEBUG */

    rf = _mm512_fmadd_pd(fh, fh, rfhp);
#ifndef NDEBUG
    (void)VDprintf(stderr, "rf\'", rf);
#endif /* !NDEBUG */
    re = _mm512_add_pd(emax, _mm512_getexp_pd(rf));
#ifndef NDEBUG
    (void)VDprintf(stderr, "re", re);
#endif /* !NDEBUG */
    rf = _mm512_getmant_pd(rf, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero);
#ifndef NDEBUG
    (void)VDprintf(stderr, "rf", rf);
#endif /* !NDEBUG */
  }
  return scalb(ef[1u], ef[0u]);
}

int main(int argc, char *argv[])
{
  alignas(VA) double r[3u];
  alignas(VA) double x[2u * VDL] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0 };
  r[2u] = dnorme(r, x, 2u * VDL);
  (void)printf("%#.17e\n", r[2u]);
  return EXIT_SUCCESS;
}
