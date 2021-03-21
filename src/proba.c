#include "proba.h"

static double dnorme_(const fnat m[static restrict 1], const double x[static restrict 1], double e[static restrict 1], double f[static restrict 1])
{
  ASSERT_VFPENV;

  *e = -HUGE_VAL;
  *f =  1.0;

  if (*m & VDL_1)
    return -1.0;
  if (!*m)
    return  0.0;

  register const VD _inf = _mm512_set1_pd(-HUGE_VAL);
  register const VD _one = _mm512_set1_pd(-1.0);
  register const VD  one = _mm512_set1_pd( 1.0);

  register VD re = _inf;
  register VD rf =  one;

  for (fnat i = 0u; i < *m; i += VDL) {
    register const VD xi = _mm512_load_pd(x + i); VDP(xi);
    register const VD fi = _mm512_getmant_pd(xi, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero); VDP(fi);
    register const VD ei = _mm512_getexp_pd(xi); VDP(ei);

    register const VD reh = VDLSB(re); VDP(reh);
    register const VD rep = _mm512_sub_pd(re, reh); VDP(rep);
    register const VD rfp = _mm512_scalef_pd(rf, reh); VDP(rfp);

    register const VD eh = _mm512_scalef_pd(ei, one); VDP(eh);
    register const VD emax = _mm512_max_pd(eh, rep); VDP(emax);
    register const VD ehp = _mm512_max_pd(_mm512_sub_pd(eh, emax), _inf); VDP(ehp);
    register const VD repp = _mm512_max_pd(_mm512_sub_pd(rep, emax), _inf); VDP(repp);
    register const VD ehpp = _mm512_scalef_pd(ehp, _one); VDP (ehpp);

    register const VD fh = _mm512_scalef_pd(fi, ehpp); VDP(fh);
    register const VD rfhp = _mm512_scalef_pd(rfp, repp); VDP(rfhp);

    rf = _mm512_fmadd_pd(fh, fh, rfhp); VDP(rf);
    re = _mm512_add_pd(emax, _mm512_getexp_pd(rf)); VDP(re);
    rf = _mm512_getmant_pd(rf, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero); VDP(rf);
  }

  return scalbln(*f, (long)*e);
}

int main(int argc, char *argv[])
{
  static const fnat m = (2u * VDL);
  alignas(VA) double x[m] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0 };
  double e = 0.0, f = 0.0, d = dnorme_(&m, x, &e, &f);
  (void)printf("%#.17e\n", d);
  return EXIT_SUCCESS;
}
