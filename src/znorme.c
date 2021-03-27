#include "znorme.h"

#include "dkvsrt.h"
#include "defops.h"

double znorme_(const fnat m[static restrict 1], const double zr[static restrict VDL], const double zi[static restrict VDL], double e[static restrict 2], double f[static restrict 2])
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

  register const VD _inf = _mm512_set1_pd(-HUGE_VAL);
  register const VD _one = _mm512_set1_pd(-1.0);
  register const VD  one = _mm512_set1_pd( 1.0);
#ifndef NDEBUG
  register const VD  inf = _mm512_set1_pd( HUGE_VAL);
#endif /* !NDEBUG */

#ifdef __AVX512DQ__
  register const __m512i ione = _mm512_set1_epi64(1);
#else /* !__AVX512DQ__ */
  register const __m256i ione = _mm256_set1_epi32(1);
#endif /* ?__AVX512DQ__ */

  register VD re = _inf;
  register VD rf =  one;

  for (fnat i = 0u; i < *m; i += VDL) {
    register const VD xi = _mm512_load_pd(zr + i); VDP(xi);
    register const VD fi = VDMANT(xi); VDP(fi);
    register const VD ei = _mm512_getexp_pd(xi); VDP(ei);
#ifndef NDEBUG
    if (MD2U(_mm512_cmplt_pd_mask(ei, inf)) != 0xFFu)
      return -4.0;
#endif /* !NDEBUG */

    register const VD reh = VDLSB(re); VDP(reh);
    register const VD rep = _mm512_sub_pd(re, reh); VDP(rep);
    register const VD rfp = _mm512_scalef_pd(rf, reh); VDP(rfp);

    register const VD eh = _mm512_scalef_pd(ei, one); VDP(eh);
    register const VD emax = _mm512_max_pd(eh, rep); VDP(emax);
    register const VD ehp = VDSUBE(eh, emax); VDP(ehp);
    register const VD repp = VDSUBE(rep, emax); VDP(repp);
    register const VD ehpp = _mm512_scalef_pd(ehp, _one); VDP (ehpp);

    register const VD fh = _mm512_scalef_pd(fi, ehpp); VDP(fh);
    register const VD rfhp = _mm512_scalef_pd(rfp, repp); VDP(rfhp);

    rf = _mm512_fmadd_pd(fh, fh, rfhp); VDP(rf);
    re = _mm512_add_pd(emax, _mm512_getexp_pd(rf)); VDP(re);
    rf = VDMANT(rf); VDP(rf);
  }

  for (fnat i = 0u; i < *m; i += VDL) {
    register const VD xi = _mm512_load_pd(zi + i); VDP(xi);
    register const VD fi = VDMANT(xi); VDP(fi);
    register const VD ei = _mm512_getexp_pd(xi); VDP(ei);
#ifndef NDEBUG
    if (MD2U(_mm512_cmplt_pd_mask(ei, inf)) != 0xFFu)
      return -4.0;
#endif /* !NDEBUG */

    register const VD reh = VDLSB(re); VDP(reh);
    register const VD rep = _mm512_sub_pd(re, reh); VDP(rep);
    register const VD rfp = _mm512_scalef_pd(rf, reh); VDP(rfp);

    register const VD eh = _mm512_scalef_pd(ei, one); VDP(eh);
    register const VD emax = _mm512_max_pd(eh, rep); VDP(emax);
    register const VD ehp = VDSUBE(eh, emax); VDP(ehp);
    register const VD repp = VDSUBE(rep, emax); VDP(repp);
    register const VD ehpp = _mm512_scalef_pd(ehp, _one); VDP (ehpp);

    register const VD fh = _mm512_scalef_pd(fi, ehpp); VDP(fh);
    register const VD rfhp = _mm512_scalef_pd(rfp, repp); VDP(rfhp);

    rf = _mm512_fmadd_pd(fh, fh, rfhp); VDP(rf);
    re = _mm512_add_pd(emax, _mm512_getexp_pd(rf)); VDP(re);
    rf = VDMANT(rf); VDP(rf);
  }

  VDKVSORT(re, rf);
  VDEFRED(re, rf);

  _mm512_mask_storeu_pd((e + 1u), 0x01u, re);
  _mm512_mask_storeu_pd((f + 1u), 0x01u, rf);

  if ((long)(e[1u]) & 1l) {
    e[0u] = scalbn((e[1u] - 1.0), -1);
    f[0u] = sqrt(scalbn(f[1u], 1));
  }
  else {
    e[0u] = scalbn(e[1u], -1);
    f[0u] = sqrt(f[1u]);
  }
  return scalb(f[0u], e[0u]);
}
