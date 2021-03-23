#include "dnorme.h"

#ifdef MDOR
#error MDOR already defined
#else /* !MDOR */
#ifdef __AVX512DQ__
#define MDOR(a,b) _kor_mask8((a),(b))
#else /* !__AVX512DQ__ */
#define MDOR(a,b) (__mmask8)_kor_mask16((a),(b))
#endif /* ?__AVX512DQ__ */
#endif /* ?MDOR */

// VDKVBITONIC and VDKVSORT have been inspired by the code of Berenger Bramas.
// Please, see https://gitlab.inria.fr/bramas/avx-512-sort and its LICENSE.
#ifdef VDKVBITONIC
#error VDKVBITONIC already defined
#else /* !VDKVBITONIC */
#define VDKVBITONIC(e,f,i,m)                                                                       \
  {                                                                                                \
    register const VD pe = _mm512_permutexvar_pd(i, e); VDP(pe);                                   \
    register const VD pf = _mm512_permutexvar_pd(i, f); VDP(pf);                                   \
    register const MD mf = _mm512_mask_cmple_pd_mask(_mm512_cmpeq_pd_mask(e, pe), f, pf); MDP(mf); \
    register const MD me = _mm512_cmplt_pd_mask(e, pe); MDP(me);                                   \
    register const MD mm = MDOR(me, mf); MDP(mm);                                                  \
    register const VD em = _mm512_mask_blend_pd(mm, pe, e); VDP(em);                               \
    register const VD ex = _mm512_mask_blend_pd(mm, e, pe); VDP(ex);                               \
    register const VD fm = _mm512_mask_blend_pd(mm, pf, f); VDP(fm);                               \
    register const VD fx = _mm512_mask_blend_pd(mm, f, pf); VDP(fx);                               \
    e = _mm512_mask_mov_pd(em, m, ex); VDP(e);                                                     \
    f = _mm512_mask_mov_pd(fm, m, fx); VDP(f);                                                     \
  }
#endif /* ?VDKVBITONIC */
#ifdef VDKVSORT
#error VDKVSORT already defined
#else /* !VDKVSORT */
#define VDKVSORT(e,f)                                         \
  {                                                           \
    register VI i = _mm512_set_epi64(6, 7, 4, 5, 2, 3, 0, 1); \
    VDKVBITONIC(e,f,i,0xAAu);                                 \
    i = _mm512_set_epi64(4, 5, 6, 7, 0, 1, 2, 3);             \
    VDKVBITONIC(e,f,i,0xCCu);                                 \
    i = _mm512_set_epi64(6, 7, 4, 5, 2, 3, 0, 1);             \
    VDKVBITONIC(e,f,i,0xAAu);                                 \
    i = _mm512_set_epi64(0, 1, 2, 3, 4, 5, 6, 7);             \
    VDKVBITONIC(e,f,i,0xF0u);                                 \
    i = _mm512_set_epi64(5, 4, 7, 6, 1, 0, 3, 2);             \
    VDKVBITONIC(e,f,i,0xCCu);                                 \
    i = _mm512_set_epi64(6, 7, 4, 5, 2, 3, 0, 1);             \
    VDKVBITONIC(e,f,i,0xAAu);                                 \
  }
#endif /* ?VDKVSORT */

#ifdef VDMANT
#error VDMANT already defined
#else /* !VDMANT */
#define VDMANT(x) _mm512_getmant_pd((x),_MM_MANT_NORM_1_2,_MM_MANT_SIGN_zero)
#endif /* ?VDMANT */

#ifdef VDSUBE
#error VDSUBE already defined
#else /* !VDSUBE */
#define VDSUBE(x,y) _mm512_max_pd(_mm512_sub_pd((x),(y)),_inf)
#endif /* ?VDSUBE */

#ifdef VDLSB
#error VDLSB already defined
#else /* !VDLSB */
#ifdef __AVX512DQ__
#define VDLSB(x) _mm512_cvtepi64_pd(_mm512_and_epi64(_mm512_cvtpd_epi64(x),ione))
#else /* !__AVX512DQ__ */
#define VDLSB(x) _mm512_cvtepi32_pd(_mm256_and_si256(_mm512_cvtpd_epi32(x),ione))
#endif /* ?__AVX512DQ__ */
#endif /* ?VDLSB */

#ifdef VDEFADD
#error VDEFADD already defined
#else /* !VDEFADD */
#define VDEFADD(e,f,ma,mb)                                                      \
  {                                                                             \
    register const VD ae = _mm512_maskz_compress_pd(ma, e); VDP(ae);            \
    register const VD be = _mm512_maskz_compress_pd(mb, e); VDP(be);            \
    register const VD af = _mm512_maskz_compress_pd(ma, f); VDP(af);            \
    register const VD bf = _mm512_maskz_compress_pd(mb, f); VDP(bf);            \
    f = _mm512_fmadd_pd(_mm512_scalef_pd(one, VDSUBE(ae, be)), af, bf); VDP(f); \
    e = _mm512_add_pd(be, _mm512_getexp_pd(f)); VDP(e);                         \
    f = VDMANT(f); VDP(f);                                                      \
  }
#endif /* ?VDEFADD */
#ifdef VDEFRED
#error VDEFRED already defined
#else /* !VDEFRED */
#define VDEFRED(e,f)          \
  {                           \
    VDEFADD(e,f,0x55u,0xAAu); \
    VDEFADD(e,f,0x05u,0x0Au); \
    VDEFADD(e,f,0x01u,0x02u); \
  }
#endif /* ?VDEFRED */

double dnorme_(const fnat m[static restrict 1], const double x[static restrict 1], double e[static restrict 1], double f[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -4.0;
  if (*m & VDL_1)
    return -1.0;
  if (IS_NOT_ALIGNED(x))
    return -2.0;
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
    register const VD xi = _mm512_load_pd(x + i); VDP(xi);
    register const VD fi = VDMANT(xi); VDP(fi);
    register const VD ei = _mm512_getexp_pd(xi); VDP(ei);
#ifndef NDEBUG
    if (MD2U(_mm512_cmplt_pd_mask(ei, inf)) != 0xFFu)
      return -3.0;
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

  _mm512_mask_storeu_pd(e, 0x01u, re);
  _mm512_mask_storeu_pd(f, 0x01u, rf);
  return scalbln(*f, (long)*e);
}
