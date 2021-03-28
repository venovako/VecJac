#ifndef DEFOPS_H
#define DEFOPS_H

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

#ifdef MDOR
#error MDOR already defined
#else /* !MDOR */
#ifdef __AVX512DQ__
#define MDOR(a,b) _kor_mask8((a),(b))
#else /* !__AVX512DQ__ */
#define MDOR(a,b) (__mmask8)_kor_mask16((a),(b))
#endif /* ?__AVX512DQ__ */
#endif /* ?MDOR */

#ifdef VDEFLE
#error VDEFLE already defined
#else /* !VDEFLE */
#define VDEFLE(e0,e1,f0,f1) MDOR(_mm512_cmplt_pd_mask(e0,e1),_mm512_mask_cmple_pd_mask(_mm512_cmpeq_pd_mask(e0,e1),f0,f1))
#endif /* ?VDEFLE */

#endif /* !DEFOPS_H */
