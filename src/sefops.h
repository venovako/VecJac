#ifndef SEFOPS_H
#define SEFOPS_H

#ifdef VSMANT
#error VSMANT already defined
#else /* !VSMANT */
#define VSMANT(x) _mm512_getmant_ps((x),_MM_MANT_NORM_1_2,_MM_MANT_SIGN_zero)
#endif /* ?VSMANT */

// _inff has to be defined to -\inftys
#ifdef VSSUBE
#error VSSUBE already defined
#else /* !VSSUBE */
#define VSSUBE(x,y) _mm512_max_ps(_mm512_sub_ps((x),(y)),_inff)
#endif /* ?VSSUBE */

// iones has to be defined to integer ones
#ifdef VSLSB
#error VSLSB already defined
#else /* !VSLSB */
#define VSLSB(x) _mm512_cvtepi32_ps(_mm512_and_si256(_mm512_cvtps_epi32(x),iones))
#endif /* ?VSLSB */

#ifdef MSOR
#error MSOR already defined
#else /* !MSOR */
#define MSOR(a,b) _kor_mask16((a),(b))
#endif /* ?MSOR */

#ifdef VSEFLE
#error VSEFLE already defined
#else /* !VSEFLE */
#define VSEFLE(e0,e1,f0,f1) MSOR(_mm512_cmplt_ps_mask(e0,e1),_mm512_mask_cmple_ps_mask(_mm512_cmpeq_ps_mask(e0,e1),f0,f1))
#endif /* ?VSEFLE */

#ifdef FLT_MAX_FIN_EXP
#error FLT_MAX_FIN_EXP already defined
#else /* !FLT_MAX_FIN_EXP */
#define FLT_MAX_FIN_EXP 127.0f
#endif /* ?FLT_MAX_FIN_EXP */

#ifdef SCNRME_SEQRED
static inline void efswpf(float e0[static restrict 1], float f0[static restrict 1], float e1[static restrict 1], float f1[static restrict 1])
{
  if ((*e0 > *e1) || ((*e0 == *e1) && (*f0 > *f1))) {
    float t = *e1;
    *e1 = *e0;
    *e0 = t;
    t = *f1;
    *f1 = *f0;
    *f0 = t;
  }
}

static inline void efsumf(float e0[static restrict 1], float f0[static restrict 1], float e1[static restrict 1], float f1[static restrict 1], int ef[static restrict 1])
{
  efswpf(e0, f0, e1, f1);
  *f1 = fmaf(scalbf(1.0f, fmaxf((*e0 - *e1), -HUGE_VALF)), *f0, *f1);
  *f1 = scalbnf(frexpf(*f1, ef), 1);
  *e1 += --*ef;
}

static inline void efredf(float e[static restrict VSL], float f[static restrict VSL])
{
  int ef
#ifndef NDEBUG
    = 0
#endif /* !NDEBUG */
    ;
  for (unsigned i = 0u; i < VSL_1; ) {
    float *const e0 = (e + i);
    float *const f0 = (f + i);
    ++i;
    float *const e1 = (e + i);
    float *const f1 = (f + i);
    efsumf(e0, f0, e1, f1, &ef);
  }
}
#else /* !SCNRME_SEQRED */
#ifdef VSEFADD
#error VSEFADD already defined
#else /* !VSEFADD */
#define VSEFADD(e,f,ma,mb)                                                       \
  {                                                                              \
    register const VS ae = _mm512_maskz_compress_ps(ma, e); VSP(ae);             \
    register const VS be = _mm512_maskz_compress_ps(mb, e); VSP(be);             \
    register const VS af = _mm512_maskz_compress_ps(ma, f); VSP(af);             \
    register const VS bf = _mm512_maskz_compress_ps(mb, f); VSP(bf);             \
    f = _mm512_fmadd_ps(_mm512_scalef_ps(onef, VSSUBE(ae, be)), af, bf); VSP(f); \
    e = _mm512_add_ps(be, _mm512_getexp_ps(f)); VSP(e);                          \
    f = VSMANT(f); VSP(f);                                                       \
  }
#endif /* ?VSEFADD */

#ifdef VSEFRED
#error VSEFRED already defined
#else /* !VSEFRED */
#define VSEFRED(e,f)              \
  {                               \
    VSEFADD(e,f,0x5555u,0xAAAAu); \
    VSEFADD(e,f,0x0055u,0x00AAu); \
    VSEFADD(e,f,0x0005u,0x000Au); \
    VSEFADD(e,f,0x0001u,0x0002u); \
  }
#endif /* ?VSEFRED */
#endif /* ?SCNRME_SEQRED */

#endif /* !SEFOPS_H */
