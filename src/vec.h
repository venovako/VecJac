#ifndef VEC_H
#define VEC_H

#include "common.h"

/*** redefine if changing the vector type ***/

#ifdef __AVX512F__
#include <immintrin.h>
#else /* !__AVX512F__ */
#error AVX-512 instructions not available
#endif /* ?__AVX512F__ */

/* vector types */

/* vector type containing floats */
#ifdef VS
#error VS already defined
#else /* !VS */
#define VS __m512
#endif /* ?VS */

/* vector type containing doubles */
#ifdef VD
#error VD already defined
#else /* !VD */
#define VD __m512d
#endif /* ?VD */

/* mask types */

/* mask type for float lanes */
#ifdef MS
#error MS already defined
#else /* !MS */
#define MS __mmask16
#endif /* ?MS */

/* mask type for double lanes */
#ifdef MD
#error MD already defined
#else /* !MD */
#define MD __mmask8
#endif /* ?MD */

/* various vectorization constants */

/* vector length in 32-bit lanes */
#ifdef VSL
#error VSL already defined
#else /* !VSL */
#define VSL 16u
#endif /* ?VSL */

/* VSL - 1 */
#ifdef VSL_1
#error VSL_1 already defined
#else /* !VSL_1 */
#define VSL_1 15u
#endif /* ?VSL_1 */

/* lg(VSL) */
#ifdef VSLlg
#error VSLlg already defined
#else /* !VSLlg */
#define VSLlg 4u
#endif /* ?VSLlg */

/* vector length in 64-bit lanes */
#ifdef VDL
#error VDL already defined
#else /* !VDL */
#define VDL 8u
#endif /* ?VDL */

/* VDL - 1 */
#ifdef VDL_1
#error VDL_1 already defined
#else /* !VDL_1 */
#define VDL_1 7u
#endif /* ?VDL_1 */

/* lg(VDL) */
#ifdef VDLlg
#error VDLlg already defined
#else /* !VDLlg */
#define VDLlg 3u
#endif /* ?VDLlg */

/* vector alignment in bytes */
#ifdef VA
#error VA already defined
#else /* !VA */
#define VA 64u
#endif /* ?VA */

/* vector instructions */

/* VS instruction name */
#ifdef VSI
#error VSI already defined
#else /* !VSI */
#define VSI(x) _mm512_##x##_ps
#endif /* ?VSI */

/* VD instruction name */
#ifdef VDI
#error VDI already defined
#else /* !VDI */
#define VDI(x) _mm512_##x##_pd
#endif /* ?VDI */

/* float <-comparison */
#ifdef VSLT
#error VSLT already defined
#else /* !VSLT */
#define VSLT(x,y) _mm512_cmplt_ps_mask((x),(y))
#endif /* ?VSLT */

/* double <-comparison */
#ifdef VDLT
#error VDLT already defined
#else /* !VDLT */
#define VDLT(x,y) _mm512_cmplt_pd_mask((x),(y))
#endif /* ?VDLT */

/* bitwise operations */

#ifdef VSAND
#error VSAND already defined
#endif /* VSAND */
#ifdef VDAND
#error VDAND already defined
#endif /* VDAND */

#ifdef VSANDNOT
#error VSANDNOT already defined
#endif /* VSANDNOT */
#ifdef VDANDNOT
#error VDANDNOT already defined
#endif /* VDANDNOT */

#ifdef VSOR
#error VSOR already defined
#endif /* VSOR */
#ifdef VDOR
#error VDOR already defined
#endif /* VDOR */

#ifdef VSXOR
#error VSXOR already defined
#endif /* VSXOR */
#ifdef VDXOR
#error VDXOR already defined
#endif /* VDXOR */

#ifdef VSLSB
#error VSLSB already defined
#else /* !VSLSB */
#define VSLSB(x) VSI(cvtepi32)(_mm256_and_si256(_mm512_cvtps_epi32(x),_mm256_set1_epi32(1)))
#endif /* ?VSLSB */
#ifdef VDLSB
#error VDLSB already defined
#endif /* VDLSB */

#ifdef __AVX512DQ__
#define VSAND(x,y) VSI(and)((x),(y))
#define VSANDNOT(x,y) VSI(andnot)((x),(y))
#define VSOR(x,y) VSI(or)((x),(y))
#define VSXOR(x,y) VSI(xor)((x),(y))

#define VDAND(x,y) VDI(and)((x),(y))
#define VDANDNOT(x,y) VDI(andnot)((x),(y))
#define VDOR(x,y) VDI(or)((x),(y))
#define VDXOR(x,y) VDI(xor)((x),(y))
#define VDLSB(x) VDI(cvtepi64)(_mm512_and_epi64(_mm512_cvtpd_epi64(x),_mm512_set1_epi64(INT64_C(1))))
#else /* AVX512F only */
#define VSAND(x,y) VSI(castsi512)(_mm512_and_epi32(_mm512_castps_si512(x),_mm512_castps_si512(y)))
#define VSANDNOT(x,y) VSI(castsi512)(_mm512_andnot_epi32(_mm512_castps_si512(x),_mm512_castps_si512(y)))
#define VSOR(x,y) VSI(castsi512)(_mm512_or_epi32(_mm512_castps_si512(x),_mm512_castps_si512(y)))
#define VSXOR(x,y) VSI(castsi512)(_mm512_xor_epi32(_mm512_castps_si512(x),_mm512_castps_si512(y)))

#define VDAND(x,y) VDI(castsi512)(_mm512_and_epi64(_mm512_castpd_si512(x),_mm512_castpd_si512(y)))
#define VDANDNOT(x,y) VDI(castsi512)(_mm512_andnot_epi64(_mm512_castpd_si512(x),_mm512_castpd_si512(y)))
#define VDOR(x,y) VDI(castsi512)(_mm512_or_epi64(_mm512_castpd_si512(x),_mm512_castpd_si512(y)))
#define VDXOR(x,y) VDI(castsi512)(_mm512_xor_epi64(_mm512_castpd_si512(x),_mm512_castpd_si512(y)))
#define VDLSB(x) VDI(cvtepi32)(_mm256_and_si256(_mm512_cvtpd_epi32(x),_mm256_set1_epi32(1)))
#endif /* ?__AVX512DQ__ */

/* mask conversion */

#ifdef MS2U
#error MS2U already defined
#else /* !MS2U */
#define MS2U(m) _cvtmask16_u32(m)
#endif /* ?MS2U */

#ifdef MD2U
#error MD2U already defined
#else /* !MD2U */
#ifdef __AVX512DQ__
#define MD2U(m) _cvtmask8_u32(m)
#else /* AVX512F only */
#define MD2U(m) _cvtmask16_u32((__mmask16)(m))
#endif /* ?__AVX512DQ__ */
#endif /* ?MD2U */

/*** end of vector definitions ***/

static inline size_t n2VS(const size_t n)
{
  return ((n + VSL_1) >> VSLlg);
}

static inline size_t n2VD(const size_t n)
{
  return ((n + VDL_1) >> VDLlg);
}

static inline size_t n2NS(const size_t n)
{
  return (n2VS(n) << VSLlg);
}

static inline size_t n2ND(const size_t n)
{
  return (n2VD(n) << VDLlg);
}

extern int VSprintf(FILE f[static 1], const char *const h, const VS v);
extern int VDprintf(FILE f[static 1], const char *const h, const VD v);

extern int MSprintf(FILE f[static 1], const char *const h, const MS m);
extern int MDprintf(FILE f[static 1], const char *const h, const MD m);

extern size_t VSread(float s[static VSL], const size_t V, FILE f[static 1]);
extern size_t VDread(double d[static VDL], const size_t V, FILE f[static 1]);

extern size_t VSwrite(const float s[static VSL], const size_t V, FILE f[static 1]);
extern size_t VDwrite(const double d[static VDL], const size_t V, FILE f[static 1]);

#endif /* !VEC_H */
