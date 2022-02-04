#ifndef VEC_H
#define VEC_H

#include "common.h"

/*** redefine if changing the vector type ***/

#ifdef __AVX512F__
#include <immintrin.h>
#else /* !__AVX512F__ */
#error AVX-512 instructions not available
#endif /* ?__AVX512F__ */

/* default MXCSR */

#ifdef STD_MXCSR
#error STD_MXCSR already defined
#else /* !STD_MXCSR */
#define STD_MXCSR 0x1F80u
#endif /* ?STD_MXCSR */

#ifdef IS_NOT_VFPENV
#error IS_NOT_VFPENV already defined
#else /* !IS_NOT_VFPENV */
#define IS_NOT_VFPENV ((_mm_getcsr() & 0xFFC0u) != (STD_MXCSR))
#endif /* ?IS_NOT_VFPENV */

/* default alignment */

/* vector alignment in bytes */
#ifdef VA
#error VA already defined
#else /* !VA */
#define VA 64u
#endif /* ?VA */

#ifdef IS_NOT_ALIGNED
#error IS_NOT_ALIGNED already defined
#else /* !IS_NOT_ALIGNED */
#define IS_NOT_ALIGNED(x) ((uintptr_t)(x) & 0x3Fu)
#endif /* ?IS_NOT_ALIGNED */

/* vector types */

/* vector type containing integers */

#ifdef VI
#error VI already defined
#else /* !VI */
#define VI __m512i
#endif /* ?VI */

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

/* VSL / 2 */
#ifdef VSL_2
#error VSL_2 already defined
#else /* !VSL_2 */
#define VSL_2 8u
#endif /* ?VSL_2 */

/* (VSL - 1) / 2 */
#ifdef VSL__2
#error VSL__2 already defined
#else /* !VSL__2 */
#define VSL__2 7u
#endif /* ?VSL__2 */

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

#ifdef VDL2
#error VDL2 already defined
#else /* !VDL2 */
#define VDL2 16u
#endif /* ?VDL2 */

/* VDL - 1 */
#ifdef VDL_1
#error VDL_1 already defined
#else /* !VDL_1 */
#define VDL_1 7u
#endif /* ?VDL_1 */

/* VDL / 2 */
#ifdef VDL_2
#error VDL_2 already defined
#else /* !VDL_2 */
#define VDL_2 4u
#endif /* ?VDL_2 */

/* (VDL - 1) / 2 */
#ifdef VDL__2
#error VDL__2 already defined
#else /* !VDL__2 */
#define VDL__2 3u
#endif /* ?VDL__2 */

/* lg(VDL) */
#ifdef VDLlg
#error VDLlg already defined
#else /* !VDLlg */
#define VDLlg 3u
#endif /* ?VDLlg */

/* vector instructions */

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

#ifdef __AVX512DQ__
#define VSAND(x,y) _mm512_and_ps((x),(y))
#define VSANDNOT(x,y) _mm512_andnot_ps((x),(y))
#define VSOR(x,y) _mm512_or_ps((x),(y))
#define VSXOR(x,y) _mm512_xor_ps((x),(y))

#define VDAND(x,y) _mm512_and_pd((x),(y))
#define VDANDNOT(x,y) _mm512_andnot_pd((x),(y))
#define VDOR(x,y) _mm512_or_pd((x),(y))
#define VDXOR(x,y) _mm512_xor_pd((x),(y))
#else /* !__AVX512DQ__ */
#define VSAND(x,y) _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(x),_mm512_castps_si512(y)))
#define VSANDNOT(x,y) _mm512_castsi512_ps(_mm512_andnot_epi32(_mm512_castps_si512(x),_mm512_castps_si512(y)))
#define VSOR(x,y) _mm512_castsi512_ps(_mm512_or_epi32(_mm512_castps_si512(x),_mm512_castps_si512(y)))
#define VSXOR(x,y) _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(x),_mm512_castps_si512(y)))

#define VDAND(x,y) _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(x),_mm512_castpd_si512(y)))
#define VDANDNOT(x,y) _mm512_castsi512_pd(_mm512_andnot_epi64(_mm512_castpd_si512(x),_mm512_castpd_si512(y)))
#define VDOR(x,y) _mm512_castsi512_pd(_mm512_or_epi64(_mm512_castpd_si512(x),_mm512_castpd_si512(y)))
#define VDXOR(x,y) _mm512_castsi512_pd(_mm512_xor_epi64(_mm512_castpd_si512(x),_mm512_castpd_si512(y)))
#endif /* ?__AVX512DQ__ */

/* mask operations */

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
#else /* !__AVX512DQ__ */
#define MD2U(m) _cvtmask16_u32(m)
#endif /* ?__AVX512DQ__ */
#endif /* ?MD2U */

/* printout */

#ifdef VSP
#error VSP already defined
#else /* !VSP */
#if (defined(PRINTOUT) && !defined(_OPENMP))
#define VSP(v) (void)VSprintf((PRINTOUT), #v, (v))
#else /* !PRINTOUT || _OPENMP */
#define VSP(v) (void)0
#endif /* ?PRINTOUT */
#endif /* ?VSP */

#ifdef VDP
#error VDP already defined
#else /* !VDP */
#if (defined(PRINTOUT) && !defined(_OPENMP))
#define VDP(v) (void)VDprintf((PRINTOUT), #v, (v))
#else /* !PRINTOUT || _OPENMP */
#define VDP(v) (void)0
#endif /* ?PRINTOUT */
#endif /* ?VDP */

#ifdef MSP
#error MSP already defined
#else /* !MSP */
#if (defined(PRINTOUT) && !defined(_OPENMP))
#define MSP(m) (void)MSprintf((PRINTOUT), #m, (m))
#else /* !PRINTOUT || _OPENMP */
#define MSP(m) (void)0
#endif /* ?PRINTOUT */
#endif /* ?MSP */

#ifdef MDP
#error MDP already defined
#else /* !MDP */
#if (defined(PRINTOUT) && !defined(_OPENMP))
#define MDP(m) (void)MDprintf((PRINTOUT), #m, (m))
#else /* !PRINTOUT || _OPENMP */
#define MDP(m) (void)0
#endif /* ?PRINTOUT */
#endif /* ?MDP */

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
