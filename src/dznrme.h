#ifndef DZNRME_H
#define DZNRME_H

#ifdef DZNRME_ITSORT
#include "vecdef.h"
#include "d8sort.h"
#endif /* DZNRME_ITSORT */
#include "defops.h"

#ifdef DZNRME_VARS0
#error DZNRME_VARS0 already defined
#else /* !DZNRME_VARS0 */
#define DZNRME_VARS0                                  \
  register const VD _inf = _mm512_set1_pd(-HUGE_VAL); \
  register const VD _one = _mm512_set1_pd(-1.0);      \
  register const VD  one = _mm512_set1_pd( 1.0);      \
  register VD re = _inf;                              \
  register VD rf =  one
#endif /* ?DZNRME_VARS0 */

#ifdef DZNRME_VARS1
#error DZNRME_VARS1 already defined
#else /* !DZNRME_VARS1 */
#ifdef __AVX512DQ__
#define DZNRME_VARS1                                  \
  register const __m512i ione = _mm512_set1_epi64(1); \
  DZNRME_VARS0
#else /* !__AVX512DQ__ */
#define DZNRME_VARS1                                  \
  register const __m256i ione = _mm256_set1_epi32(1); \
  DZNRME_VARS0
#endif /* ?__AVX512DQ__ */
#endif /* ?DZNRME_VARS1 */

#ifdef DZNRME_VARS2
#error DZNRME_VARS2 already defined
#else /* !DZNRME_VARS2 */
#ifdef DZNRME_ITSORT
#define DZNRME_VARS2                           \
  register const VD m0 = _mm512_set1_pd(-0.0); \
  DZNRME_VARS1
#else /* !DZNRME_ITSORT */
#define DZNRME_VARS2 DZNRME_VARS1
#endif /* ?DZNRME_ITSORT */
#endif /* ?DZNRME_VARS2 */

#ifdef DZNRME_VARS
#error DZNRME_VARS already defined
#else /* !DZNRME_VARS */
#ifdef NDEBUG
#define DZNRME_VARS DZNRME_VARS2
#else /* !NDEBUG */
#define DZNRME_VARS                                 \
  register const VD inf = _mm512_set1_pd(HUGE_VAL); \
  DZNRME_VARS2
#endif /* ?NDEBUG */
#endif /* ?DZNRME_VARS */

#ifdef ASSERT_FINITE
#error ASSERT_FINITE already defined
#else /* !ASSERT_FINITE */
#ifdef NDEBUG
#define ASSERT_FINITE(ec) (void)(ec)
#else /* !NDEBUG */
#define ASSERT_FINITE(ec)                           \
  if (MD2U(_mm512_cmplt_pd_mask(ei, inf)) != 0xFFu) \
    return (ec)
#endif /* ?NDEBUG */
#endif /* ?ASSERT_FINITE */

#ifdef DZNRME_XISORT
#error DZNMRE_XISORT already defined
#else /* !DZNRME_XISORT */
#ifdef DZNRME_ITSORT
#define DZNRME_XISORT      \
  xi = VDABS(xi); VDP(xi); \
  VDSORT(xi); VDP(xi)
#else /* !DZNRME_ITSORT */
#define DZNRME_XISORT
#endif /* ?DZNRME_ITSORT */
#endif /* ?DZNRME_XISORT */

#ifdef DZNRME_LOOP
#error DZNRME_LOOP already defined
#else /* !DZNRME_LOOP */
#define DZNRME_LOOP(x,ec)                                             \
  for (fnat i = 0u; i < *m; i += VDL) {                               \
    register VD xi = _mm512_load_pd((x) + i); VDP(xi);                \
    register const VD reh = VDLSB(re); VDP(reh);                      \
    register const VD rep = _mm512_sub_pd(re, reh); VDP(rep);         \
    register const VD rfp = _mm512_scalef_pd(rf, reh); VDP(rfp);      \
    DZNRME_XISORT;                                                    \
    register const VD fi = VDMANT(xi); VDP(fi);                       \
    register const VD ei = _mm512_getexp_pd(xi); VDP(ei);             \
    ASSERT_FINITE(ec);                                                \
    register const VD eh = _mm512_scalef_pd(ei, one); VDP(eh);        \
    register const VD emax = _mm512_max_pd(eh, rep); VDP(emax);       \
    register const VD ehp = VDSUBE(eh, emax); VDP(ehp);               \
    register const VD repp = VDSUBE(rep, emax); VDP(repp);            \
    register const VD ehpp = _mm512_scalef_pd(ehp, _one); VDP (ehpp); \
    register const VD fh = _mm512_scalef_pd(fi, ehpp); VDP(fh);       \
    register const VD rfhp = _mm512_scalef_pd(rfp, repp); VDP(rfhp);  \
    rf = _mm512_fmadd_pd(fh, fh, rfhp); VDP(rf);                      \
    re = _mm512_add_pd(emax, _mm512_getexp_pd(rf)); VDP(re);          \
    rf = VDMANT(rf); VDP(rf);                                         \
  }
#endif /* ?DZNRME_LOOP */

#ifdef DZNRME_RET
#error DZNRME_RET already defined
#else /* !DZNRME_RET */
#define DZNRME_RET                      \
  _mm512_mask_storeu_pd(e1, 0x01u, re); \
  _mm512_mask_storeu_pd(f1, 0x01u, rf); \
  if ((long)*e1 & 1l) {                 \
    *e0 = scalbn((*e1 - 1.0), -1);      \
    *f0 = sqrt(scalbn(*f1, 1));         \
  }                                     \
  else {                                \
    *e0 = scalbn(*e1, -1);              \
    *f0 = sqrt(*f1);                    \
  }                                     \
  return scalb(*f0, *e0)
#endif /* ?DZNRME_RET */

#endif /* !DZNRME_H */
