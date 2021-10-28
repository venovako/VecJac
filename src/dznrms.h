#ifndef DZNRMS_H
#define DZNRMS_H

#include "vecdef.h"
#include "d8sort.h"
#include "defops.h"

#ifdef DZNRMS_VARS0
#error DZNRMS_VARS0 already defined
#else /* !DZNRMS_VARS0 */
#define DZNRMS_VARS0                                   \
  register const VD  _inf = _mm512_set1_pd(-HUGE_VAL); \
  register const VD  _one = _mm512_set1_pd(-1.0);      \
  register const VD _zero = _mm512_set1_pd(-0.0);      \
  register const VD   one = _mm512_set1_pd( 1.0);      \
  register VD re = _inf;                               \
  register VD rf =  one
#endif /* ?DZNRMS_VARS0 */

#ifdef DZNRMS_VARS1
#error DZNRMS_VARS1 already defined
#else /* !DZNRMS_VARS1 */
#ifdef __AVX512DQ__
#define DZNRMS_VARS1                                  \
  register const __m512i ione = _mm512_set1_epi64(1); \
  DZNRMS_VARS0
#else /* !__AVX512DQ__ */
#define DZNRMS_VARS1                                  \
  register const __m256i ione = _mm256_set1_epi32(1); \
  DZNRMS_VARS0
#endif /* ?__AVX512DQ__ */
#endif /* ?DZNRMS_VARS1 */

#ifdef DZNRMS_VARS
#error DZNRMS_VARS already defined
#else /* !DZNRMS_VARS */
#ifdef NDEBUG
#define DZNRMS_VARS DZNRMS_VARS1
#else /* !NDEBUG */
#define DZNRMS_VARS                                 \
  register const VD inf = _mm512_set1_pd(HUGE_VAL); \
  DZNRMS_VARS1
#endif /* ?NDEBUG */
#endif /* ?DZNRMS_VARS */

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

#ifdef DZNRMS_LOOP
#error DZNRMS_LOOP already defined
#else /* !DZNRMS_LOOP */
#define DZNRMS_LOOP(x,ec)                                             \
  for (fnat i = 0u; i < *m; i += VDL) {                               \
    register VD xi = _mm512_load_pd((x) + i); VDP(xi);                \
    register const VD reh = VDLSB(re); VDP(reh);                      \
    register const VD rep = _mm512_sub_pd(re, reh); VDP(rep);         \
    register const VD rfp = _mm512_scalef_pd(rf, reh); VDP(rfp);      \
    xi = VDABS(xi); VDP(xi);                                          \
    VDSORT(xi); VDP(xi);                                              \
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
#endif /* ?DZNRMS_LOOP */

#ifdef DZNRMS_RET
#error DZNRMS_RET already defined
#else /* !DZNRMS_RET */
#define DZNRMS_RET                      \
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
#endif /* ?DZNRMS_RET */

#endif /* !DZNRMS_H */