#include "sjac2s.h"

#include "scjac2.h"

#ifdef SJAC2_PARAMS
#error SJAC2_PARAMS already defined
#else /* !SJAC2_PARAMS */
#define SJAC2_PARAMS                                 \
  SCJAC2_PARAMS;                                     \
  register const VS huge = _mm512_set1_ps(FLT_MAX);  \
  register const VS be = _mm512_set1_ps(FLT_BIG_EXP)
#endif /* ?SJAC2_PARAMS */

#ifdef SJAC2_LOOP
#error SJAC2_LOOP already defined
#else /* !SJAC2_LOOP */
#define SJAC2_LOOP                                                                                                       \
register VS a1 = _mm512_load_ps(a11 + i); VSP(a1);                                                                       \
register VS a2 = _mm512_load_ps(a22 + i); VSP(a2);                                                                       \
register VS ar = _mm512_load_ps(a21 + i); VSP(ar);                                                                       \
register const VS e1 = _mm512_sub_ps(be, _mm512_getexp_ps(a1)); VSP(e1);                                                 \
register const VS e2 = _mm512_sub_ps(be, _mm512_getexp_ps(a2)); VSP(e2);                                                 \
register const VS er = _mm512_sub_ps(be, _mm512_getexp_ps(ar)); VSP(er);                                                 \
register VS es = _mm512_min_ps(_mm512_min_ps(e1, e2), _mm512_min_ps(er, huge)); VSP(es);                                 \
ar = _mm512_scalef_ps(ar, es); VSP(ar);                                                                                  \
a1 = _mm512_scalef_ps(a1, es); VSP(a1);                                                                                  \
a2 = _mm512_scalef_ps(a2, es); VSP(a2);                                                                                  \
es = VSNEG(es); VSP(es);                                                                                                 \
register const VS aa = VSABS(ar); VSP(aa);                                                                               \
register const VS as = VSSGN(ar); VSP(as);                                                                               \
register const VS ab = _mm512_scalef_ps(aa, onef); VSP(ab);                                                              \
register const VS ad = _mm512_sub_ps(a1, a2); VSP(ad);                                                                   \
register const VS t2 = VSOR(_mm512_min_ps(_mm512_max_ps(_mm512_div_ps(ab, VSABS(ad)), zerof), shf), VSSGN(ad)); VSP(t2); \
register const VS t1 = _mm512_div_ps(t2, _mm512_add_ps(onef, _mm512_sqrt_ps(_mm512_fmadd_ps(t2, t2, onef)))); VSP(t1);   \
register const VS s2 = _mm512_fmadd_ps(t1, t1, onef); VSP(s2);                                                           \
register const VS s1 = _mm512_sqrt_ps(s2); VSP(s1);                                                                      \
register const VS co = _mm512_div_ps(onef, s1); VSP(co);                                                                 \
register const VS si = VSXOR(_mm512_div_ps(t1, s1), as); VSP(si);                                                        \
_mm512_store_ps((cs + i), co);                                                                                           \
register const VS L1 = _mm512_div_ps(_mm512_fmadd_ps(t1, _mm512_fmadd_ps(a2, t1, ab), a1), s2); VSP(L1);                 \
_mm512_store_ps((sn + i), si);                                                                                           \
register const VS L2 = _mm512_div_ps(_mm512_fmadd_ps(t1, _mm512_fmsub_ps(a1, t1, ab), a2), s2); VSP(L2);                 \
_mm512_store_ps((l1 + i), _mm512_scalef_ps(L1, es));                                                                     \
register const MS P = _mm512_cmplt_ps_mask(L1, L2); MSP(P);                                                              \
_mm512_store_ps((l2 + i), _mm512_scalef_ps(L2, es));                                                                     \
p[i >> VSLlg] = MS2U(P)
#endif /* ?SJAC2_LOOP */

fint sjac2s_(const fnat n[static restrict 1], const float a11[static restrict VSL], const float a22[static restrict VSL], const float a21[static restrict VSL], float cs[static restrict VSL], float sn[static restrict VSL], float l1[static restrict VSL], float l2[static restrict VSL], unsigned p[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -10;
  if (*n & VSL_1)
    return -1;
  if (IS_NOT_ALIGNED(a11))
    return -2;
  if (IS_NOT_ALIGNED(a22))
    return -3;
  if (IS_NOT_ALIGNED(a21))
    return -4;
  if (IS_NOT_ALIGNED(cs))
    return -5;
  if (IS_NOT_ALIGNED(sn))
    return -6;
  if (IS_NOT_ALIGNED(l1))
    return -7;
  if (IS_NOT_ALIGNED(l2))
    return -8;
#endif /* !NDEBUG */

#ifdef _OPENMP
  fint th = 0;

#pragma omp parallel for default(none) shared(n,a11,a22,a21,cs,sn,l1,l2,p) reduction(max:th)
  for (fnat i = 0u; i < *n; i += VSL) {
    SJAC2_PARAMS;
    SJAC2_LOOP;
    th = imax(th, omp_get_thread_num());
  }

  return (th + 1);
#else /* !_OPENMP */
  SJAC2_PARAMS;

  for (fnat i = 0u; i < *n; i += VSL) {
    SJAC2_LOOP;
  }

  return 0;
#endif /* ?_OPENMP */
}
