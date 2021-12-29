#include "cjac2s.h"

#include "scjac2.h"

#ifdef CJAC2_PARAMS
#error CJAC2_PARAMS already defined
#else /* !CJAC2_PARAMS */
#define CJAC2_PARAMS                                   \
  SCJAC2_PARAMS;                                       \
  register const VS huge = _mm512_set1_ps(FLT_MAX);    \
  register const VS be = _mm512_set1_ps(FLT_BIG_EXP);  \
  register const VS ftm = _mm512_set1_ps(FLT_TRUE_MIN)
#endif /* ?CJAC2_PARAMS */

#ifdef CJAC2_LOOP
#error CJAC2_LOOP already defined
#else /* !CJAC2_LOOP */
#define CJAC2_LOOP                                                                                                       \
register VS a1 = _mm512_load_ps(a11 + i); VSP(a1);                                                                       \
register VS a2 = _mm512_load_ps(a22 + i); VSP(a2);                                                                       \
register VS ar = _mm512_load_ps(a21r + i); VSP(ar);                                                                      \
register VS ai = _mm512_load_ps(a21i + i); VSP(ai);                                                                      \
register const VS e1 = _mm512_sub_ps(be, _mm512_getexp_ps(a1)); VSP(e1);                                                 \
register const VS e2 = _mm512_sub_ps(be, _mm512_getexp_ps(a2)); VSP(e2);                                                 \
register const VS er = _mm512_sub_ps(be, _mm512_getexp_ps(ar)); VSP(er);                                                 \
register const VS ei = _mm512_sub_ps(be, _mm512_getexp_ps(ai)); VSP(ei);                                                 \
register VS es = _mm512_min_ps(huge, _mm512_min_ps(_mm512_min_ps(e1, e2), _mm512_min_ps(er, ei))); VSP(es);              \
ar = _mm512_scalef_ps(ar, es); VSP(ar);                                                                                  \
ai = _mm512_scalef_ps(ai, es); VSP(ai);                                                                                  \
a1 = _mm512_scalef_ps(a1, es); VSP(a1);                                                                                  \
a2 = _mm512_scalef_ps(a2, es); VSP(a2);                                                                                  \
register VS ar_ = VSABS(ar); VSP(ar_);                                                                                   \
register VS ai_ = VSABS(ai); VSP(ai_);                                                                                   \
register const VS as = VSSGN(ar); VSP(as);                                                                               \
register const VS am = _mm512_min_ps(ar_, ai_); VSP(am);                                                                 \
register const VS aM = _mm512_max_ps(ar_, ai_); VSP(aM);                                                                 \
es = VSNEG(es); VSP(es);                                                                                                 \
register VS aa = _mm512_max_ps(_mm512_div_ps(am, aM), zerof); VSP(aa);                                                   \
aa = _mm512_mul_ps(_mm512_sqrt_ps(_mm512_fmadd_ps(aa, aa, onef)), aM); VSP(aa);                                          \
ar_ = VSOR(_mm512_min_ps(_mm512_div_ps(ar_, aa), onef), as); VSP(ar_);                                                   \
ai_ = _mm512_div_ps(ai, _mm512_max_ps(aa, ftm)); VSP(ai_);                                                               \
register const VS ab = _mm512_scalef_ps(aa, onef); VSP(ab);                                                              \
register const VS ad = _mm512_sub_ps(a1, a2); VSP(ad);                                                                   \
register const VS t2 = VSOR(_mm512_min_ps(_mm512_max_ps(_mm512_div_ps(ab, VSABS(ad)), zerof), shf), VSSGN(ad)); VSP(t2); \
register const VS t1 = _mm512_div_ps(t2, _mm512_add_ps(onef, _mm512_sqrt_ps(_mm512_fmadd_ps(t2, t2, onef)))); VSP(t1);   \
register const VS s2 = _mm512_fmadd_ps(t1, t1, onef); VSP(s2);                                                           \
register const VS s1 = _mm512_sqrt_ps(s2); VSP(s1);                                                                      \
register const VS co = _mm512_div_ps(onef, s1); VSP(co);                                                                 \
register const VS si = _mm512_div_ps(t1, s1); VSP(si);                                                                   \
_mm512_store_ps((cs + i), co);                                                                                           \
ar_ = _mm512_mul_ps(ar_, si); VSP(ar_);                                                                                  \
ai_ = _mm512_mul_ps(ai_, si); VSP(ai_);                                                                                  \
_mm512_store_ps((ca + i), ar_);                                                                                          \
register const VS L1 = _mm512_div_ps(_mm512_fmadd_ps(t1, _mm512_fmadd_ps(a2, t1, ab), a1), s2); VSP(L1);                 \
_mm512_store_ps((sa + i), ai_);                                                                                          \
register const VS L2 = _mm512_div_ps(_mm512_fmadd_ps(t1, _mm512_fmsub_ps(a1, t1, ab), a2), s2); VSP(L2);                 \
_mm512_store_ps((l1 + i), _mm512_scalef_ps(L1, es));                                                                     \
register const MS P = _mm512_cmplt_ps_mask(L1, L2); MSP(P);                                                              \
_mm512_store_ps((l2 + i), _mm512_scalef_ps(L2, es));                                                                     \
p[i >> VSLlg] = MS2U(P)
#endif /* ?CJAC2_LOOP */

fint cjac2s_(const fnat n[static restrict 1], const float a11[static restrict VSL], const float a22[static restrict VSL], const float a21r[static restrict VSL], const float a21i[static restrict VSL], float cs[static restrict VSL], float ca[static restrict VSL], float sa[static restrict VSL], float l1[static restrict VSL], float l2[static restrict VSL], unsigned p[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -12;
  if (*n & VSL_1)
    return -1;
  if (IS_NOT_ALIGNED(a11))
    return -2;
  if (IS_NOT_ALIGNED(a22))
    return -3;
  if (IS_NOT_ALIGNED(a21r))
    return -4;
  if (IS_NOT_ALIGNED(a21i))
    return -5;
  if (IS_NOT_ALIGNED(cs))
    return -6;
  if (IS_NOT_ALIGNED(ca))
    return -7;
  if (IS_NOT_ALIGNED(sa))
    return -8;
  if (IS_NOT_ALIGNED(l1))
    return -9;
  if (IS_NOT_ALIGNED(l2))
    return -10;
#endif /* !NDEBUG */

#ifdef _OPENMP
  fint th = 0;

#pragma omp parallel for default(none) shared(n,a11,a22,a21r,a21i,cs,ca,sa,l1,l2,p) reduction(max:th)
  for (fnat i = 0u; i < *n; i += VSL) {
    CJAC2_PARAMS;
    CJAC2_LOOP;
    th = imax(th, omp_get_thread_num());
  }

  return (th + 1);
#else /* !_OPENMP */
  CJAC2_PARAMS;

  for (fnat i = 0u; i < *n; i += VSL) {
    CJAC2_LOOP;
  }

  return 0;
#endif /* ?_OPENMP */
}
