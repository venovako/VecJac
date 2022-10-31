#include "cbjac2.h"

#include "scjac2.h"

#ifdef C8JAC2_PARAMS
#error C8JAC2_PARAMS already defined
#else /* !C8JAC2_PARAMS */
#define C8JAC2_PARAMS                                  \
  SCJAC2_PARAMS;                                       \
  register const VS huge = _mm512_set1_ps(FLT_MAX);    \
  register const VS be = _mm512_set1_ps(FLT_BIG_EXP);  \
  register const VS ftm = _mm512_set1_ps(FLT_TRUE_MIN)
#endif /* ?C8JAC2_PARAMS */

#ifdef C8JAC2_LOOP
#error C8JAC2_LOOP already defined
#else /* !C8JAC2_LOOP */
#define C8JAC2_LOOP                                                                                                      \
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
ar_ = _mm512_mul_ps(ar_, t1); VSP(ar_);                                                                                  \
_mm512_store_ps((c + i), co);                                                                                            \
ai_ = _mm512_mul_ps(ai_, t1); VSP(ai_);                                                                                  \
_mm512_store_ps((cat + i), ar_);                                                                                         \
register const VS L1 = _mm512_div_ps(_mm512_fmadd_ps(t1, _mm512_fmadd_ps(a2, t1, ab), a1), s2); VSP(L1);                 \
_mm512_store_ps((sat + i), ai_);                                                                                         \
register const VS L2 = _mm512_div_ps(_mm512_fmadd_ps(t1, _mm512_fmsub_ps(a1, t1, ab), a2), s2); VSP(L2);                 \
_mm512_store_ps((l1 + i), _mm512_scalef_ps(L1, es));                                                                     \
register const MS P = _mm512_cmplt_ps_mask(L1, L2); MSP(P);                                                              \
_mm512_store_ps((l2 + i), _mm512_scalef_ps(L2, es))
#endif /* ?C8JAC2_LOOP */

// return the sines instead of the tangents
#ifdef C8JACL_LOOP
#error C8JACL_LOOP already defined
#else /* !C8JACL_LOOP */
#define C8JACL_LOOP                                                                                                      \
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
ar_ = _mm512_div_ps(_mm512_mul_ps(ar_, t1), s1); VSP(ar_);                                                               \
_mm512_store_ps((c + i), co);                                                                                            \
ai_ = _mm512_div_ps(_mm512_mul_ps(ai_, t1), s1); VSP(ai_);                                                               \
_mm512_store_ps((cat + i), ar_);                                                                                         \
register const VS L1 = _mm512_div_ps(_mm512_fmadd_ps(t1, _mm512_fmadd_ps(a2, t1, ab), a1), s2); VSP(L1);                 \
_mm512_store_ps((sat + i), ai_);                                                                                         \
register const VS L2 = _mm512_div_ps(_mm512_fmadd_ps(t1, _mm512_fmsub_ps(a1, t1, ab), a2), s2); VSP(L2);                 \
_mm512_store_ps((l1 + i), _mm512_scalef_ps(L1, es));                                                                     \
register const MS P = _mm512_cmplt_ps_mask(L1, L2); MSP(P);                                                              \
_mm512_store_ps((l2 + i), _mm512_scalef_ps(L2, es))
#endif /* ?C8JACL_LOOP */

fint cbjac2_(const fint n[static restrict 1], const float a11[static restrict VSL], const float a22[static restrict VSL], const float a21r[static restrict VSL], const float a21i[static restrict VSL], float c[static restrict VSL], float cat[static restrict VSL], float sat[static restrict VSL], float l1[static restrict VSL], float l2[static restrict VSL], unsigned p[static restrict 1])
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
  if (IS_NOT_ALIGNED(c))
    return -6;
  if (IS_NOT_ALIGNED(cat))
    return -7;
  if (IS_NOT_ALIGNED(sat))
    return -8;
  if (IS_NOT_ALIGNED(l1))
    return -9;
  if (IS_NOT_ALIGNED(l2))
    return -10;
#endif /* !NDEBUG */

  if (*n < 0) {
    const fnat _n = (fnat)-*n;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(_n,a11,a22,a21r,a21i,c,cat,sat,l1,l2,p)
    for (fnat i = 0u; i < _n; i += VSL) {
      C8JAC2_PARAMS;
      C8JACL_LOOP;
      p[i >> VSLlg] = MS2U(P);
    }
    return 1;
#else /* !_OPENMP */
    C8JAC2_PARAMS;
    for (fnat i = 0u; i < _n; i += VSL) {
      C8JACL_LOOP;
      p[i >> VSLlg] = MS2U(P);
    }
    return 0;
#endif
  }
  else {
    const fnat _n = (fnat)*n;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(_n,a11,a22,a21r,a21i,c,cat,sat,l1,l2,p)
    for (fnat i = 0u; i < _n; i += VSL) {
      C8JAC2_PARAMS;
      C8JAC2_LOOP;
      p[i >> VSLlg] = MS2U(P);
    }
    return 1;
#else /* !_OPENMP */
    C8JAC2_PARAMS;
    for (fnat i = 0u; i < _n; i += VSL) {
      C8JAC2_LOOP;
      p[i >> VSLlg] = MS2U(P);
    }
    return 0;
#endif /* ?_OPENMP */
  }
}

// store the secants instead of the cosines
#ifdef C8JACI_LOOP
#error C8JACI_LOOP already defined
#else /* !C8JACI_LOOP */
#define C8JACI_LOOP                                                                                                      \
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
ar_ = _mm512_mul_ps(ar_, t1); VSP(ar_);                                                                                  \
_mm512_store_ps((c + i), s1);                                                                                            \
ai_ = _mm512_mul_ps(ai_, t1); VSP(ai_);                                                                                  \
_mm512_store_ps((cat + i), ar_);                                                                                         \
register const VS L1 = _mm512_div_ps(_mm512_fmadd_ps(t1, _mm512_fmadd_ps(a2, t1, ab), a1), s2); VSP(L1);                 \
_mm512_store_ps((sat + i), ai_);                                                                                         \
register const VS L2 = _mm512_div_ps(_mm512_fmadd_ps(t1, _mm512_fmsub_ps(a1, t1, ab), a2), s2); VSP(L2);                 \
_mm512_store_ps((l1 + i), _mm512_scalef_ps(L1, es));                                                                     \
register const MS P = _mm512_cmplt_ps_mask(L1, L2); MSP(P);                                                              \
_mm512_store_ps((l2 + i), _mm512_scalef_ps(L2, es))
#endif /* ?C8JACI_LOOP */

// for internal use only
fint cbjac2i(const fint n[static restrict 1], const float a11[static restrict VSL], const float a22[static restrict VSL], const float a21r[static restrict VSL], const float a21i[static restrict VSL], float c[static restrict VSL], float cat[static restrict VSL], float sat[static restrict VSL], float l1[static restrict VSL], float l2[static restrict VSL], unsigned p[static restrict 1])
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
  if (IS_NOT_ALIGNED(c))
    return -6;
  if (IS_NOT_ALIGNED(cat))
    return -7;
  if (IS_NOT_ALIGNED(sat))
    return -8;
  if (IS_NOT_ALIGNED(l1))
    return -9;
  if (IS_NOT_ALIGNED(l2))
    return -10;
#endif /* !NDEBUG */

  if (*n < 0) {
    const fnat _n = (fnat)-*n;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(_n,a11,a22,a21r,a21i,c,cat,sat,l1,l2,p)
    for (fnat i = 0u; i < _n; i += VSL) {
      const fnat j = (i >> VSLlg);
      if (p[j] & 0xFFFFu) {
        C8JAC2_PARAMS;
        C8JACI_LOOP;
        p[j] = ((p[j] & 0xFFFF0000u) | MS2U(P));
      }
    }
    return 1;
#else /* !_OPENMP */
    C8JAC2_PARAMS;
    for (fnat i = 0u; i < _n; i += VSL) {
      const fnat j = (i >> VSLlg);
      if (p[j] & 0xFFFFu) {
        C8JACI_LOOP;
        p[j] = ((p[j] & 0xFFFF0000u) | MS2U(P));
      }
    }
    return 0;
#endif /* ?_OPENMP */
  }
  else {
    const fnat _n = (fnat)*n;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(_n,a11,a22,a21r,a21i,c,cat,sat,l1,l2,p)
    for (fnat i = 0u; i < _n; i += VSL) {
      const fnat j = (i >> VSLlg);
      if (p[j] & 0xFFFFu) {
        C8JAC2_PARAMS;
        C8JAC2_LOOP;
        p[j] = ((p[j] & 0xFFFF0000u) | MS2U(P));
      }
    }
    return 1;
#else /* !_OPENMP */
    C8JAC2_PARAMS;
    for (fnat i = 0u; i < _n; i += VSL) {
      const fnat j = (i >> VSLlg);
      if (p[j] & 0xFFFFu) {
        C8JAC2_LOOP;
        p[j] = ((p[j] & 0xFFFF0000u) | MS2U(P));
      }
    }
    return 0;
#endif /* ?_OPENMP */
  }
}
