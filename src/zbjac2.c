#include "zbjac2.h"

#include "dzjac2.h"

#ifdef Z8JAC2_PARAMS
#error Z8JAC2_PARAMS already defined
#else /* !Z8JAC2_PARAMS */
#define Z8JAC2_PARAMS                                  \
  DZJAC2_PARAMS;                                       \
  register const VD huge = _mm512_set1_pd(DBL_MAX);    \
  register const VD be = _mm512_set1_pd(DBL_BIG_EXP);  \
  register const VD dtm = _mm512_set1_pd(DBL_TRUE_MIN)
#endif /* ?Z8JAC2_PARAMS */

#ifdef Z8JAC2_LOOP
#error Z8JAC2_LOOP already defined
#else /* !Z8JAC2_LOOP */
#define Z8JAC2_LOOP                                                                                                    \
register VD a1 = _mm512_load_pd(a11 + i); VDP(a1);                                                                     \
register VD a2 = _mm512_load_pd(a22 + i); VDP(a2);                                                                     \
register VD ar = _mm512_load_pd(a21r + i); VDP(ar);                                                                    \
register VD ai = _mm512_load_pd(a21i + i); VDP(ai);                                                                    \
register const VD e1 = _mm512_sub_pd(be, _mm512_getexp_pd(a1)); VDP(e1);                                               \
register const VD e2 = _mm512_sub_pd(be, _mm512_getexp_pd(a2)); VDP(e2);                                               \
register const VD er = _mm512_sub_pd(be, _mm512_getexp_pd(ar)); VDP(er);                                               \
register const VD ei = _mm512_sub_pd(be, _mm512_getexp_pd(ai)); VDP(ei);                                               \
register VD es = _mm512_min_pd(huge, _mm512_min_pd(_mm512_min_pd(e1, e2), _mm512_min_pd(er, ei))); VDP(es);            \
ar = _mm512_scalef_pd(ar, es); VDP(ar);                                                                                \
ai = _mm512_scalef_pd(ai, es); VDP(ai);                                                                                \
a1 = _mm512_scalef_pd(a1, es); VDP(a1);                                                                                \
a2 = _mm512_scalef_pd(a2, es); VDP(a2);                                                                                \
register VD ar_ = VDABS(ar); VDP(ar_);                                                                                 \
register VD ai_ = VDABS(ai); VDP(ai_);                                                                                 \
register const VD as = VDSGN(ar); VDP(as);                                                                             \
register const VD am = _mm512_min_pd(ar_, ai_); VDP(am);                                                               \
register const VD aM = _mm512_max_pd(ar_, ai_); VDP(aM);                                                               \
es = VDNEG(es); VDP(es);                                                                                               \
register VD aa = _mm512_max_pd(_mm512_div_pd(am, aM), zero); VDP(aa);                                                  \
aa = _mm512_mul_pd(_mm512_sqrt_pd(_mm512_fmadd_pd(aa, aa, one)), aM); VDP(aa);                                         \
ar_ = VDOR(_mm512_min_pd(_mm512_div_pd(ar_, aa), one), as); VDP(ar_);                                                  \
ai_ = _mm512_div_pd(ai, _mm512_max_pd(aa, dtm)); VDP(ai_);                                                             \
register const VD ab = _mm512_scalef_pd(aa, one); VDP(ab);                                                             \
register const VD ad = _mm512_sub_pd(a1, a2); VDP(ad);                                                                 \
register const VD t2 = VDOR(_mm512_min_pd(_mm512_max_pd(_mm512_div_pd(ab, VDABS(ad)), zero), sh), VDSGN(ad)); VDP(t2); \
register const VD t1 = _mm512_div_pd(t2, _mm512_add_pd(one, _mm512_sqrt_pd(_mm512_fmadd_pd(t2, t2, one)))); VDP(t1);   \
register const VD s2 = _mm512_fmadd_pd(t1, t1, one); VDP(s2);                                                          \
register const VD s1 = _mm512_sqrt_pd(s2); VDP(s1);                                                                    \
register const VD co = _mm512_div_pd(one, s1); VDP(co);                                                                \
ar_ = _mm512_mul_pd(ar_, t1); VDP(ar_);                                                                                \
_mm512_store_pd((c + i), co);                                                                                          \
ai_ = _mm512_mul_pd(ai_, t1); VDP(ai_);                                                                                \
_mm512_store_pd((cat + i), ar_);                                                                                       \
register const VD L1 = _mm512_div_pd(_mm512_fmadd_pd(t1, _mm512_fmadd_pd(a2, t1, ab), a1), s2); VDP(L1);               \
_mm512_store_pd((sat + i), ai_);                                                                                       \
register const VD L2 = _mm512_div_pd(_mm512_fmadd_pd(t1, _mm512_fmsub_pd(a1, t1, ab), a2), s2); VDP(L2);               \
_mm512_store_pd((l1 + i), _mm512_scalef_pd(L1, es));                                                                   \
register const MD P = _mm512_cmplt_pd_mask(L1, L2); MDP(P);                                                            \
_mm512_store_pd((l2 + i), _mm512_scalef_pd(L2, es));                                                                   \
p[i >> VDLlg] = MD2U(P)
#endif /* ?Z8JAC2_LOOP */

#ifdef Z8JACL_LOOP
#error Z8JACL_LOOP already defined
#else /* !Z8JACL_LOOP */
#define Z8JACL_LOOP                                                                                                    \
register VD a1 = _mm512_load_pd(a11 + i); VDP(a1);                                                                     \
register VD a2 = _mm512_load_pd(a22 + i); VDP(a2);                                                                     \
register VD ar = _mm512_load_pd(a21r + i); VDP(ar);                                                                    \
register VD ai = _mm512_load_pd(a21i + i); VDP(ai);                                                                    \
register const VD e1 = _mm512_sub_pd(be, _mm512_getexp_pd(a1)); VDP(e1);                                               \
register const VD e2 = _mm512_sub_pd(be, _mm512_getexp_pd(a2)); VDP(e2);                                               \
register const VD er = _mm512_sub_pd(be, _mm512_getexp_pd(ar)); VDP(er);                                               \
register const VD ei = _mm512_sub_pd(be, _mm512_getexp_pd(ai)); VDP(ei);                                               \
register VD es = _mm512_min_pd(huge, _mm512_min_pd(_mm512_min_pd(e1, e2), _mm512_min_pd(er, ei))); VDP(es);            \
ar = _mm512_scalef_pd(ar, es); VDP(ar);                                                                                \
ai = _mm512_scalef_pd(ai, es); VDP(ai);                                                                                \
a1 = _mm512_scalef_pd(a1, es); VDP(a1);                                                                                \
a2 = _mm512_scalef_pd(a2, es); VDP(a2);                                                                                \
register VD ar_ = VDABS(ar); VDP(ar_);                                                                                 \
register VD ai_ = VDABS(ai); VDP(ai_);                                                                                 \
register const VD as = VDSGN(ar); VDP(as);                                                                             \
register const VD am = _mm512_min_pd(ar_, ai_); VDP(am);                                                               \
register const VD aM = _mm512_max_pd(ar_, ai_); VDP(aM);                                                               \
es = VDNEG(es); VDP(es);                                                                                               \
register VD aa = _mm512_max_pd(_mm512_div_pd(am, aM), zero); VDP(aa);                                                  \
aa = _mm512_mul_pd(_mm512_sqrt_pd(_mm512_fmadd_pd(aa, aa, one)), aM); VDP(aa);                                         \
ar_ = VDOR(_mm512_min_pd(_mm512_div_pd(ar_, aa), one), as); VDP(ar_);                                                  \
ai_ = _mm512_div_pd(ai, _mm512_max_pd(aa, dtm)); VDP(ai_);                                                             \
register const VD ab = _mm512_scalef_pd(aa, one); VDP(ab);                                                             \
register const VD ad = _mm512_sub_pd(a1, a2); VDP(ad);                                                                 \
register const VD t2 = VDOR(_mm512_min_pd(_mm512_max_pd(_mm512_div_pd(ab, VDABS(ad)), zero), sh), VDSGN(ad)); VDP(t2); \
register const VD t1 = _mm512_div_pd(t2, _mm512_add_pd(one, _mm512_sqrt_pd(_mm512_fmadd_pd(t2, t2, one)))); VDP(t1);   \
register const VD s2 = _mm512_fmadd_pd(t1, t1, one); VDP(s2);                                                          \
register const VD s1 = _mm512_sqrt_pd(s2); VDP(s1);                                                                    \
register const VD co = _mm512_div_pd(one, s1); VDP(co);                                                                \
ar_ = _mm512_div_pd(_mm512_mul_pd(ar_, t1), s1); VDP(ar_);                                                             \
_mm512_store_pd((c + i), co);                                                                                          \
ai_ = _mm512_div_pd(_mm512_mul_pd(ai_, t1), s1); VDP(ai_);                                                             \
_mm512_store_pd((cat + i), ar_);                                                                                       \
register const VD L1 = _mm512_div_pd(_mm512_fmadd_pd(t1, _mm512_fmadd_pd(a2, t1, ab), a1), s2); VDP(L1);               \
_mm512_store_pd((sat + i), ai_);                                                                                       \
register const VD L2 = _mm512_div_pd(_mm512_fmadd_pd(t1, _mm512_fmsub_pd(a1, t1, ab), a2), s2); VDP(L2);               \
_mm512_store_pd((l1 + i), _mm512_scalef_pd(L1, es));                                                                   \
register const MD P = _mm512_cmplt_pd_mask(L1, L2); MDP(P);                                                            \
_mm512_store_pd((l2 + i), _mm512_scalef_pd(L2, es));                                                                   \
p[i >> VDLlg] = MD2U(P)
#endif /* ?Z8JACL_LOOP */

fint zbjac2_(const fint n[static restrict 1], const double a11[static restrict VDL], const double a22[static restrict VDL], const double a21r[static restrict VDL], const double a21i[static restrict VDL], double c[static restrict VDL], double cat[static restrict VDL], double sat[static restrict VDL], double l1[static restrict VDL], double l2[static restrict VDL], unsigned p[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -12;
  if (*n & VDL_1)
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
    fint th = 0;

#pragma omp parallel for default(none) shared(_n,a11,a22,a21r,a21i,c,cat,sat,l1,l2,p) reduction(max:th)
    for (fnat i = 0u; i < _n; i += VDL) {
      Z8JAC2_PARAMS;
      Z8JACL_LOOP;
      th = imax(th, omp_get_thread_num());
    }

    return (th + 1);
#else /* !_OPENMP */
    Z8JAC2_PARAMS;

    for (fnat i = 0u; i < _n; i += VDL) {
      Z8JACL_LOOP;
    }

    return 0;
#endif /* ?_OPENMP */
  }
  else {
    const fnat _n = (fnat)*n;
#ifdef _OPENMP
    fint th = 0;

#pragma omp parallel for default(none) shared(_n,a11,a22,a21r,a21i,c,cat,sat,l1,l2,p) reduction(max:th)
    for (fnat i = 0u; i < _n; i += VDL) {
      Z8JAC2_PARAMS;
      Z8JAC2_LOOP;
      th = imax(th, omp_get_thread_num());
    }

    return (th + 1);
#else /* !_OPENMP */
    Z8JAC2_PARAMS;

    for (fnat i = 0u; i < _n; i += VDL) {
      Z8JAC2_LOOP;
    }

    return 0;
#endif /* ?_OPENMP */
  }
}

// for internal use only
fint zbjac2i(const fnat n[static restrict 1], const double a11[static restrict VDL], const double a22[static restrict VDL], const double a21r[static restrict VDL], const double a21i[static restrict VDL], double c[static restrict VDL], double cat[static restrict VDL], double sat[static restrict VDL], double l1[static restrict VDL], double l2[static restrict VDL], unsigned p[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -12;
  if (*n & VDL_1)
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

#ifdef _OPENMP
  fint th = 0;

#pragma omp parallel for default(none) shared(n,a11,a22,a21r,a21i,c,cat,sat,l1,l2,p) reduction(max:th)
  for (fnat i = 0u; i < *n; i += VDL) {
    if (p[i >> VDLlg]) {
      Z8JAC2_PARAMS;
      Z8JAC2_LOOP;
    }
    th = imax(th, omp_get_thread_num());
  }

  return (th + 1);
#else /* !_OPENMP */
  Z8JAC2_PARAMS;

  for (fnat i = 0u; i < *n; i += VDL) {
    if (p[i >> VDLlg]) {
      Z8JAC2_LOOP;
    }
  }

  return 0;
#endif /* ?_OPENMP */
}
