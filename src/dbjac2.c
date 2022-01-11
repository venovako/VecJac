#include "dbjac2.h"

#include "dzjac2.h"

#ifdef D8JAC2_PARAMS
#error D8JAC2_PARAMS already defined
#else /* !D8JAC2_PARAMS */
#define D8JAC2_PARAMS                                \
  DZJAC2_PARAMS;                                     \
  register const VD huge = _mm512_set1_pd(DBL_MAX);  \
  register const VD be = _mm512_set1_pd(DBL_BIG_EXP)
#endif /* ?D8JAC2_PARAMS */

#ifdef D8JAC2_LOOP
#error D8JAC2_LOOP already defined
#else /* !D8JAC2_LOOP */
#define D8JAC2_LOOP                                                                                                    \
register VD a1 = _mm512_load_pd(a11 + i); VDP(a1);                                                                     \
register VD a2 = _mm512_load_pd(a22 + i); VDP(a2);                                                                     \
register VD ar = _mm512_load_pd(a21 + i); VDP(ar);                                                                     \
register const VD e1 = _mm512_sub_pd(be, _mm512_getexp_pd(a1)); VDP(e1);                                               \
register const VD e2 = _mm512_sub_pd(be, _mm512_getexp_pd(a2)); VDP(e2);                                               \
register const VD er = _mm512_sub_pd(be, _mm512_getexp_pd(ar)); VDP(er);                                               \
register VD es = _mm512_min_pd(_mm512_min_pd(e1, e2), _mm512_min_pd(er, huge)); VDP(es);                               \
ar = _mm512_scalef_pd(ar, es); VDP(ar);                                                                                \
a1 = _mm512_scalef_pd(a1, es); VDP(a1);                                                                                \
a2 = _mm512_scalef_pd(a2, es); VDP(a2);                                                                                \
register const VD aa = VDABS(ar); VDP(aa);                                                                             \
register const VD as = VDSGN(ar); VDP(as);                                                                             \
es = VDNEG(es); VDP(es);                                                                                               \
register const VD ab = _mm512_scalef_pd(aa, one); VDP(ab);                                                             \
register const VD ad = _mm512_sub_pd(a1, a2); VDP(ad);                                                                 \
register const VD t2 = VDOR(_mm512_min_pd(_mm512_max_pd(_mm512_div_pd(ab, VDABS(ad)), zero), sh), VDSGN(ad)); VDP(t2); \
register const VD t1 = _mm512_div_pd(t2, _mm512_add_pd(one, _mm512_sqrt_pd(_mm512_fmadd_pd(t2, t2, one)))); VDP(t1);   \
register const VD s2 = _mm512_fmadd_pd(t1, t1, one); VDP(s2);                                                          \
register const VD s1 = _mm512_sqrt_pd(s2); VDP(s1);                                                                    \
register const VD co = _mm512_div_pd(one, s1); VDP(co);                                                                \
_mm512_store_pd((at + i), VDXOR(t1, as));                                                                              \
register const VD L1 = _mm512_div_pd(_mm512_fmadd_pd(t1, _mm512_fmadd_pd(a2, t1, ab), a1), s2); VDP(L1);               \
_mm512_store_pd((c + i), co);                                                                                          \
register const VD L2 = _mm512_div_pd(_mm512_fmadd_pd(t1, _mm512_fmsub_pd(a1, t1, ab), a2), s2); VDP(L2);               \
_mm512_store_pd((l1 + i), _mm512_scalef_pd(L1, es));                                                                   \
register const MD P = _mm512_cmplt_pd_mask(L1, L2); MDP(P);                                                            \
_mm512_store_pd((l2 + i), _mm512_scalef_pd(L2, es));                                                                   \
p[i >> VDLlg] = MD2U(P)
#endif /* ?D8JAC2_LOOP */

fint dbjac2_(const fnat n[static restrict 1], const double a11[static restrict VDL], const double a22[static restrict VDL], const double a21[static restrict VDL], double c[static restrict VDL], double at[static restrict VDL], double l1[static restrict VDL], double l2[static restrict VDL], unsigned p[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -10;
  if (*n & VDL_1)
    return -1;
  if (IS_NOT_ALIGNED(a11))
    return -2;
  if (IS_NOT_ALIGNED(a22))
    return -3;
  if (IS_NOT_ALIGNED(a21))
    return -4;
  if (IS_NOT_ALIGNED(c))
    return -5;
  if (IS_NOT_ALIGNED(at))
    return -6;
  if (IS_NOT_ALIGNED(l1))
    return -7;
  if (IS_NOT_ALIGNED(l2))
    return -8;
#endif /* !NDEBUG */

#ifdef _OPENMP
  fint th = 0;

#pragma omp parallel for default(none) shared(n,a11,a22,a21,c,at,l1,l2,p) reduction(max:th)
  for (fnat i = 0u; i < *n; i += VDL) {
    D8JAC2_PARAMS;
    D8JAC2_LOOP;
    th = imax(th, omp_get_thread_num());
  }

  return (th + 1);
#else /* !_OPENMP */
  D8JAC2_PARAMS;

  for (fnat i = 0u; i < *n; i += VDL) {
    D8JAC2_LOOP;
  }

  return 0;
#endif /* ?_OPENMP */
}