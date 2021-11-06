#include "djac2f.h"

#include "dzjac2.h"

#ifdef BIG_EXP
#error BIG_EXP already defined
#else /* !BIG_EXP */
// (double)(DBL_MAX_EXP - 3)
#define BIG_EXP 1021.0
#endif /* ?BIG_EXP */

#ifdef DJAC2_PARAMS
#error DJAC2_PARAMS already defined
#else /* !DJAC2_PARAMS */
#define DJAC2_PARAMS                                 \
  DZJAC2_PARAMS;                                     \
  register const VD  huge = _mm512_set1_pd(DBL_MAX); \
  register const VD be = _mm512_set1_pd(BIG_EXP)
#endif /* ?DJAC2_PARAMS */

#ifdef DJAC2_LOOP
#error DJAC2_LOOP already defined
#else /* !DJAC2_LOOP */
#define DJAC2_LOOP                                                                                                     \
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
register const VD ab = _mm512_scalef_pd(aa, one); VDP(ab);                                                             \
register const VD ad = _mm512_sub_pd(a1, a2); VDP(ad);                                                                 \
register const VD t2 = VDOR(_mm512_min_pd(_mm512_max_pd(_mm512_div_pd(ab, VDABS(ad)), zero), sh), VDSGN(ad)); VDP(t2); \
register const VD t1 = _mm512_div_pd(t2, _mm512_add_pd(one, _mm512_sqrt_pd(_mm512_fmadd_pd(t2, t2, one)))); VDP(t1);   \
_mm512_store_pd((t + i), VDXOR(t1, as));                                                                               \
register VD C = _mm512_invsqrt_pd(_mm512_fmadd_pd(t1, t1, one)); VDP(C);                                               \
_mm512_store_pd((c + i), C);                                                                                           \
C = _mm512_mul_pd(C, C);                                                                                               \
register VD L1 = _mm512_fmadd_pd(t1, _mm512_fmadd_pd(a2, t1, ab), a1); VDP(L1);                                        \
register VD L2 = _mm512_fmadd_pd(t1, _mm512_fmsub_pd(a1, t1, ab), a2); VDP(L2);                                        \
L1 = _mm512_mul_pd(L1, C); VDP(L1);                                                                                    \
L2 = _mm512_mul_pd(L2, C); VDP(L2);                                                                                    \
register const MD P = _mm512_cmplt_pd_mask(L1, L2); MDP(P);                                                            \
p[i >> VDLlg] = MD2U(P);                                                                                               \
es = VDNEG(es);                                                                                                        \
L1 = _mm512_scalef_pd(L1, es);                                                                                         \
_mm512_store_pd((l1 + i), L1);                                                                                         \
L2 = _mm512_scalef_pd(L2, es);                                                                                         \
_mm512_store_pd((l2 + i), L2)
#endif /* ?DJAC2_LOOP */

fint djac2f_(const fnat n[static restrict 1], const double a11[static restrict VDL], const double a22[static restrict VDL], const double a21[static restrict VDL], double t[static restrict VDL], double c[static restrict VDL], double l1[static restrict VDL], double l2[static restrict VDL], unsigned p[static restrict 1])
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
  if (IS_NOT_ALIGNED(t))
    return -5;
  if (IS_NOT_ALIGNED(c))
    return -6;
  if (IS_NOT_ALIGNED(l1))
    return -7;
  if (IS_NOT_ALIGNED(l2))
    return -8;
#endif /* !NDEBUG */

#ifdef _OPENMP
  fint th = 0;

#pragma omp parallel for default(none) shared(n,a11,a22,a21,t,c,l1,l2,p) reduction(max:th)
  for (fnat i = 0u; i < *n; i += VDL) {
    DJAC2_PARAMS;
    DJAC2_LOOP;
    th = imax(th, omp_get_thread_num());
  }

  return (th + 1);
#else /* !_OPENMP */
  DJAC2_PARAMS;

  for (fnat i = 0u; i < *n; i += VDL) {
    DJAC2_LOOP;
  }

  return 0;
#endif /* ?_OPENMP */
}
