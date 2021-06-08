#include "zjac2.h"

#include "dzjac2.h"

#ifdef BIG_EXP
#error BIG_EXP already defined
#else /* !BIG_EXP */
// (double)(DBL_MAX_EXP - 4)
#define BIG_EXP 1020.0
#endif /* ?BIG_EXP */

#ifdef ZJAC2_PARAMS
#error ZJAC2_PARAMS already defined
#else /* !ZJAC2_PARAMS */
#define ZJAC2_PARAMS                                   \
  DZJAC2_PARAMS;                                       \
  register const VD be = _mm512_set1_pd(BIG_EXP);      \
  register const VD dtm = _mm512_set1_pd(DBL_TRUE_MIN)
#endif /* ?ZJAC2_PARAMS */

#ifdef ZJAC2_LOOP
#error ZJAC2_LOOP already defined
#else /* !ZJAC2_LOOP */
#define ZJAC2_LOOP                                                                                                     \
register VD a1 = _mm512_load_pd(a11 + i); VDP(a1);                                                                     \
register VD a2 = _mm512_load_pd(a22 + i); VDP(a2);                                                                     \
register VD ar = _mm512_load_pd(a21r + i); VDP(ar);                                                                    \
register VD ai = _mm512_load_pd(a21i + i); VDP(ai);                                                                    \
register const VD e1 = _mm512_sub_pd(be, _mm512_getexp_pd(a1)); VDP(e1);                                               \
register const VD e2 = _mm512_sub_pd(be, _mm512_getexp_pd(a2)); VDP(e2);                                               \
register const VD er = _mm512_sub_pd(be, _mm512_getexp_pd(ar)); VDP(er);                                               \
register const VD ei = _mm512_sub_pd(be, _mm512_getexp_pd(ai)); VDP(ei);                                               \
register const VD es = _mm512_min_pd(huge, _mm512_min_pd(_mm512_min_pd(e1, e2), _mm512_min_pd(er, ei))); VDP(es);      \
_mm512_store_pd((s + i), es);                                                                                          \
ar = _mm512_scalef_pd(ar, es); VDP(ar);                                                                                \
ai = _mm512_scalef_pd(ai, es); VDP(ai);                                                                                \
a1 = _mm512_scalef_pd(a1, es); VDP(a1);                                                                                \
a2 = _mm512_scalef_pd(a2, es); VDP(a2);                                                                                \
register const VD aa = _mm512_hypot_pd(ar, ai); VDP(aa);                                                               \
_mm512_store_pd((ca + i), VDOR(_mm512_min_pd(_mm512_div_pd(VDABS(ar), aa), one), VDSGN(ar)));                          \
_mm512_store_pd((sa + i), _mm512_div_pd(ai, _mm512_max_pd(aa, dtm)));                                                  \
register const VD ab = _mm512_scalef_pd(aa, one); VDP(ab);                                                             \
register const VD ad = _mm512_sub_pd(a1, a2); VDP(ad);                                                                 \
register const VD t2 = VDOR(_mm512_min_pd(_mm512_max_pd(_mm512_div_pd(ab, VDABS(ad)), zero), sh), VDSGN(ad)); VDP(t2); \
register const VD t1 = _mm512_div_pd(t2, _mm512_add_pd(one, _mm512_sqrt_pd(_mm512_fmadd_pd(t2, t2, one)))); VDP(t1);   \
_mm512_store_pd((t + i), t1);                                                                                          \
register const VD C = _mm512_invsqrt_pd(_mm512_fmadd_pd(t1, t1, one)); VDP(C);                                         \
_mm512_store_pd((c + i), C);                                                                                           \
register const VD L1 = _mm512_fmadd_pd(t1, _mm512_fmadd_pd(a2, t1, ab), a1); VDP(L1);                                  \
_mm512_store_pd((l1 + i), L1);                                                                                         \
register const VD L2 = _mm512_fmadd_pd(t1, _mm512_fmsub_pd(a1, t1, ab), a2); VDP(L2);                                  \
_mm512_store_pd((l2 + i), L2);                                                                                         \
register const MD P = _mm512_cmplt_pd_mask(L1, L2); MDP(P);                                                            \
p[i >> VDLlg] = MD2U(P)
#endif /* ?ZJAC2_LOOP */

int zjac2_(const fnat n[static restrict 1], const double a11[static restrict VDL], const double a22[static restrict VDL], const double a21r[static restrict VDL], const double a21i[static restrict VDL], double s[static restrict VDL], double t[static restrict VDL], double c[static restrict VDL], double ca[static restrict VDL], double sa[static restrict VDL], double l1[static restrict VDL], double l2[static restrict VDL], unsigned p[static restrict 1])
{
#ifdef _OPENMP
  int th = 0;

#pragma omp parallel for default(none) shared(n,a11,a22,a21r,a21i,s,t,c,ca,sa,l1,l2,p) reduction(max:th)
  for (fnat i = 0u; i < *n; i += VDL) {
    ZJAC2_PARAMS;
    ZJAC2_LOOP;
    th = imax(th, omp_get_thread_num());
  }

  return (th + 1);
#else /* !_OPENMP */
  ZJAC2_PARAMS;

  for (fnat i = 0u; i < *n; i += VDL) {
    ZJAC2_LOOP;
  }

  return 0;
#endif /* ?_OPENMP */
}
