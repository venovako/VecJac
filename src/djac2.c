#include "djac2.h"

#include "dzjac2.h"

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
register const VD es = _mm512_min_pd(_mm512_min_pd(e1, e2), _mm512_min_pd(er, huge)); VDP(es);                         \
_mm512_store_pd((s + i), es);                                                                                          \
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
_mm512_store_pd((c + i), _mm512_invsqrt_pd(_mm512_fmadd_pd(t1, t1, one)));                                             \
_mm512_store_pd((l1 + i), _mm512_fmadd_pd(t1, _mm512_fmadd_pd(a2, t1, ab), a1));                                       \
_mm512_store_pd((l2 + i), _mm512_fmadd_pd(t1, _mm512_fmsub_pd(a1, t1, ab), a2))
#endif /* ?DJAC2_LOOP */

int djac2_(const fnat n[static restrict 1], const double a11[static restrict VDL], const double a22[static restrict VDL], const double a21[static restrict VDL], double s[static restrict VDL], double t[static restrict VDL], double c[static restrict VDL], double l1[static restrict VDL], double l2[static restrict VDL])
{
#ifdef _OPENMP
  int th = 0;

#pragma omp parallel for default(none) shared(n,a11,a22,a21,s,t,c,l1,l2) reduction(max:th)
  for (fnat i = 0u; i < *n; i += VDL) {
    DZJAC2_PARAMS;
    DJAC2_LOOP;
    th = imax(th, omp_get_thread_num());
  }

  return (th + 1);
#else /* !_OPENMP */
  DZJAC2_PARAMS;

  for (fnat i = 0u; i < *n; i += VDL) {
    DJAC2_LOOP;
  }

  return 0;
#endif /* ?_OPENMP */
}
