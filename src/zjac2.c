#include "zjac2.h"

#include "dzjac2.h"
#include "zcmplx.h"

#ifdef ZJAC2_PARAMS
#error ZJAC2_PARAMS already defined
#else /* !ZJAC2_PARAMS */
#define ZJAC2_PARAMS                                   \
  DZJAC2_PARAMS;                                       \
  register const VD dtm = _mm512_set1_pd(DBL_TRUE_MIN)
#endif /* ?ZJAC2_PARAMS */

#ifdef ZJAC2_LOOP
#error ZJAC2_LOOP already defined
#else /* !ZJAC2_LOOP */
#define ZJAC2_LOOP                                                                                            \
register const VD a1 = _mm512_load_pd(a11 + i);                                                               \
register const VD a2 = _mm512_load_pd(a22 + i);                                                               \
register const VD ar = _mm512_load_pd(a21r + i);                                                              \
register const VD ai = _mm512_load_pd(a21i + i);                                                              \
register const VD aa = _mm512_hypot_pd(ar, ai);                                                               \
_mm512_store_pd((ca + i), VDOR(_mm512_min_pd(_mm512_div_pd(VDABS(ar), aa), one), VDSGN(ar)));                 \
_mm512_store_pd((sa + i), _mm512_div_pd(ai, _mm512_max_pd(aa, dtm)));                                         \
register const VD ab = _mm512_scalef_pd(aa, one);                                                             \
register const VD ad = _mm512_sub_pd(a1, a2);                                                                 \
register const VD t2 = VDOR(_mm512_min_pd(_mm512_max_pd(_mm512_div_pd(ab, VDABS(ad)), zero), sh), VDSGN(ad)); \
register const VD t1 = _mm512_div_pd(t2, _mm512_add_pd(one, _mm512_sqrt_pd(_mm512_fmadd_pd(t2, t2, one))));   \
_mm512_store_pd((t + i), t1);                                                                                 \
register const VD c1 = _mm512_invsqrt_pd(_mm512_fmadd_pd(t1, t1, one));                                       \
_mm512_store_pd((c + i), c1);                                                                                 \
_mm512_store_pd((l1 + i), _mm512_fmadd_pd(t1, _mm512_fmadd_pd(a2, t1, ab), a1));                              \
_mm512_store_pd((l2 + i), _mm512_fmadd_pd(t1, _mm512_fmsub_pd(a1, t1, ab), a2))
#endif /* ?ZJAC2_LOOP */

extern int zjac2_(const fnat n[static restrict 1], const double a11[static restrict VDL], const double a22[static restrict VDL], const double a21r[static restrict VDL], const double a21i[static restrict VDL], double t[static restrict VDL], double c[static restrict VDL], double ca[static restrict VDL], double sa[static restrict VDL], double l1[static restrict VDL], double l2[static restrict VDL])
{
#ifdef _OPENMP
  int th = 0;

#pragma omp parallel for default(none) shared(n,a11,a22,a21r,a21i,t,c,ca,sa,l1,l2) reduction(max:th)
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
