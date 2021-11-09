#ifndef DZNRMX_H
#define DZNRMX_H

#include "vecdef.h"

#ifdef DZNRMX_LOOP
#error DZNRMX_LOOP already defined
#else /* !DZNRMX_LOOP */
#ifdef _OPENMP
#define DZNRMX_LOOP(A,ldA)                                                     \
  for (fnat j = 0u; j < *n; ++j) {                                             \
    register const VD _zero = _mm512_set1_pd(-0.0);                            \
    register const VD inf = _mm512_set1_pd(HUGE_VAL);                          \
    register VD x = _mm512_setzero_pd();                                       \
    const double *const Aj = (A) + j * (size_t)(*(ldA));                       \
    for (fnat i = 0u; i < *m; i += VDL)                                        \
      x = _mm512_max_pd(_mm512_min_pd(VDABS(_mm512_load_pd(Aj + i)), inf), x); \
    y = fmax(_mm512_reduce_max_pd(x), y);                                      \
  }
#else /* !_OPENMP */
#define DZNRMX_LOOP(A,ldA)                                                             \
  for (fnat j = 0u; j < *n; ++j) {                                                     \
    const double *const Aj = (A) + j * (size_t)(*(ldA));                               \
    for (fnat i = 0u; i < *m; i += VDL) {                                              \
      x = _mm512_max_pd(_mm512_min_pd(VDABS(_mm512_load_pd(Aj + i)), inf), x); VDP(x); \
    }                                                                                  \
  }
#endif /* ?_OPENMP */
#endif /* ?DZNRMX_LOOP */

#endif /* !DZNRMX_H */
