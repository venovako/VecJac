#ifndef SCNRMX_H
#define SCNRMX_H

#include "vecdef.h"

#ifdef SCNRMX_LOOP
#error SCNRMX_LOOP already defined
#else /* !SCNRMX_LOOP */
#ifdef _OPENMP
#define SCNRMX_LOOP(A,ldA)                                                      \
  for (fnat j = 0u; j < *n; ++j) {                                              \
    register const VS _zerof = _mm512_set1_ps(-0.0f);                           \
    register const VS inff = _mm512_set1_ps(HUGE_VALF);                         \
    register VS x = _mm512_setzero_ps();                                        \
    const float *const Aj = (A) + j * (size_t)(*(ldA));                         \
    for (fnat i = 0u; i < *m; i += VSL)                                         \
      x = _mm512_max_ps(_mm512_min_ps(VSABS(_mm512_load_ps(Aj + i)), inff), x); \
    y = fmaxf(_mm512_reduce_max_ps(x), y);                                      \
  }
#else /* !_OPENMP */
#define SCNRMX_LOOP(A,ldA)                                                              \
  for (fnat j = 0u; j < *n; ++j) {                                                      \
    const float *const Aj = (A) + j * (size_t)(*(ldA));                                 \
    for (fnat i = 0u; i < *m; i += VDL) {                                               \
      x = _mm512_max_ps(_mm512_min_ps(VSABS(_mm512_load_ps(Aj + i)), inff), x); VSP(x); \
    }                                                                                   \
  }
#endif /* ?_OPENMP */
#endif /* ?SCNRMX_LOOP */

#endif /* !SCNRMX_H */
