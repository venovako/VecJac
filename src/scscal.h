#ifndef SCSCAL_H
#define SCSCAL_H

#ifdef SCSCAL_LOOP
#error SCSCAL_LOOP already defined
#else /* !SCSCAL_LOOP */
#ifdef _OPENMP
#define SCSCAL_LOOP(A,ldA)                                            \
  for (fnat j = 0u; j < *n; ++j) {                                    \
    register const VS s = _mm512_set1_ps(e_);                         \
    float *const Aj = (A) + j * (size_t)(*(ldA));                     \
    for (fnat i = 0u; i < *m; i += VSL) {                             \
      float *const Aij = Aj + i;                                      \
      _mm512_store_ps(Aij, _mm512_scalef_ps(_mm512_load_ps(Aij), s)); \
    }                                                                 \
    t = imax(t, omp_get_thread_num());                                \
  }
#else /* !_OPENMP */
#define SCSCAL_LOOP(A,ldA)                                            \
  for (fnat j = 0u; j < *n; ++j) {                                    \
    float *const Aj = (A) + j * (size_t)(*(ldA));                     \
    for (fnat i = 0u; i < *m; i += VSL) {                             \
      float *const Aij = Aj + i;                                      \
      _mm512_store_ps(Aij, _mm512_scalef_ps(_mm512_load_ps(Aij), s)); \
    }                                                                 \
  }
#endif /* ?_OPENMP */
#endif /* ?SCSCAL_LOOP */

#endif /* !SCSCAL_H */
