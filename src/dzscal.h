#ifndef DZSCAL_H
#define DZSCAL_H

#ifdef DZSCAL_LOOP
#error DZSCAL_LOOP already defined
#else /* !DZSCAL_LOOP */
#ifdef _OPENMP
#define DZSCAL_LOOP(A,ldA)                                            \
  for (fnat j = 0u; j < *n; ++j) {                                    \
    register const VD s = _mm512_set1_pd(*e);                         \
    double *const Aj = (A) + j * (*(ldA));                            \
    for (fnat i = 0u; i < *m; i += VDL) {                             \
      double *const Aij = Aj + i;                                     \
      _mm512_store_pd(Aij, _mm512_scalef_pd(_mm512_load_pd(Aij), s)); \
    }                                                                 \
    t = imax(t, omp_get_thread_num());                                \
  }
#else /* !_OPENMP */
#define DZSCAL_LOOP(A,ldA)                                            \
  for (fnat j = 0u; j < *n; ++j) {                                    \
    double *const Aj = (A) + j * (*(ldA));                            \
    for (fnat i = 0u; i < *m; i += VDL) {                             \
      double *const Aij = Aj + i;                                     \
      _mm512_store_pd(Aij, _mm512_scalef_pd(_mm512_load_pd(Aij), s)); \
    }                                                                 \
  }
#endif /* ?_OPENMP */
#endif /* ?DZSCAL_LOOP */

#ifdef EDOMEGA
#error EDOMEGA already defined
#else /* !EDOMEGA */
#define EDOMEGA 1024
#endif /* ?EDOMEGA */

static inline fint s_k(const double M_k, const fint l)
{
  int e_k
#ifndef NDEBUG
    = 0
#endif /* !NDEBUG */
    ;
  (void)frexp(M_k, &e_k);
  return (EDOMEGA - l - e_k);
}

#endif /* !DZSCAL_H */
