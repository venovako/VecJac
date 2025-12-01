#include "pjs.h"

extern int pvn_cjs_init(void *const js, const int id, const int n);
extern int pvn_cjs_next(void *const js, int *const arr);
extern int pvn_cjs_free(void *const js);

unsigned *pjs(const long id, const unsigned n, unsigned stp[static restrict 1])
{
  long js[4] = { 0l, 0l, 0l, 0l };
  unsigned *st = (unsigned*)NULL;
  int *arr = (int*)NULL;

  *stp = 0u;
  if (n <= 1u)
    goto theEnd;

  switch (id) {
  case PJS_ME:
    *stp = n - 1u;
    break;
  case PJS_MM:
    *stp = n;
    break;
  default:
    goto theEnd;
  }

  if (pvn_cjs_init((void*)js, (int)id, (int)n) != (int)(*stp))
    goto theEnd;
  if (!(st = (unsigned*)malloc(n * (*stp * sizeof(unsigned)))))
    goto theEnd;
  if (!(arr = (int*)malloc(((id == PJS_MM) ? 2u : 1u) * (n * sizeof(int)))))
    goto theEnd;

  const int n_2 = (int)(n >> 1u);
  for (unsigned s = 0u; s < *stp; ++s) {
    if (pvn_cjs_next((void*)js, arr) != n_2)
      goto theEnd;
    unsigned *const r = st + n * (size_t)s;
    if (id == PJS_ME) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,arr)
#endif /* _OPENMP */
      for (unsigned i = 0u; i < n; ++i)
        r[i] = (unsigned)(arr[i]);
    }
    else {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,arr)
#endif /* _OPENMP */
      for (unsigned i = 0u; i < n; ++i)
        r[i] = (unsigned)(((const unsigned long*)arr)[i]);
    }
  }

  free(arr);
  (void)pvn_cjs_free((void*)js);
  return st;

 theEnd:
  free(arr);
  free(st);
  (void)pvn_cjs_free((void*)js);
  return (unsigned*)NULL;
}
