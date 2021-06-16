#include "pjs.h"

#include "jstrat.h"

unsigned *pjs(const long id, const unsigned n, unsigned stp[static restrict 1])
{
  jstrat_common js;
  (void)memset(&js, 0, sizeof(js));
  unsigned *st = (unsigned*)NULL;
  integer *arr = (integer*)NULL;

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

  if (jstrat_init(&js, id, (integer)n) != (integer)(*stp))
    goto theEnd;
  if (!(st = (unsigned*)malloc(n * (*stp * sizeof(unsigned)))))
    goto theEnd;
  if (!(arr = (integer*)malloc(n * sizeof(integer))))
    goto theEnd;

  const long n_2 = (long)(n >> 1u);
  for (unsigned s = 0u; s < *stp; ++s) {
    if (labs(jstrat_next(&js, arr)) != n_2)
      goto theEnd;
    unsigned *const r = st + n * (size_t)s;
    for (unsigned i = 0u; i < n; ++i)
      r[i] = (unsigned)(arr[i]);
  }

  free(arr);
  jstrat_free(&js);
  return st;

 theEnd:
  free(arr);
  free(st);
  jstrat_free(&js);
  return (unsigned*)NULL;
}
