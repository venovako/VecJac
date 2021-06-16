#include "pjs.h"

int main(int argc, char *argv[])
{
  if (argc != 3) {
    (void)fprintf(stderr, "%s ID N\n", *argv);
    return EXIT_FAILURE;
  }

  const long id = atol(argv[1u]);
  const unsigned n = (unsigned)atol(argv[2u]);

  const char *fmt = (const char*)NULL;
  if (n <= 10u)
    fmt = "{%1u,%1u}";
  else if (n <= 100u)
    fmt = "{%2u,%2u}";
  else if (n <= 1000u)
    fmt = "{%3u,%3u}";
  else if (n <= 10000u)
    fmt = "{%4u,%4u}";
  else if (n <= 100000u)
    fmt = "{%5u,%5u}";
  else // n > 100000
    fmt = "{%u,%u}";

  unsigned stp = 0u;
  unsigned *const st = pjs(id, n, &stp);

  if (st) {
    for (unsigned s = 0u; s < stp; ++s) {
      unsigned *const r = st + n * (size_t)s;
      for (unsigned p = 0u; p < n; p += 2u) {
        if (p)
          (void)printf(",");
        (void)printf(fmt, r[p], r[p + 1u]);
      }
      (void)printf("\n");
    }
  }
  else
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
