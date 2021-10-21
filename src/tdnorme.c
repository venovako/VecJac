#include "wdp.h"
#include "rnd.h"
#include "timer.h"

static int dcmp(const double x[static 1], const double y[static 1])
{
  return ((*x < *y) ? -1 : ((*y < *x) ? 1 : 0));
}

int main(int argc, char *argv[])
{
  if (4 != argc) {
    (void)fprintf(stderr, "%s 2^{EXP} {AUB} {nIT}\n", *argv);
    return EXIT_FAILURE;
  }

  const size_t n = (((size_t)1u) << atoz(argv[1u]));
  if (n < VDL) {
    (void)fprintf(stderr, "%zu = 2^{EXP} < %u\n", n, VDL);
    return EXIT_FAILURE;
  }
  const double aub = atof(argv[2u]);
  if (aub <= 0.0) {
    (void)fprintf(stderr, "{AUB} <= 0\n");
    return EXIT_FAILURE;
  }
  const size_t nit = atoz(argv[3u]);
  if (!nit) {
    (void)fprintf(stderr, "{nIT} == 0\n");
    return EXIT_FAILURE;
  }

  double *x = (double*)NULL;
  if (posix_memalign((void**)&x, (VDL * sizeof(double)), (n * sizeof(double))))
    return EXIT_FAILURE;

  unsigned rd[2u] = { 0u, 0u };
  const uint64_t hz = tsc_get_freq_hz_(rd);
#ifndef NDEBUG
  (void)fprintf(stderr, "TSC frequency: %lu Hz\n", hz);
  (void)fflush(stderr);
#endif /* !NDEBUG */
  uint64_t b = 0u, e = 0u;
  char s[31u] = { '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0' };
  (void)fprintf(stdout, "\"IT\",\"DPs\",\"N2sd\",\"NEsd\",\"WDP\",\"DPre\",\"N2re\",\"NEre\"\n");
  (void)fflush(stdout);

  for (size_t it = 0u; it < nit; ++it) {
    (void)fprintf(stdout, "%3zu,", it);
    (void)fflush(stdout);

    gendfrand_(&n, &aub, x);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,x)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i) {
      x[i] = fabs(x[i]);
    }

    b = rdtsc_beg(rd);
    const wide dp = wddp(n, x);
    e = rdtsc_end(rd);
    const long double dpl = tsc_lap(hz, b, e);
    (void)fprintf(stdout, "%s,", xtoa(s, dpl));
    (void)fflush(stdout);

    b = rdtsc_beg(rd);
    const wide n2 = wdn2(n, x);
    e = rdtsc_end(rd);
    const long double n2l = tsc_lap(hz, b, e);
    (void)fprintf(stdout, "%s,", xtoa(s, (n2l / dpl)));
    (void)fflush(stdout);

    b = rdtsc_beg(rd);
    const wide ne = wdne(n, x);
    e = rdtsc_end(rd);
    const long double nel = tsc_lap(hz, b, e);
    (void)fprintf(stdout, "%s,", xtoa(s, (nel / dpl)));
    (void)fflush(stdout);

    qsort(x, n, sizeof(double), (int (*)(const void*, const void*))dcmp);
    const wide sq = wdsq(n, x);
    (void)fprintf(stdout, "%s,", xtoa(s, (long double)sq));
    (void)fprintf(stdout, "%s,", xtoa(s, qdnre(dp, sq)));
    (void)fprintf(stdout, "%s,", xtoa(s, qdnre(n2, sq)));
    (void)fprintf(stdout, "%s\n", xtoa(s, qdnre(ne, sq)));
    (void)fflush(stdout);
  }

  free(x);
  return EXIT_SUCCESS;
}
