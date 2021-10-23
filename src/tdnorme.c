#include "wdp.h"
#include "rnd.h"
#include "timer.h"

#ifdef TDNORME_SEQSRT
static int dcmp(const double x[static 1], const double y[static 1])
{
  return ((*x < *y) ? -1 : ((*y < *x) ? 1 : 0));
}
#else /* !TDNORME_SEQSRT */
#include "psort.h"
#endif /* ?TDNORME_SEQSRT */

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
  char s[26u] = { '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0' };
  (void)fprintf(stdout, "\"IT\",\"DPs\",\"NFsd\",\"N2sd\",\"XPsd\",\"NEsd\",\"WDP\",\"DPre\",\"NFre\",\"N2re\",\"XPre\",\"NEre\"\n");
  (void)fflush(stdout);

  for (size_t it = 0u; it < nit; ++it) {
    (void)fprintf(stdout, "%3zu,", it);
    (void)fflush(stdout);

    gendfrand_(&n, &aub, x);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,x)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; i += VDL) {
      double *const xi = x + i;
      _mm512_store_pd(xi, _mm512_abs_pd(_mm512_load_pd(xi)));
    }

    b = rdtsc_beg(rd);
    const double nf = dnf(n, x);
    e = rdtsc_end(rd);
    const double nfl = (double)tsc_lap(hz, b, e);

    b = rdtsc_beg(rd);
    const double dp = ddp(n, x);
    e = rdtsc_end(rd);
    const double dpl = (double)tsc_lap(hz, b, e);
    (void)fprintf(stdout, "%s,", dtoa(s, dpl));
    (void)fprintf(stdout, "%s,", dtoa(s, (nfl / dpl)));
    (void)fflush(stdout);

    b = rdtsc_beg(rd);
    const double n2 = dn2(n, x);
    e = rdtsc_end(rd);
    const double n2l = (double)tsc_lap(hz, b, e);
    (void)fprintf(stdout, "%s,", dtoa(s, (n2l / dpl)));
    (void)fflush(stdout);

    b = rdtsc_beg(rd);
    const double xp = xdp(n, x);
    e = rdtsc_end(rd);
    const double xpl = (double)tsc_lap(hz, b, e);
    (void)fprintf(stdout, "%s,", dtoa(s, (xpl / dpl)));
    (void)fflush(stdout);

    b = rdtsc_beg(rd);
    const double ne = dne(n, x);
    e = rdtsc_end(rd);
    const double nel = (double)tsc_lap(hz, b, e);
    (void)fprintf(stdout, "%s,", dtoa(s, (nel / dpl)));
    (void)fflush(stdout);

#ifdef TDNORME_SEQSRT
    qsort(x, n, sizeof(double), (int (*)(const void*, const void*))dcmp);
#else /* !TDNORME_SEQSRT */
    dpsort(n, x);
#endif /* ?TDNORME_SEQSRT */

    const double sq = wsq(n, x);
    (void)fprintf(stdout, "%s,", dtoa(s, sq));
    (void)fprintf(stdout, "%s,", dtoa(s, dre(dp, sq)));
    (void)fprintf(stdout, "%s,", dtoa(s, dre(nf, sq)));
    (void)fprintf(stdout, "%s,", dtoa(s, dre(n2, sq)));
    (void)fprintf(stdout, "%s,", dtoa(s, dre(xp, sq)));
    (void)fprintf(stdout, "%s\n", dtoa(s, dre(ne, sq)));
    (void)fflush(stdout);
  }

  free(x);
  return EXIT_SUCCESS;
}
