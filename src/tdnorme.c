#include "wdp.h"
#include "rnd.h"
#include "timer.h"

#ifdef TDNORME_PARSRT
#include "psort.h"
#else /* !TDNORME_PARSRT */
static int dcmp(const double x[static 1], const double y[static 1])
{
  return ((*x < *y) ? -1 : ((*y < *x) ? 1 : 0));
}
#endif /* ?TDNORME_PARSRT */

int main(int argc, char *argv[])
{
  (void)fprintf(stderr, "CBWR: %d\n", set_cbwr());
  (void)fflush(stderr);

  if (4 != argc) {
    (void)fprintf(stderr, "%s 2^{EXP} 2^{AUB} {nIT}\n", *argv);
    return EXIT_FAILURE;
  }

  const size_t n = (((size_t)1u) << atoz(argv[1u]));
  if (n < VDL) {
    (void)fprintf(stderr, "%zu = 2^{EXP} < %u\n", n, VDL);
    return EXIT_FAILURE;
  }
  const double aub = scalbn(1.0, atoi(argv[2u]));
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
  (void)fprintf(stderr, "TSC frequency: %llu+(%u/%u) Hz.\n", (unsigned long long)hz, rd[0u], rd[1u]);
  (void)fflush(stderr);
#endif /* !NDEBUG */
  uint64_t b = 0u, e = 0u;
  char s[26u] = { '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0' };
  (void)fprintf(stdout, "\"IT\",\"NBs\",\"DPs\",\"N2s\",\"NFs\",\"NEs\",\"NSs\",\"NCs\",\"NDs\",\"WDP\",\"NBre\",\"DPre\",\"N2re\",\"NFre\",\"NEre\",\"NSre\",\"NCre\",\"NDre\"\n");
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

    // call the reference BLAS for warmup
    b = rdtsc_beg(rd);
    const double nb = dnb(n, x);
    e = rdtsc_end(rd);
    (void)fprintf(stdout, "%s,", dtoa(s, (double)tsc_lap(hz, b, e)));
    (void)fflush(stdout);

    // call the MKL: ddot
    b = rdtsc_beg(rd);
    const double dp = ddp(n, x);
    e = rdtsc_end(rd);
    (void)fprintf(stdout, "%s,", dtoa(s, (double)tsc_lap(hz, b, e)));
    (void)fflush(stdout);

    // call the MKL: dnrm2
    b = rdtsc_beg(rd);
    const double n2 = dn2(n, x);
    e = rdtsc_end(rd);
    (void)fprintf(stdout, "%s,", dtoa(s, (double)tsc_lap(hz, b, e)));
    (void)fflush(stdout);

    // use long double: dnormf
    b = rdtsc_beg(rd);
    const double nf = dnf(n, x);
    e = rdtsc_end(rd);
    (void)fprintf(stdout, "%s,", dtoa(s, (double)tsc_lap(hz, b, e)));
    (void)fflush(stdout);

    // dnorme
    b = rdtsc_beg(rd);
    const double ne = dne(n, x);
    e = rdtsc_end(rd);
    (void)fprintf(stdout, "%s,", dtoa(s, (double)tsc_lap(hz, b, e)));
    (void)fflush(stdout);

    // dnorms
    b = rdtsc_beg(rd);
    const double ns = dns(n, x);
    e = rdtsc_end(rd);
    (void)fprintf(stdout, "%s,", dtoa(s, (double)tsc_lap(hz, b, e)));
    (void)fflush(stdout);

    // simple overflowing dnorme
    b = rdtsc_beg(rd);
    const double nc = dnc(n, x);
    e = rdtsc_end(rd);
    (void)fprintf(stdout, "%s,", dtoa(s, (double)tsc_lap(hz, b, e)));
    (void)fflush(stdout);

    // simple overflowing dnorms
    b = rdtsc_beg(rd);
    const double nd = dnd(n, x);
    e = rdtsc_end(rd);
    (void)fprintf(stdout, "%s,", dtoa(s, (double)tsc_lap(hz, b, e)));
    (void)fflush(stdout);

#ifdef TDNORME_PARSRT
    dpsort(n, x);
#else /* !TDNORME_PARSRT */
    qsort(x, n, sizeof(double), (int (*)(const void*, const void*))dcmp);
#endif /* ?TDNORME_PARSRT */

    // use wide precision
    const double sq = wsq(n, x);
    (void)fprintf(stdout, "%s,", dtoa(s, sq));

    (void)fprintf(stdout, "%s,", dtoa(s, dre(nb, sq)));
    (void)fprintf(stdout, "%s,", dtoa(s, dre(dp, sq)));
    (void)fprintf(stdout, "%s,", dtoa(s, dre(n2, sq)));
    (void)fprintf(stdout, "%s,", dtoa(s, dre(nf, sq)));
    (void)fprintf(stdout, "%s,", dtoa(s, dre(ne, sq)));
    (void)fprintf(stdout, "%s,", dtoa(s, dre(ns, sq)));
    (void)fprintf(stdout, "%s,", dtoa(s, dre(nc, sq)));
    (void)fprintf(stdout, "%s\n", dtoa(s, dre(nd, sq)));
    (void)fflush(stdout);
  }

  free(x);
  return EXIT_SUCCESS;
}
