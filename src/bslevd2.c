#include "laev2.h"
#include "wnrme.h"
#include "rnd.h"
#include "timer.h"

int main(int argc, char *argv[])
{
  (void)set_cbwr();

  if (4 != argc) {
    (void)fprintf(stderr, "%s filename 2^{batch_size} #batches\n", *argv);
    return EXIT_FAILURE;
  }

  const size_t n = ((size_t)1u << atoz(argv[2u]));
  int th = 0;
#ifdef _OPENMP
  th = omp_get_max_threads();
  if (n % th) {
    (void)fprintf(stderr, "batch_size has to be a multiple of %d.\n", th);
    return EXIT_FAILURE;
  }
#endif /* _OPENMP */
  const size_t b = atoz(argv[3u]);
  if (!b)
    return EXIT_SUCCESS;
  const size_t
    nl = strlen(argv[1u]),
    nl1 = (nl + 1u);
  char *const fn = calloc((nl + 3u), sizeof(char));
  assert(fn);
  strcpy(fn, argv[1u])[nl] = '.';
  int fm = O_RDONLY;
#ifdef _LARGEFILE64_SOURCE
  fm |= O_LARGEFILE;
#endif /* _LARGEFILE64_SOURCE */

  fn[nl1] = 'k';
  const int fk = open(fn, fm);
  if (-1 >= fk) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }
  fn[nl1] = 'l';
  const int fl = open(fn, fm);
  if (-1 >= fl) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }
  fn[nl1] = 'f';
  const int ff = open(fn, fm);
  if (-1 >= ff) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }
  fn[nl1] = 'g';
  const int fg = open(fn, fm);
  if (-1 >= fg) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }
  fn[nl1] = 'h';
  const int fh = open(fn, fm);
  if (-1 >= fh) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }

  const size_t nt = n * sizeof(float);
  float
    *const a11 = (float*)aligned_alloc(sizeof(float), nt),
    *const a22 = (float*)aligned_alloc(sizeof(float), nt),
    *const a21 = (float*)aligned_alloc(sizeof(float), nt),
    *const cs1 = (float*)aligned_alloc(sizeof(float), nt),
    *const sn1 = (float*)aligned_alloc(sizeof(float), nt),
    *const l1 = (float*)aligned_alloc(sizeof(float), nt),
    *const l2 = (float*)aligned_alloc(sizeof(float), nt);
  assert(a11);
  assert(a22);
  assert(a21);
  assert(cs1);
  assert(sn1);
  assert(l1);
  assert(l2);

  unsigned rd[2u] = { 0u, 0u };
  const uint64_t hz = tsc_get_freq_hz_(rd);
  (void)fprintf(stderr, "TSC frequency: %llu+(%u/%u) Hz.\n", (unsigned long long)hz, rd[0u], rd[1u]);
  (void)fflush(stderr);

  (void)fprintf(stdout, "\"B\",\"Ts\",\"ORT\",\"REN\",\"RLN\",\"RLX\",\"RLM\"\n");
  (void)fflush(stdout);
  const char *bf = (const char*)NULL;
  if (b <= 10u)
    bf = "%1zu";
  else if (b <= 100u)
    bf = "%2zu";
  else if (b <= 1000u)
    bf = "%3zu";
  else // b > 1000
    bf = "%zu";
  const size_t n_t = n / imax(th, 1);
  const size_t cnt = n_t * sizeof(float);
  char a[31u] = { '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0' };

  for (size_t j = 0u; j < b; ++j) {
    (void)fprintf(stdout, bf, j);
    (void)fflush(stdout);
    const size_t jn = j * n;
#ifdef _OPENMP
#pragma omp parallel default(none) shared(ff,fg,fh,a11,a22,a21,n,n_t,cnt,jn)
#endif /* _OPENMP */
    {
      const int mt =
#ifdef _OPENMP
        omp_get_thread_num()
#else /* !_OPENMP */
        0
#endif /* ?_OPENMP */
        ;
      const size_t tnt = mt * n_t;
      const off_t off = (jn + tnt) * sizeof(float);
      if ((ssize_t)cnt != pread(ff, (a11 + tnt), cnt, off))
        exit(EXIT_FAILURE);
      if ((ssize_t)cnt != pread(fg, (a22 + tnt), cnt, off))
        exit(EXIT_FAILURE);
      if ((ssize_t)cnt != pread(fh, (a21 + tnt), cnt, off))
        exit(EXIT_FAILURE);
    }
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);

    uint64_t be[2u] = { UINT64_C(0), UINT64_C(0) };
    be[0u] = rdtsc_beg(rd);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,a11,a22,a21,l1,l2,cs1,sn1)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i)
      LAPACK_S(laev2)((a11 + i), (a22 + i), (a21 + i), (l1 + i), (l2 + i), (cs1 + i), (sn1 + i));
    be[1u] = rdtsc_end(rd);
    (void)fprintf(stdout, "%15.9Lf,", tsc_lap(hz, be[0u], be[1u]));
    (void)fflush(stdout);
    wide o = W_ZERO, r = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,a11,a22,a21,cs1,sn1,l1,l2) reduction(max:o,r)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i) {
      wide AE = W_ZERO, AN = W_ZERO;
      o = fmaxw(o, worr(cs1[i], sn1[i]));
      r = fmaxw(r, wrer(a11[i], a22[i], a21[i], cs1[i], sn1[i], l1[i], l2[i], &AE, &AN));
    }
    (void)fprintf(stdout, "%s,", xtoa(a, (long double)o));
    (void)fprintf(stdout, "%s", xtoa(a, (long double)r));
    (void)fflush(stdout);
#ifdef _OPENMP
#pragma omp parallel default(none) shared(fk,fl,cs1,sn1,n,n_t,cnt,jn)
#endif /* _OPENMP */
    {
      const int mt =
#ifdef _OPENMP
        omp_get_thread_num()
#else /* !_OPENMP */
        0
#endif /* ?_OPENMP */
        ;
      const size_t tnt = mt * n_t;
      const off_t off = (jn + tnt) * sizeof(float);
      if ((ssize_t)cnt != pread(fk, (cs1 + tnt), cnt, off))
        exit(EXIT_FAILURE);
      if ((ssize_t)cnt != pread(fl, (sn1 + tnt), cnt, off))
        exit(EXIT_FAILURE);
    }
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);
    wide x = W_ZERO, m = W_ZERO;
    r = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,l1,l2,cs1,sn1) reduction(max:r,x,m)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i) {
      wide AE = W_ZERO, AN = W_ZERO;
      const wide RE = wlam(l1[i], l2[i], cs1[i], sn1[i], &AE, &AN);
      r = fmaxw(r, RE);
      x = fmaxw(x, AE);
      m = fmaxw(m, AN);
    }
    (void)fprintf(stdout, "%s,", xtoa(a, (long double)r));
    (void)fprintf(stdout, "%s,", xtoa(a, (long double)x));
    (void)fprintf(stdout, "%s\n", xtoa(a, (long double)m));
    (void)fflush(stdout);
  }

  (void)close(fh);
  (void)close(fg);
  (void)close(ff);
  (void)close(fl);
  (void)close(fk);

  free(l2);
  free(l1);
  free(sn1);
  free(cs1);
  free(a21);
  free(a22);
  free(a11);

  free(fn);
  return EXIT_SUCCESS;
}
