#include "laev2.h"
#include "wnrme.h"
#include "rnd.h"
#include "timer.h"

int main(int argc, char *argv[])
{
  if (4 != argc) {
    (void)fprintf(stderr, "%s filename batch_size #batches\n", *argv);
    return EXIT_FAILURE;
  }

  const size_t n = atoz(argv[2u]);
  if (!n)
    return EXIT_SUCCESS;
  fint th = 0;
#ifdef _OPENMP
  th = omp_get_max_threads();
  if (n % th) {
    (void)fprintf(stderr, "batch_size has to be a multiple of %d.\n", (int)th);
    return EXIT_FAILURE;
  }
#endif /* _OPENMP */
  const size_t b = atoz(argv[3u]);
  if (!b)
    return EXIT_SUCCESS;
  if (b > 10000u) {
    (void)fprintf(stderr, "#batches is limited to at most 10000.\n");
    return EXIT_FAILURE;
  }
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
  fn[nl1] = 'r';
  const int fr = open(fn, fm);
  if (-1 >= fr) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }
  fn[nl1] = 'j';
  const int fj = open(fn, fm);
  if (-1 >= fj) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }

  const size_t nt = n * sizeof(double);
  double
    *const a11 = (double*)aligned_alloc(sizeof(double), nt),
    *const a22 = (double*)aligned_alloc(sizeof(double), nt),
    *const a21r = (double*)aligned_alloc(sizeof(double), nt),
    *const a21i = (double*)aligned_alloc(sizeof(double), nt),
    *const cs1 = (double*)aligned_alloc(sizeof(double), nt),
    *const snr = (double*)aligned_alloc(sizeof(double), nt),
    *const sni = (double*)aligned_alloc(sizeof(double), nt),
    *const l1 = (double*)aligned_alloc(sizeof(double), nt),
    *const l2 = (double*)aligned_alloc(sizeof(double), nt),
    *const t = (double*)aligned_alloc(sizeof(double), nt);
  assert(a11);
  assert(a22);
  assert(a21r);
  assert(a21i);
  assert(cs1);
  assert(snr);
  assert(sni);
  assert(l1);
  assert(l2);
  assert(t);

  const size_t nw = n * sizeof(wide);
  wide
    *const RE = (wide*)aligned_alloc(sizeof(wide), nw),
    *const AE = (wide*)aligned_alloc(sizeof(wide), nw),
    *const AN = (wide*)aligned_alloc(sizeof(wide), nw);
  assert(RE);
  assert(AE);
  assert(AN);

  unsigned rd[2u] = { 0u, 0u };
  uint64_t hz = tsc_get_freq_hz_(rd), be[2u] = { UINT64_C(0), UINT64_C(0) };
#ifndef NDEBUG
  (void)fprintf(stderr, "TSC frequency: %llu+(%u/%u) Hz.\n", (unsigned long long)hz, rd[0u], rd[1u]);
#endif /* !NDEBUG */

  (void)fprintf(stdout, "\"B\",\"Ts\",\"REN\",\"RLN\",\"RLX\",\"RLM\"\n");
  const char *bf = (const char*)NULL;
  if (b <= 10u)
    bf = "%zu";
  else if (b <= 100u)
    bf = "%2zu";
  else if (b <= 1000u)
    bf = "%3zu";
  else // b > 1000
    bf = "%4zu";
  const size_t n_t = n / imax(th, 1);
  const size_t cnt = n_t * sizeof(double);
  char a[31u] = { '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0' };
  th = 0;

  for (size_t j = 0u; j < b; ++j) {
    (void)fprintf(stdout, bf, j);
    (void)fflush(stdout);
    const size_t jn = j * n;
#ifdef _OPENMP
#pragma omp parallel default(none) shared(ff,fg,fr,fj,a11,a22,a21r,a21i,n,n_t,cnt,jn)
#endif /* _OPENMP */
    {
      const fint mt =
#ifdef _OPENMP
        omp_get_thread_num()
#else /* !_OPENMP */
        0
#endif /* ?_OPENMP */
        ;
      const size_t tnt = mt * n_t;
      const off_t off = (jn + tnt) * sizeof(double);
      if ((ssize_t)cnt != pread(ff, (a11 + tnt), cnt, off))
        exit(EXIT_FAILURE);
      if ((ssize_t)cnt != pread(fg, (a22 + tnt), cnt, off))
        exit(EXIT_FAILURE);
      if ((ssize_t)cnt != pread(fr, (a21r + tnt), cnt, off))
        exit(EXIT_FAILURE);
      if ((ssize_t)cnt != pread(fj, (a21i + tnt), cnt, off))
        exit(EXIT_FAILURE);
    }
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);
    be[0u] = rdtsc_beg(rd);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,a11,a22,a21r,a21i,l1,l2,cs1,snr,sni) reduction(max:th)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i) {
      zlevd2(a11[i], a22[i], a21r[i], a21i[i], (l1 + i), (l2 + i), (cs1 + i), (snr + i), (sni + i));
#ifdef _OPENMP
      th = imax(th, (omp_get_thread_num() + 1));
#endif /* _OPENMP */
    }
    be[1u] = rdtsc_end(rd);
    (void)fprintf(stdout, "%15.9Lf", tsc_lap(hz, be[0u], be[1u]));
    (void)fflush(stdout);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,cs1,snr,sni,t)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i)
      t[i] = zlevd2_pp(cs1[i], (snr + i), (sni + i));
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);
    wide r = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,a11,a22,a21r,a21i,t,cs1,snr,sni,l1,l2,AE,AN,RE) reduction(max:r)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i) {
      RE[i] = wrecf(a11[i], a22[i], a21r[i], a21i[i], t[i], cs1[i], snr[i], sni[i], l1[i], l2[i], (AE + i), (AN + i));
      r = fmaxw(r, RE[i]);
    }
    (void)fprintf(stdout, "%30s", xtoa(a, (long double)r));
    (void)fflush(stdout);
#ifdef _OPENMP
#pragma omp parallel default(none) shared(fk,fl,t,cs1,n,n_t,cnt,jn)
#endif /* _OPENMP */
    {
      const fint mt =
#ifdef _OPENMP
        omp_get_thread_num()
#else /* !_OPENMP */
        0
#endif /* ?_OPENMP */
        ;
      const size_t tnt = mt * n_t;
      const off_t off = (jn + tnt) * sizeof(double);
      if ((ssize_t)cnt != pread(fk, (t + tnt), cnt, off))
        exit(EXIT_FAILURE);
      if ((ssize_t)cnt != pread(fl, (cs1 + tnt), cnt, off))
        exit(EXIT_FAILURE);
    }
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);
    wide x = W_ZERO, m = W_ZERO;
    r = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,l1,l2,t,cs1,AE,AN,RE) reduction(max:r,x,m)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i) {
      RE[i] = wlam(l1[i], l2[i], t[i], cs1[i], (AE + i), (AN + i));
      r = fmaxw(r, RE[i]);
      x = fmaxw(x, AE[i]);
      m = fmaxw(m, AN[i]);
    }
    (void)fprintf(stdout, "%30s,", xtoa(a, (long double)r));
    (void)fprintf(stdout, "%30s,", xtoa(a, (long double)x));
    (void)fprintf(stdout, "%30s\n", xtoa(a, (long double)m));
    (void)fflush(stdout);
  }
  (void)fprintf(stderr, "max(#threads) = %u\n", (unsigned)th);

  (void)close(fj);
  (void)close(fr);
  (void)close(fg);
  (void)close(ff);
  (void)close(fl);
  (void)close(fk);

  free(AN);
  free(AE);
  free(RE);

  free(t);
  free(l2);
  free(l1);
  free(sni);
  free(snr);
  free(cs1);
  free(a21i);
  free(a21r);
  free(a22);
  free(a11);

  free(fn);
  return EXIT_SUCCESS;
}
