#include "cjac2f.h"
#include "wnrme.h"
#include "rnd.h"
#include "timer.h"

int main(int argc, char *argv[])
{
  if (4 != argc) {
    (void)fprintf(stderr, "%s filename 2^{batch_size} #batches\n", *argv);
    return EXIT_FAILURE;
  }

  const size_t n = ((size_t)1u << atoz(argv[2u]));
  if (n % VSL) {
    (void)fprintf(stderr, "batch_size has to be a multiple of %u.\n", VSL);
    return EXIT_FAILURE;
  }
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

  const size_t nt = n * sizeof(float);
  float
    *const a11 = (float*)aligned_alloc(VA, nt),
    *const a22 = (float*)aligned_alloc(VA, nt),
    *const a21r = (float*)aligned_alloc(VA, nt),
    *const a21i = (float*)aligned_alloc(VA, nt),
    *const t = (float*)aligned_alloc(VA, nt),
    *const c = (float*)aligned_alloc(VA, nt),
    *const ca = (float*)aligned_alloc(VA, nt),
    *const sa = (float*)aligned_alloc(VA, nt),
    *const l1 = (float*)aligned_alloc(VA, nt),
    *const l2 = (float*)aligned_alloc(VA, nt);
  assert(a11);
  assert(a22);
  assert(a21r);
  assert(a21i);
  assert(t);
  assert(c);
  assert(ca);
  assert(sa);
  assert(l1);
  assert(l2);

  unsigned *const p = (unsigned*)malloc((n >> VSLlg) * sizeof(unsigned));
  assert(p);

  unsigned rd[2u] = { 0u, 0u };
  uint64_t hz = tsc_get_freq_hz_(rd), be[2u] = { UINT64_C(0), UINT64_C(0) };
  (void)fprintf(stderr, "TSC frequency: %llu+(%u/%u) Hz.\n", (unsigned long long)hz, rd[0u], rd[1u]);
  (void)fflush(stderr);

  (void)fprintf(stdout, "\"B\",\"Ts\",\"REN\",\"RLN\",\"RLX\",\"RLM\"\n");
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
      const off_t off = (jn + tnt) * sizeof(float);
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
    th = imax(th, cjac2f_((const fnat*)&n, a11, a22, a21r, a21i, t, c, ca, sa, l1, l2, p));
    be[1u] = rdtsc_end(rd);
    (void)fprintf(stdout, "%15.9Lf,", tsc_lap(hz, be[0u], be[1u]));
    (void)fflush(stdout);
    wide r = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,a11,a22,a21r,a21i,t,c,ca,sa,l1,l2) reduction(max:r)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i) {
      wide AE = W_ZERO, AN = W_ZERO;
      const wide RE = wrecf(a11[i], a22[i], a21r[i], a21i[i], t[i], c[i], ca[i], sa[i], l1[i], l2[i], &AE, &AN);
      r = fmaxw(r, RE);
    }
    (void)fprintf(stdout, "%s", xtoa(a, (long double)r));
    (void)fflush(stdout);
#ifdef _OPENMP
#pragma omp parallel default(none) shared(fk,fl,t,c,n,n_t,cnt,jn)
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
      const off_t off = (jn + tnt) * sizeof(float);
      if ((ssize_t)cnt != pread(fk, (t + tnt), cnt, off))
        exit(EXIT_FAILURE);
      if ((ssize_t)cnt != pread(fl, (c + tnt), cnt, off))
        exit(EXIT_FAILURE);
    }
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);
    wide x = W_ZERO, m = W_ZERO;
    r = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,l1,l2,t,c) reduction(max:r,x,m)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i) {
      wide AE = W_ZERO, AN = W_ZERO;
      const wide RE = wlam(l1[i], l2[i], t[i], c[i], &AE, &AN);
      r = fmaxw(r, RE);
      x = fmaxw(x, AE);
      m = fmaxw(m, AN);
    }
    (void)fprintf(stdout, "%s,", xtoa(a, (long double)r));
    (void)fprintf(stdout, "%s,", xtoa(a, (long double)x));
    (void)fprintf(stdout, "%s\n", xtoa(a, (long double)m));
    (void)fflush(stdout);
  }
  (void)fprintf(stderr, "max(#threads) = %u\n", (unsigned)th);
  (void)fflush(stderr);

  (void)close(fj);
  (void)close(fr);
  (void)close(fg);
  (void)close(ff);
  (void)close(fl);
  (void)close(fk);

  free(p);

  free(l2);
  free(l1);
  free(sa);
  free(ca);
  free(c);
  free(t);
  free(a21i);
  free(a21r);
  free(a22);
  free(a11);

  free(fn);
  return EXIT_SUCCESS;
}
