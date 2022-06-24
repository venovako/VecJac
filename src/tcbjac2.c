#include "cbjac2.h"
#include "wnrme.h"
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
  char *const fn = (char*)calloc((nl + 3u), sizeof(char));
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
    *const c = (float*)aligned_alloc(VA, nt),
    *const cat = (float*)aligned_alloc(VA, nt),
    *const sat = (float*)aligned_alloc(VA, nt),
    *const l1 = (float*)aligned_alloc(VA, nt),
    *const l2 = (float*)aligned_alloc(VA, nt);
  assert(a11);
  assert(a22);
  assert(a21r);
  assert(a21i);
  assert(c);
  assert(cat);
  assert(sat);
  assert(l1);
  assert(l2);

  wide *const w = (wide*)aligned_alloc(sizeof(wide), (n * sizeof(wide)));
  assert(w);

  unsigned *const p = (unsigned*)malloc((n >> VSLlg) * sizeof(unsigned));
  assert(p);

  unsigned rd[2u] = { 0u, 0u };
  uint64_t hz = tsc_get_freq_hz_(rd), be[2u] = { UINT64_C(0), UINT64_C(0) };
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
#pragma omp parallel default(none) shared(ff,fg,fr,fj,a11,a22,a21r,a21i,n,n_t,cnt,jn)
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
      if ((ssize_t)cnt != pread(fr, (a21r + tnt), cnt, off))
        exit(EXIT_FAILURE);
      if ((ssize_t)cnt != pread(fj, (a21i + tnt), cnt, off))
        exit(EXIT_FAILURE);
    }
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);
    be[0u] = rdtsc_beg(rd);
    const fint _n = -(fint)n;
    (void)cbjac2_(&_n, a11, a22, a21r, a21i, c, cat, sat, l1, l2, p);
    be[1u] = rdtsc_end(rd);
    (void)fprintf(stdout, "%15.9Lf,", tsc_lap(hz, be[0u], be[1u]));
    (void)fflush(stdout);
    wide o = W_ZERO, r = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,a11,a22,a21r,a21i,c,cat,sat,l1,l2,w) reduction(max:o,r)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i) {
      const wide CS = (wide)(c[i]);
      const wide SNR = (wide)(cat[i]);
      const wide SNI = (wide)(sat[i]);
      wide AE = W_ZERO, AN = W_ZERO;
      o = fmaxw(o, (w[i] = worc(CS, SNR, SNI)));
      r = fmaxw(r, wrec(a11[i], a22[i], a21r[i], a21i[i], CS, SNR, SNI, l1[i], l2[i], &AE, &AN));
    }
    (void)fprintf(stdout, "%s,", xtoa(a, (long double)o));
    (void)fprintf(stdout, "%s", xtoa(a, (long double)r));
    (void)fflush(stdout);
    size_t ix = n;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,o,w) reduction(min:ix)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i)
      if (w[i] == o)
        ix = i;
    (void)fprintf(stderr, "%zu,%zu,%s;", j, ix, xtoa(a, (long double)o));
    (void)fprintf(stderr, "%s,", xtoa(a, a11[ix]));
    (void)fprintf(stderr, "%s,", xtoa(a, a22[ix]));
    (void)fprintf(stderr, "(%s,", xtoa(a, a21r[ix]));
    (void)fprintf(stderr, "%s);", xtoa(a, a21i[ix]));
    (void)fprintf(stderr, "%s,", xtoa(a, c[ix]));
    (void)fprintf(stderr, "(%s,", xtoa(a, cat[ix]));
    (void)fprintf(stderr, "%s);", xtoa(a, sat[ix]));
    (void)fprintf(stderr, "%s,", xtoa(a, l1[ix]));
    (void)fprintf(stderr, "%s\n", xtoa(a, l2[ix]));
    (void)fflush(stderr);
#ifdef _OPENMP
#pragma omp parallel default(none) shared(fk,fl,cat,sat,n,n_t,cnt,jn)
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
      if ((ssize_t)cnt != pread(fk, (cat + tnt), cnt, off))
        exit(EXIT_FAILURE);
      if ((ssize_t)cnt != pread(fl, (sat + tnt), cnt, off))
        exit(EXIT_FAILURE);
    }
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);
    wide x = W_ZERO, m = W_ZERO;
    r = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,l1,l2,cat,sat) reduction(max:r,x,m)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i) {
      wide AE = W_ZERO, AN = W_ZERO;
      const wide RE = wlam(l1[i], l2[i], cat[i], sat[i], &AE, &AN);
      r = fmaxw(r, RE);
      x = fmaxw(x, AE);
      m = fmaxw(m, AN);
    }
    (void)fprintf(stdout, "%s,", xtoa(a, (long double)r));
    (void)fprintf(stdout, "%s,", xtoa(a, (long double)x));
    (void)fprintf(stdout, "%s\n", xtoa(a, (long double)m));
    (void)fflush(stdout);
  }

  (void)close(fj);
  (void)close(fr);
  (void)close(fg);
  (void)close(ff);
  (void)close(fl);
  (void)close(fk);

  free(p);
  free(w);
  free(l2);
  free(l1);
  free(sat);
  free(cat);
  free(c);
  free(a21i);
  free(a21r);
  free(a22);
  free(a11);

  free(fn);
  return EXIT_SUCCESS;
}
