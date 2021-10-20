#include "zjac2.h"
#include "dzjac2.h"
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
  if (n % VDL) {
    (void)fprintf(stderr, "batch_size has to be a multiple of %u.\n", VDL);
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
  const int fm = O_RDONLY | O_LARGEFILE;

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
    *const a11 = (double*)aligned_alloc(VA, nt),
    *const a22 = (double*)aligned_alloc(VA, nt),
    *const a21r = (double*)aligned_alloc(VA, nt),
    *const a21i = (double*)aligned_alloc(VA, nt),
    *const s = (double*)aligned_alloc(VA, nt),
    *const t = (double*)aligned_alloc(VA, nt),
    *const c = (double*)aligned_alloc(VA, nt),
    *const ca = (double*)aligned_alloc(VA, nt),
    *const sa = (double*)aligned_alloc(VA, nt),
    *const l1 = (double*)aligned_alloc(VA, nt),
    *const l2 = (double*)aligned_alloc(VA, nt);
  assert(a11);
  assert(a22);
  assert(a21r);
  assert(a21i);
  assert(s);
  assert(t);
  assert(c);
  assert(ca);
  assert(sa);
  assert(l1);
  assert(l2);

  unsigned *const p = (unsigned*)malloc((n >> VDLlg) * sizeof(unsigned));
  assert(p);

  const size_t nw = n * sizeof(wide);
  wide
    *const RE = (wide*)aligned_alloc(sizeof(wide), nw),
    *const AE = (wide*)aligned_alloc(sizeof(wide), nw),
    *const AN = (wide*)aligned_alloc(sizeof(wide), nw),
    *const L1 = (wide*)aligned_alloc(sizeof(wide), nw),
    *const L2 = (wide*)aligned_alloc(sizeof(wide), nw);
  assert(RE);
  assert(AE);
  assert(AN);
  assert(L1);
  assert(L2);

  unsigned rd[2u] = { 0u, 0u };
  uint64_t hz = tsc_get_freq_hz_(rd), be[2u] = { UINT64_C(0), UINT64_C(0) };
  (void)fprintf(stderr, "TSC frequency: %lu+(%u/%u) Hz.\n", hz, rd[0u], rd[1u]);

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
    th = imax(th, zjac2_((const fnat*)&n, a11, a22, a21r, a21i, s, t, c, ca, sa, l1, l2, p));
    be[1u] = rdtsc_end(rd);
    (void)fprintf(stdout, "%15.9Lf", tsc_lap(hz, be[0u], be[1u]));
    (void)fflush(stdout);
    (void)dzjac2_pp((fnat)n, s, c, l1, l2, p, (double*)L1, (double*)L2);
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);
    wide r = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,a11,a22,a21r,a21i,s,t,c,ca,sa,l1,l2,AE,AN,L1,L2,RE) reduction(max:r)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i) {
      RE[i] = wrec(a11[i], a22[i], a21r[i], a21i[i], s[i], t[i], c[i], ca[i], sa[i], l1[i], l2[i], (AE + i), (AN + i), (L1 + i), (L2 + i));
      r = fmaxw(r, RE[i]);
    }
    (void)fprintf(stdout, "%30s", xtoa(a, (long double)r));
    (void)fflush(stdout);
#ifdef _OPENMP
#pragma omp parallel default(none) shared(fk,fl,l1,l2,n,n_t,cnt,jn)
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
      if ((ssize_t)cnt != pread(fk, (l1 + tnt), cnt, off))
        exit(EXIT_FAILURE);
      if ((ssize_t)cnt != pread(fl, (l2 + tnt), cnt, off))
        exit(EXIT_FAILURE);
    }
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);
    wide x = W_ZERO, m = W_ZERO;
    r = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,L1,L2,l1,l2,AE,AN,RE) reduction(max:r,x,m)
#endif /* _OPENMP */
    for (size_t i = 0u; i < n; ++i) {
      RE[i] = wlam(L1[i], L2[i], l1[i], l2[i], (AE + i), (AN + i));
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

  free(L2);
  free(L1);
  free(AN);
  free(AE);
  free(RE);

  free(p);

  free(l2);
  free(l1);
  free(sa);
  free(ca);
  free(c);
  free(t);
  free(s);
  free(a21i);
  free(a21r);
  free(a22);
  free(a11);

  free(fn);
  return EXIT_SUCCESS;
}
