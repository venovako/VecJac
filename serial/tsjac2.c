#include "serial.h"
#include "timer.h"
#include "wnrme.h"

int main(int argc, char *argv[])
{
  if (5 != argc) {
    (void)fprintf(stderr, "%s {C|c|D|d|S|s|Z|z} filename lg(batch_size) #batches\n", *argv);
    return EXIT_FAILURE;
  }

  int kind = toupper(*(argv[1u]));
  (void)fprintf(stderr, "%c ", (char)kind);
  switch ((char)kind)
  {
  case 'C':
    kind = (int)sizeof(float);
    break;
  case 'D':
    kind = -(int)sizeof(double);
    break;
  case 'S':
    kind = -(int)sizeof(float);
    break;
  case 'Z':
    kind = (int)sizeof(double);
    break;
  default:
    (void)fprintf(stderr, "the first argument is invalid.\n");
    return EXIT_FAILURE;
  }

  const size_t n = ((size_t)1u << atoz(argv[3u]));
  (void)fprintf(stderr, "%llu.\n", (unsigned long long)n);
  (void)fflush(stderr);

  const size_t b = atoz(argv[4u]);
  if (!b)
    return EXIT_SUCCESS;

  const size_t
    nl = strlen(argv[2u]),
    nl1 = (nl + 1u);
#ifdef _LARGEFILE64_SOURCE
#ifdef _GNU_SOURCE
  const int fm = O_RDONLY | O_LARGEFILE;
#else /* Darwin */
  const int fm = O_RDONLY;
#endif /* ?_GNU_SOURCE */
  char fn[nl + 3u];
#else /* !_LARGEFILE64_SOURCE */
  const int fm = _O_RDONLY | _O_BINARY;
  char *const fn = (char*)_alloca(nl + 3u);
#endif /* ?_LARGEFILE64_SOURCE */
  (void)strcpy((strcpy(fn, argv[2u]) + nl), ". ");

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
  const int fh = ((kind > 0) ? 0 : open(fn, fm));
  if (-1 >= fh) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }
  fn[nl1] = 'r';
  const int fr = ((kind > 0) ? open(fn, fm) : 0);
  if (-1 >= fr) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }
  fn[nl1] = 'j';
  const int fj = ((kind > 0) ? open(fn, fm) : 0);
  if (-1 >= fj) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }

  wide *const w = (wide*)malloc(n * sizeof(wide));
  if (!w)
    return EXIT_FAILURE;
  const size_t cnt = n * abs(kind);

  void
    *const ha11 = malloc(cnt),
    *const ha22 = malloc(cnt),
    *const ha21r = malloc(cnt),
    *const ha21i = ((kind > 0) ? malloc(cnt) : NULL),
    *const hc = malloc(cnt),
    *const hcat = malloc(cnt),
    *const hsat = ((kind > 0) ? malloc(cnt) : NULL),
    *const hl1 = malloc(cnt),
    *const hl2 = malloc(cnt);

  if (!ha11)
    return EXIT_FAILURE;
  if (!ha22)
    return EXIT_FAILURE;
  if (!ha21r)
    return EXIT_FAILURE;
  if ((kind > 0) && !ha21i)
    return EXIT_FAILURE;
  if (!hc)
    return EXIT_FAILURE;
  if (!hcat)
    return EXIT_FAILURE;
  if ((kind > 0) && !hsat)
    return EXIT_FAILURE;
  if (!hl1)
    return EXIT_FAILURE;
  if (!hl2)
    return EXIT_FAILURE;

  int aux = 0;
  uint64_t be[2u] = { UINT64_C(0), UINT64_C(0) };
  const uint64_t hz = tsc_get_freq_hz_(&aux);

  (void)fprintf(stdout, "\"B\",\"Ts\",\"ORT\",\"REN\",\"RLN\",\"RLX\",\"RLM\"\n");
  (void)fflush(stdout);
  const char *bf = (const char*)NULL;
  if (b <= 10u)
    bf = "%1llu";
  else if (b <= 100u)
    bf = "%2llu";
  else if (b <= 1000u)
    bf = "%3llu";
  else // b > 1000
    bf = "%llu";
  char a[26u] = { '\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0' };

  for (size_t j = 0u; j < b; ++j) {
    (void)fprintf(stdout, bf, j);
    (void)fflush(stdout);
    if (cnt != (size_t)read(ff, ha11, cnt))
      return EXIT_FAILURE;
    if (cnt != (size_t)read(fg, ha22, cnt))
      return EXIT_FAILURE;
    if (kind > 0) {
      if (cnt != (size_t)read(fr, ha21r, cnt))
        return EXIT_FAILURE;
      if (cnt != (size_t)read(fj, ha21i, cnt))
        return EXIT_FAILURE;
    }
    else {
      if (cnt != (size_t)read(fh, ha21r, cnt))
        return EXIT_FAILURE;
    }
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);
    be[0u] = rdtsc_beg(&aux);
    if (kind > 0) {
      if (kind == (int)sizeof(float)) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,ha11,ha22,ha21r,ha21i,hc,hcat,hsat,hl1,hl2) private(aux)
#endif /* _OPENMP */
        for (size_t i = 0u; i < n; ++i) {
          aux = csjac2(((const float*)ha11)[i], ((const float*)ha22)[i], ((const float*)ha21r)[i], ((const float*)ha21i)[i], ((float*)hc + i), ((float*)hcat + i), ((float*)hsat + i), ((float*)hl1 + i), ((float*)hl2 + i));
        }
      }
      else {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,ha11,ha22,ha21r,ha21i,hc,hcat,hsat,hl1,hl2) private(aux)
#endif /* _OPENMP */
        for (size_t i = 0u; i < n; ++i) {
          aux = zsjac2(((const double*)ha11)[i], ((const double*)ha22)[i], ((const double*)ha21r)[i], ((const double*)ha21i)[i], ((double*)hc + i), ((double*)hcat + i), ((double*)hsat + i), ((double*)hl1 + i), ((double*)hl2 + i));
        }
      }
    }
    else if (kind == -(int)sizeof(float)) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,ha11,ha22,ha21r,hc,hcat,hl1,hl2) private(aux)
#endif /* _OPENMP */
      for (size_t i = 0u; i < n; ++i) {
        aux = ssjac2(((const float*)ha11)[i], ((const float*)ha22)[i], ((const float*)ha21r)[i], ((float*)hc + i), ((float*)hcat + i), ((float*)hl1 + i), ((float*)hl2 + i));
      }
    }
    else {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,ha11,ha22,ha21r,hc,hcat,hl1,hl2) private(aux)
#endif /* _OPENMP */
      for (size_t i = 0u; i < n; ++i) {
        aux = dsjac2(((const double*)ha11)[i], ((const double*)ha22)[i], ((const double*)ha21r)[i], ((double*)hc + i), ((double*)hcat + i), ((double*)hl1 + i), ((double*)hl2 + i));
      }
    }
    be[1u] = rdtsc_end(&aux);
    (void)fprintf(stdout, "%15.9Lf,", tsc_lap(hz, be[0u], be[1u]));
    (void)fflush(stdout);
    wide o = W_ZERO, r = W_ZERO;
    if (kind > 0) {
      if (kind == (int)sizeof(float)) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,ha11,ha22,ha21r,ha21i,hc,hcat,hsat,hl1,hl2,w) reduction(max:o,r)
#endif /* _OPENMP */
        for (size_t i = 0u; i < n; ++i) {
          const wide CS = ((const float*)hc)[i];
          const wide SNR = ((const float*)hcat)[i];
          const wide SNI = ((const float*)hsat)[i];
          wide AE = W_ZERO, AN = W_ZERO;
          o = fmaxw(o, (w[i] = worc(CS, SNR, SNI)));
          r = fmaxw(r, wrec(((const float*)ha11)[i], ((const float*)ha22)[i], ((const float*)ha21r)[i], ((const float*)ha21i)[i], CS, SNR, SNI, ((const float*)hl1)[i], ((const float*)hl2)[i], &AE, &AN));
        }
      }
      else {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,ha11,ha22,ha21r,ha21i,hc,hcat,hsat,hl1,hl2,w) reduction(max:o,r)
#endif /* _OPENMP */
        for (size_t i = 0u; i < n; ++i) {
          const wide CS = ((const double*)hc)[i];
          const wide SNR = ((const double*)hcat)[i];
          const wide SNI = ((const double*)hsat)[i];
          wide AE = W_ZERO, AN = W_ZERO;
          o = fmaxw(o, (w[i] = worc(CS, SNR, SNI)));
          r = fmaxw(r, wrec(((const double*)ha11)[i], ((const double*)ha22)[i], ((const double*)ha21r)[i], ((const double*)ha21i)[i], CS, SNR, SNI, ((const double*)hl1)[i], ((const double*)hl2)[i], &AE, &AN));
        }
      }
    }
    else if (kind == -(int)sizeof(float)) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,ha11,ha22,ha21r,hc,hcat,hl1,hl2) reduction(max:o,r)
#endif /* _OPENMP */
      for (size_t i = 0u; i < n; ++i) {
        const wide CS = ((const float*)hc)[i];
        const wide SN = ((const float*)hcat)[i];
        wide AE = W_ZERO, AN = W_ZERO;
        o = fmaxw(o, worr(CS, SN));
        r = fmaxw(r, wrer(((const float*)ha11)[i], ((const float*)ha22)[i], ((const float*)ha21r)[i], CS, SN, ((const float*)hl1)[i], ((const float*)hl2)[i], &AE, &AN));
      }
    }
    else {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,ha11,ha22,ha21r,hc,hcat,hl1,hl2) reduction(max:o,r)
#endif /* _OPENMP */
      for (size_t i = 0u; i < n; ++i) {
        const wide CS = ((const double*)hc)[i];
        const wide SN = ((const double*)hcat)[i];
        wide AE = W_ZERO, AN = W_ZERO;
        o = fmaxw(o, worr(CS, SN));
        r = fmaxw(r, wrer(((const double*)ha11)[i], ((const double*)ha22)[i], ((const double*)ha21r)[i], CS, SN, ((const double*)hl1)[i], ((const double*)hl2)[i], &AE, &AN));
      }
    }
    (void)fprintf(stdout, "%s,", dtoa(a, (double)o));
    (void)fprintf(stdout, "%s", dtoa(a, (double)r));
    (void)fflush(stdout);
    size_t ix = n;
    if (kind > 0) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,o,w) reduction(min:ix)
#endif /* _OPENMP */
      for (size_t i = 0u; i < n; ++i)
        if (w[i] == o)
          ix = i;
      (void)fprintf(stderr, bf, j);
      (void)fprintf(stderr, ",%11llu,%s;", (unsigned long long)ix, dtoa(a, (double)o));
      if (kind == (int)sizeof(float)) {
        (void)fprintf(stderr, "%s,", dtoa(a, ((const float*)ha11)[ix]));
        (void)fprintf(stderr, "%s,", dtoa(a, ((const float*)ha22)[ix]));
        (void)fprintf(stderr, "(%s,", dtoa(a, ((const float*)ha21r)[ix]));
        (void)fprintf(stderr, "%s);", dtoa(a, ((const float*)ha21i)[ix]));
        (void)fprintf(stderr, "%s,", dtoa(a, ((const float*)hc)[ix]));
        (void)fprintf(stderr, "(%s,", dtoa(a, ((const float*)hcat)[ix]));
        (void)fprintf(stderr, "%s);", dtoa(a, ((const float*)hsat)[ix]));
        (void)fprintf(stderr, "%s,", dtoa(a, ((const float*)hl1)[ix]));
        (void)fprintf(stderr, "%s\n", dtoa(a, ((const float*)hl2)[ix]));
      }
      else {
        (void)fprintf(stderr, "%s,", dtoa(a, ((const double*)ha11)[ix]));
        (void)fprintf(stderr, "%s,", dtoa(a, ((const double*)ha22)[ix]));
        (void)fprintf(stderr, "(%s,", dtoa(a, ((const double*)ha21r)[ix]));
        (void)fprintf(stderr, "%s);", dtoa(a, ((const double*)ha21i)[ix]));
        (void)fprintf(stderr, "%s,", dtoa(a, ((const double*)hc)[ix]));
        (void)fprintf(stderr, "(%s,", dtoa(a, ((const double*)hcat)[ix]));
        (void)fprintf(stderr, "%s);", dtoa(a, ((const double*)hsat)[ix]));
        (void)fprintf(stderr, "%s,", dtoa(a, ((const double*)hl1)[ix]));
        (void)fprintf(stderr, "%s\n", dtoa(a, ((const double*)hl2)[ix]));
      }
      (void)fflush(stderr);
    }
    if (cnt != (size_t)read(fk, hc, cnt))
      return EXIT_FAILURE;
    if (cnt != (size_t)read(fl, hcat, cnt))
      return EXIT_FAILURE;
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);
    wide x = W_ZERO, m = W_ZERO;
    r = W_ZERO;
    if (abs(kind) == (int)sizeof(float)) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,hl1,hl2,hc,hcat) reduction(max:r,x,m)
#endif /* _OPENMP */
      for (size_t i = 0u; i < n; ++i) {
        wide AE = W_ZERO, AN = W_ZERO;
        const wide RE = wlam(((const float*)hl1)[i], ((const float*)hl2)[i], ((const float*)hc)[i], ((const float*)hcat)[i], &AE, &AN);
        r = fmaxw(r, RE);
        x = fmaxw(x, AE);
        m = fmaxw(m, AN);
      }
    }
    else {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,hl1,hl2,hc,hcat) reduction(max:r,x,m)
#endif /* _OPENMP */
      for (size_t i = 0u; i < n; ++i) {
        wide AE = W_ZERO, AN = W_ZERO;
        const wide RE = wlam(((const double*)hl1)[i], ((const double*)hl2)[i], ((const double*)hc)[i], ((const double*)hcat)[i], &AE, &AN);
        r = fmaxw(r, RE);
        x = fmaxw(x, AE);
        m = fmaxw(m, AN);
      }
    }
    (void)fprintf(stdout, "%s,", dtoa(a, (double)r));
    (void)fprintf(stdout, "%s,", dtoa(a, (double)x));
    (void)fprintf(stdout, "%s\n", dtoa(a, (double)m));
    (void)fflush(stdout);
  }

  free(hl2);
  free(hl1);
  free(hsat);
  free(hcat);
  free(hc);
  free(ha21i);
  free(ha21r);
  free(ha22);
  free(ha11);
  free(w);

  if (kind > 0) {
    (void)close(fj);
    (void)close(fr);
  }
  else /* real */
    (void)close(fh);
  (void)close(fg);
  (void)close(ff);
  (void)close(fl);
  (void)close(fk);

  return EXIT_SUCCESS;
}
