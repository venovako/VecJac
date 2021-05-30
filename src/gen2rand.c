#include "rnd.h"
#include "timer.h"

static const char fm[3] = { 'w', 'b', '\0' };

int main(int argc, char *argv[])
{
  if (argc != 4) {
    (void)fprintf(stderr, "%s { S | D | C | Z } n filename\n", *argv);
    (void)fprintf(stderr, "Outputs the following binary files: filename.\n");
    (void)fprintf(stderr, "k: \\lambda_1\n");
    (void)fprintf(stderr, "l: \\lambda_2\n");
    (void)fprintf(stderr, "f: A_{11}\n");
    (void)fprintf(stderr, "g: A_{22}\n");
    (void)fprintf(stderr, "r: \\Re{A_{12}}\n");
    (void)fprintf(stderr, "j: \\Im{A_{12}}\n");
    return EXIT_FAILURE;
  }
  const size_t n = atoz(argv[2]);
  if (!n) {
    (void)fprintf(stderr, "%s == 0\n", argv[2]);
    return EXIT_FAILURE;
  }
  const size_t nl = strlen(argv[3]);
  char *const fn = (char*)calloc((nl + 3u), sizeof(char));
  if (!fn)
    return EXIT_FAILURE;
  strcpy(fn, argv[3])[nl] = '.';
  fn[nl + 1u] = 'k';
  FILE *const fk = fopen(fn, fm);
  if (!fk) {
    perror(fn);
    return EXIT_FAILURE;
  }
  fn[nl + 1u] = 'l';
  FILE *const fl = fopen(fn, fm);
  if (!fl) {
    perror(fn);
    return EXIT_FAILURE;
  }
  fn[nl + 1u] = 'f';
  FILE *const ff = fopen(fn, fm);
  if (!ff) {
    perror(fn);
    return EXIT_FAILURE;
  }
  fn[nl + 1u] = 'g';
  FILE *const fg = fopen(fn, fm);
  if (!fg) {
    perror(fn);
    return EXIT_FAILURE;
  }
  fn[nl + 1u] = 'r';
  FILE *const fr = fopen(fn, fm);
  if (!fr) {
    perror(fn);
    return EXIT_FAILURE;
  }
  fn[nl + 1u] = 'j';
  FILE *const fj = fopen(fn, fm);
  if (!fj) {
    perror(fn);
    return EXIT_FAILURE;
  }
  unsigned rd[2u] = { 0u, 0u };
  uint64_t hz = tsc_get_freq_hz_(rd), be[2u] = { UINT64_C(0), UINT64_C(0) };
  (void)fprintf(stdout, "TSC frequency: %lu+(%u/%u) Hz.\n", hz, rd[0u], rd[1u]);
  const char t = (char)toupper(argv[1][0]);
  if (!t)
    return EXIT_FAILURE;
  else if (t == 'S') {
    const size_t ns = n * sizeof(float);
    float *const l1 = (float*)malloc(ns);
    if (!l1)
      return EXIT_FAILURE;
    float *const l2 = (float*)malloc(ns);
    if (!l2)
      return EXIT_FAILURE;
    float *const f = (float*)malloc(ns);
    if (!f)
      return EXIT_FAILURE;
    float *const g = (float*)malloc(ns);
    if (!g)
      return EXIT_FAILURE;
    float *const h = (float*)malloc(ns);
    if (!h)
      return EXIT_FAILURE;
    be[0u] = rdtsc_beg(rd);
    ssym2rand_(&n, l1, l2, f, g, h);
    be[1u] = rdtsc_end(rd);
    if (n != fwrite(h, sizeof(float), n, fr)) {
      perror("fwrite(r)");
      return EXIT_FAILURE;
    }
    if (n != fwrite(memset(h, 0, ns), sizeof(float), n, fj)) {
      perror("fwrite(j)");
      return EXIT_FAILURE;
    }
    free(h);
    if (n != fwrite(g, sizeof(float), n, fg)) {
      perror("fwrite(g)");
      return EXIT_FAILURE;
    }
    free(g);
    if (n != fwrite(f, sizeof(float), n, ff)) {
      perror("fwrite(f)");
      return EXIT_FAILURE;
    }
    free(f);
    if (n != fwrite(l2, sizeof(float), n, fl)) {
      perror("fwrite(l)");
      return EXIT_FAILURE;
    }
    free(l2);
    if (n != fwrite(l1, sizeof(float), n, fk)) {
      perror("fwrite(k)");
      return EXIT_FAILURE;
    }
    free(l1);
  }
  else if (t == 'D') {
    const size_t ns = n * sizeof(double);
    double *const l1 = (double*)malloc(ns);
    if (!l1)
      return EXIT_FAILURE;
    double *const l2 = (double*)malloc(ns);
    if (!l2)
      return EXIT_FAILURE;
    double *const f = (double*)malloc(ns);
    if (!f)
      return EXIT_FAILURE;
    double *const g = (double*)malloc(ns);
    if (!g)
      return EXIT_FAILURE;
    double *const h = (double*)malloc(ns);
    if (!h)
      return EXIT_FAILURE;
    be[0u] = rdtsc_beg(rd);
    dsym2rand_(&n, l1, l2, f, g, h);
    be[1u] = rdtsc_end(rd);
    if (n != fwrite(h, sizeof(double), n, fr)) {
      perror("fwrite(r)");
      return EXIT_FAILURE;
    }
    if (n != fwrite(memset(h, 0, ns), sizeof(double), n, fj)) {
      perror("fwrite(j)");
      return EXIT_FAILURE;
    }
    free(h);
    if (n != fwrite(g, sizeof(double), n, fg)) {
      perror("fwrite(g)");
      return EXIT_FAILURE;
    }
    free(g);
    if (n != fwrite(f, sizeof(double), n, ff)) {
      perror("fwrite(f)");
      return EXIT_FAILURE;
    }
    free(f);
    if (n != fwrite(l2, sizeof(double), n, fl)) {
      perror("fwrite(l)");
      return EXIT_FAILURE;
    }
    free(l2);
    if (n != fwrite(l1, sizeof(double), n, fk)) {
      perror("fwrite(k)");
      return EXIT_FAILURE;
    }
    free(l1);
  }
  else if (t == 'C') {
    const size_t ns = n * sizeof(float);
    float *const l1 = (float*)malloc(ns);
    if (!l1)
      return EXIT_FAILURE;
    float *const l2 = (float*)malloc(ns);
    if (!l2)
      return EXIT_FAILURE;
    float *const f = (float*)malloc(ns);
    if (!f)
      return EXIT_FAILURE;
    float *const g = (float*)malloc(ns);
    if (!g)
      return EXIT_FAILURE;
    float *const hr = (float*)malloc(ns);
    if (!hr)
      return EXIT_FAILURE;
    float *const hi = (float*)malloc(ns);
    if (!hi)
      return EXIT_FAILURE;
    be[0u] = rdtsc_beg(rd);
    cher2rand_(&n, l1, l2, f, g, hi, hr);
    be[1u] = rdtsc_end(rd);
    if (n != fwrite(hi, sizeof(float), n, fj)) {
      perror("fwrite(j)");
      return EXIT_FAILURE;
    }
    free(hi);
    if (n != fwrite(hr, sizeof(float), n, fr)) {
      perror("fwrite(r)");
      return EXIT_FAILURE;
    }
    free(hr);
    if (n != fwrite(g, sizeof(float), n, fg)) {
      perror("fwrite(g)");
      return EXIT_FAILURE;
    }
    free(g);
    if (n != fwrite(f, sizeof(float), n, ff)) {
      perror("fwrite(f)");
      return EXIT_FAILURE;
    }
    free(f);
    if (n != fwrite(l2, sizeof(float), n, fl)) {
      perror("fwrite(l)");
      return EXIT_FAILURE;
    }
    free(l2);
    if (n != fwrite(l1, sizeof(float), n, fk)) {
      perror("fwrite(k)");
      return EXIT_FAILURE;
    }
    free(l1);
  }
  else if (t == 'Z') {
    const size_t ns = n * sizeof(double);
    double *const l1 = (double*)malloc(ns);
    if (!l1)
      return EXIT_FAILURE;
    double *const l2 = (double*)malloc(ns);
    if (!l2)
      return EXIT_FAILURE;
    double *const f = (double*)malloc(ns);
    if (!f)
      return EXIT_FAILURE;
    double *const g = (double*)malloc(ns);
    if (!g)
      return EXIT_FAILURE;
    double *const hr = (double*)malloc(ns);
    if (!hr)
      return EXIT_FAILURE;
    double *const hi = (double*)malloc(ns);
    if (!hi)
      return EXIT_FAILURE;
    be[0u] = rdtsc_beg(rd);
    zher2rand_(&n, l1, l2, f, g, hi, hr);
    be[1u] = rdtsc_end(rd);
    if (n != fwrite(hi, sizeof(double), n, fj)) {
      perror("fwrite(j)");
      return EXIT_FAILURE;
    }
    free(hi);
    if (n != fwrite(hr, sizeof(double), n, fr)) {
      perror("fwrite(r)");
      return EXIT_FAILURE;
    }
    free(hr);
    if (n != fwrite(g, sizeof(double), n, fg)) {
      perror("fwrite(g)");
      return EXIT_FAILURE;
    }
    free(g);
    if (n != fwrite(f, sizeof(double), n, ff)) {
      perror("fwrite(f)");
      return EXIT_FAILURE;
    }
    free(f);
    if (n != fwrite(l2, sizeof(double), n, fl)) {
      perror("fwrite(l)");
      return EXIT_FAILURE;
    }
    free(l2);
    if (n != fwrite(l1, sizeof(double), n, fk)) {
      perror("fwrite(k)");
      return EXIT_FAILURE;
    }
    free(l1);
  }
  else {
    (void)fprintf(stderr, "%s does not start with S, D, C, or Z\n", argv[1]);
    return EXIT_FAILURE;
  }
  (void)fprintf(stdout, "The generation took %Lf s.\n", tsc_lap(hz, be[0u], be[1u]));
  if (fclose(fj))
    perror("j");
  if (fclose(fr))
    perror("r");
  if (fclose(fg))
    perror("g");
  if (fclose(ff))
    perror("f");
  if (fclose(fl))
    perror("l");
  if (fclose(fk))
    perror("k");
  return EXIT_SUCCESS;
}
