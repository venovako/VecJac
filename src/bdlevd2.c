#include "rnd.h"
#include "timer.h"
#include "laev2.h"
#include "wnrme.h"

int main(int argc, char *argv[])
{
  if (4 != argc) {
    (void)fprintf(stderr, "%s filename batch_size #batches\n", *argv);
    return EXIT_FAILURE;
  }

  const size_t n = atoz(argv[2u]);
  if (!n)
    return EXIT_SUCCESS;
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
  const char fm[3u] = { 'r', 'b', '\0' };

  fn[nl1] = 'f';
  FILE *const ff = fopen(fn, fm);
  if (!ff) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }
  fn[nl1] = 'g';
  FILE *const fg = fopen(fn, fm);
  if (!fg) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }
  fn[nl1] = 'r';
  FILE *const fr = fopen(fn, fm);
  if (!fr) {
    (void)fprintf(stderr, "Cannot open %s for reading!\n", fn);
    return EXIT_FAILURE;
  }

  const size_t nt = n * sizeof(double);
  double
    *const a11 = (double*)aligned_alloc(sizeof(double), nt),
    *const a22 = (double*)aligned_alloc(sizeof(double), nt),
    *const a21 = (double*)aligned_alloc(sizeof(double), nt),
    *const cs1 = (double*)aligned_alloc(sizeof(double), nt),
    *const sn1 = (double*)aligned_alloc(sizeof(double), nt),
    *const l1 = (double*)aligned_alloc(sizeof(double), nt),
    *const l2 = (double*)aligned_alloc(sizeof(double), nt);
  assert(a11);
  assert(a22);
  assert(a21);
  assert(cs1);
  assert(sn1);
  assert(l1);
  assert(l2);

  unsigned rd[2u] = { 0u, 0u };
  uint64_t hz = tsc_get_freq_hz_(rd), be[2u] = { UINT64_C(0), UINT64_C(0) };
  (void)fprintf(stderr, "TSC frequency: %lu+(%u/%u) Hz.\n", hz, rd[0u], rd[1u]);

  (void)fprintf(stdout, "\"B\",\"Ts\"\n");
  const char *bf = (const char*)NULL;
  if (b <= 10u)
    bf = "%zu";
  else if (b <= 100u)
    bf = "%2zu";
  else if (b <= 1000u)
    bf = "%3zu";
  else // b > 1000
    bf = "%4zu";
  int th = 0;

  for (size_t j = 0u; j < b; ++j) {
    (void)fprintf(stdout, bf, j);
    (void)fflush(stdout);
    if (n != fread(a11, sizeof(double), n, ff))
      return EXIT_FAILURE;
    if (n != fread(a22, sizeof(double), n, fg))
      return EXIT_FAILURE;
    if (n != fread(a21, sizeof(double), n, fr))
      return EXIT_FAILURE;
    (void)fprintf(stdout, ",");
    (void)fflush(stdout);
    be[0u] = rdtsc_beg(rd);
#pragma omp parallel for default(none) shared(n,a11,a22,a21,l1,l2,cs1,sn1) reduction(max:th)
    for (size_t i = 0u; i < n; ++i) {
      dlevd2(a11[i], a22[i], a21[i], (l1 + i), (l2 + i), (cs1 + i), (sn1 + i));
      th = imax(th, omp_get_num_threads());
    }
    be[1u] = rdtsc_end(rd);
    (void)fprintf(stdout, "%15.9Lf\n", tsc_lap(hz, be[0u], be[1u]));
    (void)fflush(stdout);
  }
  (void)fprintf(stderr, "max(#threads) = %u\n", (unsigned)th);

  (void)fclose(fr);
  (void)fclose(fg);
  (void)fclose(ff);

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
