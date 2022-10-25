// Estimates the number of floating-point operations per second.
// One scalar FMA is counted as a single FLOP, as well as SCALEF.
// IT, the number of iterations, should be set sufficiently high; e.g., to 100000.
#include "ddpscl.h"
#include "djrotf.h"
#include "rnd.h"
#include "timer.h"

int main(int argc, char *argv[])
{
  if (argc != 3) {
    (void)fprintf(stderr, "%s N IT\n", *argv);
    return EXIT_FAILURE;
  }
  unsigned rd[2u] = { 0u, 0u };
  const uint64_t hz = tsc_get_freq_hz_(rd);
  if (!hz)
    return EXIT_FAILURE;
  const size_t n = atoz(argv[1u]);
  if (!n)
    return EXIT_SUCCESS;
  const size_t it = atoz(argv[2u]);
  if (!it)
    return EXIT_SUCCESS;
  if (set_cbwr() < 0)
    return EXIT_FAILURE;
  const size_t sz = (n * sizeof(double));
  double *const x = (double*)aligned_alloc(64u, sz);
  if (!x)
    return EXIT_FAILURE;
  double *const y = (double*)aligned_alloc(64u, sz);
  if (!y)
    return EXIT_FAILURE;
  const fint incx = 1ll, incy = 1ll;
  // warmup
  (void)(BLAS_D(dot)((const fint*)&n, x, &incx, y, &incy));
  const double aub = sqrt(DBL_MAX / (16u * n));
  uint64_t b = 0u, e = 0u, l = 0u;
  for (size_t i = 0u; i < it; ++i) {
    gendfrand_(&n, &aub, x);
    gendfrand_(&n, &aub, y);
    b = rdtsc_beg(rd);
    (void)(BLAS_D(dot)((const fint*)&n, x, &incx, y, &incy)); // n FMAs
    e = rdtsc_end(rd);
    if (e > b)
      l += (e - b);
  }
  long double t = tsc_lap(hz, 0u, l);
  long double f = ((n * it) / t);
  char a[31] = { '\0' };
  (void)fprintf(stdout, "%15.9Lf,%s,", t, xtoa(a, f));
  (void)fflush(stdout);
  const double rt3 = sqrt(3.0);
  const double cp6 = (0.5 * rt3);
  const double sp6 = 0.5;
  // warmup
  BLAS_D(rot)((const fint*)&n, x, &incx, y, &incy, &cp6, &sp6);
  l = 0u;
  for (size_t i = 0u; i < it; ++i) {
    gendfrand_(&n, &aub, x);
    gendfrand_(&n, &aub, y);
    b = rdtsc_beg(rd);
    BLAS_D(rot)((const fint*)&n, x, &incx, y, &incy, &cp6, &sp6); // 2n FMAs + 2n MULs
    e = rdtsc_end(rd);
    if (e > b)
      l += (e - b);
  }
  t = tsc_lap(hz, 0u, l);
  f = ((4u * n * it) / t);
  (void)fprintf(stdout, "%15.9Lf,%s,", t, xtoa(a, f));
  (void)fflush(stdout);
  double e2[2u] = { 1.0, 1.0 };
  double f2[2u] = { 1.5, 1.5 };
  // warmup
  (void)ddpscl_((const fnat*)&n, x, y, e2, f2);
  l = 0u;
  for (size_t i = 0u; i < it; ++i) {
    gendfrand_(&n, &aub, x);
    gendfrand_(&n, &aub, y);
    b = rdtsc_beg(rd);
    (void)ddpscl_((const fnat*)&n, x, y, e2, f2);
    e = rdtsc_end(rd);
    if (e > b)
      l += (e - b);
  }
  t = tsc_lap(hz, 0u, l);
#ifdef USE_2SUM
  f = (((18u * n + 35u) * it) / t);
#else /* !USE_2SUM */
  f = (((3u * n + 11u) * it) / t);
#endif /* ?USE_2SUM */
  (void)fprintf(stdout, "%15.9Lf,%s,", t, xtoa(a, f));
  (void)fflush(stdout);
  const double tp6 = (1.0 / rt3);
  // warmup
  (void)djrotf_((const fint*)&n, x, y, &cp6, &tp6);
  l = 0u;
  for (size_t i = 0u; i < it; ++i) {
    gendfrand_(&n, &aub, x);
    gendfrand_(&n, &aub, y);
    b = rdtsc_beg(rd);
    (void)djrotf_((const fint*)&n, x, y, &cp6, &tp6); // 2n FMAs + 2n MULs
    e = rdtsc_end(rd);
    if (e > b)
      l += (e - b);
  }
  t = tsc_lap(hz, 0u, l);
  f = ((4u * n * it) / t);
  (void)fprintf(stdout, "%15.9Lf,%s\n", t, xtoa(a, f));
  (void)fflush(stdout);
  free(y);
  free(x);
  return EXIT_SUCCESS;
}
