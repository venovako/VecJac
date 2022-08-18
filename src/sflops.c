// Estimates the number of floating-point operations per second.
// One scalar FMA is counted as a single FLOP, as well as SCALEF.
// IT, the number of iterations, should be set sufficiently high; e.g., to 200000.
#include "sdpscl.h"
#include "sjrotf.h"
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
  const size_t sz = (n * sizeof(float));
  float *const x = (float*)aligned_alloc(64u, sz);
  if (!x)
    return EXIT_FAILURE;
  float *const y = (float*)aligned_alloc(64u, sz);
  if (!y)
    return EXIT_FAILURE;
  const fint incx = 1ll, incy = 1ll;
  // warmup
  (void)(BLAS_S(dot)((const fint*)&n, x, &incx, y, &incy));
  const float aub = sqrtf(FLT_MAX / (16u * n));
  uint64_t b = 0u, e = 0u, l = 0u;
  for (size_t i = 0u; i < it; ++i) {
    gensfrand_(&n, &aub, x);
    gensfrand_(&n, &aub, y);
    b = rdtsc_beg(rd);
    (void)(BLAS_S(dot)((const fint*)&n, x, &incx, y, &incy)); // n FMAs
    e = rdtsc_end(rd);
    if (e >= b)
      l += (e - b);
  }
  long double t = tsc_lap(hz, 0u, l);
  long double f = ((n * it) / t);
  char a[31] = { '\0' };
  (void)fprintf(stdout, "%15.9Lf,%s,", t, xtoa(a, f));
  (void)fflush(stdout);
  const float rt3 = sqrtf(3.0f);
  const float cp6 = (0.5f * rt3);
  const float sp6 = 0.5f;
  // warmup
  BLAS_S(rot)((const fint*)&n, x, &incx, y, &incy, &cp6, &sp6);
  l = 0u;
  for (size_t i = 0u; i < it; ++i) {
    gensfrand_(&n, &aub, x);
    gensfrand_(&n, &aub, y);
    b = rdtsc_beg(rd);
    BLAS_S(rot)((const fint*)&n, x, &incx, y, &incy, &cp6, &sp6); // 2n FMAs + 2n MULs
    e = rdtsc_end(rd);
    if (e >= b)
      l += (e - b);
  }
  t = tsc_lap(hz, 0u, l);
  f = ((4u * n * it) / t);
  (void)fprintf(stdout, "%15.9Lf,%s,", t, xtoa(a, f));
  (void)fflush(stdout);
  float e2[2u] = { 1.0f, 1.0f };
  float f2[2u] = { 1.5f, 1.5f };
  // warmup
  (void)sdpscl_((const fnat*)&n, x, y, e2, f2);
  l = 0u;
  for (size_t i = 0u; i < it; ++i) {
    gensfrand_(&n, &aub, x);
    gensfrand_(&n, &aub, y);
    b = rdtsc_beg(rd);
    (void)sdpscl_((const fnat*)&n, x, y, e2, f2); // n FMAs + 2n SCALEFs + O(1)
    e = rdtsc_end(rd);
    if (e >= b)
      l += (e - b);
  }
  t = tsc_lap(hz, 0u, l);
  f = (((3u * n + 19u) * it) / t);
  (void)fprintf(stdout, "%15.9Lf,%s,", t, xtoa(a, f));
  (void)fflush(stdout);
  const float tp6 = (1.0f / rt3);
  // warmup
  (void)sjrotf_((const fint*)&n, x, y, &cp6, &tp6);
  l = 0u;
  for (size_t i = 0u; i < it; ++i) {
    gensfrand_(&n, &aub, x);
    gensfrand_(&n, &aub, y);
    b = rdtsc_beg(rd);
    (void)sjrotf_((const fint*)&n, x, y, &cp6, &tp6); // 2n FMAs + 2n MULs
    e = rdtsc_end(rd);
    if (e >= b)
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
