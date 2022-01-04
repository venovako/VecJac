#include "mtxio.h"
#include "timer.h"

int main(int argc, char *argv[])
{
  (void)set_cbwr();

  if (argc != 4)
    goto err;

  const fint m = (fint)atoz(argv[1u]);
  if (!m)
    goto err;
  const fint n = (fint)atoz(argv[2u]);
  if (!n)
    goto err;
  const char *const bn = argv[3u];
  if (!*bn)
    goto err;

  const fint ldG = m;
  double complex *const G = (double complex*)aligned_alloc(alignof(double complex), (ldG * (n * sizeof(double complex))));
  if (!G)
    return EXIT_FAILURE;
  const fint ldV = n;
  double complex *const V = (double complex*)aligned_alloc(alignof(double complex), (ldV * (n * sizeof(double complex))));
  if (!V)
    return EXIT_FAILURE;
  double *const sva = (double*)calloc((size_t)n, sizeof(double));
  if (!sva)
    return EXIT_FAILURE;
  const fint mv = n;
  fint lwork = -1;
  fint lrwork = -1;
  fint info = 0;
  double complex cwork1[1u] = { CMPLX(0.0, 0.0) };
  double rwork1[1u] = { 1.0 };
  LAPACK_Z(gesvj)("G", "U", "V", &m, &n, G, &ldG, sva, &mv, V, &ldV, cwork1, &lwork, rwork1, &lrwork, &info);
  if (info)
    return EXIT_FAILURE;
  lwork = (fint)(creal(*cwork1));
  if (lwork < (m + n))
    return EXIT_FAILURE;
  lrwork = (fint)*rwork1;
  double complex *const cwork = calloc((size_t)lwork, sizeof(double complex));
  if (!cwork)
    return EXIT_FAILURE;
  double *const rwork = calloc((size_t)lrwork, sizeof(double));
  if (!rwork)
    return EXIT_FAILURE;

  const int gd = open_ro_(bn, "G");
  if (gd < 0)
    return EXIT_FAILURE;

  if (zread2_(&m, &n, G, &ldG, &gd))
    return EXIT_FAILURE;

  if (close(gd))
    return EXIT_FAILURE;

  unsigned rd[2u] = { 0u, 0u };
  const uint64_t hz = tsc_get_freq_hz_(rd);
  const uint64_t b = rdtsc_beg(rd);
  LAPACK_Z(gesvj)("G", "U", "V", &m, &n, G, &ldG, sva, &mv, V, &ldV, cwork, &lwork, rwork, &lrwork, &info);
  while (info > 0) {
    const fint swp = (fint)(rwork[3u]);
    LAPACK_Z(gesvj)("G", "U", "A", &m, &n, G, &ldG, sva, &mv, V, &ldV, cwork, &lwork, rwork, &lrwork, &info);
    rwork[3u] += swp;
  }
  const uint64_t e = rdtsc_end(rd);

  (void)fprintf(stdout, "%4lld,%4lld,%1lld,%15.9Lf,%#.17e,%4lld,%4lld,%2lld,%#.17e,%#.17e\n", m, n, info, tsc_lap(hz, b, e), rwork[0u], (fint)(rwork[1u]), (fint)(rwork[2u]), (fint)(rwork[3u]), rwork[4u], rwork[5u]);
  (void)fflush(stdout);
  *rwork1 = *rwork;
  free(rwork);

  wide *const ws = (wide*)memset(cwork, 0, (n * sizeof(wide)));
  if (*rwork1 == 1.0)
    for (fint i = 0; i < n; ++i)
      ws[i] = sva[i];
  else { // SCALE .NE. ONE
    const wide s = (wide)*rwork1;
    for (fint i = 0; i < n; ++i)
      ws[i] = sva[i] * s;
  }
  free(sva);

  *(size_t*)cwork1 = (n * sizeof(wide));
  const int sd = open_wo_(bn, "SL");
  if (sd < 0)
    return EXIT_FAILURE;
  if (resizef_(&sd, (const size_t*)cwork1))
    return EXIT_FAILURE;
  *(fnat*)cwork1 = 1u;
  if (zwrite2_(&n, (const fnat*)cwork1, cwork, &n, &sd))
    return EXIT_FAILURE;
  if (close(sd))
    return EXIT_FAILURE;
  free(cwork);

  *(size_t*)cwork1 = (ldV * (n * sizeof(double)));
  const int vd = open_wo_(bn, "VL");
  if (vd < 0)
    return EXIT_FAILURE;
  if (resizef_(&vd, (const size_t*)cwork1))
    return EXIT_FAILURE;
  if (zwrite2_(&n, &n, V, &ldV, &vd))
    return EXIT_FAILURE;
  if (close(vd))
    return EXIT_FAILURE;
  free(V);

  *(size_t*)cwork1 = (ldG * (n * sizeof(double)));
  const int ud = open_wo_(bn, "UL");
  if (ud < 0)
    return EXIT_FAILURE;
  if (resizef_(&ud, (const size_t*)cwork1))
    return EXIT_FAILURE;
  if (zwrite2_(&m, &n, G, &ldG, &ud))
    return EXIT_FAILURE;
  if (close(ud))
    return EXIT_FAILURE;
  free(G);

  return EXIT_SUCCESS;

 err:
  (void)fprintf(stderr, "%s M N BaseName\n", *argv);
  return EXIT_FAILURE;
}
