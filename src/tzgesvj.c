#include "mtxio.h"
#include "timer.h"
#include "vec.h"

#ifndef USE_MKL
extern void zgesvj_(const char *const, const char *const, const char *const, const fint *const, const fint *const, double complex *const, const fint *const, double *const, const fint *const, double complex *const, const fint *const, double complex *const, const fint *const, double *const, const fint *const, fint *const);
#endif /* !USE_MKL */
extern void znssvj_(const char *const, const char *const, const char *const, const fint *const, const fint *const, double complex *const, const fint *const, double *const, const fint *const, double complex *const, const fint *const, double complex *const, const fint *const, double *const, const fint *const, const fint *const, fint *const);

int main(int argc, char *argv[])
{
  (void)set_cbwr();

  if ((argc < 4) || (argc > 5))
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
  const fint nsweep = ((argc == 5) ? (fint)atoz(argv[4u]) : (fint)0);

  const fint ldG = m;
  double complex *const G = (double complex*)aligned_alloc(VA, (ldG * (n * sizeof(double complex))));
  if (!G)
    return EXIT_FAILURE;
  const fint ldV = n;
  double complex *const V = (double complex*)aligned_alloc(VA, (ldV * (n * sizeof(double complex))));
  if (!V)
    return EXIT_FAILURE;
  double *const sva = (double*)aligned_alloc(VA, (n * sizeof(double)));
  if (!sva)
    return EXIT_FAILURE;
  const fint mv = n;
  fint lwork = -1;
  fint lrwork = -1;
  fint info = 0;
  double complex cwork1[1u] = { 0.0 };
  double rwork1[1u] = { 1.0 };
  LAPACK_Z(gesvj)("G", "U", "V", &m, &n, G, &ldG, sva, &mv, V, &ldV, cwork1, &lwork, rwork1, &lrwork, &info);
  if (info)
    return EXIT_FAILURE;
  lwork = (fint)(creal(*cwork1));
  if (lwork < (m + n))
    return EXIT_FAILURE;
  lrwork = (fint)*rwork1;
  double complex *const cwork = (double complex*)aligned_alloc(VA, (lwork * sizeof(double complex)));
  if (!cwork)
    return EXIT_FAILURE;
  double *const rwork = (double*)aligned_alloc(VA, (lrwork * sizeof(double)));
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
  if (nsweep)
    znssvj_("G", "U", "V", &m, &n, G, &ldG, sva, &mv, V, &ldV, cwork, &lwork, rwork, &lrwork, &nsweep, &info);
  else
    LAPACK_Z(gesvj)("G", "U", "V", &m, &n, G, &ldG, sva, &mv, V, &ldV, cwork, &lwork, rwork, &lrwork, &info);
  const uint64_t e = rdtsc_end(rd);

  (void)fprintf(stdout, "\"%s\",%4lld,%4lld,%1lld,%15.9Lf,%#.17e,%4lld,%4lld,%3lld,%#.17e,%#.17e\n", bn, (long long)m, (long long)n, (long long)info, tsc_lap(hz, b, e), rwork[0u], (long long)(rwork[1u]), (long long)(rwork[2u]), (long long)(rwork[3u]), rwork[4u], rwork[5u]);
  (void)fflush(stdout);
  *rwork1 = *rwork;
  free(rwork);

  *(size_t*)cwork1 = (n * sizeof(wide));
  wide *const ws = (wide*)memset(cwork, 0, *(const size_t*)cwork1);
  if (*rwork1 == 1.0)
    for (fint i = 0; i < n; ++i)
      ws[i] = sva[i];
  else { // SCALE .NE. ONE
    const wide s = (wide)*rwork1;
    for (fint i = 0; i < n; ++i)
      ws[i] = sva[i] * s;
  }
  free(sva);

  const int sd = open_wo_(bn, "SL");
  if (sd < 0)
    return EXIT_FAILURE;
  if (resizef_(&sd, (const size_t*)cwork1))
    return EXIT_FAILURE;
  if (wwrite1_(&n, ws, &sd))
    return EXIT_FAILURE;
  if (close(sd))
    return EXIT_FAILURE;
  free(cwork);

  *(size_t*)cwork1 = (n * (n * sizeof(double complex)));
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

  *(size_t*)cwork1 = (m * (n * sizeof(double complex)));
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
  (void)fprintf(stderr, "%s M N BaseName [NSWEEP]\n", *argv);
  return EXIT_FAILURE;
}
