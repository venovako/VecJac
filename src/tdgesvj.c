#include "mtxio.h"
#include "timer.h"
#include "vec.h"

#ifndef USE_MKL
extern void dgesvj_(const char *const, const char *const, const char *const, const fint *const, const fint *const, double *const, const fint *const, double *const, const fint *const, double *const, const fint *const, double *const, const fint *const, fint *const);
#endif /* !USE_MKL */
extern void dnssvj_(const char *const, const char *const, const char *const, const fint *const, const fint *const, double *const, const fint *const, double *const, const fint *const, double *const, const fint *const, double *const, const fint *const, const fint *const, fint *const);

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
  double *const G = (double*)aligned_alloc(VA, (ldG * (n * sizeof(double))));
  if (!G)
    return EXIT_FAILURE;
  const fint ldV = n;
  double *const V = (double*)aligned_alloc(VA, (ldV * (n * sizeof(double))));
  if (!V)
    return EXIT_FAILURE;
  double *const sva = (double*)aligned_alloc(VA, (n * sizeof(double)));
  if (!sva)
    return EXIT_FAILURE;
  fint mv = n;
  const fint n2 = (n << 1);
  const fint lwork = ((n2 < (fint)6) ? (fint)6 : n2);
  double *const work = (double*)aligned_alloc(VA, (lwork * sizeof(double)));
  if (!work)
    return EXIT_FAILURE;
  fint info = 0;

  const int gd = open_ro_(bn, "G");
  if (gd < 0)
    return EXIT_FAILURE;

  if (dread2_(&m, &n, G, &ldG, &gd))
    return EXIT_FAILURE;

  if (close(gd))
    return EXIT_FAILURE;

  unsigned rd[2u] = { 0u, 0u };
  const uint64_t hz = tsc_get_freq_hz_(rd);
  const uint64_t b = rdtsc_beg(rd);
  if (nsweep)
    dnssvj_("G", "U", "V", &m, &n, G, &ldG, sva, &mv, V, &ldV, work, &lwork, &nsweep, &info);
  else
    LAPACK_D(gesvj)("G", "U", "V", &m, &n, G, &ldG, sva, &mv, V, &ldV, work, &lwork, &info);
  const uint64_t e = rdtsc_end(rd);

  (void)fprintf(stdout, "\"%s\",%4lld,%4lld,%1lld,%15.9Lf,%#.17e,%4lld,%4lld,%3lld,%#.17e,%#.17e\n", bn, m, n, info, tsc_lap(hz, b, e), work[0u], (fint)(work[1u]), (fint)(work[2u]), (fint)(work[3u]), work[4u], work[5u]);
  (void)fflush(stdout);

  const double ds = *work;
  size_t l = (n * sizeof(wide));
  wide *const ws = (wide*)memset(work, 0, l);
  if (ds == 1.0)
    for (fint i = 0; i < n; ++i)
      ws[i] = sva[i];
  else { // SCALE .NE. ONE
    const wide s = ds;
    for (fint i = 0; i < n; ++i)
      ws[i] = sva[i] * s;
  }
  free(sva);

  const int sd = open_wo_(bn, "SL");
  if (sd < 0)
    return EXIT_FAILURE;
  if (resizef_(&sd, &l))
    return EXIT_FAILURE;
  mv = 1;
  if (dwrite2_(&n2, &mv, work, &n2, &sd))
    return EXIT_FAILURE;
  if (close(sd))
    return EXIT_FAILURE;
  free(work);

  l = (n * (n * sizeof(double)));
  const int vd = open_wo_(bn, "VL");
  if (vd < 0)
    return EXIT_FAILURE;
  if (resizef_(&vd, &l))
    return EXIT_FAILURE;
  if (dwrite2_(&n, &n, V, &ldV, &vd))
    return EXIT_FAILURE;
  if (close(vd))
    return EXIT_FAILURE;
  free(V);

  l = (m * (n * sizeof(double)));
  const int ud = open_wo_(bn, "UL");
  if (ud < 0)
    return EXIT_FAILURE;
  if (resizef_(&ud, &l))
    return EXIT_FAILURE;
  if (dwrite2_(&m, &n, G, &ldG, &ud))
    return EXIT_FAILURE;
  if (close(ud))
    return EXIT_FAILURE;
  free(G);

  return EXIT_SUCCESS;

 err:
  (void)fprintf(stderr, "%s M N BaseName [NSWEEP]\n", *argv);
  return EXIT_FAILURE;
}
