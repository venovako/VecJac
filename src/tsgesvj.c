#include "mtxio.h"
#include "timer.h"
#include "vec.h"

#ifndef USE_MKL
extern void sgesvj_(const char *const, const char *const, const char *const, const fint *const, const fint *const, float *const, const fint *const, float *const, const fint *const, float *const, const fint *const, float *const, const fint *const, fint *const);
#endif /* !USE_MKL */
extern void snssvj_(const char *const, const char *const, const char *const, const fint *const, const fint *const, float *const, const fint *const, float *const, const fint *const, float *const, const fint *const, float *const, const fint *const, const fint *const, fint *const);

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
  float *const G = (float*)aligned_alloc(VA, (ldG * (n * sizeof(float))));
  if (!G)
    return EXIT_FAILURE;
  const fint ldV = n;
  float *const V = (float*)aligned_alloc(VA, (ldV * (n * sizeof(float))));
  if (!V)
    return EXIT_FAILURE;
  float *const sva = (float*)aligned_alloc(VA, (n * sizeof(float)));
  if (!sva)
    return EXIT_FAILURE;
  const fint mv = n;
  const fint lwork = imax(6, (n << 2));
  float *const work = (float*)aligned_alloc(VA, (lwork * sizeof(float)));
  if (!work)
    return EXIT_FAILURE;
  fint info = 0;

  const int gd = open_ro_(bn, "G");
  if (gd < 0)
    return EXIT_FAILURE;

  if (sread2_((const fnat*)&m, (const fnat*)&n, G, (const fnat*)&ldG, &gd))
    return EXIT_FAILURE;

  if (close(gd))
    return EXIT_FAILURE;

  unsigned rd[2u] = { 0u, 0u };
  const uint64_t hz = tsc_get_freq_hz_(rd);
  const uint64_t b = rdtsc_beg(rd);
  if (nsweep)
    snssvj_("G", "U", "V", &m, &n, G, &ldG, sva, &mv, V, &ldV, work, &lwork, &nsweep, &info);
  else
    LAPACK_S(gesvj)("G", "U", "V", &m, &n, G, &ldG, sva, &mv, V, &ldV, work, &lwork, &info);
  const uint64_t e = rdtsc_end(rd);

  (void)fprintf(stdout, "\"%s\",%4lld,%4lld,%1lld,%15.9Lf,%#.9e,%4lld,%4lld,%3lld,%#.9e,%#.9e\n", bn, (long long)m, (long long)n, (long long)info, tsc_lap(hz, b, e), work[0u], (long long)(work[1u]), (long long)(work[2u]), (long long)(work[3u]), work[4u], work[5u]);
  (void)fflush(stdout);

  const float ds = *work;
  size_t l = (n * sizeof(wide));
  wide *const ws = (wide*)memset(work, 0, l);
  if (ds == 1.0f)
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
  if (wwrite1_((const fnat*)&n, ws, &sd))
    return EXIT_FAILURE;
  if (close(sd))
    return EXIT_FAILURE;
  free(work);

  l = (n * (n * sizeof(float)));
  const int vd = open_wo_(bn, "VL");
  if (vd < 0)
    return EXIT_FAILURE;
  if (resizef_(&vd, &l))
    return EXIT_FAILURE;
  if (swrite2_((const fnat*)&n, (const fnat*)&n, V, (const fnat*)&ldV, &vd))
    return EXIT_FAILURE;
  if (close(vd))
    return EXIT_FAILURE;
  free(V);

  l = (m * (n * sizeof(float)));
  const int ud = open_wo_(bn, "UL");
  if (ud < 0)
    return EXIT_FAILURE;
  if (resizef_(&ud, &l))
    return EXIT_FAILURE;
  if (swrite2_((const fnat*)&m, (const fnat*)&n, G, (const fnat*)&ldG, &ud))
    return EXIT_FAILURE;
  if (close(ud))
    return EXIT_FAILURE;
  free(G);

  return EXIT_SUCCESS;

 err:
  (void)fprintf(stderr, "%s M N BaseName [NSWEEP]\n", *argv);
  return EXIT_FAILURE;
}
