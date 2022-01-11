#include "aalloc.h"
#include "mtxio.h"
#include "pjs.h"
#include "timer.h"
#include "dvjsvd.h"

int main(int argc, char *argv[])
{
  (void)set_cbwr();

  if (argc != 5) {
    (void)fprintf(stderr, "%s J M N BaseName\n", *argv);
    return 1;
  }

  const fnat m = (fnat)atoz(argv[2u]);
  if (!m)
    return 3;
  fnat ldG = m;

  const fnat n = (fnat)atoz(argv[3u]);
  if (!n)
    return 4;
  if (n > m)
    return 4;
  if (n & 1u)
    return 4;
  if ((n >> 1u) & VDL_1)
    return 4;
  fnat ldV = n;

  const long j = atol(argv[1u]);
  switch (j) {
  case PJS_ME:
  case PJS_MM:
    break;
  default:
    return 2;
  }
  unsigned stp = 0u;
  const unsigned *const js = pjs(j, (unsigned)n, &stp);
  if (!js)
    return 2;

  const char *const bn = argv[4u];
  if (!*bn)
    return 5;
  const int gd = open_ro_(bn, "G");
  if (gd < 0)
    return 5;

  double *G = (double*)NULL;
  if (dalloc2_(&m, &n, &G, &ldG) < 0)
    return 6;

  if (dread2_(&m, &n, G, &ldG, &gd))
    return 5;
  if (close(gd))
    return 5;

  double *V = (double*)NULL;
  if (dalloc2_(&n, &n, &V, &ldV) < 0)
    return 7;

  double *const w = (double*)aligned_alloc(VA, (6u * (n * sizeof(double))));
  if (!w)
    return 8;
  double *const eS = w;
  double *const fS = eS + n;
  double *const work = fS + n;
  wide *const ws = (wide*)work;
#ifdef JTRACE
  (void)sprintf((char*)work, "%s.%ld", bn,
#ifdef _OPENMP
                j
#else /* !_OPENMP */
                (j - 1l)
#endif /* ?_OPENMP */
                );
#endif /* JTRACE */
  unsigned *const iwork = (unsigned*)aligned_alloc(VA, ((n >> VDLlg) * sizeof(unsigned)));
  if (!iwork)
    return 9;

  const unsigned swp = 999u;
  unsigned rd[2u] = { 0u, 0u };
  const uint64_t hz = tsc_get_freq_hz_(rd);
  const uint64_t b = rdtsc_beg(rd);
  const fint o = dvjsvd_(&m, &n, G, &ldG, V, &ldV, eS, fS, js, &stp, &swp, work, iwork);
  const uint64_t e = rdtsc_end(rd);
  (void)fprintf(stdout, "\"%s\",%1ld,%4llu,%4llu,%15.9Lf,%3lld\n", bn,
#ifdef _OPENMP
                j
#else /* !_OPENMP */
                (j - 1l)
#endif /* ?_OPENMP */
                , m, n, tsc_lap(hz, b, e), o);
  (void)fflush(stdout);
  free(iwork);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(eS,fS,ws,n)
#endif /* _OPENMP */
  for (fnat i = 0u; i < n; ++i) {
    if (eS[i] != 0.0)
      ws[i] = scalbw(fS[i], eS[i]);
    else // 2^0 == 1
      ws[i] = fS[i];
  }

  const int sd = open_wo_(bn,
#ifdef _OPENMP
                          ((j == PJS_ME) ? "S2" : "S4")
#else /* !_OPENMP */
                          ((j == PJS_ME) ? "S1" : "S3")
#endif /* ?_OPENMP */
                          );
  if (sd < 0)
    return 10;
  *(size_t*)js = (n * sizeof(wide));
  if (resizef_(&sd, (const size_t*)js))
    return 10;
  const fnat n2 = (n << 1u), n1 = 1u;
  if (dwrite2_(&n2, &n1, work, &n2, &sd))
    return 10;
  if (close(sd))
    return 10;
  free(w);

  const int vd = open_wo_(bn,
#ifdef _OPENMP
                          ((j == PJS_ME) ? "V2" : "V4")
#else /* !_OPENMP */
                          ((j == PJS_ME) ? "V1" : "V3")
#endif /* ?_OPENMP */
                          );
  if (vd < 0)
    return 11;
  *(size_t*)js = (n * (n * sizeof(double)));
  if (resizef_(&vd, (const size_t*)js))
    return 11;
  if (dwrite2_(&n, &n, V, &ldV, &vd))
    return 11;
  if (close(vd))
    return 11;
  free(V);

  const int ud = open_wo_(bn,
#ifdef _OPENMP
                          ((j == PJS_ME) ? "U2" : "U4")
#else /* !_OPENMP */
                          ((j == PJS_ME) ? "U1" : "U3")
#endif /* ?_OPENMP */
                          );
  if (ud < 0)
    return 12;
  *(size_t*)js = (m * (n * sizeof(double)));
  if (resizef_(&ud, (const size_t*)js))
    return 12;
  if (dwrite2_(&m, &n, G, &ldG, &ud))
    return 12;
  if (close(ud))
    return 12;
  free(G);

  free(js);
  return EXIT_SUCCESS;
}
