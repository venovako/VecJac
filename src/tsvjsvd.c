#include "aalloc.h"
#include "mtxio.h"
#include "pjs.h"
#include "timer.h"
#include "svjsvd.h"

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
  if ((n >> 1u) & VSL_1)
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

  float *G = (float*)NULL;
  if (salloc2_(&m, &n, &G, &ldG) < 0)
    return 6;

  if (sread2_(&m, &n, G, &ldG, &gd))
    return 5;
  if (close(gd))
    return 5;

  float *V = (float*)NULL;
  if (salloc2_(&n, &n, &V, &ldV) < 0)
    return 7;

  float *const w = (float*)aligned_alloc(VA, (7u * (n * sizeof(float))));
  if (!w)
    return 8;
  float *const eS = w;
  float *const fS = eS + n;
  float *const work = fS + n;
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
  unsigned *const iwork = (unsigned*)aligned_alloc(VA, ((n >> VSLlg) * sizeof(unsigned)));
  if (!iwork)
    return 9;
  *iwork = 0u;

  const unsigned swp = 999u;
  unsigned rd[2u] = { 0u, 0u };
  const uint64_t hz = tsc_get_freq_hz_(rd);
  const uint64_t b = rdtsc_beg(rd);
  const fint o = svjsvd_(&m, &n, G, &ldG, V, &ldV, eS, fS, js, &stp, &swp, work, iwork);
  const uint64_t e = rdtsc_end(rd);
  (void)fprintf(stdout, "\"%s\",%1ld,%4llu,%4llu,%15.9Lf,%3lld\n", bn,
#ifdef _OPENMP
                j
#else /* !_OPENMP */
                (j - 1l)
#endif /* ?_OPENMP */
                , (unsigned long long)m, (unsigned long long)n, tsc_lap(hz, b, e), (long long)o);
  (void)fflush(stdout);
  free(iwork);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(eS,fS,ws,n)
#endif /* _OPENMP */
  for (fnat i = 0u; i < n; ++i) {
    if (eS[i] != 0.0f)
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
  if (wwrite1_(&n, ws, &sd))
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
  *(size_t*)js = (n * (n * sizeof(float)));
  if (resizef_(&vd, (const size_t*)js))
    return 11;
  if (swrite2_(&n, &n, V, &ldV, &vd))
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
  *(size_t*)js = (m * (n * sizeof(float)));
  if (resizef_(&ud, (const size_t*)js))
    return 12;
  if (swrite2_(&m, &n, G, &ldG, &ud))
    return 12;
  if (close(ud))
    return 12;
  free(G);

  free((void*)js);
  return EXIT_SUCCESS;
}
