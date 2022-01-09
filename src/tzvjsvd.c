#include "aalloc.h"
#include "mtxio.h"
#include "pjs.h"
#include "timer.h"
#include "zmerge.h"
#include "zsplit.h"
#include "zvjsvd.h"

int main(int argc, char *argv[])
{
  (void)set_cbwr();

  if (argc != 4) {
    (void)fprintf(stderr, "%s J M N BaseName\n", *argv);
    return 1;
  }

  const fnat m = (fnat)atoz(argv[2u]);
  if (!m)
    return 3;
  fnat ldG = m, ldGr = m, ldGi = m;

  const fnat n = (fnat)atoz(argv[3u]);
  if (!n)
    return 4;
  if (n > m)
    return 4;
  if (n & 1u)
    return 4;
  if ((n >> 1u) & VDL_1)
    return 4;
  fnat ldV = n, ldVr = n, ldVi = n;

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

  double complex *G = (double complex*)NULL;
  double *Gr = (double*)NULL;
  double *Gi = (double*)NULL;
  if (zalloc2_(&m, &n, &G, &ldG, &Gr, &ldGr, &Gi, &ldGi) < 0)
    return 6;

  if (zread2_(&m, &n, G, &ldG, &gd))
    return 5;
  if (close(gd))
    return 5;

  (void)fprintf(stdout, "%4llu,%4llu,", m, n);
  (void)fflush(stdout);

  unsigned rd[2u] = { 0u, 0u };
  const uint64_t hz = tsc_get_freq_hz_(rd);

  uint64_t b = rdtsc_beg(rd);
  if (zsplit_(&m, &n, G, &ldG, Gr, &ldGr, Gi, &ldGi) < 0)
    return 6;
  uint64_t e = rdtsc_end(rd);
  (void)fprintf(stdout, "%15.9Lf,", tsc_lap(hz, b, e));
  (void)fflush(stdout);

  double complex *V = (double complex*)NULL;
  double *Vr = (double*)NULL;
  double *Vi = (double*)NULL;
  if (zalloc2_(&n, &n, &V, &ldV, &Vr, &ldVr, &Vi, &ldVi) < 0)
    return 7;

  double *const w = (double*)aligned_alloc(VA, (7u * (n * sizeof(double))));
  if (!w)
    return 8;
  double *const eS = w;
  double *const fS = eS + n;
  double *const work = fS + n;
  wide *const ws = (wide*)work;
  unsigned *const iwork = (unsigned*)aligned_alloc(VA, ((n >> VDLlg) * sizeof(unsigned)));
  if (!iwork)
    return 9;

  const unsigned swp = 999u;
  b = rdtsc_beg(rd);
  const fint o = zvjsvd_(&m, &n, Gr, &ldGr, Gi, &ldGi, Vr, &ldVr, Vi, &ldVi, eS, fS, js, &stp, &swp, work, iwork);
  e = rdtsc_end(rd);
  (void)fprintf(stdout, "%15.9Lf,%3lld,", tsc_lap(hz, b, e), o);
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

  const int sd = open_wo_(bn, ((j == PJS_ME) ? "S2" : "S4"));
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

  b = rdtsc_beg(rd);
  if (zmerge_(&n, &n, Vr, &ldVr, Vi, &ldVi, V, &ldV) < 0)
    return 11;
  e = rdtsc_end(rd);
  (void)fprintf(stdout, "%15.9Lf,", tsc_lap(hz, b, e));
  (void)fflush(stdout);
  free(Vi);
  free(Vr);

  const int vd = open_wo_(bn, ((j == PJS_ME) ? "V2" : "V4"));
  if (vd < 0)
    return 12;
  *(size_t*)js = (n * (n * sizeof(double complex)));
  if (resizef_(&vd, (const size_t*)js))
    return 12;
  if (zwrite2_(&n, &n, V, &ldV, &vd))
    return 12;
  if (close(vd))
    return 12;
  free(V);

  b = rdtsc_beg(rd);
  if (zmerge_(&m, &n, Gr, &ldGr, Gi, &ldGi, G, &ldG) < 0)
    return 13;
  e = rdtsc_end(rd);
  (void)fprintf(stdout, "%15.9Lf\n", tsc_lap(hz, b, e));
  (void)fflush(stdout);
  free(Gi);
  free(Gr);

  const int ud = open_wo_(bn, ((j == PJS_ME) ? "U2" : "U4"));
  if (ud < 0)
    return 14;
  *(size_t*)js = (m * (n * sizeof(double complex)));
  if (resizef_(&ud, (const size_t*)js))
    return 14;
  if (zwrite2_(&m, &n, G, &ldG, &ud))
    return 14;
  if (close(ud))
    return 14;
  free(G);

  free(js);
  return EXIT_SUCCESS;
}
