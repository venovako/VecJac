#include "dvjsvd.h"

#include "dscale.h"

fint dvjsvd_(const fnat m[static restrict 1], const fnat n[static restrict 1], double G[static restrict VDL], const fnat ldG[static restrict 1], double V[static restrict VDL], const fnat ldV[static restrict 1], double eS[static restrict 1], double fS[static restrict 1], const unsigned js[static restrict 1], const unsigned stp[static restrict 1], const unsigned swp[static restrict 1], double work[static restrict VDL], unsigned iwork[static restrict 1])
{
  if (IS_NOT_VFPENV)
    return -14;
  if (!*n)
    return 0;
  if (*m < *n)
    return -1;
  if (*n & VDL_1)
    return -2;
  if (IS_NOT_ALIGNED(G))
    return -3;
  if (*ldG < *m)
    return -4;
  if (*ldG & VDL_1)
    return -4;
  if (IS_NOT_ALIGNED(V))
    return -5;
  if (*ldV < *n)
    return -6;
  if (*ldV & VDL_1)
    return -6;
  if (IS_NOT_ALIGNED(work))
    return -12;

  const fnat l = 2;
  double M = 0.0;
  fint e = 0;

  double *const a11 = work;
  double *const a22 = a11 + *n;
  double *const a21 = a22 + *n;
  double *const s = a21 + *n;
  double *const t = s + *n;
  double *const c = t + *n;
  double *const l1 = c + *n;
  double *const l2 = l1 + *n;
  unsigned *const p = iwork;

  unsigned sw = 0u;
  while (sw < *swp) {
    size_t swt = 0u;
    for (unsigned st = 0u; st < *stp; ++st) {
      if (dlscal_(m, n, G, ldG, &l, &M, &e) < 0)
        return -15;
      if (!isfinite(M))
        return -3;
      const unsigned *const r = js + st * (size_t)(*n);
      size_t stt = 0u;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,l) reduction(+:stt)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const unsigned _p = r[pq];
        const unsigned _q = r[pq + 1u];
      }
      swt += stt;
    }
    if (!swt)
      break;
    ++sw;
  }

  // TODO

  return (fint)sw;
}
