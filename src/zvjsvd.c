#include "zvjsvd.h"

#include "zscale.h"

fint zvjsvd_(const fnat m[static restrict 1], const fnat n[static restrict 1], double Gr[static restrict VDL], const fnat ldGr[static restrict 1], double Gi[static restrict VDL], const fnat ldGi[static restrict 1], double Vr[static restrict VDL], const fnat ldVr[static restrict 1], double Vi[static restrict VDL], const fnat ldVi[static restrict 1], double eS[static restrict 1], double fS[static restrict 1], const unsigned js[static restrict 1], const unsigned stp[static restrict 1], const unsigned swp[static restrict 1], double work[static restrict VDL], unsigned iwork[static restrict 1])
{
  if (IS_NOT_VFPENV)
    return -18;
  if (!*n)
    return 0;
  if (*m < *n)
    return -1;
  if (*m & VDL_1)
    return -1;
  if (*n & VDL_1)
    return -2;
  if (IS_NOT_ALIGNED(Gr))
    return -3;
  if (*ldGr < *m)
    return -4;
  if (*ldGr & VDL_1)
    return -4;
  if (IS_NOT_ALIGNED(Gi))
    return -5;
  if (*ldGi < *m)
    return -6;
  if (*ldGi & VDL_1)
    return -6;
  if (IS_NOT_ALIGNED(Vr))
    return -7;
  if (*ldVr < *n)
    return -8;
  if (*ldVr & VDL_1)
    return -8;
  if (IS_NOT_ALIGNED(Vi))
    return -9;
  if (*ldVi < *n)
    return -10;
  if (*ldVi & VDL_1)
    return -10;
  if (IS_NOT_ALIGNED(work))
    return -16;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,Vr,ldVr,Vi,ldVi)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    register const VD z = _mm512_setzero_pd();
    double *const Vrj = Vr + j * (size_t)(*ldVr);
    double *const Vij = Vi + j * (size_t)(*ldVi);
    for (fnat i = 0u; i < *n; i += VDL) {
      _mm512_store_pd((Vrj + i), z);
      _mm512_store_pd((Vij + i), z);
    }
    Vrj[j] = 1.0;
  }

  const fnat l = 2;
  double M = 0.0;
  fint e = 0;

  double *const a11 = work;
  double *const a22 = a11 + *n;
  double *const a21r = a22 + *n;
  double *const a21i = a21r + *n;
  double *const s = a21i + *n;
  double *const t = s + *n;
  double *const c = t + *n;
  double *const ca = c + *n;
  double *const sa = ca + *n;
  double *const l1 = sa + *n;
  double *const l2 = l1 + *n;
  unsigned *const p = iwork;

  unsigned sw = 0u;
  while (sw < *swp) {
    size_t swt = 0u;
    for (unsigned st = 0u; st < *stp; ++st) {
      if (zlscal_(m, n, Gr, ldGr, Gi, ldGi, &l, &M, &e) < 0)
        return -19;
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
