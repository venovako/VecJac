#include "zvjsvd.h"

#include "zscale.h"
#include "znorme.h"
#include "zdpscl.h"

fint zvjsvd_(const fnat m[static restrict 1], const fnat n[static restrict 1], double Gr[static restrict VDL], const fnat ldGr[static restrict 1], double Gi[static restrict VDL], const fnat ldGi[static restrict 1], double Vr[static restrict VDL], const fnat ldVr[static restrict 1], double Vi[static restrict VDL], const fnat ldVi[static restrict 1], double eS[static restrict 1], double fS[static restrict 1], const unsigned js[static restrict 1], const unsigned stp[static restrict 1], const unsigned swp[static restrict 1], double work[static restrict VDL], unsigned iwork[static restrict 1])
{
  const fnat n_2 = (*n >> 1u);
  if (IS_NOT_VFPENV)
    return -18;
  if (!*n)
    return 0;
  if (*m < *n)
    return -1;
  if (*m & VDL_1)
    return -1;
  if (n_2 & VDL_1)
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

  double *const a11 = work;
  double *const a22 = a11 + n_2;
  double *const a21r = a22 + n_2;
  double *const a21i = a21r + n_2;
  double *const s = a21i + n_2;
  double *const t = s + n_2;
  double *const c = t + n_2;
  double *const ca = c + n_2;
  double *const sa = ca + n_2;
  double *const l1 = sa + n_2;
  double *const l2 = l1 + n_2;
  unsigned *const p = iwork;

  // see LAPACK's ZGESVJ
  const double tol = sqrt((double)(*m)) * scalbn(DBL_EPSILON, -1);
  const fnat l = 2;
  fint e = 0;
  double M = 0.0;
  unsigned sw = 0u;

  while (sw < *swp) {
    size_t swt = 0u;
    for (unsigned st = 0u; st < *stp; ++st) {
      if ((zlscal_(m, n, Gr, ldGr, Gi, ldGi, &l, &M, &e) < 0) || !isfinite(M))
        return -19;
      const unsigned *const r = js + st * (size_t)(*n);
      size_t stt = 0u;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,Gr,ldGr,Gi,ldGi,eS,fS,a11,a22,a21r,a21i,s,t,c,ca,sa,l1,l2,p,tol) reduction(+:stt)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat pq_ = pq + 1u;
        const fnat _pq = (pq >> 1u);
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        double *const Grp = Gr + _p * (*ldGr);
        double *const Gip = Gi + _p * (*ldGi);
        if ((a11[_p] = znorme_(m, Grp, Gip, (eS + _p), (fS + _p), (s + _p), (c + _p))) < 0.0)
          exit(EXIT_FAILURE);
        double *const Grq = Gr + _q * (*ldGr);
        double *const Giq = Gi + _q * (*ldGi);
        if ((a11[_q] = znorme_(m, Grq, Giq, (eS + _q), (fS + _q), (s + _q), (c + _q))) < 0.0)
          exit(EXIT_FAILURE);
        s[pq] = eS[_p];
        c[pq] = fS[_p];
        s[pq_] = eS[_q];
        c[pq_] = fS[_q];
        const double complex z = zdpscl_(m, Grp, Gip, Grq, Giq, (s + pq), (c + pq));
        if (cabs(z) > tol) {
          // do not increment if the rotation turns out to be identity
          ++stt;
        }
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
