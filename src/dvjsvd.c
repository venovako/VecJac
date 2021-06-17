#include "dvjsvd.h"

#include "dscale.h"
#include "dnorme.h"
#include "ddpscl.h"

fint dvjsvd_(const fnat m[static restrict 1], const fnat n[static restrict 1], double G[static restrict VDL], const fnat ldG[static restrict 1], double V[static restrict VDL], const fnat ldV[static restrict 1], double eS[static restrict 1], double fS[static restrict 1], const unsigned js[static restrict 1], const unsigned stp[static restrict 1], const unsigned swp[static restrict 1], double work[static restrict VDL], unsigned iwork[static restrict 1])
{
  const fnat n_2 = (*n >> 1u);
  if (IS_NOT_VFPENV)
    return -14;
  if (!*n)
    return 0;
  if (*m < *n)
    return -1;
  if (*m & VDL_1)
    return -1;
  if (*n & 1u)
    return -2;
  if (n_2 & VDL_1)
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

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,V,ldV)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    register const VD z = _mm512_setzero_pd();
    double *const Vj = V + j * (size_t)(*ldV);
    for (fnat i = 0u; i < *n; i += VDL)
      _mm512_store_pd((Vj + i), z);
    Vj[j] = 1.0;
  }

  double *const a11 = work;
  double *const a22 = a11 + n_2;
  double *const a21 = a22 + n_2;
  double *const s = a21 + n_2;
  double *const t = s + n_2;
  double *const c = t + n_2;
  double *const l1 = c + n_2;
  double *const l2 = l1 + n_2;
  unsigned *const p = iwork;

  // see LAPACK's DGESVJ
  const double tol = sqrt((double)(*m)) * scalbn(DBL_EPSILON, -1);
  const fnat l = 2;
  fint e = 0;
  double M = 0.0;
  unsigned sw = 0u;

  while (sw < *swp) {
    size_t swt = 0u;
    for (unsigned st = 0u; st < *stp; ++st) {
      if ((dlscal_(m, n, G, ldG, &l, &M, &e) < 0) || !isfinite(M))
        return -15;
      const unsigned *const r = js + st * (size_t)(*n);
      size_t stt = 0u;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,G,ldG,eS,fS,a11,a22,a21,s,t,c,l1,l2,p,tol) reduction(+:stt)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat pq_ = pq + 1u;
        const fnat _pq = (pq >> 1u);
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        double *const Gp = G + _p * (*ldG);
        double *const Gq = G + _q * (*ldG);
        if (dnorme_(m, Gp, (eS + _p), (fS + _p), (s + _p), (c + _p)) < 0.0)
          exit(EXIT_FAILURE);
        if (dnorme_(m, Gq, (eS + _q), (fS + _q), (s + _q), (c + _q)) < 0.0)
          exit(EXIT_FAILURE);
        s[pq] = eS[_p];
        c[pq] = fS[_p];
        s[pq_] = eS[_q];
        c[pq_] = fS[_q];
        a21[_pq] = ddpscl_(m, Gp, Gq, (s + pq), (c + pq));
        // repack data
        a11[_pq] = fS[_p];
        a22[_pq] = fS[_q];
        l1[_pq] = eS[_p];
        l2[_pq] = eS[_q];
        if (fabs(a21[_pq]) > tol)
          ++stt;
      }
      if (stt) {
        stt = 0u;
        // TODO
        swt += stt;
      }
    }
    if (!swt)
      break;
    ++sw;
  }

  // TODO: normalize U

  return (fint)sw;
}
