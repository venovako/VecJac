#include "dvjsvd.h"

#include "dnormx.h"
#include "dscale.h"
#include "dnorm2.h"
#include "dznrm2.h"
#include "ddpscl.h"
#include "dbjac2.h"
#include "djrot.h"
#include "vecdef.h"
#include "defops.h"

#ifdef DBL_MAX_ROT_EXP
#error DBL_MAX_ROT_EXP already defined
#else /* !DBL_MAX_ROT_EXP */
#define DBL_MAX_ROT_EXP 1022
#endif /* ?DBL_MAX_ROT_EXP */

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
  if (*n >= 0x10000000u) // 2^28
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

#ifdef JTRACE
  FILE *const jtr = fopen((const char*)work, "w");
  if (!jtr)
    return -13;
  (void)fprintf(jtr, "M=");
  (void)fflush(jtr);
#endif /* JTRACE */

  double M = dnormx_(m, n, G, ldG);
  if (!(M <= DBL_MAX))
    return -15;
  if (copysign(1.0, M) == -1.0)
    return -16;

#ifdef JTRACE
  (void)fprintf(jtr, "%#.17e\n", M);
  (void)fflush(jtr);
#endif /* JTRACE */

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,V,ldV,eS,fS)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    register const VD z = _mm512_setzero_pd();
    double *const Vj = V + j * (size_t)(*ldV);
    for (fnat i = 0u; i < *n; i += VDL)
      _mm512_store_pd((Vj + i), z);
    fS[j] = Vj[j] = 1.0;
    eS[j] = -HUGE_VAL;
  }

  double *const a11 = work;
  double *const a22 = a11 + n_2;
  double *const a21 = a22 + n_2;
  double *const c = a22 + n_2;
  double *const at = c + n_2;
  double *const l1 = at + n_2;
  double *const l2 = l1 + n_2;
  double *const w = l2 + n_2;
  unsigned *const p = iwork;
  unsigned *const pc = p + (n_2 >> VDLlg);

  if (M == 0.0)
    return 0;
  const double M_m = (DBL_MAX / (*m << 1u));
  double es = 0.0, fs = 0.0;
  dbl2ef(M_m, &es, &fs);
  const int DBL_MAX_NRM_EXP = (int)es;
  dbl2ef(M, &es, &fs);
  int eM = (int)es;
  int sR = DBL_MAX_ROT_EXP - eM;
  int sN = DBL_MAX_NRM_EXP - eM - 1;
#ifdef JTRACE
  (void)fprintf(jtr, "eM=%d, sR=%d, sN=%d, M=", eM, sR, sN);
  (void)fflush(jtr);
#endif /* JTRACE */
  if (sN) {
    *(fint*)&es = sN;
    if (dscale_(m, n, G, ldG, (const fint*)&es) < 0)
      return -17;
    M = scalbn(M, sN);
  }
  int sT = sN;
#ifdef JTRACE
  (void)fprintf(jtr, "%#.17e\n", M);
  (void)fflush(jtr);
#endif /* JTRACE */

  // see LAPACK's DGESVJ
  const double tol = sqrt((double)(*m)) * scalbn(DBL_EPSILON, -1);
  unsigned sw = 0u;

  while (sw < *swp) {
#ifdef JTRACE
    (void)fprintf(jtr, "Sweep %u\n", sw);
    (void)fflush(jtr);
#endif /* JTRACE */
    size_t swt = 0u;
    for (unsigned st = 0u; st < *stp; ++st) {
      // rescale according to M if necessary and update M
      dbl2ef(M, &es, &fs);
      eM = (int)es;
      sR = DBL_MAX_ROT_EXP - eM;
      sN = DBL_MAX_NRM_EXP - eM - 1;
      if (sR < 0) {
#ifdef JTRACE
        (void)fprintf(jtr, "step=%u, eM=%d, sR=%d, sN=%d, M=", st, eM, sR, sN);
        (void)fflush(jtr);
#endif /* JTRACE */
        *(fint*)&es = sN;
        if (dscale_(m, n, G, ldG, (const fint*)&es) < 0)
          return -18;
        M = scalbn(M, sN);
        sT += sN;
#ifdef JTRACE
        (void)fprintf(jtr, "%#.17e\n", M);
        (void)fflush(jtr);
#endif /* JTRACE */
      }
      // compute the norms, overflow-aware
      const unsigned *const r = js + st * (size_t)(*n);
      double nM = -0.0;
      bool overflow = false;
      do {
        nM = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,G,ldG,eS,fS,c,at,l1,l2) reduction(max:nM)
#endif /* _OPENMP */
        for (fnat pq = 0u; pq < *n; pq += 2u) {
          const fnat _pq = (pq >> 1u);
          if (nM > DBL_MAX) {
            l1[_pq] = NAN;
            l2[_pq] = NAN;
            continue;
          }
          const fnat pq_ = pq + 1u;
          const size_t _p = r[pq];
          const size_t _q = r[pq_];
          double *const Gp = G + _p * (*ldG);
          nM = fmax(nM, fmin((l1[_pq] = dnorm2_(m, Gp, (eS + _p), (fS + _p), (c + _pq), (at + _pq))), HUGE_VAL));
          if (nM > DBL_MAX) {
            l2[_pq] = NAN;
            continue;
          }
          double *const Gq = G + _q * (*ldG);
          nM = fmax(nM, fmin((l2[_pq] = dnorm2_(m, Gq, (eS + _q), (fS + _q), (c + _pq), (at + _pq))), HUGE_VAL));
        }
        if (overflow = (nM > DBL_MAX)) {
#ifdef JTRACE
          (void)fprintf(jtr, "step=%u, M=", st);
          (void)fflush(jtr);
#endif /* JTRACE */
          *(fint*)&es = sN;
          if (dscale_(m, n, G, ldG, (const fint*)&es) < 0)
            return -19;
          M = scalbn(M, sN);
          sT += sN;
#ifdef JTRACE
          (void)fprintf(jtr, "%#.17e\n", M);
          (void)fflush(jtr);
#endif /* JTRACE */
        }
      } while (overflow);
      // scaled dot-products
      nM = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,G,ldG,eS,fS,w) reduction(min:nM)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat _pq = (pq >> 1u);
        if (nM < 0.0) {
          w[_pq] = NAN;
          continue;
        }
        const fnat pq_ = pq + 1u;
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        // pack the norms
        const double e2[2u] = { eS[_p], eS[_q] };
        const double f2[2u] = { fS[_p], fS[_q] };
        double *const Gp = G + _p * (*ldG);
        double *const Gq = G + _q * (*ldG);
        w[_pq] = ddpscl_(m, Gp, Gq, e2, f2);
        if (!(isfinite(w[_pq])))
          nM = fmin(nM, -20.0);
      }
      if (nM < 0.0)
        return (fint)nM;
      // repack data
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,eS,fS,c,at,l1,l2)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat pq_ = pq + 1u;
        const fnat _pq = (pq >> 1u);
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        c[_pq] = eS[_p];
        at[_pq] = eS[_q];
        l1[_pq] = fS[_p];
        l2[_pq] = fS[_q];
      }
      fnat stt = 0u, k = 0u;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n_2,a11,a22,a21,c,at,w,p,pc,tol,k) reduction(+:stt)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; i += VDL) {
        const fnat j = (i >> VDLlg);
        // convergence check
        register VD _a21 = _mm512_load_pd(w + i);
        register const VD _zero = _mm512_set1_pd(-0.0);
        register const VD _tol = _mm512_set1_pd(tol);
        register const VD _a21_ = VDABS(_a21);
        pc[j] = MD2U(_mm512_cmple_pd_mask(_tol, _a21_));
        if (!(p[j] = _mm_popcnt_u32(pc[j])))
          continue;
        stt += p[j];
        // Grammian pre-scaling into the double precision range
        register const VD f1 = _mm512_load_pd(l1 + i);
        register const VD f2 = _mm512_load_pd(l2 + i);
        register const VD e1 = _mm512_load_pd(c + i);
        register const VD e2 = _mm512_load_pd(at + i);
        register VD f12 = _mm512_div_pd(f1, f2);
        register VD e12 = _mm512_sub_pd(e1, e2);
        register VD f21 = _mm512_div_pd(f2, f1);
        register VD e21 = _mm512_sub_pd(e2, e1);
        e12 = _mm512_add_pd(e12, _mm512_getexp_pd(f12));
        f12 = VDMANT(f12);
        e21 = _mm512_add_pd(e21, _mm512_getexp_pd(f21));
        f21 = VDMANT(f21);
        register const MD c12 = VDEFLE(e12,e21,f12,f21);
        register const VD E = _mm512_mask_blend_pd(c12, e12, e21);
        register const VD d = _mm512_min_pd(_mm512_sub_pd(_mm512_set1_pd(DBL_MAX_FIN_EXP), E), _mm512_setzero_pd());
        e12 = _mm512_add_pd(e12, d);
        e21 = _mm512_add_pd(e21, d);
        register const VD _a11 = _mm512_scalef_pd(f12, e12);
        register const VD _a22 = _mm512_scalef_pd(f21, e21);
        _a21 = _mm512_scalef_pd(_a21, d);
        // pack the data and record the translation in pc
        fnat kk;
#ifdef _OPENMP
#pragma omp atomic capture seq_cst
#endif /* _OPENMP */
        kk = k++;
        // lower 8 bits: mask, upper 24 bits: j (assumes n < 2^28)
        pc[kk] |= (j << VDL);
        kk <<= VDLlg;
        _mm512_store_pd((a11 + kk), _a11);
        _mm512_store_pd((a22 + kk), _a22);
        _mm512_store_pd((a21 + kk), _a21);
      }
      if (!stt)
        continue;
      swt += stt;
      const fnat kk = (k << VDLlg);
      if (dbjac2_(&kk, a11, a22, a21, c, at, l1, l2, p) < 0)
        return -21;
      fnat np = 0u; // number of swaps
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(a11,a22,a21,p,pc,r,k) reduction(+:np)
#endif /* _OPENMP */
      for (fnat i = 0u; i < k; ++i) {
        const fnat i_ = (i << VDLlg);
        const fnat _pq = ((pc[i] >> VDL) << VDLlg);
        for (unsigned b = (pc[i] & 0xFFu), x = p[i], l_ = 0u; b; (b >>= 1u), (x >>= 1u), ++l_) {
          const fnat k_ = (i_ + l_);
          const fnat pq = ((_pq + l_) << 1u);
          const size_t _p = r[pq];
          const size_t _q = r[pq + 1u];
          *(size_t*)(a11 + k_) = _p;
          *(size_t*)(a22 + k_) = _q;
          if (x & 1u) {
            a21[k_] = ((b & 1u) ? -2.0 : -1.0);
            ++np;
          }
          else // no swap
            a21[k_] = ((b & 1u) ? 2.0 : 1.0);
        }
      }
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,G,ldG,V,ldV,a11,a22,a21,c,at,w,kk) reduction(max:M)
#endif /* _OPENMP */
      for (fnat i = 0u; i < kk; ++i) {
        if (M > DBL_MAX) {
          w[i] = NAN;
          continue;
        }
        const size_t _p = *(const size_t*)(a11 + i);
        const size_t _q = *(const size_t*)(a22 + i);
        double _at, _c;
        fint _m, _n;
        bool triv = false;
        if (a21[i] == -2.0) {
          _m = -(fint)*m;
          _n = -(fint)*n;
          _c = c[i];
          _at = at[i];
        }
        else if (a21[i] == -1.0) {
          _m = -(fint)*m;
          _n = -(fint)*n;
          _c = 1.0;
          _at = 0.0;
        }
        else if (a21[i] == 2.0) {
          _m = (fint)*m;
          _n = (fint)*n;
          _c = c[i];
          _at = at[i];
        }
        else // no-op
          triv = true;
        if (triv)
          M = fmax(M, 0.0);
        else {
          w[i] = djrot_(&_m, (G + _p * (*ldG)), (G + _q * (*ldG)), &_c, &_at);
          M = fmax(M, (!(w[i] >= 0.0) ? HUGE_VAL : w[i]));
          if (M > DBL_MAX) {
            w[i] = NAN;
            continue;
          }
          w[i] = djrot_(&_n, (V + _p * (*ldV)), (V + _q * (*ldV)), &_c, &_at);
          // V should not overflow but check anyway
          if (!(w[i] >= 0.0) || !(w[i] <= DBL_MAX))
            M = HUGE_VAL;
        }
      }
      if (M > DBL_MAX)
        return -22;
    }
    if (!swt)
      break;
    ++sw;
  }

  // TODO: normalize U and extract S

#ifdef JTRACE
  (void)fclose(jtr);
#endif /* JTRACE */
  return (fint)sw;
}
