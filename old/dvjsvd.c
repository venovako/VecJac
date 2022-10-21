#include "dvjsvd.h"

#include "dnormx.h"
#include "dscale.h"
#include "dnorm2.h"
#include "dznrm2.h"
#include "ddpscl.h"
#include "dbjac2.h"
#include "djrotf.h"
#include "djrot.h"
#include "dswp.h"
#include "vecdef.h"
#include "defops.h"

#ifdef JTRACE
#include "timer.h"
#endif /* JTRACE */

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

  if (M == 0.0)
    return 0;
  const double M_m = (DBL_MAX / (*m << 1u));
  double es = 0.0, fs = 0.0;
  dbl2ef(M_m, &es, &fs);
  const int DBL_MAX_NRM_EXP = (int)es;
  dbl2ef(M, &es, &fs);
  int eM = (int)es;
  int sR = DBL_MAX_ROT_EXP - eM;
  fint sN = DBL_MAX_NRM_EXP - eM - 1;
#ifdef JTRACE
  (void)fprintf(jtr, "eM=%d, sR=%d, sN=%d, M=", eM, sR, (int)sN);
  (void)fflush(jtr);
#endif /* JTRACE */
  if (sN) {
    if (dscale_(m, n, G, ldG, &sN) < 0)
      return -17;
    M = scalbn(M, (int)sN);
  }
  int sT = (int)sN;
#ifdef JTRACE
  (void)fprintf(jtr, "%#.17e\n", M);
  (void)fflush(jtr);
#endif /* JTRACE */

  const fnat n_16 = (n_2 >> VDLlg);

  double *const a11 = work;
  double *const a22 = a11 + n_2;
  double *const a21 = a22 + n_2;
  double *const c = a21 + n_2;
  double *const at = c + n_2;
  double *const l1 = at + n_2;
  double *const l2 = l1 + n_2;
  double *const w = l2 + n_2;
  unsigned *const p = iwork;
  unsigned *const pc = p + n_16;

  if (*swp) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(l1,n)
#endif /* _OPENMP */
    for (fnat i = 0u; i < *n; ++i)
      l1[i] = 1.0;
  }

  // see LAPACK's DGESVJ
  const double tol = sqrt((double)(*m)) * scalbn(DBL_EPSILON, -1);
  unsigned sw = 0u;

#ifdef JTRACE
  unsigned rd[2u] = { 0u, 0u };
  const uint64_t hz = tsc_get_freq_hz_(rd);

  long double Tn = 0.0L, Tp = 0.0L, Ta = 0.0L, Te = 0.0L, Tr = 0.0L;
  uint64_t T = UINT64_C(0);
#endif /* JTRACE */

  while (sw < *swp) {
    size_t swt = 0u;
    for (unsigned st = 0u; st < *stp; ++st) {
      // rescale according to M if necessary and update M
      dbl2ef(M, &es, &fs);
      eM = (int)es;
      sR = DBL_MAX_ROT_EXP - eM;
      sN = DBL_MAX_NRM_EXP - eM - 1;
      if (sR < 0) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u, eM=%d, sR=%d, sN=%d, M=", sw, st, eM, sR, (int)sN);
        (void)fflush(jtr);
#endif /* JTRACE */
        if (dscale_(m, n, G, ldG, &sN) < 0)
          return -18;
        M = scalbn(M, (int)sN);
        sT += (int)sN;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(l1,n)
#endif /* _OPENMP */
        for (fnat i = 0u; i < *n; ++i)
          l1[i] = 1.0;
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
#ifdef JTRACE
        T = rdtsc_beg(rd);
#endif /* JTRACE */
        nM = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,G,ldG,eS,fS,a11,a22,c,at,l1) reduction(max:nM)
#endif /* _OPENMP */
        for (fnat pq = 0u; pq < *n; pq += 2u) {
          const fnat _pq = (pq >> 1u);
          if (!(nM <= DBL_MAX)) {
            a11[_pq] = NAN;
            a22[_pq] = NAN;
            continue;
          }
          const fnat pq_ = pq + 1u;
          const size_t _p = r[pq];
          const size_t _q = r[pq_];
          if (l1[_p] == 1.0) {
            double *const Gp = G + _p * (*ldG);
            nM = fmax(nM, fmin((a11[_pq] = dnorm2_(m, Gp, (eS + _p), (fS + _p), (c + _pq), (at + _pq))), HUGE_VAL));
            if (!(nM <= DBL_MAX)) {
              a22[_pq] = NAN;
              continue;
            }
          }
          if (l1[_q] == 1.0) {
            double *const Gq = G + _q * (*ldG);
            nM = fmax(nM, fmin((a22[_pq] = dnorm2_(m, Gq, (eS + _q), (fS + _q), (c + _pq), (at + _pq))), HUGE_VAL));
          }
        }
#ifdef JTRACE
        Tn += tsc_lap(hz, T, rdtsc_end(rd));
#endif /* JTRACE */
        if (overflow = !(nM <= DBL_MAX)) {
#ifdef JTRACE
          (void)fprintf(jtr, "sweep=%u, step=%u, M=", sw, st);
          (void)fflush(jtr);
#endif /* JTRACE */
          if (dscale_(m, n, G, ldG, &sN) < 0)
            return -19;
          M = scalbn(M, (int)sN);
          sT += (int)sN;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(l1,n)
#endif /* _OPENMP */
          for (fnat i = 0u; i < *n; ++i)
            l1[i] = 1.0;
#ifdef JTRACE
          (void)fprintf(jtr, "%#.17e\n", M);
          (void)fflush(jtr);
#endif /* JTRACE */
        }
      } while (overflow);
      // scaled dot-products
#ifdef JTRACE
      T = rdtsc_beg(rd);
#endif /* JTRACE */
      nM = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,G,ldG,eS,fS,w) reduction(min:nM)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat _pq = (pq >> 1u);
        if (!(nM >= 0.0)) {
          w[_pq] = NAN;
          continue;
        }
        const fnat pq_ = pq + 1u;
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        // pack the norms
        const double e2[2u] = { eS[_q], eS[_p] };
        const double f2[2u] = { fS[_q], fS[_p] };
        double *const Gp = G + _p * (*ldG);
        double *const Gq = G + _q * (*ldG);
        w[_pq] = ddpscl_(m, Gq, Gp, e2, f2);
        if (!(isfinite(w[_pq])))
          nM = fmin(nM, -20.0);
      }
#ifdef JTRACE
      Tp += tsc_lap(hz, T, rdtsc_end(rd));
#endif /* JTRACE */
      if (!(nM >= 0.0)) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u\n", sw, st);
        (void)fflush(jtr);
#endif /* JTRACE */
        return (fint)nM;
      }
      // repack data
#ifdef JTRACE
      T = rdtsc_beg(rd);
#endif /* JTRACE */
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
      fnat stt = 0u;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n_2,a11,a22,a21,c,at,l1,l2,w,p,pc,tol) reduction(+:stt)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; i += VDL) {
        const fnat j = (i >> VDLlg);
        // convergence check
        register VD _a21 = _mm512_load_pd(w + i);
        register const VD _zero = _mm512_set1_pd(-0.0);
        register const VD zero = _mm512_setzero_pd();
        register const VD _tol = _mm512_set1_pd(tol);
        register const VD _a21_ = VDABS(_a21);
        pc[j] = MD2U(_mm512_cmple_pd_mask(_tol, _a21_));
        if (p[j] = _mm_popcnt_u32(pc[j])) {
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
          register const VD mxe = _mm512_set1_pd(DBL_MAX_FIN_EXP);
          register const VD E = _mm512_mask_blend_pd(c12, e12, e21);
          register const VD d = _mm512_min_pd(_mm512_sub_pd(mxe, E), zero);
          e12 = _mm512_add_pd(e12, d);
          e21 = _mm512_add_pd(e21, d);
          register const VD _a11 = _mm512_scalef_pd(f12, e12);
          register const VD _a22 = _mm512_scalef_pd(f21, e21);
          _a21 = _mm512_scalef_pd(_a21, d);
          _mm512_store_pd((a11 + i), _a11);
          _mm512_store_pd((a22 + i), _a22);
          _mm512_store_pd((a21 + i), _a21);
        }
      }
      swt += stt;
#ifdef JTRACE
      Ta += tsc_lap(hz, T, rdtsc_end(rd));
      T = rdtsc_beg(rd);
#endif /* JTRACE */
      const fint _n_2 =
#ifdef USE_SECANTS
        -(fint)n_2
#else /* !USE_SECANTS */
        (fint)n_2
#endif /* ?USE_SECANTS */
        ;
      if (dbjac2i(&_n_2, a11, a22, a21, c, at, l1, l2, p) < 0)
        return -21;
#ifdef JTRACE
      Te += tsc_lap(hz, T, rdtsc_end(rd));
      T = rdtsc_beg(rd);
#endif /* JTRACE */
      fnat np = 0u; // number of swaps
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(a11,a22,a21,eS,fS,p,pc,r,n_2) reduction(+:np)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; i += VDL) {
        const fnat j = (i >> VDLlg);
        unsigned trans = (pc[j] & 0xFFu);
        unsigned perm = (p[j] & 0xFFu);
        for (fnat k = 0u; k < VDL; ++k) {
          const fnat l = (i + k);
          const fnat pq = (l << 1u);
          const uint64_t _p = r[pq];
          const uint64_t _q = r[pq + 1u];
          *(uint64_t*)(a11 + l) = _p;
          *(uint64_t*)(a22 + l) = _q;
          if (trans & 1u) {
            if (perm & 1u) {
              a21[l] = -2.0;
              ++np;
            }
            else // no swap
              a21[l] = 2.0;
          }
          else if (efcmp((eS + _p), (fS + _p), (eS + _q), (fS + _q)) < 0) {
            a21[l] = eS[_p];
            eS[_p] = eS[_q];
            eS[_q] = a21[l];
            a21[l] = fS[_p];
            fS[_p] = fS[_q];
            fS[_q] = a21[l];
            a21[l] = -1.0;
            ++np;
          }
          else // no swap
            a21[l] = 1.0;
          trans >>= 1u;
          perm >>= 1u;
        }
      }
      nM = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,G,ldG,V,ldV,a11,a22,a21,c,at,l1,w,eS,fS,n_2) reduction(max:nM)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; ++i) {
        const size_t _p = *(const uint64_t*)(a11 + i);
        const size_t _q = *(const uint64_t*)(a22 + i);
        l1[_q] = l1[_p] = 0.0;
        if (!(nM <= DBL_MAX)) {
          w[i] = NAN;
          continue;
        }
        double _at, _c;
        fint _m, _n;
        if (a21[i] == -2.0) {
          _m = -(fint)*m;
          _n = -(fint)*n;
          _c = c[i];
          _at = at[i];
        }
        else if (a21[i] == -1.0) {
          double *const Gp = G + _p * (*ldG);
          double *const Gq = G + _q * (*ldG);
          if (_m = dswp_(m, Gp, Gq)) {
            w[i] = (double)_m;
            nM = HUGE_VAL;
            continue;
          }
          double *const Vp = V + _p * (*ldV);
          double *const Vq = V + _q * (*ldV);
          if (_n = dswp_(n, Vp, Vq)) {
            w[i] = (double)_n;
            nM = HUGE_VAL;
            continue;
          }
          nM = fmax(nM, (w[i] = 0.0));
          continue;
        }
        else if (a21[i] == 1.0) {
          nM = fmax(nM, (w[i] = 0.0));
          continue;
        }
        else if (a21[i] == 2.0) {
          _m = (fint)*m;
          _n = (fint)*n;
          _c = c[i];
          _at = at[i];
        }
        else { // should never happen
          w[i] = NAN;
          nM = HUGE_VAL;
          continue;
        }
        w[i] = djrot_(&_m, (G + _p * (*ldG)), (G + _q * (*ldG)), &_c, &_at);
        if (!(w[i] >= 0.0) || !(w[i] <= DBL_MAX)) {
          nM = w[i] = HUGE_VAL;
          continue;
        }
        else // no overflow
          nM = fmax(nM, w[i]);
        if (_m = djrotf_(&_n, (V + _p * (*ldV)), (V + _q * (*ldV)), &_c, &_at)) {
          w[i] = (double)_m;
          nM = HUGE_VAL;
          continue;
        }
        l1[_q] = l1[_p] = 1.0;
      }
      M = fmax(M, nM);
#ifdef JTRACE
      Tr += tsc_lap(hz, T, rdtsc_end(rd));
#endif /* JTRACE */
      if (!(M <= DBL_MAX)) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u\n", sw, st);
        (void)fflush(jtr);
#endif /* JTRACE */
        return -22;
      }
    }
    if (!swt)
      break;
    ++sw;
  }

  if (sw < *swp) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,G,ldG,eS,fS,sT)
#endif /* _OPENMP */
    for (fnat j = 0u; j < *n; ++j) {
      double *const Gj = G + j * (size_t)(*ldG);
      register const VD _f = _mm512_set1_pd(fS[j]);
      register const VD _s = _mm512_set1_pd(-(eS[j]));
      for (fnat i = 0u; i < *m; i += VDL) {
        double *const Gij = Gj + i;
        _mm512_store_pd(Gij, _mm512_scalef_pd(_mm512_div_pd(_mm512_load_pd(Gij), _f), _s));
      }
      eS[j] -= sT;
    }
  }

#ifdef JTRACE
  (void)fprintf(jtr, "sT=%d, M=%#.17e\n", sT, M);
  (void)fprintf(jtr, "Tn=%15.9Lf, Tp=%15.9Lf, Ta=%15.9Lf, Te=%15.9Lf, Tr=%15.9Lf\n", Tn, Tp, Ta, Te, Tr);
  (void)fclose(jtr);
#endif /* JTRACE */
  return (fint)sw;
}
