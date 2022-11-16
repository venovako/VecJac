#include "dvjsvd.h"

#include "dnormx.h"
#include "dscale.h"
#include "dnorm2.h"
#include "dznrm2.h"
#include "ddpscl.h"
#include "dgsscl.h"
#include "dbjac2.h"
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

#ifdef SWAP_EFS
#error SWAP_EFS already defined
#else /* !SWAP_EFS */
#define SWAP_EFS(t) \
  t = eS[_p];       \
  eS[_p] = eS[_q];  \
  eS[_q] = t;       \
  t = fS[_p];       \
  fS[_p] = fS[_q];  \
  fS[_q] = t
#endif /* ?SWAP_EFS */

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
    return -__LINE__;
  (void)fprintf(jtr, "M=");
  (void)fflush(jtr);
#endif /* JTRACE */

  double MG = dnormx_(m, n, G, ldG);
  if (!(MG <= DBL_MAX))
    return -__LINE__;
  if (copysign(1.0, MG) == -1.0)
    return -__LINE__;
#ifdef JTRACE
  (void)fprintf(jtr, "%#.17e\n", MG);
  (void)fflush(jtr);
#endif /* JTRACE */

  double MV = 1.0;
  if (*iwork) {
    MV = dnormx_(n, n, V, ldV);
    if (!(MV <= DBL_MAX))
      return -__LINE__;
    // a dirty hack to avoid accumulation on a zero matrix
    if (MV <= 0.0)
      return -13;
  }
  else { // V = I
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
  }

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,eS,fS)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    eS[j] = -HUGE_VAL;
    fS[j] = 1.0;
  }

  if (MG == 0.0)
    return 0;
  const double M_m = (DBL_MAX / (*m << 1u));
  double es = 0.0, fs = 0.0;
  dbl2ef(M_m, &es, &fs);
  const int DBL_MAX_NRM_EXP = (int)es;
  dbl2ef(MG, &es, &fs);
  int eM = (int)es;
  int sR = DBL_MAX_ROT_EXP - eM;
  fint sN = DBL_MAX_NRM_EXP - eM - 1;
#ifdef JTRACE
  (void)fprintf(jtr, "eM=%d, sR=%d, sN=%d, M=", eM, sR, (int)sN);
  (void)fflush(jtr);
#endif /* JTRACE */
  if (sN) {
    if (dscale_(m, n, G, ldG, &sN) < 0)
      return -__LINE__;
    MG = scalbn(MG, (int)sN);
  }
  int sT = (int)sN;
#ifdef JTRACE
  (void)fprintf(jtr, "%#.17e\n", MG);
  (void)fflush(jtr);
#endif /* JTRACE */

  dbl2ef(MV, &es, &fs);
  eM = (int)es;
  sR = DBL_MAX_ROT_EXP - eM;
  sN = sR;
  if (sN) {
    if (dscale_(n, n, V, ldV, &sN) < 0)
      return -__LINE__;
    MV = scalbn(MV, sR);
  }
  int sV = sR;

  const fnat n_16 = (n_2 >> VDLlg);

  double *const a11 = work;
  double *const a22 = a11 + n_2;
  double *const a21 = a22 + n_2;
  double *const c = a21 + n_2;
  double *const at = c + n_2;
  double *const l1 = at + n_2;
  double *const l2 = l1 + n_2;
  double *const w0 = l2 + n_2;
  double *const w1 = w0 + n_2;
  double *const w2 = w1 + n_2;
  unsigned *const p = iwork;
  unsigned *const pc = p + n_16;

  if (*swp) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(w1,n)
#endif /* _OPENMP */
    for (fnat i = 0u; i < *n; ++i)
      w1[i] = 1.0;
  }

  // see LAPACK's DGESVJ
  const double tol = sqrt((double)(*m)) * scalbn(DBL_EPSILON, -1);
  const double gst = scalb(tol, DBL_MAX_FIN_EXP);
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
      // rescale G according to MG if necessary and update MG
      dbl2ef(MG, &es, &fs);
      eM = (int)es;
      sR = DBL_MAX_ROT_EXP - eM;
      sN = DBL_MAX_NRM_EXP - eM - 1;
      if (sR < 0) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u, eM=%d, sR=%d, sN=%d, M=", sw, st, eM, sR, (int)sN);
        (void)fflush(jtr);
#endif /* JTRACE */
        if (dscale_(m, n, G, ldG, &sN) < 0)
          return -__LINE__;
        sR = (int)sN;
        MG = scalbn(MG, sR);
        sT += sR;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(w1,n)
#endif /* _OPENMP */
        for (fnat i = 0u; i < *n; ++i)
          w1[i] = 1.0;
#ifdef JTRACE
        (void)fprintf(jtr, "%#.17e\n", MG);
        (void)fflush(jtr);
#endif /* JTRACE */
      }
      // rescale V according to MV if necessary and update MV
      dbl2ef(MV, &es, &fs);
      eM = (int)es;
      sR = DBL_MAX_ROT_EXP - eM;
      sN = sR;
      if (sR < 0) {
        if (dscale_(n, n, V, ldV, &sN) < 0)
          return -__LINE__;
        MV = scalbn(MV, sR);
        sV += sR;
      }
      // compute the norms, overflow-aware
      const unsigned *const r = js + st * (size_t)(*n);
      double nMG = -0.0;
      bool overflow = false;
      do {
#ifdef JTRACE
        T = rdtsc_beg(rd);
#endif /* JTRACE */
        nMG = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,G,ldG,eS,fS,a11,a22,c,at,l1,l2,w1) reduction(max:nMG)
#endif /* _OPENMP */
        for (fnat pq = 0u; pq < *n; pq += 2u) {
          const fnat _pq = (pq >> 1u);
          if (!(nMG <= DBL_MAX)) {
            a11[_pq] = NAN;
            a22[_pq] = NAN;
            continue;
          }
          const fnat pq_ = pq + 1u;
          const size_t _p = r[pq];
          const size_t _q = r[pq_];
          if (w1[_p] == 1.0) {
            double *const Gp = G + _p * (size_t)(*ldG);
            nMG = fmax(nMG, fmin((a11[_pq] = dnorm2_(m, Gp, (eS + _p), (fS + _p), (c + _pq), (at + _pq))), HUGE_VAL));
            if (!(nMG <= DBL_MAX)) {
              a22[_pq] = NAN;
              continue;
            }
          }
          else
            a11[_pq] = scalb(fS[_p], eS[_p]);
          if (w1[_q] == 1.0) {
            double *const Gq = G + _q * (size_t)(*ldG);
            nMG = fmax(nMG, fmin((a22[_pq] = dnorm2_(m, Gq, (eS + _q), (fS + _q), (l1 + _pq), (l2 + _pq))), HUGE_VAL));
          }
          else
            a22[_pq] = scalb(fS[_q], eS[_q]);
        }
#ifdef JTRACE
        Tn += tsc_lap(hz, T, rdtsc_end(rd));
#endif /* JTRACE */
        if (overflow = !(nMG <= DBL_MAX)) {
#ifdef JTRACE
          (void)fprintf(jtr, "sweep=%u, step=%u, M=", sw, st);
          (void)fflush(jtr);
#endif /* JTRACE */
          if (dscale_(m, n, G, ldG, &sN) < 0)
            return -__LINE__;
          sR = (int)sN;
          MG = scalbn(MG, sR);
          sT += sR;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(w1,n)
#endif /* _OPENMP */
          for (fnat i = 0u; i < *n; ++i)
            w1[i] = 1.0;
#ifdef JTRACE
          (void)fprintf(jtr, "%#.17e\n", MG);
          (void)fflush(jtr);
#endif /* JTRACE */
        }
      } while (overflow);
      // scaled dot-products
#ifdef JTRACE
      T = rdtsc_beg(rd);
#endif /* JTRACE */
      nMG = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,G,ldG,eS,fS,a21) reduction(min:nMG)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat _pq = (pq >> 1u);
        if (!(nMG >= 0.0)) {
          a21[_pq] = NAN;
          continue;
        }
        const fnat pq_ = pq + 1u;
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        // pack the norms
        const double e2[2u] = { eS[_q], eS[_p] };
        const double f2[2u] = { fS[_q], fS[_p] };
        double *const Gp = G + _p * (size_t)(*ldG);
        double *const Gq = G + _q * (size_t)(*ldG);
        a21[_pq] = ddpscl_(m, Gq, Gp, e2, f2);
        nMG = fmin(nMG, (isfinite(a21[_pq]) ? 0.0 : (double)-__LINE__));
      }
#ifdef JTRACE
      Tp += tsc_lap(hz, T, rdtsc_end(rd));
#endif /* JTRACE */
      if (!(nMG >= 0.0)) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u\n", sw, st);
        (void)fflush(jtr);
#endif /* JTRACE */
        return (fint)nMG;
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
#pragma omp parallel for default(none) shared(n_2,a21,c,at,l1,l2,w0,w1,w2,p,pc,tol,gst) reduction(+:stt)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; i += VDL) {
        const fnat j = (i >> VDLlg);
        // convergence check
        register VD _a21 = _mm512_load_pd(a21 + i);
        register const VD _zero = _mm512_set1_pd(-0.0);
        register const VD zero = _mm512_setzero_pd();
        register const VD _tol = _mm512_set1_pd(tol);
        register const VD _a21_ = VDABS(_a21);
        pc[j] = MD2U(_mm512_cmple_pd_mask(_tol, _a21_));
        if (p[j] = _mm_popcnt_u32(pc[j])) {
          stt += p[j];
          register const VD f1 = _mm512_load_pd(l1 + i);
          register const VD f2 = _mm512_load_pd(l2 + i);
          register const VD e1 = _mm512_load_pd(c + i);
          register const VD e2 = _mm512_load_pd(at + i);
          register const VD _gst = _mm512_set1_pd(gst);
          // might not yet be sorted, so check both cases
          register const MD ngs = VDEFLT(e1,e2,f1,f2);
          register VD maf = _mm512_mask_blend_pd(ngs, f2, f1);
          register VD mae = _mm512_mask_blend_pd(ngs, e2, e1);
          register const VD Maf = _mm512_mask_blend_pd(ngs, f1, f2);
          register const VD Mae = _mm512_mask_blend_pd(ngs, e1, e2);
          maf = _mm512_mul_pd(_gst, maf);
          mae = _mm512_add_pd(mae, _mm512_getexp_pd(maf));
          maf = VDMANT(maf);
          register const MD cgs = VDEFLT(mae,Mae,maf,Maf);
          const unsigned gsp = (MD2U(MDANDN(ngs,cgs)) << VDL);
          if (gsp) {
            p[j] |= gsp;
            stt -= _mm_popcnt_u32(gsp);
          }
          const unsigned gsn = (MD2U(MDAND(ngs,cgs)) << VDL);
          if (gsn) {
            pc[j] |= gsn;
            stt -= _mm_popcnt_u32(gsn);
          }
          // Grammian pre-scaling into the double precision range
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
          _mm512_store_pd((w0 + i), _a11);
          _mm512_store_pd((w1 + i), _a22);
          _mm512_store_pd((w2 + i), _a21);
        }
      }
      swt += stt;
#ifdef JTRACE
      Ta += tsc_lap(hz, T, rdtsc_end(rd));
      T = rdtsc_beg(rd);
#endif /* JTRACE */
#ifdef USE_SECANTS
      const fint _n_2 = -(fint)n_2;
#else /* !USE_SECANTS */
      const fint _n_2 = (fint)n_2;
#endif /* ?USE_SECANTS */
      if (dbjac2i(&_n_2, w0, w1, w2, c, at, l1, l2, p) < 0)
        return -__LINE__;
#ifdef JTRACE
      Te += tsc_lap(hz, T, rdtsc_end(rd));
      T = rdtsc_beg(rd);
#endif /* JTRACE */
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(eS,fS,l1,l2,w0,p,pc,r,n_2)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; i += VDL) {
        const fnat j = (i >> VDLlg);
        unsigned gsp = ((p[j] & 0xFFFFFF00u) >> VDL);
        unsigned gsn = ((pc[j] & 0xFFFFFF00u) >> VDL);
        unsigned trans = (pc[j] & 0xFFu);
        unsigned perm = (p[j] & 0xFFu);
        for (fnat k = 0u; k < VDL; ++k) {
          const fnat l = (i + k);
          const fnat pq = (l << 1u);
          const uint64_t _p = r[pq];
          const uint64_t _q = r[pq + 1u];
          *(uint64_t*)(l1 + l) = _p;
          *(uint64_t*)(l2 + l) = _q;
          if (trans & 1u) {
            if (gsp & 1u)
              w0[l] = 3.0;
            else if (gsn & 1u)
              w0[l] = -3.0;
            else if (perm & 1u)
              w0[l] = -2.0;
            else // no swap
              w0[l] = 2.0;
          }
          else if (efcmp((eS + _p), (fS + _p), (eS + _q), (fS + _q)) < 0) {
            SWAP_EFS(w0[l]);
            w0[l] = -1.0;
          }
          else // no swap
            w0[l] = 1.0;
          gsp >>= 1u;
          gsn >>= 1u;
          trans >>= 1u;
          perm >>= 1u;
        }
      }
      nMG = 0.0;
      double nMV = -0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,G,ldG,V,ldV,eS,fS,a21,c,at,l1,l2,w0,w1,n_2) reduction(max:nMG,nMV)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; ++i) {
        const size_t _p = *(const uint64_t*)(l1 + i);
        const size_t _q = *(const uint64_t*)(l2 + i);
        w1[_q] = w1[_p] = 0.0;
        if (!(nMG <= DBL_MAX)) {
          nMG = HUGE_VAL;
          continue;
        }
        if (!(nMV <= DBL_MAX)) {
          nMV = HUGE_VAL;
          continue;
        }
        double *const G_p = G + _p * (size_t)(*ldG);
        double *const G_q = G + _q * (size_t)(*ldG);
        double *const V_p = V + _p * (size_t)(*ldV);
        double *const V_q = V + _q * (size_t)(*ldV);
        if (w0[i] == -3.0) {
          const fint _m = -(fint)*m;
          const double e2[2u] = { eS[_p], eS[_q] };
          const double f2[2u] = { fS[_p], fS[_q] };
          double tG = dgsscl_(&_m, (a21 + i), G_p, G_q, e2, f2);
          if (!isfinite(tG)) {
            nMG = HUGE_VAL;
            continue;
          }
          nMG = fmax(nMG, tG);
          if (dswp_(n, V_p, V_q)) {
            nMV = HUGE_VAL;
            continue;
          }
          nMV = fmax(nMV, 0.0);
          SWAP_EFS(tG);
          w1[_q] = 1.0;
          continue;
        }
        else if (w0[i] == -2.0) {
          const fint _m = -(fint)*m;
          const fint _n = -(fint)*n;
          const double *const _c = (c + i);
          const double *const _at = (at + i);
          const double tG = djrot_(&_m, G_p, G_q, _c, _at);
          if (!isfinite(tG)) {
            nMG = HUGE_VAL;
            continue;
          }
          nMG = fmax(nMG, tG);
          const double tV = djrot_(&_n, V_p, V_q, _c, _at);
          if (!isfinite(tV)) {
            nMV = HUGE_VAL;
            continue;
          }
          nMV = fmax(nMV, tV);
        }
        else if (w0[i] == -1.0) {
          if (dswp_(m, G_p, G_q)) {
            nMG = HUGE_VAL;
            continue;
          }
          nMG = fmax(nMG, 0.0);
          if (dswp_(n, V_p, V_q)) {
            nMV = HUGE_VAL;
            continue;
          }
          nMV = fmax(nMV, 0.0);
          continue;
        }
        else if (w0[i] == 1.0) {
          nMG = fmax(nMG, 0.0);
          nMV = fmax(nMV, 0.0);
          continue;
        }
        else if (w0[i] == 2.0) {
          const fint _m = (fint)*m;
          const fint _n = (fint)*n;
          const double *const _c = (c + i);
          const double *const _at = (at + i);
          const double tG = djrot_(&_m, G_p, G_q, _c, _at);
          if (!isfinite(tG)) {
            nMG = HUGE_VAL;
            continue;
          }
          nMG = fmax(nMG, tG);
          const double tV = djrot_(&_n, V_p, V_q, _c, _at);
          if (!isfinite(tV)) {
            nMV = HUGE_VAL;
            continue;
          }
          nMV = fmax(nMV, tV);
        }
        else if (w0[i] == 3.0) {
          const fint _m = (fint)*m;
          const double e2[2u] = { eS[_p], eS[_q] };
          const double f2[2u] = { fS[_p], fS[_q] };
          double tG = dgsscl_(&_m, (a21 + i), G_p, G_q, e2, f2);
          if (!isfinite(tG)) {
            nMG = HUGE_VAL;
            continue;
          }
          nMG = fmax(nMG, tG);
          nMV = fmax(nMV, 0.0);
          w1[_q] = 1.0;
          continue;
        }
        else { // should never happen
          nMG = HUGE_VAL;
          nMV = HUGE_VAL;
          continue;
        }
        w1[_q] = w1[_p] = 1.0;
      }
      if (!(nMG <= DBL_MAX)) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u\n", sw, st);
        (void)fflush(jtr);
#endif /* JTRACE */
        return -__LINE__;
      }
      MG = fmax(MG, nMG);
      if (!(nMV <= DBL_MAX)) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u\n", sw, st);
        (void)fflush(jtr);
#endif /* JTRACE */
        return -__LINE__;
      }
      MV = fmax(MV, nMV);
#ifdef JTRACE
      Tr += tsc_lap(hz, T, rdtsc_end(rd));
#endif /* JTRACE */
    }
#ifdef JTRACE
    (void)fprintf(jtr, "%3u %zu\n", sw, swt);
    (void)fflush(jtr);
#endif /* JTRACE */
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
    if (sV) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,V,ldV,sV)
#endif /* _OPENMP */
      for (fnat j = 0u; j < *n; ++j) {
        double *const Vj = V + j * (size_t)(*ldV);
        register const VD _s = _mm512_set1_pd((double)-sV);
        for (fnat i = 0u; i < *n; i += VDL) {
          double *const Vij = Vj + i;
          _mm512_store_pd(Vij, _mm512_scalef_pd(_mm512_load_pd(Vij), _s));
        }
      }
    }
  }

#ifdef JTRACE
  (void)fprintf(jtr, "sT=%d, M=%#.17e\n", sT, MG);
  (void)fprintf(jtr, "Tn=%15.9Lf, Tp=%15.9Lf, Ta=%15.9Lf, Te=%15.9Lf, Tr=%15.9Lf\n", Tn, Tp, Ta, Te, Tr);
  (void)fclose(jtr);
#endif /* JTRACE */
  return (fint)sw;
}
