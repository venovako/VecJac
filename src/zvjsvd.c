#include "zvjsvd.h"

#include "znormx.h"
#include "zscale.h"
#include "znorm2.h"
#include "dznrm2.h"
#include "zdpscl.h"
#include "zgsscl.h"
#include "zbjac2.h"
#include "zjrot.h"
#include "dswp.h"
#include "vecdef.h"
#include "defops.h"

#ifdef JTRACE
#include "timer.h"
#endif /* JTRACE */

#ifdef DBL_MAX_ROT_EXP
#error DBL_MAX_ROT_EXP already defined
#else /* !DBL_MAX_ROT_EXP */
#define DBL_MAX_ROT_EXP 1021
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
  if (*n & 1u)
    return -2;
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

#ifdef JTRACE
  FILE *const jtr = fopen((const char*)work, "w");
  if (!jtr)
    return -__LINE__;
  (void)fprintf(jtr, "M=");
  (void)fflush(jtr);
#endif /* JTRACE */

  double MG = znormx_(m, n, Gr, ldGr, Gi, ldGi);
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
    MV = znormx_(n, n, Vr, ldVr, Vi, ldVi);
    if (!(MV <= DBL_MAX))
      return -__LINE__;
    // a dirty hack to avoid accumulation on a zero matrix
    if (MV <= 0.0)
      return -17;
  }
  else { // V = I
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
  const double M_m = (DBL_MAX / ((*m << 2u) * M_SQRT2));
  double es = 0.0, fs = 0.0;
  dbl2ef(M_m, &es, &fs);
  const int DBL_MAX_NRM_EXP = (int)es;
  dbl2ef(MG, &es, &fs);
  int eM = (int)es;
  int sR = DBL_MAX_ROT_EXP - eM - 1;
  fint sN = DBL_MAX_NRM_EXP - eM - 1;
#ifdef JTRACE
  (void)fprintf(jtr, "eM=%d, sR=%d, sN=%d, M=", eM, sR, (int)sN);
  (void)fflush(jtr);
#endif /* JTRACE */
  if (sN) {
    if (zscale_(m, n, Gr, ldGr, Gi, ldGi, &sN) < 0)
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
  sR = DBL_MAX_ROT_EXP - eM - 1;
  sN = sR;
  if (sN) {
    if (zscale_(n, n, Vr, ldVr, Vi, ldVi, &sN) < 0)
      return -__LINE__;
    MV = scalbn(MV, sR);
  }
  int sV = sR;

  const fnat n_16 = (n_2 >> VDLlg);

  double *const a11 = work;
  double *const a22 = a11 + n_2;
  double *const a21r = a22 + n_2;
  double *const a21i = a21r + n_2;
  double *const c = a21i + n_2;
  double *const cat = c + n_2;
  double *const sat = cat + n_2;
  double *const l1 = sat + n_2;
  double *const l2 = l1 + n_2;
  double *const w0 = l2 + n_2;
  double *const w1 = w0 + n_2;
  double *const w2 = w1 + n_2;
  double *const w3 = w2 + n_2;
  double *const w4 = w3 + n_2;
  unsigned *const p = iwork;
  unsigned *const pc = p + n_16;

  if (*swp) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(w1,n)
#endif /* _OPENMP */
    for (fnat i = 0u; i < *n; ++i)
      w1[i] = 1.0;
  }

  // see LAPACK's ZGESVJ
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
      sR = DBL_MAX_ROT_EXP - eM - 1;
      sN = DBL_MAX_NRM_EXP - eM - 1;
      if (sR < 0) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u, eM=%d, sR=%d, sN=%d, M=", sw, st, eM, sR, (int)sN);
        (void)fflush(jtr);
#endif /* JTRACE */
        if (zscale_(m, n, Gr, ldGr, Gi, ldGi, &sN) < 0)
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
      sR = DBL_MAX_ROT_EXP - eM - 1;
      sN = sR;
      if (sN) {
        if (zscale_(n, n, Vr, ldVr, Vi, ldVi, &sN) < 0)
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
#pragma omp parallel for default(none) shared(n,r,m,Gr,ldGr,Gi,ldGi,eS,fS,a11,a22,cat,sat,l1,l2,w1) reduction(max:nMG)
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
            double *const Grp = Gr + _p * (size_t)(*ldGr);
            double *const Gip = Gi + _p * (size_t)(*ldGi);
            nMG = fmax(nMG, fmin((a11[_pq] = znorm2_(m, Grp, Gip, (eS + _p), (fS + _p), (cat + _pq), (sat + _pq))), HUGE_VAL));
            if (!(nMG <= DBL_MAX)) {
              a22[_pq] = NAN;
              continue;
            }
          }
          else
            a11[_pq] = scalb(fS[_p], eS[_p]);
          if (w1[_q] == 1.0) {
            double *const Grq = Gr + _q * (size_t)(*ldGr);
            double *const Giq = Gi + _q * (size_t)(*ldGi);
            nMG = fmax(nMG, fmin((a22[_pq] = znorm2_(m, Grq, Giq, (eS + _q), (fS + _q), (l1 + _pq), (l2 + _pq))), HUGE_VAL));
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
          if (zscale_(m, n, Gr, ldGr, Gi, ldGi, &sN) < 0)
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
#pragma omp parallel for default(none) shared(n,r,m,Gr,ldGr,Gi,ldGi,eS,fS,a21r,a21i) reduction(min:nMG)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat _pq = (pq >> 1u);
        if (!(nMG >= 0.0)) {
          a21r[_pq] = NAN;
          a21i[_pq] = NAN;
          continue;
        }
        const fnat pq_ = pq + 1u;
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        // pack the norms
        const double e2[2u] = { eS[_q], eS[_p] };
        const double f2[2u] = { fS[_q], fS[_p] };
        double *const Grp = Gr + _p * (size_t)(*ldGr);
        double *const Gip = Gi + _p * (size_t)(*ldGi);
        double *const Grq = Gr + _q * (size_t)(*ldGr);
        double *const Giq = Gi + _q * (size_t)(*ldGi);
        const double complex z = zdpscl_(m, Grq, Giq, Grp, Gip, e2, f2);
        a21r[_pq] = creal(z);
        nMG = fmin(nMG, (isfinite(a21r[_pq]) ? 0.0 : (double)-__LINE__));
        a21i[_pq] = cimag(z);
        nMG = fmin(nMG, (isfinite(a21i[_pq]) ? 0.0 : (double)-__LINE__));
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
#pragma omp parallel for default(none) shared(n,r,eS,fS,cat,sat,l1,l2)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat pq_ = pq + 1u;
        const fnat _pq = (pq >> 1u);
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        cat[_pq] = eS[_p];
        sat[_pq] = eS[_q];
        l1[_pq] = fS[_p];
        l2[_pq] = fS[_q];
      }
      fnat stt = 0u;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n_2,a21r,a21i,cat,sat,l1,l2,w1,w2,w3,w4,p,pc,tol,gst) reduction(+:stt)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; i += VDL) {
        const fnat j = (i >> VDLlg);
        // convergence check
        register VD _a21r = _mm512_load_pd(a21r + i);
        register VD _a21i = _mm512_load_pd(a21i + i);
        register const VD _zero = _mm512_set1_pd(-0.0);
        register const VD zero = _mm512_setzero_pd();
        register const VD one = _mm512_set1_pd(1.0);
        register const VD _tol = _mm512_set1_pd(tol);
        register VD _a21_ /* = _mm512_hypot_pd(_a21r, _a21i) */;
        VDHYPOT(_a21_, _a21r, _a21i);
        pc[j] = MD2U(_mm512_cmple_pd_mask(_tol, _a21_));
        if (p[j] = _mm_popcnt_u32(pc[j])) {
          stt += p[j];
          register const VD f1 = _mm512_load_pd(l1 + i);
          register const VD f2 = _mm512_load_pd(l2 + i);
          register const VD e1 = _mm512_load_pd(cat + i);
          register const VD e2 = _mm512_load_pd(sat + i);
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
          _a21r = _mm512_scalef_pd(_a21r, d);
          _a21i = _mm512_scalef_pd(_a21i, d);
          _mm512_store_pd((w1 + i), _a11);
          _mm512_store_pd((w2 + i), _a22);
          _mm512_store_pd((w3 + i), _a21r);
          _mm512_store_pd((w4 + i), _a21i);
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
      if (zbjac2i(&_n_2, w1, w2, w3, w4, c, cat, sat, l1, l2, p) < 0)
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
#pragma omp parallel for default(none) shared(m,n,Gr,ldGr,Gi,ldGi,Vr,ldVr,Vi,ldVi,eS,fS,a21r,a21i,c,cat,sat,l1,l2,w0,w1,n_2) reduction(max:nMG,nMV)
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
        double *const Gr_p = Gr + _p * (size_t)(*ldGr);
        double *const Gr_q = Gr + _q * (size_t)(*ldGr);
        double *const Gi_p = Gi + _p * (size_t)(*ldGi);
        double *const Gi_q = Gi + _q * (size_t)(*ldGi);
        double *const Vr_p = Vr + _p * (size_t)(*ldVr);
        double *const Vr_q = Vr + _q * (size_t)(*ldVr);
        double *const Vi_p = Vi + _p * (size_t)(*ldVi);
        double *const Vi_q = Vi + _q * (size_t)(*ldVi);
        if (w0[i] == -3.0) {
          const fint _m = -(fint)*m;
          const double e2[2u] = { eS[_p], eS[_q] };
          const double f2[2u] = { fS[_p], fS[_q] };
          double tG = zgsscl_(&_m, (a21r + i), (a21i + i), Gr_p, Gi_p, Gr_q, Gi_q, e2, f2);
          if (!isfinite(tG)) {
            nMG = HUGE_VAL;
            continue;
          }
          nMG = fmax(nMG, tG);
          if (dswp_(n, Vr_p, Vr_q)) {
            nMV = HUGE_VAL;
            continue;
          }
          if (dswp_(n, Vi_p, Vi_q)) {
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
          const double *const _cat = (cat + i);
          const double *const _sat = (sat + i);
          const double tG = zjrot_(&_m, Gr_p, Gi_p, Gr_q, Gi_q, _c, _cat, _sat);
          if (!isfinite(tG)) {
            nMG = HUGE_VAL;
            continue;
          }
          nMG = fmax(nMG, tG);
          const double tV = zjrot_(&_n, Vr_p, Vi_p, Vr_q, Vi_q, _c, _cat, _sat);
          if (!isfinite(tV)) {
            nMV = HUGE_VAL;
            continue;
          }
          nMV = fmax(nMV, tV);
        }
        else if (w0[i] == -1.0) {
          if (dswp_(m, Gr_p, Gr_q)) {
            nMG = HUGE_VAL;
            continue;
          }
          if (dswp_(m, Gi_p, Gi_q)) {
            nMG = HUGE_VAL;
            continue;
          }
          nMG = fmax(nMG, 0.0);
          if (dswp_(n, Vr_p, Vr_q)) {
            nMV = HUGE_VAL;
            continue;
          }
          if (dswp_(n, Vi_p, Vi_q)) {
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
          const double *const _cat = (cat + i);
          const double *const _sat = (sat + i);
          const double tG = zjrot_(&_m, Gr_p, Gi_p, Gr_q, Gi_q, _c, _cat, _sat);
          if (!isfinite(tG)) {
            nMG = HUGE_VAL;
            continue;
          }
          nMG = fmax(nMG, tG);
          const double tV = zjrot_(&_n, Vr_p, Vi_p, Vr_q, Vi_q, _c, _cat, _sat);
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
          double tG = zgsscl_(&_m, (a21r + i), (a21i + i), Gr_p, Gi_p, Gr_q, Gi_q, e2, f2);
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
    if (!swt)
      break;
    ++sw;
  }

  if (sw < *swp) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,Gr,ldGr,Gi,ldGi,eS,fS,sT)
#endif /* _OPENMP */
    for (fnat j = 0u; j < *n; ++j) {
      double *const Gr_j = Gr + j * (size_t)(*ldGr);
      double *const Gi_j = Gi + j * (size_t)(*ldGi);
      register const VD _f = _mm512_set1_pd(fS[j]);
      register const VD _s = _mm512_set1_pd(-(eS[j]));
      for (fnat i = 0u; i < *m; i += VDL) {
        double *const Gr_ij = Gr_j + i;
        double *const Gi_ij = Gi_j + i;
        _mm512_store_pd(Gr_ij, _mm512_scalef_pd(_mm512_div_pd(_mm512_load_pd(Gr_ij), _f), _s));
        _mm512_store_pd(Gi_ij, _mm512_scalef_pd(_mm512_div_pd(_mm512_load_pd(Gi_ij), _f), _s));
      }
      eS[j] -= sT;
    }
    if (sV) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,Vr,ldVr,Vi,ldVi,sV)
#endif /* _OPENMP */
      for (fnat j = 0u; j < *n; ++j) {
        double *const Vr_j = Vr + j * (size_t)(*ldVr);
        double *const Vi_j = Vi + j * (size_t)(*ldVi);
        register const VD _s = _mm512_set1_pd((double)-sV);
        for (fnat i = 0u; i < *n; i += VDL) {
          double *const Vr_ij = Vr_j + i;
          double *const Vi_ij = Vi_j + i;
          _mm512_store_pd(Vr_ij, _mm512_scalef_pd(_mm512_load_pd(Vr_ij), _s));
          _mm512_store_pd(Vi_ij, _mm512_scalef_pd(_mm512_load_pd(Vi_ij), _s));
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
