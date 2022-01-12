#include "zvjsvd.h"

#include "znormx.h"
#include "zscale.h"
#include "znorm2.h"
#include "dznrm2.h"
#include "zdpscl.h"
#include "zbjac2.h"
#include "zjrotf.h"
#include "zjrot.h"
#include "dswp.h"
#include "vecdef.h"
#include "defops.h"

#ifdef DBL_MAX_ROT_EXP
#error DBL_MAX_ROT_EXP already defined
#else /* !DBL_MAX_ROT_EXP */
#define DBL_MAX_ROT_EXP 1021
#endif /* ?DBL_MAX_ROT_EXP */

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
    return -13;
  (void)fprintf(jtr, "M=");
  (void)fflush(jtr);
#endif /* JTRACE */

  double M = znormx_(m, n, Gr, ldGr, Gi, ldGi);
  if (!(M <= DBL_MAX))
    return -19;
  if (copysign(1.0, M) == -1.0)
    return -20;

#ifdef JTRACE
  (void)fprintf(jtr, "%#.17e\n", M);
  (void)fflush(jtr);
#endif /* JTRACE */

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,Vr,ldVr,Vi,ldVi,eS,fS)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    register const VD z = _mm512_setzero_pd();
    double *const Vrj = Vr + j * (size_t)(*ldVr);
    double *const Vij = Vi + j * (size_t)(*ldVi);
    for (fnat i = 0u; i < *n; i += VDL) {
      _mm512_store_pd((Vrj + i), z);
      _mm512_store_pd((Vij + i), z);
    }
    fS[j] = Vrj[j] = 1.0;
    eS[j] = -HUGE_VAL;
  }

  double *const a11 = work;
  double *const a22 = a11 + n_2;
  double *const a21r = a22 + n_2;
  double *const a21i = a21r + n_2;
  double *const c = a21i + n_2;
  double *const cat = c + n_2;
  double *const sat = cat + n_2;
  double *const l1 = sat + n_2;
  double *const l2 = l1 + n_2;
  double *const w = l2 + n_2;
  unsigned *const p = iwork;
  unsigned *const pc = p + (n_2 >> VDLlg);

  if (M == 0.0)
    return 0;
  const double M_m = (DBL_MAX / ((*m << 2u) * M_SQRT2));
  double es = 0.0, fs = 0.0;
  dbl2ef(M_m, &es, &fs);
  const int DBL_MAX_NRM_EXP = (int)es;
  dbl2ef(M, &es, &fs);
  int eM = (int)es;
  int sR = DBL_MAX_ROT_EXP - eM - 1;
  int sN = DBL_MAX_NRM_EXP - eM - 1;
#ifdef JTRACE
  (void)fprintf(jtr, "eM=%d, sR=%d, sN=%d, M=", eM, sR, sN);
  (void)fflush(jtr);
#endif /* JTRACE */
  if (sN) {
    *(fint*)&es = sN;
    if (zscale_(m, n, Gr, ldGr, Gi, ldGi, (const fint*)&es) < 0)
      return -21;
    M = scalbn(M, sN);
  }
  int sT = sN;
#ifdef JTRACE
  (void)fprintf(jtr, "%#.17e\n", M);
  (void)fflush(jtr);
#endif /* JTRACE */

  // see LAPACK's ZGESVJ
  const double tol = sqrt((double)(*m)) * scalbn(DBL_EPSILON, -1);
  unsigned sw = 0u;

  while (sw < *swp) {
    size_t swt = 0u;
    for (unsigned st = 0u; st < *stp; ++st) {
      // rescale according to M if necessary and update M
      dbl2ef(M, &es, &fs);
      eM = (int)es;
      sR = DBL_MAX_ROT_EXP - eM - 1;
      sN = DBL_MAX_NRM_EXP - eM - 1;
      if (sR < 0) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u, eM=%d, sR=%d, sN=%d, M=", sw, st, eM, sR, sN);
        (void)fflush(jtr);
#endif /* JTRACE */
        *(fint*)&es = sN;
        if (zscale_(m, n, Gr, ldGr, Gi, ldGi, (const fint*)&es) < 0)
          return -22;
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
#pragma omp parallel for default(none) shared(n,r,m,Gr,ldGr,Gi,ldGi,eS,fS,cat,sat,l1,l2) reduction(max:nM)
#endif /* _OPENMP */
        for (fnat pq = 0u; pq < *n; pq += 2u) {
          const fnat _pq = (pq >> 1u);
          if (!(nM <= DBL_MAX)) {
            l1[_pq] = NAN;
            l2[_pq] = NAN;
            continue;
          }
          const fnat pq_ = pq + 1u;
          const size_t _p = r[pq];
          const size_t _q = r[pq_];
          double *const Grp = Gr + _p * (*ldGr);
          double *const Gip = Gi + _p * (*ldGi);
          nM = fmax(nM, fmin((l1[_pq] = znorm2_(m, Grp, Gip, (eS + _p), (fS + _p), (cat + _pq), (sat + _pq))), HUGE_VAL));
          if (!(nM <= DBL_MAX)) {
            l2[_pq] = NAN;
            continue;
          }
          double *const Grq = Gr + _q * (*ldGr);
          double *const Giq = Gi + _q * (*ldGi);
          nM = fmax(nM, fmin((l2[_pq] = znorm2_(m, Grq, Giq, (eS + _q), (fS + _q), (cat + _pq), (sat + _pq))), HUGE_VAL));
        }
        if (overflow = !(nM <= DBL_MAX)) {
#ifdef JTRACE
          (void)fprintf(jtr, "sweep=%u, step=%u, M=", sw, st);
          (void)fflush(jtr);
#endif /* JTRACE */
          *(fint*)&es = sN;
          if (zscale_(m, n, Gr, ldGr, Gi, ldGi, (const fint*)&es) < 0)
            return -23;
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
#pragma omp parallel for default(none) shared(n,r,m,Gr,ldGr,Gi,ldGi,eS,fS,l1,l2) reduction(min:nM)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat _pq = (pq >> 1u);
        if (!(nM >= 0.0)) {
          l1[_pq] = NAN;
          l2[_pq] = NAN;
          continue;
        }
        const fnat pq_ = pq + 1u;
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        // pack the norms
        const double e2[2u] = { eS[_p], eS[_q] };
        const double f2[2u] = { fS[_p], fS[_q] };
        double *const Grp = Gr + _p * (*ldGr);
        double *const Gip = Gi + _p * (*ldGi);
        double *const Grq = Gr + _q * (*ldGr);
        double *const Giq = Gi + _q * (*ldGi);
        const double complex z = zdpscl_(m, Grq, Giq, Grp, Gip, e2, f2);
        l1[_pq] = creal(z);
        l2[_pq] = cimag(z);
        if (!(isfinite(l2[_pq])))
          nM = fmin(nM, -25.0);
        if (!(isfinite(l1[_pq])))
          nM = fmin(nM, -24.0);
      }
      if (!(nM >= 0.0)) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u\n", sw, st);
        (void)fflush(jtr);
#endif /* JTRACE */
        return (fint)nM;
      }
      // repack data
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,eS,fS,c,cat,sat,w)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat pq_ = pq + 1u;
        const fnat _pq = (pq >> 1u);
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        c[_pq] = eS[_p];
        w[_pq] = eS[_q];
        cat[_pq] = fS[_p];
        sat[_pq] = fS[_q];
      }
      fnat stt = 0u;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n_2,a11,a22,a21r,a21i,c,cat,sat,l1,l2,w,p,pc,tol) reduction(+:stt)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; i += VDL) {
        const fnat j = (i >> VDLlg);
        // convergence check
        register VD _a21r = _mm512_load_pd(l1 + i);
        register VD _a21i = _mm512_load_pd(l2 + i);
        register const VD _zero = _mm512_set1_pd(-0.0);
        register const VD zero = _mm512_setzero_pd();
        register const VD one = _mm512_set1_pd(1.0);
        register const VD _tol = _mm512_set1_pd(tol);
        register VD _a21_ /* = _mm512_hypot_pd(_a21r, _a21i) */;
        VDHYPOT(_a21_, _a21r, _a21i);
        pc[j] = MD2U(_mm512_cmple_pd_mask(_tol, _a21_));
        stt += (p[j] = _mm_popcnt_u32(pc[j]));
        // Grammian pre-scaling into the double precision range
        register const VD f1 = _mm512_load_pd(cat + i);
        register const VD f2 = _mm512_load_pd(sat + i);
        register const VD e1 = _mm512_load_pd(c + i);
        register const VD e2 = _mm512_load_pd(w + i);
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
        _mm512_store_pd((a11 + i), _a11);
        _mm512_store_pd((a22 + i), _a22);
        _mm512_store_pd((a21r + i), _a21r);
        _mm512_store_pd((a21i + i), _a21i);
      }
      if (stt) {
        swt += stt;
#ifdef JTRACE
        //(void)fprintf(jtr, "sweep=%u, step=%u, trans=%llu\n", sw, st, stt);
        //(void)fflush(jtr);
#endif /* JTRACE */
      }
      if (zbjac2_(&n_2, a11, a22, a21r, a21i, c, cat, sat, l1, l2, p) < 0)
        return -26;
      fnat np = 0u; // number of swaps
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(a11,a22,a21r,eS,fS,p,pc,r,n_2) reduction(+:np)
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
              a21r[l] = -2.0;
              ++np;
            }
            else // no swap
              a21r[l] = 2.0;
          }
          else if (efcmp((eS + _p), (fS + _p), (eS + _q), (fS + _q)) < 0) {
            a21r[l] = eS[_p];
            eS[_p] = eS[_q];
            eS[_q] = a21r[l];
            a21r[l] = fS[_p];
            fS[_p] = fS[_q];
            fS[_q] = a21r[l];
            a21r[l] = -1.0;
            ++np;
          }
          else // no swap
            a21r[l] = 1.0;
          trans >>= 1u;
          perm >>= 1u;
        }
      }
#ifdef JTRACE
      //if (np) {
      //  (void)fprintf(jtr, "sweep=%u, step=%u, swaps=%llu\n", sw, st, np);
      //  (void)fflush(jtr);
      //}
#endif /* JTRACE */
      nM = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,Gr,ldGr,Gi,ldGi,Vr,ldVr,Vi,ldVi,a11,a22,a21r,a21i,c,cat,sat,n_2) reduction(max:nM)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; ++i) {
        if (!(nM <= DBL_MAX)) {
          a21i[i] = NAN;
          continue;
        }
        const size_t _p = *(const uint64_t*)(a11 + i);
        const size_t _q = *(const uint64_t*)(a22 + i);
        double _c, _cat, _sat;
        fint _m, _n;
        if (a21r[i] == -2.0) {
          _m = -(fint)*m;
          _n = -(fint)*n;
          _c = c[i];
          _cat = cat[i];
          _sat = sat[i];
        }
        else if (a21r[i] == -1.0) {
          double *const Gr_p = Gr + _p * (*ldGr);
          double *const Gr_q = Gr + _q * (*ldGr);
          if (_m = dswp_(m, Gr_p, Gr_q)) {
            a21i[i] = _m;
            nM = HUGE_VAL;
            continue;
          }
          double *const Gi_p = Gi + _p * (*ldGi);
          double *const Gi_q = Gi + _q * (*ldGi);
          if (_m = dswp_(m, Gi_p, Gi_q)) {
            a21i[i] = _m;
            nM = HUGE_VAL;
            continue;
          }
          double *const Vr_p = Vr + _p * (*ldVr);
          double *const Vr_q = Vr + _q * (*ldVr);
          if (_n = dswp_(n, Vr_p, Vr_q)) {
            a21i[i] = _n;
            nM = HUGE_VAL;
            continue;
          }
          double *const Vi_p = Vi + _p * (*ldVi);
          double *const Vi_q = Vi + _q * (*ldVi);
          if (_n = dswp_(n, Vi_p, Vi_q)) {
            a21i[i] = _n;
            nM = HUGE_VAL;
            continue;
          }
          nM = fmax(nM, (a21i[i] = 0.0));
          continue;
        }
        else if (a21r[i] == 1.0) {
          nM = fmax(nM, (a21i[i] = 0.0));
          continue;
        }
        else if (a21r[i] == 2.0) {
          _m = (fint)*m;
          _n = (fint)*n;
          _c = c[i];
          _cat = cat[i];
          _sat = sat[i];
        }
        else { // should never happen
          a21i[i] = NAN;
          nM = HUGE_VAL;
          continue;
        }
        a21i[i] = zjrot_(&_m, (Gr + _p * (*ldGr)), (Gi + _p * (*ldGi)), (Gr + _q * (*ldGr)), (Gi + _q * (*ldGi)), &_c, &_cat, &_sat);
        if (!(a21i[i] >= 0.0) || !(a21i[i] <= DBL_MAX)) {
          nM = a21i[i] = HUGE_VAL;
          continue;
        }
        else // no overflow
          nM = fmax(nM, a21i[i]);
        if (_m = zjrot_(&_n, (Vr + _p * (*ldVr)), (Vi + _p * (*ldVi)), (Vr + _q * (*ldVr)), (Vi + _q * (*ldVi)), &_c, &_cat, &_sat)) {
          a21i[i] = _m;
          nM = HUGE_VAL;
        }
      }
      M = fmax(M, nM);
      if (!(M <= DBL_MAX)) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u\n", sw, st);
        (void)fflush(jtr);
#endif /* JTRACE */
        return -27;
      }
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
  }

#ifdef JTRACE
  (void)fprintf(jtr, "sT=%d, M=%#.17e\n", sT, M);
  (void)fclose(jtr);
#endif /* JTRACE */
  return (fint)sw;
}
