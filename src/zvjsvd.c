#ifdef TEST_JAC2
#error TEST_JAC2 has to be undefined
#endif /* TEST_JAC2 */

#include "zvjsvd.h"

#include "znormx.h"
#include "zscale.h"
#include "znorm2.h"
#include "dznrm2.h"
#include "zdpscl.h"
#include "zjac2.h"
#include "zjrot.h"
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
  if (*n >= 0x10000000u) // 2^28
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

  double M = znormx_(m, n, Gr, ldGr, Gi, ldGi);
  if (!(M <= DBL_MAX))
    return -19;
  if (copysign(1.0, M) == -1.0)
    return -20;

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
  double *const t = a21i + n_2;
  double *const c = t + n_2;
  double *const ca = c + n_2;
  double *const sa = ca + n_2;
  double *const l1 = sa + n_2;
  double *const l2 = l1 + n_2;
  unsigned *const p = iwork;
  unsigned *const pc = p + (n_2 >> VDLlg);

  if (M == 0.0)
    return 0;
  const double M_m = (DBL_MAX / ((*m << 2u) * M_SQRT2));
  dbl2ef(M_m, t, c);
  const int DBL_MAX_NRM_EXP = (int)*t;
  dbl2ef(M, t, c);
  int eM = (int)*t;
  int sR = DBL_MAX_ROT_EXP - eM - 1;
  int sN = DBL_MAX_NRM_EXP - eM - 1;
  if (sN) {
    *(fint*)t = sN;
    if (zscale_(m, n, Gr, ldGr, Gi, ldGi, (const fint*)t) < 0)
      return -21;
    M = scalbn(M, sN);
  }
  int sT = sN;

  // see LAPACK's ZGESVJ
  const double tol = sqrt((double)(*m)) * scalbn(DBL_EPSILON, -1);
  unsigned sw = 0u;

  while (sw < *swp) {
    size_t swt = 0u;
    for (unsigned st = 0u; st < *stp; ++st) {
      // rescale according to M if necessary and update M
      dbl2ef(M, t, c);
      eM = (int)*t;
      sR = DBL_MAX_ROT_EXP - eM - 1;
      sN = DBL_MAX_NRM_EXP - eM - 1;
      if (sR < 0) {
        (void)fprintf(stderr, "Transformations in danger in sweep %u, step %u; rescaling by 2^%d.\n", sw, st, sN);
        (void)fflush(stderr);
        *(fint*)t = sN;
        if (zscale_(m, n, Gr, ldGr, Gi, ldGi, (const fint*)t) < 0)
          return -22;
        M = scalbn(M, sN);
        sT += sN;
      }
      // compute the norms, overflow-aware
      const unsigned *const r = js + st * (size_t)(*n);
      bool overflow = false;
      double nM = 0.0;
      do {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,Gr,ldGr,Gi,ldGi,eS,fS,t,c,l1,l2) reduction(max:nM)
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
          double *const Grp = Gr + _p * (*ldGr);
          double *const Gip = Gi + _p * (*ldGi);
          nM = fmax(nM, fmin((l1[_pq] = znorm2_(m, Grp, Gip, (eS + _p), (fS + _p), (t + _pq), (c + _pq))), HUGE_VAL));
          if (nM > DBL_MAX) {
            l2[_pq] = NAN;
            continue;
          }
          double *const Grq = Gr + _q * (*ldGr);
          double *const Giq = Gi + _q * (*ldGi);
          nM = fmax(nM, fmin((l2[_pq] = znorm2_(m, Grq, Giq, (eS + _q), (fS + _q), (t + _pq), (c + _pq))), HUGE_VAL));
        }
        if (overflow = (nM > DBL_MAX)) {
          (void)fprintf(stderr, "Frobenius norm overflow in sweep %u, step %u; rescaling by 2^%d.\n", sw, st, sN);
          (void)fflush(stderr);
          *(fint*)t = sN;
          if (zscale_(m, n, Gr, ldGr, Gi, ldGi, (const fint*)t) < 0)
            return -23;
          M = scalbn(M, sN);
          sT += sN;
        }
      } while (overflow);
      // scaled dot-products
      nM = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,Gr,ldGr,Gi,ldGi,eS,fS,l1,l2) reduction(min:nM)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat _pq = (pq >> 1u);
        if (nM < 0.0) {
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
        const double complex z = zdpscl_(m, Grp, Gip, Grq, Giq, e2, f2);
        l1[_pq] = creal(z);
        l2[_pq] = cimag(z);
        if (!(isfinite(l2[_pq])))
          nM = fmin(nM, -25.0);
        if (!(isfinite(l1[_pq])))
          nM = fmin(nM, -24.0);
      }
      if (nM < 0.0)
        return (fint)nM;
      // repack data
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,eS,fS,t,c,ca,sa)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat pq_ = pq + 1u;
        const fnat _pq = (pq >> 1u);
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        t[_pq] = eS[_p];
        c[_pq] = eS[_q];
        ca[_pq] = fS[_p];
        sa[_pq] = fS[_q];
      }
      fnat stt = 0u, k = 0u;
      // TODO
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n_2,a11,a22,a21r,a21i,t,c,ca,sa,l1,l2,p,pc,tol,k) reduction(+:stt)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; i += VDL) {
        const fnat j = (i >> VDLlg);
        // convergence check
        register VD _a21r = _mm512_load_pd(l1 + i);
        register VD _a21i = _mm512_load_pd(l2 + i);
        register const VD _a21_ = _mm512_hypot_pd(_a21r, _a21i);
        pc[j] = MD2U(_mm512_cmple_pd_mask(_mm512_set1_pd(tol), _a21_));
        if (!(p[j] = _mm_popcnt_u32(pc[j])))
          continue;
        stt += p[j];
        // Grammian pre-scaling into the double precision range
        register const VD f1 = _mm512_load_pd(ca + i);
        register const VD f2 = _mm512_load_pd(sa + i);
        register const VD e1 = _mm512_load_pd(t + i);
        register const VD e2 = _mm512_load_pd(c + i);
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
        _a21r = _mm512_scalef_pd(_a21r, d);
        _a21i = _mm512_scalef_pd(_a21i, d);
        // pack the data and record the translation in pc
        fnat kk;
#pragma omp atomic capture seq_cst
        kk = k++;
        // lower 8 bits: mask, upper 24 bits: j (assumes n < 2^28)
        pc[kk] |= (j << VDL);
        kk <<= VDLlg;
        _mm512_store_pd((a11 + kk), _a11);
        _mm512_store_pd((a22 + kk), _a22);
        _mm512_store_pd((a21r + kk), _a21r);
        _mm512_store_pd((a21i + kk), _a21i);
      }
      if (!stt)
        continue;
      swt += stt;
      const fnat kk = (k << VDLlg);
      if (zjac2_(&kk, a11, a22, a21r, a21i, t, c, ca, sa, p) < 0)
        return -26;
      fnat np = 0u; // number of swaps
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(a11,a22,a21r,p,pc,r,k) reduction(+:np)
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
            a21r[k_] = ((b & 1u) ? -2.0 : -1.0);
            ++np;
          }
          else // no swap
            a21r[k_] = ((b & 1u) ? 2.0 : 1.0);
        }
      }
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,Gr,ldGr,Gi,ldGi,Vr,ldVr,Vi,ldVi,a11,a22,a21r,a21i,t,c,ca,sa,kk) reduction(max:M)
#endif /* _OPENMP */
      for (fnat i = 0u; i < kk; ++i) {
        if (M > DBL_MAX) {
          a21i[i] = NAN;
          continue;
        }
        const size_t _p = *(const size_t*)(a11 + i);
        const size_t _q = *(const size_t*)(a22 + i);
        double _t, _c, _ca, _sa;
        fint _m, _n;
        bool triv = false;
        if (a21r[i] == -2.0) {
          _m = -(fint)*m;
          _n = -(fint)*n;
          _t = t[i];
          _c = c[i];
          _ca = ca[i];
          _sa = sa[i];
        }
        else if (a21r[i] == -1.0) {
          _m = -(fint)*m;
          _n = -(fint)*n;
          _t = 0.0;
          _c = 1.0;
          _ca = 1.0;
          _sa = 0.0;
        }
        else if (a21r[i] == 2.0) {
          _m = (fint)*m;
          _n = (fint)*n;
          _t = t[i];
          _c = c[i];
          _ca = ca[i];
          _sa = sa[i];
        }
        else // no-op
          triv = true;
        if (triv)
          M = fmax(M, 0.0);
        else {
          a21i[i] = zjrot_(&_m, (Gr + _p * (*ldGr)), (Gi + _p * (*ldGi)), (Gr + _q * (*ldGr)), (Gi + _q * (*ldGi)), &_t, &_c, &_ca, &_sa);
          M = fmax(M, (!(a21i[i] >= 0.0) ? HUGE_VAL : a21i[i]));
          if (M > DBL_MAX) {
            a21i[i] = NAN;
            continue;
          }
          a21i[i] = zjrot_(&_n, (Vr + _p * (*ldVr)), (Vi + _p * (*ldVi)), (Vr + _q * (*ldVr)), (Vi + _q * (*ldVi)), &_t, &_c, &_ca, &_sa);
          // V should not overflow but check anyway
          if (!(a21i[i] >= 0.0) || !(a21i[i] <= DBL_MAX))
            M = HUGE_VAL;
        }
      }
      if (M > DBL_MAX)
        return -27;
    }
    if (!swt)
      break;
    ++sw;
  }

  // TODO: normalize U and extract S

  return (fint)sw;
}
