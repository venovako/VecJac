#include "dvjsvd.h"

#include "dscale.h"
#include "dnorm2.h"
#include "ddpscl.h"
#include "djac2.h"
#include "djrot.h"
#include "vecdef.h"
#include "defops.h"

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
  unsigned *const pc = p + (n_2 >> VDLlg);

  double M = 0.0;
  // TODO: initial dnormx (to M) and dscale, return -15 if fail (infs and/or nans)

  // see LAPACK's DGESVJ
  const double tol = sqrt((double)(*m)) * scalbn(DBL_EPSILON, -1);
  unsigned sw = 0u;

  while (sw < *swp) {
    size_t swt = 0u;
    for (unsigned st = 0u; st < *stp; ++st) {
      // TODO: rescale according to M if necessary and update M
      const unsigned *const r = js + st * (size_t)(*n);
      double pe = 0.0, qe = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,G,ldG,eS,fS,a11,a22,a21,s,c,l1,l2,tol) reduction(min:pe,qe)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat pq_ = pq + 1u;
        const fnat _pq = (pq >> 1u);
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        double *const Gp = G + _p * (*ldG);
        double *const Gq = G + _q * (*ldG);
        if ((pe = fmin(pe, dnorm2_(m, Gp, (eS + _p), (fS + _p), (s + _p), (c + _p)))) < 0.0)
          continue;
        if ((qe = fmin(qe, dnorm2_(m, Gq, (eS + _q), (fS + _q), (s + _q), (c + _q)))) < 0.0)
          continue;
        // pack the norms
        s[pq] = eS[_p];
        c[pq] = fS[_p];
        s[pq_] = eS[_q];
        c[pq_] = fS[_q];
        const double d = ddpscl_(m, Gp, Gq, (s + pq), (c + pq));
#ifndef NDEBUG
        if (!isfinite(d)) {
          pe = -1.0;
          if (isnan(d))
            qe = -1.0;
        }
        if ((pe < 0.0) || (qe < 0.0))
          continue;
#endif /* !NDEBUG */
        // ensure that the compiler does not rearrange accesses to t and l1
        _mm_mfence();
        // repack data
        s[_pq] = d;
        t[_pq] = fS[_p];
        c[_pq] = fS[_q];
        l1[_pq] = eS[_p];
        l2[_pq] = eS[_q];
      }
      if (pe < 0.0)
        return -16;
      if (qe < 0.0)
        return -17;
      fnat stt = 0u, k = 0u;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n_2,a11,a22,a21,s,t,c,l1,l2,p,pc,tol,k) reduction(+:stt)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; i += VDL) {
        const fnat j = (i >> VDLlg);
        // convergence check
        register VD _a21 = _mm512_load_pd(s + i);
        register const VD _zero = _mm512_set1_pd(-0.0);
        register const VD _a21_ = VDABS(_a21);
        pc[j] = MD2U(_mm512_cmple_pd_mask(_mm512_set1_pd(tol), _a21_));
        if (!(p[j] = _mm_popcnt_u32(pc[j])))
          continue;
        stt += p[j];
        // Grammian pre-scaling into the double precision range
        register const VD f1 = _mm512_load_pd(t + i);
        register const VD f2 = _mm512_load_pd(c + i);
        register const VD e1 = _mm512_load_pd(l1 + i);
        register const VD e2 = _mm512_load_pd(l2 + i);
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
#pragma omp atomic capture seq_cst
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
      if (djac2_(&kk, a11, a22, a21, t, c, l1, l2, p) < 0)
        return -18;
      fnat np = 0u; // number of swaps
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(a11,a22,s,p,pc,r,k) reduction(+:np)
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
            s[k_] = ((b & 1u) ? -2.0 : -1.0);
            ++np;
          }
          else
            s[k_] = ((b & 1u) ? 2.0 : 1.0);
        }
      }
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,G,ldG,V,ldV,a11,a22,s,t,c,kk) reduction(max:M)
#endif /* _OPENMP */
      for (fnat i = 0u; i < kk; ++i) {
        const size_t _p = *(const size_t*)(a11 + i);
        const size_t _q = *(const size_t*)(a22 + i);
        double _t, _c;
        fint _m, _n;
        bool triv = false;
        if (s[i] == -2.0) {
          _m = -(fint)*m;
          _n = -(fint)*n;
          _t = t[i];
          _c = c[i];
        }
        else if (s[i] == -1.0) {
          _m = -(fint)*m;
          _n = -(fint)*n;
          _t = 0.0;
          _c = 1.0;
        }
        else if (s[i] == 2.0) {
          _m = (fint)*m;
          _n = (fint)*n;
          _t = t[i];
          _c = c[i];
        }
        else // no-op
          triv = true;
        if (triv)
          M = fmax(M, 0.0);
        else {
          s[i] = djrot_(&_m, (G + _p * (*ldG)), (G + _q * (*ldG)), &_t, &_c);
          M = fmax(M, ((s[i] < 0.0) ? HUGE_VAL : s[i]));
          s[i] = djrot_(&_n, (V + _p * (*ldV)), (V + _q * (*ldV)), &_t, &_c);
          // V should not overflow but check anyway
          if ((s[i] < 0.0) || (s[i] > DBL_MAX))
            M = HUGE_VAL;
        }
      }
      if (M > DBL_MAX)
        return -19;
    }
    if (!swt)
      break;
    ++sw;
  }

  // TODO: normalize U and extract S

  return (fint)sw;
}
