#include "cvjsvd.h"

#include "cnormx.h"
#include "cscale.h"
#include "cnorm2.h"
#include "scnrm2.h"
#include "cdpscl.h"
#include "cgsscl.h"
#include "cbjac2.h"
#include "cjrot.h"
#include "cjrotf.h"
#include "sswp.h"
#include "vecdef.h"
#include "sefops.h"

#ifdef JTRACE
#include "timer.h"
#endif /* JTRACE */

#ifdef FLT_MAX_ROT_EXP
#error FLT_MAX_ROT_EXP already defined
#else /* !FLT_MAX_ROT_EXP */
#define FLT_MAX_ROT_EXP 125
#endif /* ?FLT_MAX_ROT_EXP */

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

typedef float (*vrot_t)(const fint *const, float *const, float *const, float *const, float *const, const float *const, const float *const, const float *const);

fint cvjsvd_(const fnat m[static restrict 1], const fnat n[static restrict 1], float Gr[static restrict VSL], const fnat ldGr[static restrict 1], float Gi[static restrict VSL], const fnat ldGi[static restrict 1], float Vr[static restrict VSL], const fnat ldVr[static restrict 1], float Vi[static restrict VSL], const fnat ldVi[static restrict 1], float eS[static restrict 1], float fS[static restrict 1], const unsigned js[static restrict 1], const unsigned stp[static restrict 1], const unsigned swp[static restrict 1], float work[static restrict VSL], unsigned iwork[static restrict 1])
{
  const fnat n_2 = (*n >> 1u);
  if (IS_NOT_VFPENV)
    return -18;
  if (!*n)
    return 0;
  if (*m < *n)
    return -1;
  if (*m & VSL_1)
    return -1;
  if (*n & 1u)
    return -2;
  if (n_2 & VSL_1)
    return -2;
  if (IS_NOT_ALIGNED(Gr))
    return -3;
  if (*ldGr < *m)
    return -4;
  if (*ldGr & VSL_1)
    return -4;
  if (IS_NOT_ALIGNED(Gi))
    return -5;
  if (*ldGi < *m)
    return -6;
  if (*ldGi & VSL_1)
    return -6;
  if (IS_NOT_ALIGNED(Vr))
    return -7;
  if (*ldVr < *n)
    return -8;
  if (*ldVr & VSL_1)
    return -8;
  if (IS_NOT_ALIGNED(Vi))
    return -9;
  if (*ldVi < *n)
    return -10;
  if (*ldVi & VSL_1)
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

  float MG = cnormx_(m, n, Gr, ldGr, Gi, ldGi);
  if (!(MG <= FLT_MAX))
    return -__LINE__;
  if (copysignf(1.0f, MG) == -1.0f)
    return -__LINE__;
#ifdef JTRACE
  (void)fprintf(jtr, "%#.9e\n", MG);
  (void)fflush(jtr);
#endif /* JTRACE */

  float MV = 1.0f;
  if (*iwork) {
    MV = cnormx_(n, n, Vr, ldVr, Vi, ldVi);
    if (!(MV <= FLT_MAX))
      return -__LINE__;
    // a dirty hack to avoid accumulation on a zero matrix
    if (MV <= 0.0f)
      return -17;
  }
  else { // V = I
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,Vr,ldVr,Vi,ldVi)
#endif /* _OPENMP */
    for (fnat j = 0u; j < *n; ++j) {
      register const VS z = _mm512_setzero_ps();
      float *const Vrj = Vr + j * (size_t)(*ldVr);
      float *const Vij = Vi + j * (size_t)(*ldVi);
      for (fnat i = 0u; i < *n; i += VSL) {
        _mm512_store_ps((Vrj + i), z);
        _mm512_store_ps((Vij + i), z);
      }
      Vrj[j] = 1.0f;
    }
  }

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,eS,fS)
#endif /* _OPENMP */
  for (fnat j = 0u; j < *n; ++j) {
    eS[j] = -HUGE_VALF;
    fS[j] = 1.0f;
  }

  if (MG == 0.0f)
    return 0;
  const float M_m = (FLT_MAX / ((*m << 2u) * (float)(M_SQRT2)));
  float es = 0.0f, fs = 0.0f;
  flt2ef(M_m, &es, &fs);
  const int FLT_MAX_NRM_EXP = (int)es;
  flt2ef(MG, &es, &fs);
  int eM = (int)es;
  int sR = FLT_MAX_ROT_EXP - eM - 1;
  fint sN = FLT_MAX_NRM_EXP - eM - 1;
#ifdef JTRACE
  (void)fprintf(jtr, "eM=%d, sR=%d, sN=%d, M=", eM, sR, (int)sN);
  (void)fflush(jtr);
#endif /* JTRACE */
  if (sN) {
    if (cscale_(m, n, Gr, ldGr, Gi, ldGi, &sN) < 0)
      return -__LINE__;
    MG = scalbnf(MG, (int)sN);
  }
  int sT = (int)sN;
#ifdef JTRACE
  (void)fprintf(jtr, "%#.9e\n", MG);
  (void)fflush(jtr);
#endif /* JTRACE */

  if (*iwork) {
    flt2ef(MV, &es, &fs);
    eM = (int)es;
    sR = FLT_MAX_ROT_EXP - eM - 1;
    sN = sR;
    if (sN) {
      if (cscale_(n, n, Vr, ldVr, Vi, ldVi, &sN) < 0)
        return -__LINE__;
      MV = scalbnf(MV, sR);
    }
  }
  else // V = I
    sR = 0;
  int sV = sR;
  const vrot_t vrot = (*iwork ? cjrot_ : cjrotf_);

  const fnat n_32 = (n_2 >> VSLlg);

  float *const a11 = work;
  float *const a22 = a11 + n_2;
  float *const a21r = a22 + n_2;
  float *const a21i = a21r + n_2;
  float *const c = a21i + n_2;
  float *const cat = c + n_2;
  float *const sat = cat + n_2;
  float *const l1 = sat + n_2;
  float *const l2 = l1 + n_2;
  float *const w0 = l2 + n_2;
  float *const w1 = w0 + n_2;
  float *const w2 = w1 + n_2;
  float *const w3 = w2 + n_2;
  float *const w4 = w3 + n_2;
  unsigned *const p = iwork;
  unsigned *const pc = p + n_32;

  if (*swp) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(w1,n)
#endif /* _OPENMP */
    for (fnat i = 0u; i < *n; ++i)
      w1[i] = 1.0f;
  }

  // see LAPACK's CGESVJ
  const float tol = sqrtf((float)(*m)) * scalbnf(FLT_EPSILON, -1);
  const float gst = scalbf(tol, FLT_MAX_FIN_EXP);
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
      flt2ef(MG, &es, &fs);
      eM = (int)es;
      sR = FLT_MAX_ROT_EXP - eM - 1;
      sN = FLT_MAX_NRM_EXP - eM - 1;
      if (sR < 0) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u, eM=%d, sR=%d, sN=%d, M=", sw, st, eM, sR, (int)sN);
        (void)fflush(jtr);
#endif /* JTRACE */
        if (cscale_(m, n, Gr, ldGr, Gi, ldGi, &sN) < 0)
          return -__LINE__;
        sR = (int)sN;
        MG = scalbnf(MG, sR);
        sT += sR;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(w1,n)
#endif /* _OPENMP */
        for (fnat i = 0u; i < *n; ++i)
          w1[i] = 1.0f;
#ifdef JTRACE
        (void)fprintf(jtr, "%#.9e\n", MG);
        (void)fflush(jtr);
#endif /* JTRACE */
      }
      // rescale V according to MV if necessary and update MV
      flt2ef(MV, &es, &fs);
      eM = (int)es;
      sR = FLT_MAX_ROT_EXP - eM - 1;
      sN = sR;
      if (sR < 0) {
        if (cscale_(n, n, Vr, ldVr, Vi, ldVi, &sN) < 0)
          return -__LINE__;
        MV = scalbnf(MV, sR);
        sV += sR;
      }
      // compute the norms, overflow-aware
      const unsigned *const r = js + st * (size_t)(*n);
      float nMG = -0.0f;
      bool overflow = false;
      do {
#ifdef JTRACE
        T = rdtsc_beg(rd);
#endif /* JTRACE */
        nMG = 0.0f;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,Gr,ldGr,Gi,ldGi,eS,fS,a11,a22,cat,sat,l1,l2,w1) reduction(max:nMG)
#endif /* _OPENMP */
        for (fnat pq = 0u; pq < *n; pq += 2u) {
          const fnat _pq = (pq >> 1u);
          if (!(nMG <= FLT_MAX)) {
            a11[_pq] = NAN;
            a22[_pq] = NAN;
            continue;
          }
          const fnat pq_ = pq + 1u;
          const size_t _p = r[pq];
          const size_t _q = r[pq_];
          if (w1[_p] == 1.0f) {
            float *const Grp = Gr + _p * (size_t)(*ldGr);
            float *const Gip = Gi + _p * (size_t)(*ldGi);
            nMG = fmaxf(nMG, fminf((a11[_pq] = cnorm2_(m, Grp, Gip, (eS + _p), (fS + _p), (cat + _pq), (sat + _pq))), HUGE_VALF));
            if (!(nMG <= FLT_MAX)) {
              a22[_pq] = NAN;
              continue;
            }
          }
          else
            a11[_pq] = scalbf(fS[_p], eS[_p]);
          if (w1[_q] == 1.0f) {
            float *const Grq = Gr + _q * (size_t)(*ldGr);
            float *const Giq = Gi + _q * (size_t)(*ldGi);
            nMG = fmaxf(nMG, fminf((a22[_pq] = cnorm2_(m, Grq, Giq, (eS + _q), (fS + _q), (l1 + _pq), (l2 + _pq))), HUGE_VALF));
          }
          else
            a22[_pq] = scalbf(fS[_q], eS[_q]);
        }
#ifdef JTRACE
        Tn += tsc_lap(hz, T, rdtsc_end(rd));
#endif /* JTRACE */
        if ((overflow = !(nMG <= FLT_MAX))) {
#ifdef JTRACE
          (void)fprintf(jtr, "sweep=%u, step=%u, M=", sw, st);
          (void)fflush(jtr);
#endif /* JTRACE */
          if (cscale_(m, n, Gr, ldGr, Gi, ldGi, &sN) < 0)
            return -__LINE__;
          sR = (int)sN;
          MG = scalbnf(MG, sR);
          sT += sR;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(w1,n)
#endif /* _OPENMP */
          for (fnat i = 0u; i < *n; ++i)
            w1[i] = 1.0f;
#ifdef JTRACE
          (void)fprintf(jtr, "%#.9e\n", MG);
          (void)fflush(jtr);
#endif /* JTRACE */
        }
      } while (overflow);
      // scaled dot-products
#ifdef JTRACE
      T = rdtsc_beg(rd);
#endif /* JTRACE */
      nMG = 0.0f;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r,m,Gr,ldGr,Gi,ldGi,eS,fS,a21r,a21i) reduction(min:nMG)
#endif /* _OPENMP */
      for (fnat pq = 0u; pq < *n; pq += 2u) {
        const fnat _pq = (pq >> 1u);
        if (!(nMG >= 0.0f)) {
          a21r[_pq] = NAN;
          a21i[_pq] = NAN;
          continue;
        }
        const fnat pq_ = pq + 1u;
        const size_t _p = r[pq];
        const size_t _q = r[pq_];
        // pack the norms
        const float e2[2u] = { eS[_q], eS[_p] };
        const float f2[2u] = { fS[_q], fS[_p] };
        float *const Grp = Gr + _p * (size_t)(*ldGr);
        float *const Gip = Gi + _p * (size_t)(*ldGi);
        float *const Grq = Gr + _q * (size_t)(*ldGr);
        float *const Giq = Gi + _q * (size_t)(*ldGi);
        const float complex z = cdpscl_(m, Grq, Giq, Grp, Gip, e2, f2);
        a21r[_pq] = crealf(z);
        nMG = fminf(nMG, (isfinite(a21r[_pq]) ? 0.0f : (float)-__LINE__));
        a21i[_pq] = cimagf(z);
        nMG = fminf(nMG, (isfinite(a21i[_pq]) ? 0.0f : (float)-__LINE__));
      }
#ifdef JTRACE
      Tp += tsc_lap(hz, T, rdtsc_end(rd));
#endif /* JTRACE */
      if (!(nMG >= 0.0f)) {
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
      for (fnat i = 0u; i < n_2; i += VSL) {
        const fnat j = (i >> VSLlg);
        // convergence check
        register VS _a21r = _mm512_load_ps(a21r + i);
        register VS _a21i = _mm512_load_ps(a21i + i);
        register const VS _zerof = _mm512_set1_ps(-0.0f);
        register const VS zerof = _mm512_setzero_ps();
        register const VS onef = _mm512_set1_ps(1.0f);
        register const VS _tol = _mm512_set1_ps(tol);
        register VS _a21_ /* = _mm512_hypot_ps(_a21r, _a21i) */;
        VSHYPOT(_a21_, _a21r, _a21i);
        pc[j] = MS2U(_mm512_cmple_ps_mask(_tol, _a21_));
        if ((p[j] = _mm_popcnt_u32(pc[j]))) {
          stt += p[j];
          register const VS f1 = _mm512_load_ps(l1 + i);
          register const VS f2 = _mm512_load_ps(l2 + i);
          register const VS e1 = _mm512_load_ps(cat + i);
          register const VS e2 = _mm512_load_ps(sat + i);
          register const VS _gst = _mm512_set1_ps(gst);
          // might not yet be sorted, so check both cases
          register const MS ngs = VSEFLT(e1,e2,f1,f2);
          register VS maf = _mm512_mask_blend_ps(ngs, f2, f1);
          register VS mae = _mm512_mask_blend_ps(ngs, e2, e1);
          register const VS Maf = _mm512_mask_blend_ps(ngs, f1, f2);
          register const VS Mae = _mm512_mask_blend_ps(ngs, e1, e2);
          maf = _mm512_mul_ps(_gst, maf);
          mae = _mm512_add_ps(mae, _mm512_getexp_ps(maf));
          maf = VSMANT(maf);
          register const MS cgs = VSEFLT(mae,Mae,maf,Maf);
          const unsigned gsp = (MS2U(_kandn_mask16(ngs,cgs)) << VSL);
          p[j] |= gsp;
          const unsigned gsn = (MS2U(_kand_mask16(ngs,cgs)) << VSL);
          pc[j] |= gsn;
          // Grammian pre-scaling into the single precision range
          register VS f12 = _mm512_div_ps(f1, f2);
          register VS e12 = _mm512_sub_ps(e1, e2);
          register VS f21 = _mm512_div_ps(f2, f1);
          register VS e21 = _mm512_sub_ps(e2, e1);
          e12 = _mm512_add_ps(e12, _mm512_getexp_ps(f12));
          f12 = VSMANT(f12);
          e21 = _mm512_add_ps(e21, _mm512_getexp_ps(f21));
          f21 = VSMANT(f21);
          register const MS c12 = VSEFLE(e12,e21,f12,f21);
          register const VS mxe = _mm512_set1_ps(FLT_MAX_FIN_EXP);
          register const VS E = _mm512_mask_blend_ps(c12, e12, e21);
          register const VS d = _mm512_min_ps(_mm512_sub_ps(mxe, E), zerof);
          e12 = _mm512_add_ps(e12, d);
          e21 = _mm512_add_ps(e21, d);
          register const VS _a11 = _mm512_scalef_ps(f12, e12);
          register const VS _a22 = _mm512_scalef_ps(f21, e21);
          _a21r = _mm512_scalef_ps(_a21r, d);
          _a21i = _mm512_scalef_ps(_a21i, d);
          _mm512_store_ps((w1 + i), _a11);
          _mm512_store_ps((w2 + i), _a22);
          _mm512_store_ps((w3 + i), _a21r);
          _mm512_store_ps((w4 + i), _a21i);
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
      if (cbjac2i(&_n_2, w1, w2, w3, w4, c, cat, sat, l1, l2, p) < 0)
        return -__LINE__;
#ifdef JTRACE
      Te += tsc_lap(hz, T, rdtsc_end(rd));
      T = rdtsc_beg(rd);
#endif /* JTRACE */
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(eS,fS,l1,l2,w0,p,pc,r,n_2)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; i += VSL) {
        const fnat j = (i >> VSLlg);
        unsigned gsp = ((p[j] & 0xFFFF0000u) >> VSL);
        unsigned gsn = ((pc[j] & 0xFFFF0000u) >> VSL);
        unsigned trans = (pc[j] & 0xFFFFu);
        unsigned perm = (p[j] & 0xFFFFu);
        for (fnat k = 0u; k < VSL; ++k) {
          const fnat l = (i + k);
          const fnat pq = (l << 1u);
          const unsigned _p = r[pq];
          const unsigned _q = r[pq + 1u];
          *(unsigned*)(l1 + l) = _p;
          *(unsigned*)(l2 + l) = _q;
          if (trans & 1u) {
            if (gsp & 1u)
              w0[l] = 3.0f;
            else if (gsn & 1u)
              w0[l] = -3.0f;
            else if (perm & 1u)
              w0[l] = -2.0f;
            else // no swap
              w0[l] = 2.0f;
          }
          else if (efcmpf((eS + _p), (fS + _p), (eS + _q), (fS + _q)) < 0) {
            SWAP_EFS(w0[l]);
            w0[l] = -1.0f;
          }
          else // no swap
            w0[l] = 1.0f;
          gsp >>= 1u;
          gsn >>= 1u;
          trans >>= 1u;
          perm >>= 1u;
        }
      }
      nMG = 0.0f;
      float nMV = -0.0f;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(m,n,Gr,ldGr,Gi,ldGi,Vr,ldVr,Vi,ldVi,eS,fS,a21r,a21i,c,cat,sat,l1,l2,w0,w1,n_2,vrot) reduction(max:nMG,nMV)
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_2; ++i) {
        const unsigned _p = *(const unsigned*)(l1 + i);
        const unsigned _q = *(const unsigned*)(l2 + i);
        w1[_q] = w1[_p] = 0.0f;
        if (!(nMG <= FLT_MAX)) {
          nMG = HUGE_VALF;
          continue;
        }
        if (!(nMV <= FLT_MAX)) {
          nMV = HUGE_VALF;
          continue;
        }
        float *const Gr_p = Gr + _p * (size_t)(*ldGr);
        float *const Gr_q = Gr + _q * (size_t)(*ldGr);
        float *const Gi_p = Gi + _p * (size_t)(*ldGi);
        float *const Gi_q = Gi + _q * (size_t)(*ldGi);
        float *const Vr_p = Vr + _p * (size_t)(*ldVr);
        float *const Vr_q = Vr + _q * (size_t)(*ldVr);
        float *const Vi_p = Vi + _p * (size_t)(*ldVi);
        float *const Vi_q = Vi + _q * (size_t)(*ldVi);
        if (w0[i] == -3.0f) {
          const fint _m = -(fint)*m;
          const float e2[2u] = { eS[_p], eS[_q] };
          const float f2[2u] = { fS[_p], fS[_q] };
          float tG = cgsscl_(&_m, (a21r + i), (a21i + i), Gr_p, Gi_p, Gr_q, Gi_q, e2, f2);
          if (!isfinite(tG)) {
            nMG = HUGE_VALF;
            continue;
          }
          nMG = fmaxf(nMG, tG);
          if (sswp_(n, Vr_p, Vr_q)) {
            nMV = HUGE_VALF;
            continue;
          }
          if (sswp_(n, Vi_p, Vi_q)) {
            nMV = HUGE_VALF;
            continue;
          }
          nMV = fmaxf(nMV, 0.0f);
          SWAP_EFS(tG);
          w1[_q] = 1.0f;
          continue;
        }
        else if (w0[i] == -2.0f) {
          const fint _m = -(fint)*m;
          const fint _n = -(fint)*n;
          const float *const _c = (c + i);
          const float *const _cat = (cat + i);
          const float *const _sat = (sat + i);
          const float tG = cjrot_(&_m, Gr_p, Gi_p, Gr_q, Gi_q, _c, _cat, _sat);
          if (!isfinite(tG)) {
            nMG = HUGE_VALF;
            continue;
          }
          nMG = fmaxf(nMG, tG);
          const float tV = vrot(&_n, Vr_p, Vi_p, Vr_q, Vi_q, _c, _cat, _sat);
          if (!isfinite(tV)) {
            nMV = HUGE_VALF;
            continue;
          }
          nMV = fmaxf(nMV, tV);
        }
        else if (w0[i] == -1.0f) {
          if (sswp_(m, Gr_p, Gr_q)) {
            nMG = HUGE_VALF;
            continue;
          }
          if (sswp_(m, Gi_p, Gi_q)) {
            nMG = HUGE_VALF;
            continue;
          }
          nMG = fmaxf(nMG, 0.0f);
          if (sswp_(n, Vr_p, Vr_q)) {
            nMV = HUGE_VALF;
            continue;
          }
          if (sswp_(n, Vi_p, Vi_q)) {
            nMV = HUGE_VALF;
            continue;
          }
          nMV = fmaxf(nMV, 0.0f);
          continue;
        }
        else if (w0[i] == 1.0f) {
          nMG = fmaxf(nMG, 0.0f);
          nMV = fmaxf(nMV, 0.0f);
          continue;
        }
        else if (w0[i] == 2.0f) {
          const fint _m = (fint)*m;
          const fint _n = (fint)*n;
          const float *const _c = (c + i);
          const float *const _cat = (cat + i);
          const float *const _sat = (sat + i);
          const float tG = cjrot_(&_m, Gr_p, Gi_p, Gr_q, Gi_q, _c, _cat, _sat);
          if (!isfinite(tG)) {
            nMG = HUGE_VALF;
            continue;
          }
          nMG = fmaxf(nMG, tG);
          const float tV = vrot(&_n, Vr_p, Vi_p, Vr_q, Vi_q, _c, _cat, _sat);
          if (!isfinite(tV)) {
            nMV = HUGE_VALF;
            continue;
          }
          nMV = fmaxf(nMV, tV);
        }
        else if (w0[i] == 3.0f) {
          const fint _m = (fint)*m;
          const float e2[2u] = { eS[_p], eS[_q] };
          const float f2[2u] = { fS[_p], fS[_q] };
          float tG = cgsscl_(&_m, (a21r + i), (a21i + i), Gr_p, Gi_p, Gr_q, Gi_q, e2, f2);
          if (!isfinite(tG)) {
            nMG = HUGE_VALF;
            continue;
          }
          nMG = fmaxf(nMG, tG);
          nMV = fmaxf(nMV, 0.0f);
          w1[_q] = 1.0f;
          continue;
        }
        else { // should never happen
          nMG = HUGE_VALF;
          nMV = HUGE_VALF;
          continue;
        }
        w1[_q] = w1[_p] = 1.0f;
      }
      if (!(nMG <= FLT_MAX)) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u\n", sw, st);
        (void)fflush(jtr);
#endif /* JTRACE */
        return -__LINE__;
      }
      MG = fmaxf(MG, nMG);
      if (!(nMV <= FLT_MAX)) {
#ifdef JTRACE
        (void)fprintf(jtr, "sweep=%u, step=%u\n", sw, st);
        (void)fflush(jtr);
#endif /* JTRACE */
        return -__LINE__;
      }
      MV = fmaxf(MV, nMV);
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
#pragma omp parallel for default(none) shared(m,n,Gr,ldGr,Gi,ldGi,eS,fS,sT)
#endif /* _OPENMP */
    for (fnat j = 0u; j < *n; ++j) {
      float *const Gr_j = Gr + j * (size_t)(*ldGr);
      float *const Gi_j = Gi + j * (size_t)(*ldGi);
      register const VS _f = _mm512_set1_ps(fS[j]);
      register const VS _s = _mm512_set1_ps(-(eS[j]));
      for (fnat i = 0u; i < *m; i += VSL) {
        float *const Gr_ij = Gr_j + i;
        float *const Gi_ij = Gi_j + i;
        _mm512_store_ps(Gr_ij, _mm512_scalef_ps(_mm512_div_ps(_mm512_load_ps(Gr_ij), _f), _s));
        _mm512_store_ps(Gi_ij, _mm512_scalef_ps(_mm512_div_ps(_mm512_load_ps(Gi_ij), _f), _s));
      }
      eS[j] -= sT;
    }
    if (sV) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,Vr,ldVr,Vi,ldVi,sV)
#endif /* _OPENMP */
      for (fnat j = 0u; j < *n; ++j) {
        float *const Vr_j = Vr + j * (size_t)(*ldVr);
        float *const Vi_j = Vi + j * (size_t)(*ldVi);
        register const VS _s = _mm512_set1_ps((float)-sV);
        for (fnat i = 0u; i < *n; i += VSL) {
          float *const Vr_ij = Vr_j + i;
          float *const Vi_ij = Vi_j + i;
          _mm512_store_ps(Vr_ij, _mm512_scalef_ps(_mm512_load_ps(Vr_ij), _s));
          _mm512_store_ps(Vi_ij, _mm512_scalef_ps(_mm512_load_ps(Vi_ij), _s));
        }
      }
    }
  }

#ifdef JTRACE
  (void)fprintf(jtr, "sT=%d, M=%#.9e\n", sT, MG);
  (void)fprintf(jtr, "Tn=%15.9Lf, Tp=%15.9Lf, Ta=%15.9Lf, Te=%15.9Lf, Tr=%15.9Lf\n", Tn, Tp, Ta, Te, Tr);
  (void)fclose(jtr);
#endif /* JTRACE */
  return (fint)sw;
}
