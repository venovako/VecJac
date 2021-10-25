#include "wdp.h"
#include "dnorm2.h"
#include "dnorme.h"
#include "dnormf.h"
#include "dnorms.h"
#include "d8sort.h"

double dnb(const fnat n, const double x[static restrict 1])
{
  const fint incx = 1;
  return dnorm2_((const fint*)&n, x, &incx);
}

#ifndef USE_MKL
extern double BLAS_D(dot)(const fint n[static 1], const double x[static 1], const fint incx[static 1], const double y[static 1], const fint incy[static 1]);
#endif /* !USE_MKL */

double ddp(const fnat n, const double x[static restrict 1])
{
  const fint incx = 1;
  return sqrt(BLAS_D(dot)((const fint*)&n, x, &incx, x, &incx));
}

#ifndef USE_MKL
extern double BLAS_D(nrm2)(const fint n[static 1], const double x[static 1], const fint incx[static 1]);
#endif /* !USE_MKL */

double dn2(const fnat n, const double x[static restrict 1])
{
  const fint incx = 1;
  return BLAS_D(nrm2)((const fint*)&n, x, &incx);
}

double dnf(const fnat n, const double x[static restrict 1])
{
  double e0, f0, e1, f1;
  return dnormf_(&n, x, &e0, &f0, &e1, &f1);
}

double dne(const fnat n, const double x[static restrict VDL])
{
  double e0, f0, e1, f1;
  return dnorme_(&n, x, &e0, &f0, &e1, &f1);
}

double dns(const fnat n, const double x[static restrict VDL])
{
  double e0, f0, e1, f1;
  return dnorms_(&n, x, &e0, &f0, &e1, &f1);
}

double dnc(const fnat n, const double x[static restrict VDL])
{
#ifdef DZNRME_SEQRED
  alignas(VA) double sq[VDL]
#ifndef NDEBUG
    = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
#endif /* !NDEBUG */
    ;
#endif /* DZNRME_SEQRED */
  register __m512d vsq = _mm512_setzero_pd();
  for (fnat i = 0u; i < n; i += VDL) {
    register __m512d xi = _mm512_load_pd(x + i);
    vsq = _mm512_fmadd_pd(xi, xi, vsq);
  }
  VDSORT(vsq);
#ifdef DZNRME_SEQRED
  _mm512_store_pd(sq, vsq);
  sq[0u] += sq[1u];
  sq[2u] += sq[3u];
  sq[4u] += sq[5u];
  sq[6u] += sq[7u];
  sq[0u] += sq[2u];
  sq[4u] += sq[6u];
  sq[0u] += sq[4u];
  return sqrt(*sq);
#else /* !DZNRME_SEQRED */
  vsq = _mm512_permutexvar_pd(_mm512_set_epi64(7, 6, 3, 2, 5, 4, 1, 0), vsq);
  register const __m256d sq4 = _mm256_hadd_pd(_mm512_extractf64x4_pd(vsq, 0), _mm512_extractf64x4_pd(vsq, 1));
  register const __m128d sq2 = _mm_hadd_pd(_mm256_extractf128_pd(sq4, 0), _mm256_extractf128_pd(sq4, 1));
  register const __m128 sqs = _mm_castpd_ps(sq2);
  register const __m128d sqd = _mm_castps_pd(_mm_movehl_ps(sqs, sqs));
  return sqrt(_mm_cvtsd_f64(sq2) + _mm_cvtsd_f64(sqd));
#endif /* ?DZNRME_SEQRED */
}

double dnd(const fnat n, const double x[static restrict VDL])
{
#ifdef DZNRMS_SEQRED
  alignas(VA) double sq[VDL]
#ifndef NDEBUG
    = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
#endif /* !NDEBUG */
    ;
#endif /* DZNRMS_SEQRED */
  register __m512d vsq = _mm512_setzero_pd();
  for (fnat i = 0u; i < n; i += VDL) {
    register __m512d xi = _mm512_load_pd(x + i);
    VDSORT(xi);
    vsq = _mm512_fmadd_pd(xi, xi, vsq);
  }
#ifdef DZNRMS_SEQRED
  _mm512_store_pd(sq, vsq);
  sq[0u] += sq[1u];
  sq[2u] += sq[3u];
  sq[4u] += sq[5u];
  sq[6u] += sq[7u];
  sq[0u] += sq[2u];
  sq[4u] += sq[6u];
  sq[0u] += sq[4u];
  return sqrt(*sq);
#else /* !DZNRMS_SEQRED */
  vsq = _mm512_permutexvar_pd(_mm512_set_epi64(7, 6, 3, 2, 5, 4, 1, 0), vsq);
  register const __m256d sq4 = _mm256_hadd_pd(_mm512_extractf64x4_pd(vsq, 0), _mm512_extractf64x4_pd(vsq, 1));
  register const __m128d sq2 = _mm_hadd_pd(_mm256_extractf128_pd(sq4, 0), _mm256_extractf128_pd(sq4, 1));
  register const __m128 sqs = _mm_castpd_ps(sq2);
  register const __m128d sqd = _mm_castps_pd(_mm_movehl_ps(sqs, sqs));
  return sqrt(_mm_cvtsd_f64(sq2) + _mm_cvtsd_f64(sqd));
#endif /* ?DZNRMS_SEQRED */
}

double wsq(const fnat n, const double x[static restrict 1])
{
  wide sq = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,x) reduction(+:sq)
#endif /* !_OPENMP */
  for (fnat i = 0u; i < n; ++i) {
    const wide w = x[i];
    sq = fmaw(w, w, sq);
  }
  return (double)sqrtw(sq);
}

double dre(const double c, const double e)
{
  const double a = fabs(c - e);
  return ((e == 0.0) ? ((a == 0.0) ? 0.0 : (a / e)) : (a / e));
}
