#include "zjrot.h"

#include "djrot.h"
#include "vecdef.h"

int zjrot_(const fint n[static restrict 1], double xr[static restrict VDL], double xi[static restrict VDL], double yr[static restrict VDL], double yi[static restrict VDL], const double t[static restrict 1], const double c[static restrict 1], const double ca[static restrict 1], const double sa[static restrict 1])
{
  if (IS_NOT_VFPENV)
    return -10;
  if (*n & VDL_1)
    return -1;
  if (IS_NOT_ALIGNED(xr))
    return -2;
  if (IS_NOT_ALIGNED(xi))
    return -3;
  if (IS_NOT_ALIGNED(yr))
    return -4;
  if (IS_NOT_ALIGNED(yi))
    return -5;

  if (!*n)
    return 0;
  if ((*sa == 0.0) && (fabs(*ca) == 1.0)) {
    // real rotation
    const double t_ = ((*ca == 1.0) ? *t : -*t);
    const fint r = djrot_(n, xr, yr, &t_, c);
    const fint i = djrot_(n, xi, yi, &t_, c);
    return ((r <= i) ? r : i);
  }
  if (*t == 0.0) {
    if (*c == 1.0)
      return 0;
    // should never happen
    return -7;
  }

  alignas(16) double cstc[4];
  _mm_store_pd(cstc, _mm_mul_pd(_mm_set_pd(*sa, *ca), _mm_set1_pd(*t)));
  cstc[2] = -(cstc[0]);
  cstc[3] = *c;

  if (*n < 0) { // permute
    const fnat n_ = (fnat)-*n;
#ifdef _OPENMP
#pragma omp parallel default(none) shared(n_,cstc,xr,xi,yr,yi)
#endif /* _OPENMP */
    {
      register const VD cta = _mm512_set1_pd(cstc[0]);
      register const VD sta = _mm512_set1_pd(cstc[1]);
      register const VD cna = _mm512_set1_pd(cstc[2]);
      register const VD c_ = _mm512_set1_pd(cstc[3]);
#ifdef _OPENMP
#pragma omp for
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xri = xr + i;
        double *const xii = xi + i;
        double *const yri = yr + i;
        double *const yii = yi + i;
        register const VD x_r = _mm512_load_pd(xri);
        register const VD x_i = _mm512_load_pd(xii);
        register const VD y_r = _mm512_load_pd(yri);
        register const VD y_i = _mm512_load_pd(yii);
        register VD _r, _i;
        VZFMA(_r,_i,y_r,y_i,cta,sta,x_r,x_i);
        _mm512_store_pd(yri, _mm512_mul_pd(_r, c_));
        _mm512_store_pd(yii, _mm512_mul_pd(_i, c_));
        VZFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        _mm512_store_pd(xri, _mm512_mul_pd(_r, c_));
        _mm512_store_pd(xii, _mm512_mul_pd(_i, c_));
      }
    }
  }
  else { // no permute
    const fnat n_ = (fnat)*n;
#ifdef _OPENMP
#pragma omp parallel default(none) shared(n_,cstc,xr,xi,yr,yi)
#endif /* _OPENMP */
    {
      register const VD cta = _mm512_set1_pd(cstc[0]);
      register const VD sta = _mm512_set1_pd(cstc[1]);
      register const VD cna = _mm512_set1_pd(cstc[2]);
      register const VD c_ = _mm512_set1_pd(cstc[3]);
#ifdef _OPENMP
#pragma omp for
#endif /* _OPENMP */
      for (fnat i = 0u; i < n_; i += VDL) {
        double *const xri = xr + i;
        double *const xii = xi + i;
        double *const yri = yr + i;
        double *const yii = yi + i;
        register const VD x_r = _mm512_load_pd(xri);
        register const VD x_i = _mm512_load_pd(xii);
        register const VD y_r = _mm512_load_pd(yri);
        register const VD y_i = _mm512_load_pd(yii);
        register VD _r, _i;
        VZFMA(_r,_i,y_r,y_i,cta,sta,x_r,x_i);
        _mm512_store_pd(xri, _mm512_mul_pd(_r, c_));
        _mm512_store_pd(xii, _mm512_mul_pd(_i, c_));
        VZFMA(_r,_i,x_r,x_i,cna,sta,y_r,y_i);
        _mm512_store_pd(yri, _mm512_mul_pd(_r, c_));
        _mm512_store_pd(yii, _mm512_mul_pd(_i, c_));
      }
    }
  }

  return 0;
}
