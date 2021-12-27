#ifndef VECDEF_H
#define VECDEF_H

#include "vec.h"

// assumes _zerof = _mm512_set1_ps(-0.0f)
#ifdef VSABS
#error VSABS already defined
#else /* !VSABS */
#define VSABS(x) VSANDNOT(_zerof,(x))
#endif /* ?VSABS */
#ifdef VSNEG
#error VSNEG already defined
#else /* !VSNEG */
#define VSNEG(x) VSXOR((x),_zerof)
#endif /* ?VSNEG */
#ifdef VSSGN
#error VSSGN already defined
#else /* !VSSGN */
#define VSSGN(x) VSAND((x),_zerof)
#endif /* ?VSSGN */
#ifdef VCFMA
#error VCFMA already defined
#else /* !VCFMA */
#define VCFMA(dr,di,ar,ai,br,bi,cr,ci)                  \
  dr=_mm512_fmadd_ps(ar,br,_mm512_fnmadd_ps(ai,bi,cr)); \
  di=_mm512_fmadd_ps(ar,bi,_mm512_fmadd_ps(ai,br,ci))
#endif /* ?VCFMA */
#ifdef VSHYPOT
#error VSHYPOT already defined
#else /* !VSHYPOT */
#define VSHYPOT(z,x,y)                  \
  {                                     \
    const VS x_ = VSABS(x);             \
    const VS y_ = VSABS(y);             \
    const VS m_ = _mm512_min_ps(x_,y_); \
    const VS M_ = _mm512_max_ps(x_,y_); \
    z = _mm512_div_ps(m_,M_);           \
    z = _mm512_max_ps(z,zerof);         \
    z = _mm512_fmadd_ps(z,z,onef);      \
    z = _mm512_sqrt_ps(z);              \
    z = _mm512_mul_ps(z,M_);            \
  }
#endif /* ?VSHYPOT */

// assumes _zero = _mm512_set1_pd(-0.0)
#ifdef VDABS
#error VDABS already defined
#else /* !VDABS */
#define VDABS(x) VDANDNOT(_zero,(x))
#endif /* ?VDABS */
#ifdef VDNEG
#error VDNEG already defined
#else /* !VDNEG */
#define VDNEG(x) VDXOR((x),_zero)
#endif /* ?VDNEG */
#ifdef VDSGN
#error VDSGN already defined
#else /* !VDSGN */
#define VDSGN(x) VDAND((x),_zero)
#endif /* ?VDSGN */
#ifdef VZFMA
#error VZFMA already defined
#else /* !VZFMA */
#define VZFMA(dr,di,ar,ai,br,bi,cr,ci)                  \
  dr=_mm512_fmadd_pd(ar,br,_mm512_fnmadd_pd(ai,bi,cr)); \
  di=_mm512_fmadd_pd(ar,bi,_mm512_fmadd_pd(ai,br,ci))
#endif /* ?VZFMA */
#ifdef VDHYPOT
#error VDHYPOT already defined
#else /* !VDHYPOT */
#define VDHYPOT(z,x,y)                  \
  {                                     \
    const VD x_ = VDABS(x);             \
    const VD y_ = VDABS(y);             \
    const VD m_ = _mm512_min_pd(x_,y_); \
    const VD M_ = _mm512_max_pd(x_,y_); \
    z = _mm512_div_pd(m_,M_);           \
    z = _mm512_max_pd(z,zero);          \
    z = _mm512_fmadd_pd(z,z,one);       \
    z = _mm512_sqrt_pd(z);              \
    z = _mm512_mul_pd(z,M_);            \
  }
#endif /* ?VDHYPOT */

#endif /* !VECDEF_H */
