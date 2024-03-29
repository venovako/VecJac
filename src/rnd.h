#ifndef RND_H
#define RND_H

#include "vec.h"

static inline uint16_t uwrand()
{
  uint16_t w;
  while (!_rdrand16_step(&w)) /**/;
  return w;
}

static inline int16_t iwrand()
{
  int16_t w;
  while (!_rdrand16_step((uint16_t*)&w)) /**/;
  return w;
}

static inline uint32_t udrand()
{
  uint32_t d;
  while (!_rdrand32_step(&d)) /**/;
  return d;
}

static inline int32_t idrand()
{
  int32_t d;
  while (!_rdrand32_step((uint32_t*)&d)) /**/;
  return d;
}

static inline uint64_t uqrand()
{
  uint64_t q;
  while (!_rdrand64_step(&q)) /**/;
  return q;
}

static inline int64_t iqrand()
{
  int64_t q;
  while (!_rdrand64_step((uint64_t*)&q)) /**/;
  return q;
}

static inline float s_rand()
{
  float f;
  while (!_rdrand32_step((uint32_t*)&f)) /**/;
  return f;
}

static inline float sfrand(const float aub)
{
  float f;
  do {
    f = s_rand();
  } while (!(fabsf(f) <= aub));
  return f;
}

static inline double d_rand()
{
  double d;
  while (!_rdrand64_step((uint64_t*)&d)) /**/;
  return d;
}

static inline double dfrand(const double aub)
{
  double d;
  do {
    d = d_rand();
  } while (!(fabs(d) <= aub));
  return d;
}

/* the result is in [0,1) */
static inline wide w1rand()
{
  return scalbw(uqrand(), -64);
}

/* the result is in [-1,1) */
static inline wide w2rand()
{
  return scalbw(iqrand(), -63);
}

/* returns a tangent in [-1,1) and the corresponding secant squared */
static inline void wo2rand(wide t[static restrict 1], wide c[static restrict 1])
{
  *t = w2rand();
  *c = fmaw(*t, *t, W_ONE);
}

/* returns a tangent in [0,1) and the corresponding secant squared */
/* r + I*j = cos(alpha) + I*sin(alpha); alpha in [-pi,pi) */
static inline bool wu2rand(wide t[static restrict 1], wide c[static restrict 1], wide r[static restrict 1], wide j[static restrict 1])
{
  wo2rand(t, c);
  // reduce the precision in the hope of getting 1-cos(alpha)^2 exact
#ifdef USE_EXTENDED
  *r = (float)w2rand();
#else /* !USE_EXTENDED */
  *r = (double)w2rand();
#endif /* ?USE_EXTENDED */
  *j = sqrtw(fmaw(-*r, *r, W_ONE));
  // absorb the negative sign of tan(phi) into e^(alpha*I)
  if (copysignw(W_ONE, *t) != W_ONE) {
    *t = -*t;
    *r = -*r;
    *j = -*j;
    return true;
  }
  return false;
}

extern fint gen_rand_(const size_t n[static restrict 1], const size_t s[static restrict 1], void *restrict r);

extern void gensfrand_(const size_t n[static restrict 1], const float aub[static restrict 1], float r[static restrict 1]);
extern void gendfrand_(const size_t n[static restrict 1], const double aub[static restrict 1], double r[static restrict 1]);

/* [ c -s ] [ l1 0 ] [  c s ] = [ f h ] */
/* [ s  c ] [ 0 l2 ] [ -s c ]   [ h g ] */
extern void ssym2rand_(const size_t n[static restrict 1], float l1[static restrict 1], float l2[static restrict 1], float f[static restrict 1], float g[static restrict 1], float h[static restrict 1]);
extern void dsym2rand_(const size_t n[static restrict 1], double l1[static restrict 1], double l2[static restrict 1], double f[static restrict 1], double g[static restrict 1], double h[static restrict 1]);

/*                                                                       _   */
/* [          c -exp(-a*I)s ] [ l1 0 ] [           c exp(-a*I)*s ] = [ f h ] */
/* [ exp(a*I)*s           c ] [ 0 l2 ] [ -exp(a*I)*s           c ]   [ h g ] */
extern void cher2rand_(const size_t n[static restrict 1], float l1[static restrict 1], float l2[static restrict 1], float f[static restrict 1], float g[static restrict 1], float hr[static restrict 1], float hi[static restrict 1], float hs[static restrict 1]);
extern void zher2rand_(const size_t n[static restrict 1], double l1[static restrict 1], double l2[static restrict 1], double f[static restrict 1], double g[static restrict 1], double hr[static restrict 1], double hi[static restrict 1], double hd[static restrict 1]);

#endif /* !RND_H */
