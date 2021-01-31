#ifndef RND_H
#define RND_H

#include "common.h"

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
static inline void wo2rand(wide t[static 1], wide c[static 1])
{
  *t = w2rand();
  *c = fmaw(*t, *t, W_ONE);
}

/* returns a tangent in [-1,1) and the corresponding secant squared */
/* & r + I*j = cos(alpha) + I*sin(alpha); alpha in [-pi,pi) */
static inline void wu2rand(wide t[static 1], wide c[static 1], wide r[static 1], wide j[static 1])
{
  wo2rand(t, c);
  sincosw((w2rand() * W_PI), j, r);
}

extern int gen_rand(const size_t n, const size_t s, void *r);

extern void gensfrand(const size_t n, const float aub, float r[static 1]);
extern void gendfrand(const size_t n, const double aub, double r[static 1]);

/* [ c -s ] [ l1 0 ] [  c s ] = [ f h ] */
/* [ s  c ] [ 0 l2 ] [ -s c ]   [ h g ] */
extern void ssym2rand(const size_t n, float l1[static 1], float l2[static 1], float f[static 1], float g[static 1], float h[static 1]);
extern void dsym2rand(const size_t n, double l1[static 1], double l2[static 1], double f[static 1], double g[static 1], double h[static 1]);

/* [          c -exp(-a*I)s ] [ l1 0 ] [           c exp(-a*I)*s ] [ f       h ] */
/* [ exp(a*I)*s           c ] [ 0 l2 ] [ -exp(a*I)*s           c ] [ conj(h) g ] */
extern void cher2rand(const size_t n, float l1[static 1], float l2[static 1], float f[static 1], float g[static 1], float hr[static 1], float hi[static 1]);
extern void zher2rand(const size_t n, double l1[static 1], double l2[static 1], double f[static 1], double g[static 1], double hr[static 1], double hi[static 1]);

#endif /* !RND_H */
