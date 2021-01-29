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

static inline float sfrand()
{
  float f;
  do {
    f = s_rand();
  } while (!isfinite(f));
  return f;
}

static inline double d_rand()
{
  double d;
  while (!_rdrand64_step((uint64_t*)&d)) /**/;
  return d;
}

static inline double dfrand()
{
  double d;
  do {
    d = d_rand();
  } while (!isfinite(d));
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

/* returns a tangent in [-1,1) and the corresponding cosine */
static inline void wo2rand(wide t[static 1], wide c[static 1])
{
  *t = w2rand();
  *c = invsqrtw(fmaw(*t, *t, W_ONE));
}

/* returns a tangent in [-1,1) and the corresponding cosine */
/* & r + I*i = cos(alpha) + I*sin(alpha); alpha in [-pi,pi) */
static inline void wu2rand(wide t[static 1], wide c[static 1], wide r[static 1], wide i[static 1])
{
  wo2rand(t, c);
  const wide a = w2rand() * W_PI;
  sincosw(a, r, i);
}

extern int gen_rand(const size_t n, const size_t s, void *r);

extern void gensrand(const size_t n, float r[static 1]);
extern void gendrand(const size_t n, double r[static 1]);

extern void geno2rand(const size_t n, wide t[static 1], wide c[static 1]);
extern void genu2rand(const size_t n, wide t[static 1], wide c[static 1], wide r[static 1], wide i[static 1]);

#endif /* !RND_H */
