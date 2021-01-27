#ifndef RND_H
#define RND_H

#include "common.h"

static inline uint16_t uwrand()
{
  uint16_t w;
  while (!_rdrand16_step(&w)) /**/;
  return w;
}

static inline uint32_t udrand()
{
  uint32_t d;
  while (!_rdrand32_step(&d)) /**/;
  return d;
}

static inline uint64_t uqrand()
{
  uint64_t q;
  while (!_rdrand64_step(&q)) /**/;
  return q;
}

static inline float sfrand()
{
  float f;
  while (!_rdrand32_step((uint32_t*)&f) || !isfinite(f)) /**/;
  return f;
}

static inline double dfrand()
{
  double d;
  while (!_rdrand64_step((uint64_t*)&d) || !isfinite(d)) /**/;
  return d;
}

/* the result is in [0,1) */
static inline wide w1rand()
{
  return scalbw(uqrand(), -64);
}

/* returns a tangent in [0,1) and the corresponding cosine */
static inline void wrot2rand(wide t[static 1], wide c[static 1])
{
  *t = w1rand();
  *c = invsqrtw(fmaw(*t, *t, W_ONE));
}

extern int gen_rand(const size_t n, const size_t s, void *r);

extern void gensrand(const size_t n, float r[static 1]);
extern void gendrand(const size_t n, double r[static 1]);
extern void genwrand(const size_t n, wide r[static 1]);

extern void gen2rand(const size_t n, wide t[static 1], wide c[static 1]);

#endif /* !RND_H */
