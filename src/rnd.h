#ifndef RND_H
#define RND_H

#include "common.h"

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

extern int gen_rand(const size_t n, const size_t s, void *r);

extern void gensrand(const size_t n, float r[static 1]);
extern void gendrand(const size_t n, double r[static 1]);

#endif /* !RND_H */
