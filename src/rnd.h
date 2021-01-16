#ifndef RND_H
#define RND_H

#include "common.h"

extern int gen_rand(const size_t n, const size_t s, void *r);

extern void gensrand(const size_t n, float r[static 1]);
extern void gendrand(const size_t n, double r[static 1]);

#endif /* !RND_H */
