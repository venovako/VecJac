#include "rnd.h"

static void gen_rand8(const size_t n, uint8_t r[static 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < n; i += (size_t)2u) {
    const uint16_t w = uwrand();
    r[i] = (uint8_t)w;
    const size_t j = i + (size_t)1u;
    if (j < n)
      r[j] = (uint8_t)(w >> 8u);
  }
}

static void gen_rand16(const size_t n, uint16_t r[static 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < n; ++i)
    r[i] = uwrand();
}

static void gen_rand32(const size_t n, uint32_t r[static 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < n; ++i)
    r[i] = udrand();
}

static void gen_rand64(const size_t n, uint64_t r[static 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < n; ++i)
    r[i] = uqrand();
}

int gen_rand(const size_t n, const size_t s, void *r)
{
  if (!n)
    return 0;
  if (!s)
    return -2;
  if (!r)
    return -3;

  const uintptr_t p = (uintptr_t)r;
  unsigned a = 1u;
  do {
    if (p & (uintptr_t)a)
      break;
    a <<= 1u;
  } while (a <= 4u);

  for (size_t b = (n * s); a && b; a >>= 1u) {
    const size_t e = b / a;
    if (e) {
      switch (a) {
      case 8u:
        gen_rand64(e, (uint64_t*)r);
        r = ((uint64_t*)r + e);
        break;
      case 4u:
        gen_rand32(e, (uint32_t*)r);
        r = ((uint32_t*)r + e);
        break;
      case 2u:
        gen_rand16(e, (uint16_t*)r);
        r = ((uint16_t*)r + e);
        break;
      case 1u:
        gen_rand8(e, (uint8_t*)r);
        r = ((uint8_t*)r + e);
        break;
      default: /* should never happen */
        return -1;
      }
      b -= (e * a);
    }
  }

  return 0;
}

void gensrand(const size_t n, float r[static 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < n; ++i)
    r[i] = sfrand();
}

void gendrand(const size_t n, double r[static 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < n; ++i)
    r[i] = dfrand();
}
