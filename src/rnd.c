#include "rnd.h"

static const double daub = DBL_MAX / 4;
static const double zaub = DBL_MAX / 8;
static const float  saub = FLT_MAX / 4;
static const float  caub = FLT_MAX / 8;

static void gen_rand8(const size_t n, uint8_t r[static restrict 1])
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

static void gen_rand16(const size_t n, uint16_t r[static restrict 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < n; ++i)
    r[i] = uwrand();
}

static void gen_rand32(const size_t n, uint32_t r[static restrict 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < n; ++i)
    r[i] = udrand();
}

static void gen_rand64(const size_t n, uint64_t r[static restrict 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,r)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < n; ++i)
    r[i] = uqrand();
}

fint gen_rand_(const size_t n[static restrict 1], const size_t s[static restrict 1], void *restrict r)
{
  if (!*n)
    return 0;
  if (!*s)
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

  for (size_t b = ((*n) * (*s)); a && b; a >>= 1u) {
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

void gensfrand_(const size_t n[static restrict 1], const float aub[static restrict 1], float r[static restrict 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,aub,r)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < *n; ++i)
    r[i] = sfrand(*aub);
}

void gendfrand_(const size_t n[static restrict 1], const double aub[static restrict 1], double r[static restrict 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,aub,r)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < *n; ++i)
    r[i] = dfrand(*aub);
}

/* f = c^2 * (l1 + l2 * t^2) */
/* g = c^2 * (l1 * t^2 + l2) */
/* h = c^2 * exp(-alpha * I) * t * (l1 - l2) */

void ssym2rand_(const size_t n[static restrict 1], float l1[static restrict 1], float l2[static restrict 1], float f[static restrict 1], float g[static restrict 1], float h[static restrict 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,saub,l1,l2,f,g,h)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < *n; ++i) {
    wide w1, w2, t, c;
    w1 = l1[i] = sfrand(saub);
    w2 = l2[i] = sfrand(saub);
    wo2rand(&t, &c);
    h[i] = (float)((t * (w1 - w2)) / c);
    t *= t;
    f[i] = (float)(fmaw(w2, t, w1) / c);
    g[i] = (float)(fmaw(w1, t, w2) / c);
  }
}

void dsym2rand_(const size_t n[static restrict 1], double l1[static restrict 1], double l2[static restrict 1], double f[static restrict 1], double g[static restrict 1], double h[static restrict 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,daub,l1,l2,f,g,h)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < *n; ++i) {
    wide w1, w2, t, c;
    w1 = l1[i] = dfrand(daub);
    w2 = l2[i] = dfrand(daub);
    wo2rand(&t, &c);
    h[i] = (double)((t * (w1 - w2)) / c);
    t *= t;
    f[i] = (double)(fmaw(w2, t, w1) / c);
    g[i] = (double)(fmaw(w1, t, w2) / c);
  }
}

void cher2rand_(const size_t n[static restrict 1], float l1[static restrict 1], float l2[static restrict 1], float f[static restrict 1], float g[static restrict 1], float hr[static restrict 1], float hi[static restrict 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,caub,l1,l2,f,g,hr,hi)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < *n; ++i) {
    wide w1, w2, t, c, r, j, h;
    w1 = l1[i] = sfrand(caub);
    w2 = l2[i] = sfrand(caub);
    wu2rand(&t, &c, &r, &j);
    h = (t * (w1 - w2)) / c;
    hr[i] = (float)(r * h);
    hi[i] = (float)(j * h);
    t *= t;
    f[i] = (float)(fmaw(w2, t, w1) / c);
    g[i] = (float)(fmaw(w1, t, w2) / c);
  }
}

void zher2rand_(const size_t n[static restrict 1], double l1[static restrict 1], double l2[static restrict 1], double f[static restrict 1], double g[static restrict 1], double hr[static restrict 1], double hi[static restrict 1])
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,zaub,l1,l2,f,g,hr,hi)
#endif /* _OPENMP */
  for (size_t i = (size_t)0u; i < *n; ++i) {
    wide w1, w2, t, c, r, j, h;
    w1 = l1[i] = dfrand(zaub);
    w2 = l2[i] = dfrand(zaub);
    wu2rand(&t, &c, &r, &j);
    h = (t * (w1 - w2)) / c;
    hr[i] = (double)(r * h);
    hi[i] = (double)(j * h);
    t *= t;
    f[i] = (double)(fmaw(w2, t, w1) / c);
    g[i] = (double)(fmaw(w1, t, w2) / c);
  }
}
