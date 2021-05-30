#include "wnrme.h"

// fast Euclidean norms of Hermitian matrices of order two in a wider precision (assuming no over/under-flows)

wide wnrmer(const wide a11, const wide a22, const wide a21)
{
  return sqrtw(fmaw(a11, a11, fmaw(a22, a22, (W_TWO * a21 * a21))));
}

wide wnrmec(const wide a11, const wide a22, const wide a21r, const wide a21i)
{
  return sqrtw(fmaw(a11, a11, fmaw(a22, a22, (W_TWO * fmaw(a21r, a21r, (a21i * a21i))))));
}

// absolute error checkers

static wide waer(const wide a11, const wide a22, const wide a21, const wide s, const wide t, const wide c, const wide l1, const wide l2, wide L1[static restrict 1], wide L2[static restrict 1])
{
  const wide _s = -s;
  const wide c2 = c * c;
  const wide t2 = t * t;
  *L1 = c2 * scalbw(l1, _s);
  *L2 = c2 * scalbw(l2, _s);
  const wide l = *L1 - *L2;
  const wide tl = t * l;
  return wnrmer(fmaw(c2, fmaw(*L2, t2, *L1), -a11), fmaw(c2, fmaw(*L1, t2, *L2), -a22), fmaw(c2, tl, -a21));
}

static wide waec(const wide a11, const wide a22, const wide a21r, const wide a21i, const wide s, const wide t, const wide c, const wide ca, const wide sa, const wide l1, const wide l2, wide L1[static restrict 1], wide L2[static restrict 1])
{
  const wide _s = -s;
  const wide c2 = c * c;
  const wide t2 = t * t;
  *L1 = c2 * scalbw(l1, _s);
  *L2 = c2 * scalbw(l2, _s);
  const wide l = *L1 - *L2;
  const wide tl = t * l;
  return wnrmec(fmaw(c2, fmaw(*L2, t2, *L1), -a11), fmaw(c2, fmaw(*L1, t2, *L2), -a22), fmaw(c2, (ca * tl), -a21r), fmaw(c2, (sa * tl), -a21i));
}

// relative error checkers

wide wrer(const wide a11, const wide a22, const wide a21, const wide s, const wide t, const wide c, const wide l1, const wide l2, wide ae[static restrict 1], wide L1[static restrict 1], wide L2[static restrict 1])
{
  *ae = waer(a11, a22, a21, s, t, c, l1, l2, L1, L2);
  return fmaxw((*ae / wnrmer(a11, a22, a21)), W_ZERO);
}

wide wrec(const wide a11, const wide a22, const wide a21r, const wide a21i, const wide s, const wide t, const wide c, const wide ca, const wide sa, const wide l1, const wide l2, wide ae[static restrict 1], wide L1[static restrict 1], wide L2[static restrict 1])
{
  *ae = waec(a11, a22, a21r, a21i, s, t, c, ca, sa, l1, l2, L1, L2);
  return fmaxw((*ae / wnrmec(a11, a22, a21r, a21i)), W_ZERO);
}
