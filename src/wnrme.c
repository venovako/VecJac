#include "wnrme.h"

// fast Euclidean norms of 2x2 Hermitian matrices in a wider precision (assuming no over/under-flows)

wide wnrmer(const wide a11, const wide a22, const wide a21)
{
  return sqrtw(fmaw(a11, a11, fmaw(a22, a22, (W_TWO * a21 * a21))));
}

wide wnrmec(const wide a11, const wide a22, const wide a21r, const wide a21i)
{
  return sqrtw(fmaw(a11, a11, fmaw(a22, a22, (W_TWO * fmaw(a21r, a21r, (a21i * a21i))))));
}
