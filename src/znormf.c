#include "znormf.h"

#include "dznrmf.h"

double znormf_(const fnat m[static restrict 1], const double zr[static restrict 1], const double zi[static restrict 1], double e0[static restrict 1], double f0[static restrict 1], double e1[static restrict 1], double f1[static restrict 1])
{
  long double sqr = 0.0L;
  for (fnat i = 0u; i < *m; ++i) {
    const long double lr = zr[i];
    sqr += (lr * lr);
  }
  long double sqi = 0.0L;
  for (fnat i = 0u; i < *m; ++i) {
    const long double li = zi[i];
    sqi += (li * li);
  }
  long double sq = (sqr + sqi);
  ldbl2ef(sq, e1, f1);
  sq = sqrtl(sq);
  ldbl2ef(sq, e0, f0);
  return (double)sq;
}
