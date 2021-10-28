#include "dnormf.h"

#include "dznrmf.h"

double dnormf_(const fnat m[static restrict 1], const double x[static restrict 1], double e0[static restrict 1], double f0[static restrict 1], double e1[static restrict 1], double f1[static restrict 1])
{
  long double sq = 0.0L;
  for (fnat i = 0u; i < *m; ++i) {
    const long double l = x[i];
    sq += (l * l);
  }
  ldbl2ef(sq, e1, f1);
  sq = sqrtl(sq);
  ldbl2ef(sq, e0, f0);
  return (double)sq;
}
