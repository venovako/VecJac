#include "znormf.h"

#include "dznrmf.h"

double znormf_(const fnat m[static restrict 1], const double zr[static restrict 1], const double zi[static restrict 1], double e0[static restrict 1], double f0[static restrict 1], double e1[static restrict 1], double f1[static restrict 1])
{
  long double sq = 0.0L;
  for (fnat i = 0u; i < *m; ++i) {
    const long double l = zr[i];
    sq += (l * l);
  }
  for (fnat i = 0u; i < *m; ++i) {
    const long double l = zi[i];
    sq += (l * l);
  }
  if (isfinite(sq)) {
    ldbl2ef(sq, e1, f1);
    sq = sqrtl(sq);
    ldbl2ef(sq, e0, f0);
  }
  else {
    *e1 = *e0 = HUGE_VAL;
    *f1 = *f0 = (isinf(sq) ? 1.0 : 0.0);
  }
  return (double)sq;
}
