#ifndef DZNRMS_H
#define DZNRMS_H

#include "vecdef.h"
#include "d8sort.h"

#include "sleefquad.h"

#ifndef Q_ZERO
#define Q_ZERO 0.0q
#endif /* !Q_ZERO */

#ifndef Q_PMAX
#define Q_PMAX SLEEF_QUAD_MAX
#endif /* !Q_PMAX */

static inline void pquad2ef(const __float128 q[static restrict 1], double e[static restrict 1], double f[static restrict 1])
{
  if (*q <= Q_PMAX) {
    int qe = 0;
#if (defined(__ICC) || defined(__INTEL_COMPILER) || defined(__INTEL_CLANG_COMPILER) || defined(__INTEL_LLVM_COMPILER))
    const double df = (double)__frexpq(*q, &qe);
#else /* !Intel */
    const double df = (double)frexpq(*q, &qe);
#endif /* ?Intel */
    if (df == 0.0) {
      *e = -HUGE_VAL;
      *f = 1.0;
    }
    else {
      *e = (qe - 1);
      *f = scalbn(df, 1);
    }
  }
  else {
    *e = HUGE_VAL;
    *f = ((*q == *q) ? 1.0 : 0.0);
  }
}

#endif /* !DZNRMS_H */
