#ifndef DZNRM2_H
#define DZNRM2_H

static inline void dbl2ef(const double d, double e[static restrict 1], double f[static restrict 1])
{
  if (isfinite(d)) {
    int de = 0;
    const double df = frexp(d, &de);
    if (df == 0.0) {
      *e = -HUGE_VAL;
      *f = copysign(1.0, df);
    }
    else {
      *e = (de - 1);
      *f = scalbn(df, 1);
    }
  }
  else {
    *e = HUGE_VAL;
    *f = copysign((isinf(d) ? 1.0 : 0.0), d);
  }
}

#endif /* !DZNRM2_H */
