#ifndef DZNRMF_H
#define DZNRMF_H

static inline void ldbl2ef(const long double l, double e[static restrict 1], double f[static restrict 1])
{
  if (isfinite(l)) {
    int le = 0;
    const double df = (double)frexpl(l, &le);
    if (df == 0.0) {
      *e = -HUGE_VAL;
      *f = copysign(1.0, df);
    }
    else {
      *e = (le - 1);
      *f = scalbn(df, 1);
    }
  }
  else {
    *e = HUGE_VAL;
    *f = (double)copysignl((isinf(l) ? 1.0L : 0.0L), l);
  }
}

#endif /* !DZNRMF_H */
