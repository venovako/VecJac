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

static inline void sqef(const double e0[static restrict 1], const double f0[static restrict 1], double e1[static restrict 1], double f1[static restrict 1])
{
  *e1 = scalbn(*e0, 1);
  *f1 = ((*f0) * (*f0));
  if (*f1 >= 2.0) {
    ++*e1;
    *f1 = scalbn(*f1, -1);
  }
}

#endif /* !DZNRM2_H */
