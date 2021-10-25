#ifndef DZNRMF_H
#define DZNRMF_H

static inline void ldbl2ef(const long double l, double e[static restrict 1], double f[static restrict 1])
{
  int le = 0;
  const long double lf = frexpl(l, &le);
  if (lf == 0.0L) {
    *e = -HUGE_VAL;
    *f = copysign(1.0, (double)lf);
  }
  else {
    *e = (le - 1);
    *f = scalbn((double)lf, 1);
  }
}

#endif /* !DZNRMF_H */
