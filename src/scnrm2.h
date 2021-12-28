#ifndef SCNRM2_H
#define SCNRM2_H

static inline void flt2ef(const float d, float e[static restrict 1], float f[static restrict 1])
{
  if (isfinite(d)) {
    int de = 0;
    const float df = frexpf(d, &de);
    if (df == 0.0f) {
      *e = -HUGE_VALF;
      *f = copysignf(1.0f, df);
    }
    else {
      *e = (float)(de - 1);
      *f = scalbnf(df, 1);
    }
  }
  else {
    *e = HUGE_VALF;
    *f = copysignf((isinf(d) ? 1.0f : 0.0f), d);
  }
}

static inline void sqeff(const float e0[static restrict 1], const float f0[static restrict 1], float e1[static restrict 1], float f1[static restrict 1])
{
  *e1 = scalbnf(*e0, 1);
  *f1 = ((*f0) * (*f0));
  if (*f1 >= 2.0f) {
    ++*e1;
    *f1 = scalbnf(*f1, -1);
  }
}

#endif /* !SCNRM2_H */
