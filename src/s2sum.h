#ifndef S2SUM_H
#define S2SUM_H

// M\o ller's 2Sum
#ifdef TwoSum
#error TwoSum already defined
#else /* !TwoSum */
#define TwoSum(a,b,a_,b_,s,t)\
  s = _mm512_add_ps(a, b);   \
  a_ = _mm512_sub_ps(s, b);  \
  b_ = _mm512_sub_ps(s, a_); \
  a_ = _mm512_sub_ps(a, a_); \
  b_ = _mm512_sub_ps(b, b_); \
  t = _mm512_add_ps(a_, b_)
#endif /* ?TwoSum */

#endif /* !S2SUM_H */
