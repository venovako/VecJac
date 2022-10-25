#ifndef D2SUM_H
#define D2SUM_H

// M\o ller's 2Sum
#ifdef TwoSum
#error TwoSum already defined
#else /* !TwoSum */
#define TwoSum(a,b,a_,b_,s,t)\
  s = _mm512_add_pd(a, b);   \
  a_ = _mm512_sub_pd(s, b);  \
  b_ = _mm512_sub_pd(s, a_); \
  a_ = _mm512_sub_pd(a, a_); \
  b_ = _mm512_sub_pd(b, b_); \
  t = _mm512_add_pd(a_, b_)
#endif /* ?TwoSum */

#endif /* !D2SUM_H */
