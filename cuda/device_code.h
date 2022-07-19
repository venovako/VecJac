#ifndef DEVICE_CODE_H
#define DEVICE_CODE_H

#ifdef FLT_BIG_EXP
#error FLT_BIG_EXP already defined
#else /* !FLT_BIG_EXP */
// FLT_MAX_EXP - 3
#define FLT_BIG_EXP 125
#endif /* ?FLT_BIG_EXP */

#ifdef FLT_SQRT_HUGE
#error FLT_SQRT_HUGE already defined
#else /* !FLT_SQRT_HUGE */
#define FLT_SQRT_HUGE 1.844674297e+19f
#endif /* ?FLT_SQRT_HUGE */

#ifdef DBL_BIG_EXP
#error DBL_BIG_EXP already defined
#else /* !DBL_BIG_EXP */
// DBL_MAX_EXP - 3
#define DBL_BIG_EXP 1021
#endif /* ?DBL_BIG_EXP */

#ifdef DBL_SQRT_HUGE
#error DBL_SQRT_HUGE already defined
#else /* !DBL_SQRT_HUGE */
#define DBL_SQRT_HUGE 1.34078079299425956E+154
#endif /* ?DBL_SQRT_HUGE */

// in the kernels below, -_es could be returned instead of _es

__global__ void ccjac2(const float *const a11, const float *const a22, const float *const a21r, const float *const a21i, float *const cs, float *const snr, float *const sni, float *const l1, float *const l2, int *const sp)
{
  const size_t _off = (size_t)(blockIdx.x) * blockDim.x + threadIdx.x;
  const float _a11 = a11[_off];
  const float _a22 = a22[_off];
  const float _a21r = a21r[_off];
  const float _a21i = a21i[_off];
  int _e1, _e2, _er, _ei;
  (void)frexpf(fmaxf(fabsf(_a11), FLT_TRUE_MIN), &_e1);
  _e1 = FLT_BIG_EXP - _e1;
  (void)frexpf(fmaxf(fabsf(_a22), FLT_TRUE_MIN), &_e2);
  _e2 = FLT_BIG_EXP - _e2;
  (void)frexpf(fmaxf(fabsf(_a21r), FLT_TRUE_MIN), &_er);
  _er = FLT_BIG_EXP - _er;
  (void)frexpf(fmaxf(fabsf(_a21i), FLT_TRUE_MIN), &_ei);
  _ei = FLT_BIG_EXP - _ei;
  int _es = min(min(_e1, _e2), min(_er, _ei));
  const float _ar = scalbnf(_a21r, _es);
  const float _ai = scalbnf(_a21i, _es);
  float _ar_ = fabsf(_ar);
  float _ai_ = fabsf(_ai);
  const float _a1 = scalbnf(_a11, _es);
  const float _a2 = scalbnf(_a22, _es);
  const float _am = fminf(_ar_, _ai_);
  const float _aM = fmaxf(_ar_, _ai_);
  int _sp = (_es << 1);
  _es = -_es;
  float _aa = fmaxf(__fdiv_rn(_am, _aM), 0.0f);
  _aa = __fmul_rn(__fsqrt_rn(__fmaf_rn(_aa, _aa, 1.0f)), _aM);
  _ar_ = copysignf(fminf(__fdiv_rn(_ar_, _aa), 1.0f), _ar);
  _ai_ = __fdiv_rn(_ai, fmaxf(_aa, FLT_TRUE_MIN));
  const float _an = __fmul_rn(_aa, 2.0f);
  const float _ad = __fsub_rn(_a1, _a2);
  const float _t2 = copysignf(fminf(fmaxf(__fdiv_rn(_an, fabsf(_ad)), 0.0f), FLT_SQRT_HUGE), _ad);
  const float _t1 = __fdiv_rn(_t2, __fadd_rn(1.0f, __fsqrt_rn(__fmaf_rn(_t2, _t2, 1.0f))));
  const float _s2 = __fmaf_rn(_t1, _t1, 1.0f);
  const float _s1 = __fsqrt_rn(_s2);
  const float _cs = __frcp_rn(_s1);
  const float _snr = __fdiv_rn(__fmul_rn(_ar_, _t1), _s1);
  const float _sni = __fdiv_rn(__fmul_rn(_ai_, _t1), _s1);
  float _l1 = __fdiv_rn(__fmaf_rn(_t1, __fmaf_rn(_a2, _t1,  _an), _a1), _s2);
  float _l2 = __fdiv_rn(__fmaf_rn(_t1, __fmaf_rn(_a1, _t1, -_an), _a2), _s2);
  _sp |= (_l1 < _l2);
  _l1 = scalbnf(_l1, _es);
  _l2 = scalbnf(_l2, _es);
  cs[_off] = _cs;
  snr[_off] = _snr;
  sni[_off] = _sni;
  l1[_off] = _l1;
  l2[_off] = _l2;
  sp[_off] = _sp;
}

__global__ void zcjac2(const double *const a11, const double *const a22, const double *const a21r, const double *const a21i, double *const cs, double *const snr, double *const sni, double *const l1, double *const l2, int *const sp)
{
  const size_t _off = (size_t)(blockIdx.x) * blockDim.x + threadIdx.x;
  const double _a11 = a11[_off];
  const double _a22 = a22[_off];
  const double _a21r = a21r[_off];
  const double _a21i = a21i[_off];
  int _e1, _e2, _er, _ei;
  (void)frexp(fmax(fabs(_a11), DBL_TRUE_MIN), &_e1);
  _e1 = DBL_BIG_EXP - _e1;
  (void)frexp(fmax(fabs(_a22), DBL_TRUE_MIN), &_e2);
  _e2 = DBL_BIG_EXP - _e2;
  (void)frexp(fmax(fabs(_a21r), DBL_TRUE_MIN), &_er);
  _er = DBL_BIG_EXP - _er;
  (void)frexp(fmax(fabs(_a21i), DBL_TRUE_MIN), &_ei);
  _ei = DBL_BIG_EXP - _ei;
  int _es = min(min(_e1, _e2), min(_er, _ei));
  const double _ar = scalbn(_a21r, _es);
  const double _ai = scalbn(_a21i, _es);
  double _ar_ = fabs(_ar);
  double _ai_ = fabs(_ai);
  const double _a1 = scalbn(_a11, _es);
  const double _a2 = scalbn(_a22, _es);
  const double _am = fmin(_ar_, _ai_);
  const double _aM = fmax(_ar_, _ai_);
  int _sp = (_es << 1);
  _es = -_es;
  double _aa = fmax(__ddiv_rn(_am, _aM), 0.0);
  _aa = __dmul_rn(__dsqrt_rn(__fma_rn(_aa, _aa, 1.0)), _aM);
  _ar_ = copysign(fmin(__ddiv_rn(_ar_, _aa), 1.0), _ar);
  _ai_ = __ddiv_rn(_ai, fmax(_aa, DBL_TRUE_MIN));
  const double _an = __dmul_rn(_aa, 2.0);
  const double _ad = __dsub_rn(_a1, _a2);
  const double _t2 = copysign(fmin(fmax(__ddiv_rn(_an, fabs(_ad)), 0.0), DBL_SQRT_HUGE), _ad);
  const double _t1 = __ddiv_rn(_t2, __dadd_rn(1.0, __dsqrt_rn(__fma_rn(_t2, _t2, 1.0))));
  const double _s2 = __fma_rn(_t1, _t1, 1.0);
  const double _s1 = __dsqrt_rn(_s2);
  const double _cs = __drcp_rn(_s1);
  const double _snr = __ddiv_rn(__dmul_rn(_ar_, _t1), _s1);
  const double _sni = __ddiv_rn(__dmul_rn(_ai_, _t1), _s1);
  double _l1 = __ddiv_rn(__fma_rn(_t1, __fma_rn(_a2, _t1,  _an), _a1), _s2);
  double _l2 = __ddiv_rn(__fma_rn(_t1, __fma_rn(_a1, _t1, -_an), _a2), _s2);
  _sp |= (_l1 < _l2);
  _l1 = scalbn(_l1, _es);
  _l2 = scalbn(_l2, _es);
  cs[_off] = _cs;
  snr[_off] = _snr;
  sni[_off] = _sni;
  l1[_off] = _l1;
  l2[_off] = _l2;
  sp[_off] = _sp;
}

__global__ void scjac2(const float *const a11, const float *const a22, const float *const a21, float *const cs, float *const sn, float *const l1, float *const l2, int *const sp)
{
  const size_t _off = (size_t)(blockIdx.x) * blockDim.x + threadIdx.x;
  const float _a11 = a11[_off];
  const float _a22 = a22[_off];
  const float _a21 = a21[_off];
  int _e1, _e2, _er;
  (void)frexpf(fmaxf(fabsf(_a11), FLT_TRUE_MIN), &_e1);
  _e1 = FLT_BIG_EXP - _e1;
  (void)frexpf(fmaxf(fabsf(_a22), FLT_TRUE_MIN), &_e2);
  _e2 = FLT_BIG_EXP - _e2;
  (void)frexpf(fmaxf(fabsf(_a21), FLT_TRUE_MIN), &_er);
  _er = FLT_BIG_EXP - _er;
  int _es = min(min(_e1, _e2), _er);
  const float _ar = scalbnf(_a21, _es);
  const float _a1 = scalbnf(_a11, _es);
  const float _a2 = scalbnf(_a22, _es);
  const float _aa = fabsf(_ar);
  const float _as = copysignf(1.0f, _ar);
  int _sp = (_es << 1);
  _es = -_es;
  const float _an = __fmul_rn(_aa, 2.0f);
  const float _ad = __fsub_rn(_a1, _a2);
  const float _t2 = copysignf(fminf(fmaxf(__fdiv_rn(_an, fabsf(_ad)), 0.0f), FLT_SQRT_HUGE), _ad);
  const float _t1 = __fdiv_rn(_t2, __fadd_rn(1.0f, __fsqrt_rn(__fmaf_rn(_t2, _t2, 1.0f))));
  const float _s2 = __fmaf_rn(_t1, _t1, 1.0f);
  const float _s1 = __fsqrt_rn(_s2);
  const float _cs = __frcp_rn(_s1);
  const float _sn = __fdiv_rn(__fmul_rn(_t1, _as), _s1);
  float _l1 = __fdiv_rn(__fmaf_rn(_t1, __fmaf_rn(_a2, _t1,  _an), _a1), _s2);
  float _l2 = __fdiv_rn(__fmaf_rn(_t1, __fmaf_rn(_a1, _t1, -_an), _a2), _s2);
  _sp |= (_l1 < _l2);
  _l1 = scalbnf(_l1, _es);
  _l2 = scalbnf(_l2, _es);
  cs[_off] = _cs;
  sn[_off] = _sn;
  l1[_off] = _l1;
  l2[_off] = _l2;
  sp[_off] = _sp;
}

__global__ void dcjac2(const double *const a11, const double *const a22, const double *const a21, double *const cs, double *const sn, double *const l1, double *const l2, int *const sp)
{
  const size_t _off = (size_t)(blockIdx.x) * blockDim.x + threadIdx.x;
  const double _a11 = a11[_off];
  const double _a22 = a22[_off];
  const double _a21 = a21[_off];
  int _e1, _e2, _er;
  (void)frexp(fmax(fabs(_a11), DBL_TRUE_MIN), &_e1);
  _e1 = DBL_BIG_EXP - _e1;
  (void)frexp(fmax(fabs(_a22), DBL_TRUE_MIN), &_e2);
  _e2 = DBL_BIG_EXP - _e2;
  (void)frexp(fmax(fabs(_a21), DBL_TRUE_MIN), &_er);
  _er = DBL_BIG_EXP - _er;
  int _es = min(min(_e1, _e2), _er);
  const double _ar = scalbn(_a21, _es);
  const double _a1 = scalbn(_a11, _es);
  const double _a2 = scalbn(_a22, _es);
  const double _aa = fabs(_ar);
  const double _as = copysign(1.0, _ar);
  int _sp = (_es << 1);
  _es = -_es;
  const double _an = __dmul_rn(_aa, 2.0);
  const double _ad = __dsub_rn(_a1, _a2);
  const double _t2 = copysign(fmin(fmax(__ddiv_rn(_an, fabs(_ad)), 0.0), DBL_SQRT_HUGE), _ad);
  const double _t1 = __ddiv_rn(_t2, __dadd_rn(1.0, __dsqrt_rn(__fma_rn(_t2, _t2, 1.0))));
  const double _s2 = __fma_rn(_t1, _t1, 1.0);
  const double _s1 = __dsqrt_rn(_s2);
  const double _cs = __drcp_rn(_s1);
  const double _sn = __ddiv_rn(__dmul_rn(_t1, _as), _s1);
  double _l1 = __ddiv_rn(__fma_rn(_t1, __fma_rn(_a2, _t1,  _an), _a1), _s2);
  double _l2 = __ddiv_rn(__fma_rn(_t1, __fma_rn(_a1, _t1, -_an), _a2), _s2);
  _sp |= (_l1 < _l2);
  _l1 = scalbn(_l1, _es);
  _l2 = scalbn(_l2, _es);
  cs[_off] = _cs;
  sn[_off] = _sn;
  l1[_off] = _l1;
  l2[_off] = _l2;
  sp[_off] = _sp;
}

#endif /* !DEVICE_CODE_H */
