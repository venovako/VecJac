#include "dnorme.h"

#include "dznrme.h"

double dnorme_(const fnat m[static restrict 1], const double x[static restrict VDL], double e0[static restrict 1], double f0[static restrict 1], double e1[static restrict 1], double f1[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -4.0;
  if (*m & VDL_1)
    return -1.0;
  if (IS_NOT_ALIGNED(x))
    return -2.0;
#endif /* !NDEBUG */

  if (!*m) {
    *e1 = *e0 = -HUGE_VAL;
    *f1 = *f0 = 1.0;
    return -0.0;
  }

  _mm_prefetch((const char*)x, _MM_HINT_T0);
  DZNRME_VARS;

  const fnat m_ = *m - VDL;
  fnat j = 0u, i;
  for (i = j; j < m_; i = j) {
    _mm_prefetch((const char*)(x + (j += VDL)), _MM_HINT_T0);
    DZNRME_LOOP(x,-3.0);
  }
  DZNRME_LOOP(x,-3.0);

  VDKVSORT(re,rf);
  DZNRME_RED;
  DZNRME_RET;
}
