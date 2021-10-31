#include "znorme.h"

#include "dznrme.h"

double znorme_(const fnat m[static restrict 1], const double zr[static restrict VDL], const double zi[static restrict VDL], double e0[static restrict 1], double f0[static restrict 1], double e1[static restrict 1], double f1[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -5.0;
  if (*m & VDL_1)
    return -1.0;
  if (IS_NOT_ALIGNED(zr))
    return -2.0;
  if (IS_NOT_ALIGNED(zi))
    return -3.0;
#endif /* !NDEBUG */

  if (!*m) {
    *e1 = *e0 = -HUGE_VAL;
    *f1 = *f0 = 1.0;
    return 0.0;
  }

  _mm_prefetch((const char*)zr, _MM_HINT_T0);
  _mm_prefetch((const char*)zi, _MM_HINT_T0);
  DZNRME_VARS;

  const fnat m_ = *m - VDL;
  fnat j = 0u, i;
  for (i = j; j < m_; i = j) {
    j += VDL;
    _mm_prefetch((const char*)(zr + j), _MM_HINT_T0);
    _mm_prefetch((const char*)(zi + j), _MM_HINT_T0);
    DZNRME_LOOP(zr,-4.0);
    DZNRME_LOOP(zi,-4.0);
  }
  DZNRME_LOOP(zr,-4.0);
  DZNRME_LOOP(zi,-4.0);

  VDKVSORT(re,rf);
  DZNRME_RED;
  DZNRME_RET;
}
