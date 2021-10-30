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

  DZNRME_VARS;
  DZNRME_LOOP(x,-3.0);
  VDKVSORT(re,rf);
  DZNRME_RED;
  DZNRME_RET;
}
