#include "znorme.h"

#include "dznrme.h"
#include "dkvsrt.h"

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

  DZNRME_VARS;
  DZNRME_LOOP(zr,-4.0);
  DZNRME_LOOP(zi,-4.0);
  VDKVSORT(re,rf);
  VDEFRED(re,rf);
  DZNRME_RET;
}
