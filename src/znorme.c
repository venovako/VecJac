#include "znorme.h"

#include "dkvsrt.h"
#include "dznrme.h"

double znorme_(const fnat m[static restrict 1], const double zr[static restrict VDL], const double zi[static restrict VDL], double e[static restrict 2], double f[static restrict 2])
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