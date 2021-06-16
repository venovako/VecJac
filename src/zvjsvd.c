#include "zvjsvd.h"

int zvjsvd_(const fnat m[static restrict 1], const fnat n[static restrict 1], double Gr[static restrict VDL], const fnat ldGr[static restrict 1], double Gi[static restrict VDL], const fnat ldGi[static restrict 1], double Vr[static restrict VDL], const fnat ldVr[static restrict 1], double Vi[static restrict VDL], const fnat ldVi[static restrict 1], double eS[static restrict 1], double fS[static restrict 1], const unsigned js[static restrict 1], const unsigned stp[static restrict 1], unsigned swp[static restrict 1], double work[static restrict VDL], unsigned iwork[static restrict 1])
{
  if (!*n)
    return 0;
  if (*m < *n)
    return -1;
  if (*n & VDL_1)
    return -2;
  if (*ldGr < *m)
    return -4;
  if (*ldGr & VDL_1)
    return -4;
  if (*ldGi < *m)
    return -6;
  if (*ldGi & VDL_1)
    return -6;
  if (*ldVr < *n)
    return -8;
  if (*ldVr & VDL_1)
    return -8;
  if (*ldVi < *n)
    return -10;
  if (*ldVi & VDL_1)
    return -10;

  double *const a11 = work;
  double *const a22 = a11 + *n;
  double *const a21r = a22 + *n;
  double *const a21i = a21r + *n;
  double *const s = a21i + *n;
  double *const t = s + *n;
  double *const c = t + *n;
  double *const ca = c + *n;
  double *const sa = ca + *n;
  double *const l1 = sa + *n;
  double *const l2 = l1 + *n;
  unsigned *const p = iwork;

  // TODO

  return 0;
}
