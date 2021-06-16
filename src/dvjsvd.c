#include "dvjsvd.h"

int dvjsvd_(const fnat m[static restrict 1], const fnat n[static restrict 1], double G[static restrict VDL], const fnat ldG[static restrict 1], double V[static restrict VDL], const fnat ldV[static restrict 1], double eS[static restrict 1], double fS[static restrict 1], const unsigned js[static restrict 1], const unsigned stp[static restrict 1], unsigned swp[static restrict 1])
{
  if (!*n)
    return 0;
  if (*m < *n)
    return -1;
  if (*ldG < *m)
    return -4;
  if (*ldG & VDL_1)
    return -4;
  if (*ldV < *n)
    return -6;
  if (*ldV & VDL_1)
    return -6;

  // TODO

  return 0;
}
