#ifndef DVJSVD_H
#define DVJSVD_H

#include "vec.h"

// work: 8 * n elems
// iwork: n / VDL elems
extern int dvjsvd_(const fnat m[static restrict 1], const fnat n[static restrict 1], double G[static restrict VDL], const fnat ldG[static restrict 1], double V[static restrict VDL], const fnat ldV[static restrict 1], double eS[static restrict 1], double fS[static restrict 1], const unsigned js[static restrict 1], const unsigned stp[static restrict 1], unsigned swp[static restrict 1], double work[static restrict VDL], unsigned iwork[static restrict 1]);

#endif /* !DVJSVD_H */
