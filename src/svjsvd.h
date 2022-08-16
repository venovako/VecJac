#ifndef SVJSVD_H
#define SVJSVD_H

#include "vec.h"

// work: 4 * n elems
// iwork: n / VSL elems
extern fint svjsvd_(const fnat m[static restrict 1], const fnat n[static restrict 1], float G[static restrict VSL], const fnat ldG[static restrict 1], float V[static restrict VSL], const fnat ldV[static restrict 1], float eS[static restrict 1], float fS[static restrict 1], const unsigned js[static restrict 1], const unsigned stp[static restrict 1], const unsigned swp[static restrict 1], float work[static restrict VSL], unsigned iwork[static restrict 1]);

#endif /* !SVJSVD_H */
