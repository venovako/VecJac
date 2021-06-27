#ifndef ZVJSVD_H
#define ZVJSVD_H

#include "vec.h"

// work: 11 * (n/2) elems
// iwork: 2 * (n/2) / VDL elems
extern fint zvjsvd_(const fnat m[static restrict 1], const fnat n[static restrict 1], double Gr[static restrict VDL], const fnat ldGr[static restrict 1], double Gi[static restrict VDL], const fnat ldGi[static restrict 1], double Vr[static restrict VDL], const fnat ldVr[static restrict 1], double Vi[static restrict VDL], const fnat ldVi[static restrict 1], double eS[static restrict 1], double fS[static restrict 1], const unsigned js[static restrict 1], const unsigned stp[static restrict 1], const unsigned swp[static restrict 1], double work[static restrict VDL], unsigned iwork[static restrict 1]);

#endif /* !ZVJSVD_H */
