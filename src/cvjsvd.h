#ifndef CVJSVD_H
#define CVJSVD_H

#include "vec.h"

// work: 7 * n elems
// iwork: n / VSL elems
// set iwork[0] to 0 if V is to be set to I, or to something else if V is preset
extern fint cvjsvd_(const fnat m[static restrict 1], const fnat n[static restrict 1], float Gr[static restrict VSL], const fnat ldGr[static restrict 1], float Gi[static restrict VSL], const fnat ldGi[static restrict 1], float Vr[static restrict VSL], const fnat ldVr[static restrict 1], float Vi[static restrict VSL], const fnat ldVi[static restrict 1], float eS[static restrict 1], float fS[static restrict 1], const unsigned js[static restrict 1], const unsigned stp[static restrict 1], const unsigned swp[static restrict 1], float work[static restrict VSL], unsigned iwork[static restrict 1]);

#endif /* !CVJSVD_H */
