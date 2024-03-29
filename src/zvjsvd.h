#ifndef ZVJSVD_H
#define ZVJSVD_H

#include "vec.h"

// work: 7 * n elems
// if JTRACE is defined, work should start with a nul-terminated name of the trace file to be written
// iwork: n / VDL elems
// set iwork[0] to 0 if V is to be set to I, or to something else if V is preset
extern fint zvjsvd_(const fnat m[static restrict 1], const fnat n[static restrict 1], double Gr[static restrict VDL], const fnat ldGr[static restrict 1], double Gi[static restrict VDL], const fnat ldGi[static restrict 1], double Vr[static restrict VDL], const fnat ldVr[static restrict 1], double Vi[static restrict VDL], const fnat ldVi[static restrict 1], double eS[static restrict 1], double fS[static restrict 1], const unsigned js[static restrict 1], const unsigned stp[static restrict 1], const unsigned swp[static restrict 1], double work[static restrict VDL], unsigned iwork[static restrict 1]);

#endif /* !ZVJSVD_H */
