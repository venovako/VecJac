#ifndef DBJAC2_H
#define DBJAC2_H

#include "vec.h"

extern fint dbjac2_(const fint n[static restrict 1], const double a11[static restrict VDL], const double a22[static restrict VDL], const double a21[static restrict VDL], double c[static restrict VDL], double at[static restrict VDL], double l1[static restrict VDL], double l2[static restrict VDL], unsigned p[static restrict 1]);

// for internal use only
extern fint dbjac2i(const fint n[static restrict 1], const double a11[static restrict VDL], const double a22[static restrict VDL], const double a21[static restrict VDL], double c[static restrict VDL], double at[static restrict VDL], double l1[static restrict VDL], double l2[static restrict VDL], unsigned p[static restrict 1]);

#endif /* !DBJAC2_H */
