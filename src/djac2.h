#ifndef DJAC2_H
#define DJAC2_H

#include "vec.h"

extern fint djac2_(const fnat n[static restrict 1], const double a11[static restrict VDL], const double a22[static restrict VDL], const double a21[static restrict VDL], double s[static restrict VDL], double t[static restrict VDL], double c[static restrict VDL], double l1[static restrict VDL], double l2[static restrict VDL], unsigned p[static restrict 1]);

#endif /* !DJAC2_H */
