#ifndef ZJAC2_H
#define ZJAC2_H

#include "vec.h"

extern int zjac2_(const fnat n[static restrict 1], const double a11[static restrict VDL], const double a22[static restrict VDL], const double a21r[static restrict VDL], const double a21i[static restrict VDL], double s[static restrict VDL], double t[static restrict VDL], double c[static restrict VDL], double ca[static restrict VDL], double sa[static restrict VDL], double l1[static restrict VDL], double l2[static restrict VDL]);

#endif /* !ZJAC2_H */
