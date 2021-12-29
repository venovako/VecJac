#ifndef CJAC2F_H
#define CJAC2F_H

#include "vec.h"

extern fint cjac2f_(const fnat n[static restrict 1], const float a11[static restrict VSL], const float a22[static restrict VSL], const float a21r[static restrict VSL], const float a21i[static restrict VSL], float t[static restrict VSL], float c[static restrict VSL], float ca[static restrict VSL], float sa[static restrict VSL], float l1[static restrict VSL], float l2[static restrict VSL], unsigned p[static restrict 1]);

#endif /* !CJAC2F_H */
