#ifndef SBJAC2_H
#define SBJAC2_H

#include "vec.h"

extern fint sbjac2_(const fint n[static restrict 1], const float a11[static restrict VSL], const float a22[static restrict VSL], const float a21[static restrict VSL], float c[static restrict VSL], float at[static restrict VSL], float l1[static restrict VSL], float l2[static restrict VSL], unsigned p[static restrict 1]);

// for internal use only
extern fint sbjac2i(const fint n[static restrict 1], const float a11[static restrict VSL], const float a22[static restrict VSL], const float a21[static restrict VSL], float c[static restrict VSL], float at[static restrict VSL], float l1[static restrict VSL], float l2[static restrict VSL], unsigned p[static restrict 1]);

#endif /* !SBJAC2_H */
