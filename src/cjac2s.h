#ifndef CJAC2S_H
#define CJAC2S_H

#include "vec.h"

extern fint cjac2s_(const fnat n[static restrict 1], const float a11[static restrict VSL], const float a22[static restrict VSL], const float a21r[static restrict VSL], const float a21i[static restrict VSL], float cs[static restrict VSL], float ca[static restrict VSL], float sa[static restrict VSL], float l1[static restrict VSL], float l2[static restrict VSL], unsigned p[static restrict 1]);

#endif /* !CJAC2S_H */
