#ifndef CNORMX_H
#define CNORMX_H

#include "vec.h"

// Returns an approximation for ||A||_max, i.e., max{|Ar|,|Ai|}.
// Multiplied by M_SQRT2, it would be the lowest upper bound for
// ||A||_max (for a general complex A), but then it could overflow.
// M_SQRT2 should be factored in any computation with the function's result.
extern float cnormx_(const fnat m[static restrict 1], const fnat n[static restrict 1], const float Ar[static restrict VSL], const fnat ldAr[static restrict 1], const float Ai[static restrict VSL], const fnat ldAi[static restrict 1]);

#endif /* !CNORMX_H */
