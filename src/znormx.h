#ifndef ZNORMX_H
#define ZNORMX_H

#include "vec.h"

// Returns an approximation for ||A||_max, i.e., max{|Ar|,|Ai|}.
// Multiplied by M_SQRT2, it would be the lowest upper bound for
// ||A||_max (for a general complex A), but then it could overflow.
// M_SQRT2 should be factored in any computation with the function's result.
extern double znormx_(const fnat m[static restrict 1], const fnat n[static restrict 1], const double Ar[static restrict VDL], const fnat ldAr[static restrict 1], const double Ai[static restrict VDL], const fnat ldAi[static restrict 1]);

#endif /* !ZNORMX_H */
