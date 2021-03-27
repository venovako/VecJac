#ifndef DNORME_H
#define DNORME_H
#include "vec.h"
extern double dnorme_(const fnat m[static restrict 1], const double x[static restrict VDL], double e[static restrict 2], double f[static restrict 2]);
// for a complex x, properly laid out in memory as
//     [ Re(x) ]
// x = [ Im(x) ],
// set m = 2 * tilde(m) (see the paper)
#ifdef znorme_
#error znorme_ already defined
#else /* !znorme_ */
#define znorme_ dnorme_
#endif /* ?znorme_ */
#endif /* !DNORME_H */
