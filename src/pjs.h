#ifndef PJS_H
#define PJS_H

#include "common.h"

#ifdef PJS_ME
#error PJS_ME already defined
#else /* !PJS_ME */
#define PJS_ME 2l
#endif /* ?PJS_ME */

#ifdef PJS_MM
#error PJS_MM already defined
#else /* !PJS_MM */
#define PJS_MM 4l
#endif /* ?PJS_MM */

extern unsigned *pjs(const long id, const unsigned n, unsigned stp[static restrict 1]);

#endif /* !PJS_H */
