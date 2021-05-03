#ifndef VECDEF_H
#define VECDEF_H

// assumes m0 = _mm512_set1_ps(-0.0f)
#ifdef VSABS
#error VSABS already defined
#else /* !VSABS */
#define VSABS(x) VSANDNOT(m0,(x))
#endif /* ?VSABS */
#ifdef VSSGN
#error VSSGN already defined
#else /* !VSSGN */
#define VSSGN(x) VSAND((x),m0)
#endif /* ?VSSGN */

// assumes m0 = _mm512_set1_pd(-0.0)
#ifdef VDABS
#error VDABS already defined
#else /* !VDABS */
#define VDABS(x) VDANDNOT(m0,(x))
#endif /* ?VDABS */
#ifdef VDSGN
#error VDSGN already defined
#else /* !VDSGN */
#define VDSGN(x) VDAND((x),m0)
#endif /* ?VDSGN */

#endif /* !VECDEF_H */
