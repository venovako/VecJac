#ifndef D8SORT_H
#define D8SORT_H

// VDBITONIC and VDSORT have been refactored from the code of Berenger Bramas.
// Please, see https://gitlab.inria.fr/bramas/avx-512-sort and its LICENSE.

#ifdef VDBITONIC
#error VDBITONIC already defined
#else /* !VDBITONIC */
#define VDBITONIC(x,i,m)                                         \
  {                                                              \
    register const VD px = _mm512_permutexvar_pd(i, x); VDP(px); \
    register const VD xm = _mm512_min_pd(px, x); VDP(xm);        \
    register const VD xx = _mm512_max_pd(px, x); VDP(xx);        \
    x = _mm512_mask_mov_pd(xm, m, xx); VDP(x);                   \
  }
#endif /* ?VDBITONIC */

#ifdef VDSORT
#error VDSORT already defined
#else /* !VDSORT */
#define VDSORT(x)                                             \
  {                                                           \
    register VI i = _mm512_set_epi64(6, 7, 4, 5, 2, 3, 0, 1); \
    VDBITONIC(x,i,0xAAu);                                     \
    i = _mm512_set_epi64(4, 5, 6, 7, 0, 1, 2, 3);             \
    VDBITONIC(x,i,0xCCu);                                     \
    i = _mm512_set_epi64(6, 7, 4, 5, 2, 3, 0, 1);             \
    VDBITONIC(x,i,0xAAu);                                     \
    i = _mm512_set_epi64(0, 1, 2, 3, 4, 5, 6, 7);             \
    VDBITONIC(x,i,0xF0u);                                     \
    i = _mm512_set_epi64(5, 4, 7, 6, 1, 0, 3, 2);             \
    VDBITONIC(x,i,0xCCu);                                     \
    i = _mm512_set_epi64(6, 7, 4, 5, 2, 3, 0, 1);             \
    VDBITONIC(x,i,0xAAu);                                     \
  }
#endif /* ?VDSORT */

#endif /* !D8SORT_H */
