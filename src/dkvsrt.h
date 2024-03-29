#ifndef DKVSRT_H
#define DKVSRT_H

#include "defops.h"

// VDKVBITONIC and VDKVSORT have been inspired by the code of Berenger Bramas.
// Please, see https://gitlab.inria.fr/bramas/avx-512-sort and its LICENSE.

#ifdef VDKVBITONIC
#error VDKVBITONIC already defined
#else /* !VDKVBITONIC */
#define VDKVBITONIC(e,f,i,m)                                         \
  {                                                                  \
    register const VD pe = _mm512_permutexvar_pd(i, e); VDP(pe);     \
    register const VD pf = _mm512_permutexvar_pd(i, f); VDP(pf);     \
    register const MD mm = VDEFLE(e,pe,f,pf); MDP(mm);               \
    register const VD em = _mm512_mask_blend_pd(mm, pe, e); VDP(em); \
    register const VD ex = _mm512_mask_blend_pd(mm, e, pe); VDP(ex); \
    register const VD fm = _mm512_mask_blend_pd(mm, pf, f); VDP(fm); \
    register const VD fx = _mm512_mask_blend_pd(mm, f, pf); VDP(fx); \
    e = _mm512_mask_mov_pd(em, m, ex); VDP(e);                       \
    f = _mm512_mask_mov_pd(fm, m, fx); VDP(f);                       \
  }
#endif /* ?VDKVBITONIC */

#ifdef VDKVSORT
#error VDKVSORT already defined
#else /* !VDKVSORT */
#define VDKVSORT(e,f)                                         \
  {                                                           \
    register VI i = _mm512_set_epi64(6, 7, 4, 5, 2, 3, 0, 1); \
    VDKVBITONIC(e,f,i,0xAAu);                                 \
    i = _mm512_set_epi64(4, 5, 6, 7, 0, 1, 2, 3);             \
    VDKVBITONIC(e,f,i,0xCCu);                                 \
    i = _mm512_set_epi64(6, 7, 4, 5, 2, 3, 0, 1);             \
    VDKVBITONIC(e,f,i,0xAAu);                                 \
    i = _mm512_set_epi64(0, 1, 2, 3, 4, 5, 6, 7);             \
    VDKVBITONIC(e,f,i,0xF0u);                                 \
    i = _mm512_set_epi64(5, 4, 7, 6, 1, 0, 3, 2);             \
    VDKVBITONIC(e,f,i,0xCCu);                                 \
    i = _mm512_set_epi64(6, 7, 4, 5, 2, 3, 0, 1);             \
    VDKVBITONIC(e,f,i,0xAAu);                                 \
  }
#endif /* ?VDKVSORT */

#endif /* !DKVSRT_H */
