#ifndef BV_H
#define BV_H

#include "common.h"

#ifndef BV_IX_T
#define BV_IX_T uint32_t
#endif /* !BV_IX_T */

extern uint64_t *bv_new(const BV_IX_T sz);

static inline void bv_set(uint64_t bv[static restrict 1], const BV_IX_T ix)
{
  bv[ix >> 6u] |= (UINT64_C(0x01) << (ix & 0x3Fu));
}

static inline uint64_t bv_get(const uint64_t bv[static restrict 1], const BV_IX_T ix)
{
  return (bv[ix >> 6u] & (UINT64_C(0x01) << (ix & 0x3Fu)));
}

static inline void bv_clr(uint64_t bv[static restrict 1], const BV_IX_T ix)
{
  bv[ix >> 6u] &= ~(UINT64_C(0x01) << (ix & 0x3Fu));
}

#endif /* !BV_H */
