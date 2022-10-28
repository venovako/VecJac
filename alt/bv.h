#ifndef BV_H
#define BV_H

#include "common.h"

extern uint32_t *bv32_new(const fnat sz);
extern uint32_t *bv32_reset(uint32_t bv[static restrict 1], const fnat sz, const int c);

static inline void bv32_set(uint32_t bv[static restrict 1], const fnat ix)
{
  bv[ix >> 5u] |= (UINT32_C(0x01) << (ix & 0x1Fu));
}

static inline uint32_t bv32_get(const uint32_t bv[static restrict 1], const fnat ix)
{
  return (bv[ix >> 5u] & (UINT32_C(0x01) << (ix & 0x1Fu)));
}

static inline void bv32_clr(uint32_t bv[static restrict 1], const fnat ix)
{
  bv[ix >> 5u] &= ~(UINT32_C(0x01) << (ix & 0x1Fu));
}

extern uint64_t *bv64_new(const fnat sz);
extern uint64_t *bv64_reset(uint64_t bv[static restrict 1], const fnat sz, const int c);

static inline void bv64_set(uint64_t bv[static restrict 1], const fnat ix)
{
  bv[ix >> 6u] |= (UINT64_C(0x01) << (ix & 0x3Fu));
}

static inline uint64_t bv64_get(const uint64_t bv[static restrict 1], const fnat ix)
{
  return (bv[ix >> 6u] & (UINT64_C(0x01) << (ix & 0x3Fu)));
}

static inline void bv64_clr(uint64_t bv[static restrict 1], const fnat ix)
{
  bv[ix >> 6u] &= ~(UINT64_C(0x01) << (ix & 0x3Fu));
}

#endif /* !BV64_H */
