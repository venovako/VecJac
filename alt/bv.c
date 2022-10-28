#include "bv.h"

uint32_t *bv32_new(const fnat sz)
{
  return (uint32_t*)(sz ? calloc(((sz + 31u) >> 5u), sizeof(uint32_t)) : NULL);
}

uint32_t *bv32_reset(uint32_t bv[static restrict 1], const fnat sz, const int c)
{
  return (uint32_t*)memset(bv, c, ((sz + 31u) >> 5u) * sizeof(uint32_t));
}

uint64_t *bv64_new(const fnat sz)
{
  return (uint64_t*)(sz ? calloc(((sz + 63u) >> 6u), sizeof(uint64_t)) : NULL);
}

uint64_t *bv64_reset(uint64_t bv[static restrict 1], const fnat sz, const int c)
{
  return (uint64_t*)memset(bv, c, ((sz + 63u) >> 6u) * sizeof(uint64_t));
}
