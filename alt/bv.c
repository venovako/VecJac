#include "bv.h"

uint64_t *bv_new(const BV_IX_T sz)
{
  return (uint64_t*)(sz ? calloc(((sz + 63u) >> 6u), sizeof(uint64_t)) : NULL);
}
