#ifndef TIMER_H
#define TIMER_H

#include "common.h"

static inline uint64_t rdtsc_beg(unsigned aux[static 1])
{
  _mm_mfence();
  return __rdtscp(aux);
}

static inline uint64_t rdtsc_end(unsigned aux[static 1])
{
  const uint64_t tsc = __rdtscp(aux);
  _mm_lfence();
  return tsc;
}

static inline uint64_t tsc_get_freq_hz(unsigned rem_den[static 2])
{
  rem_den[0u] = 0u;
  rem_den[1u] = 0u;
#ifdef TSC_FREQ_HZ
#if (TSC_FREQ_HZ == 0ull)
  unsigned eax = 0u, ebx = 0u, ecx = 0u, edx = 0u;
  __cpuid(0x15u, eax, ebx, ecx, edx);
  if (eax) {
    const uint64_t num = (uint64_t)ebx * ecx;
    const uint64_t den = (uint64_t)eax;
    rem_den[0u] = (unsigned)(num % den);
    rem_den[1u] = eax;
    return (num / den);
  }
  return UINT64_C(0);
#else /* TSC_FREQ_HZ > 0 */
  rem_den[1u] = 1u;
  return (TSC_FREQ_HZ);
#endif /* ?TSC_FREQ_HZ */
#else /* !TSC_FREQ_HZ */
  return UINT64_C(0);
#endif /* TSC_FREQ_HZ */
}

static inline long double tsc_lap(const uint64_t freq_hz, const uint64_t beg, const uint64_t end)
{
  if (freq_hz) {
    if (end >= beg) {
      const uint64_t lap = end - beg;
      return ((lap >= freq_hz) ? (lap / (long double)freq_hz) : ((long double)lap / freq_hz));
    }
    return -0.0L;
  }
  return (long double)INFINITY;
}

#endif /* !TIMER_H */
