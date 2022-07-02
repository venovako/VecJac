#ifndef TIMER_H
#define TIMER_H

#include "common.h"

static inline long double tsc_lap(const uint64_t freq_hz, const uint64_t beg, const uint64_t end)
{
  if (freq_hz) {
    if (end >= beg) {
      const uint64_t lap = end - beg;
      return ((lap >= freq_hz) ? (lap / (long double)freq_hz) : ((long double)lap / freq_hz));
    }
    return -0.0L;
  }
  return HUGE_VALL;
}

static inline uint64_t rdtsc_beg(int *const aux)
{
#ifdef _LARGEFILE64_SOURCE
  struct timespec t;
  return (uint64_t)((*(int*)aux = clock_gettime(CLOCK_MONOTONIC_RAW, &t)) ? INT64_C(-1) : (int64_t)(t.tv_sec * INT64_C(1000000000) + t.tv_nsec));
#else /* !_LARGEFILE64_SOURCE */
  LARGE_INTEGER t;
  return (uint64_t)((*(BOOL*)aux = QueryPerformanceCounter(&t)) ? (int64_t)(t.QuadPart) : INT64_C(-1));
#endif /* ?_LARGEFILE64_SOURCE */
}

static inline uint64_t rdtsc_end(int *const aux)
{
#ifdef _LARGEFILE64_SOURCE
  struct timespec t;
  return (uint64_t)((*(int*)aux = clock_gettime(CLOCK_MONOTONIC_RAW, &t)) ? INT64_C(-1) : (int64_t)(t.tv_sec * INT64_C(1000000000) + t.tv_nsec));
#else /* !_LARGEFILE64_SOURCE */
  LARGE_INTEGER t;
  return (uint64_t)((*(BOOL*)aux = QueryPerformanceCounter(&t)) ? (int64_t)(t.QuadPart) : INT64_C(-1));  
#endif /* ?_LARGEFILE64_SOURCE */
}

static inline uint64_t tsc_get_freq_hz_(int *const aux)
{
#ifdef _LARGEFILE64_SOURCE
  *aux = 0;
  return UINT64_C(1000000000);
#else /* !_LARGEFILE64_SOURCE */
  LARGE_INTEGER t;
  return ((*(BOOL*)aux = QueryPerformanceFrequency(&t)) ? (uint64_t)(t.QuadPart) : UINT64_C(0));
#endif /* ?_LARGEFILE64_SOURCE */
}

#endif /* !TIMER_H */
