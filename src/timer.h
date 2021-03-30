#ifndef TIMER_H
#define TIMER_H

#include "common.h"

static inline uint64_t rdtsc_beg(unsigned aux[static restrict 1])
{
  _mm_mfence();
  return __rdtscp(aux);
}

static inline uint64_t rdtsc_end(unsigned aux[static restrict 1])
{
  const uint64_t tsc = __rdtscp(aux);
  _mm_lfence();
  return tsc;
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
  return HUGE_VALL;
}

static inline int64_t t2ns(const struct timespec tp[static restrict 1])
{
  return (int64_t)(tp->tv_sec * INT64_C(1000000000) + tp->tv_nsec);
}

static inline int64_t t2us(const struct timeval tp[static restrict 1])
{
  return (int64_t)(tp->tv_sec * INT64_C(1000000) + tp->tv_usec);
}

extern uint64_t tsc_get_freq_hz_(unsigned rem_den[static restrict 2]);

extern int64_t get_thread_ns_();
extern int64_t get_sys_us_();

#endif /* !TIMER_H */
