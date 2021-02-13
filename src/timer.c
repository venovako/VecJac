#include "timer.h"

uint64_t tsc_get_freq_hz(unsigned rem_den[static 2])
{
  const uint64_t hz = (uint64_t)atoz(getenv("TSC_FREQ_HZ"));
  if (hz) {
    rem_den[0u] = 0u;
    rem_den[1u] = 1u;
    return hz;
  }
  int abcd[4] = { 0, 0, 0, 0 };
  __cpuid(abcd, 0x15);
  if (abcd[0u]) {
    const uint64_t num = (uint64_t)(abcd[1u]) * (uint64_t)(abcd[2u]);
    const uint64_t den = (uint64_t)(abcd[0u]);
    rem_den[0u] = (unsigned)(num % den);
    rem_den[1u] = (unsigned)den;
    return (num / den);
  }
  rem_den[0u] = 0u;
  rem_den[1u] = 0u;
  return UINT64_C(0);
}

int64_t get_thread_ns()
{
  struct timespec t;
  return (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t) ? INT64_C(-1) : t2ns(&t));
}

int64_t get_sys_us()
{
  struct timeval t;
  return (gettimeofday(&t, NULL) ? INT64_C(-1) : t2us(&t));
}
