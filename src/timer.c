#include "timer.h"

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
