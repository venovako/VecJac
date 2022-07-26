#include "psort.h"

#include <oneapi/tbb/parallel_sort.h>

void dpsort(const unsigned long long n, double *const x) throw()
{
  oneapi::tbb::parallel_sort(x, (x + n));
}
