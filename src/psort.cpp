#include "psort.h"

#include <tbb/parallel_sort.h>

void dpsort(const unsigned long long n, double *const x) throw()
{
  tbb::parallel_sort(x, (x + n));
}
