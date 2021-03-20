#include "proba.h"

static double dnormE(double r[static restrict 2], const double x[static restrict 1], const size_t m)
{
  ASSERT_VFPENV;
  r[0u] = -INFINITY;
  r[1u] = 1.0;
  if (m & (size_t)VDL_1)
    return -0.0;
  if (!m)
    goto end;

 end:
  return scalb(r[1u], r[0u]);
}

int main(int argc, char *argv[])
{
  alignas(VA) double r[3u];
  alignas(VA) double x[VDL];
  r[2u] = dnormE(r, x, VDL);
  (void)printf("%#.17e\n", r[2u]);
  return EXIT_SUCCESS;
}
