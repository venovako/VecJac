#include <mathimf.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#define N 16
static double A[N * N];

int main(void)
{
  for (int j = 0; j < N/2; ++j)
    for (int i = 0; i < N; ++i)
      A[j * N + i] = scalbn(1.0, DBL_MAX_EXP - 1 - j - abs(i - j));
  for (int j = N/2; j < N; ++j)
    for (int i = 0; i < N; ++i)
      A[j * N + i] = scalbn(1.0, -51 - j - abs(i - j));
  FILE *const f = fopen("dgsex.G", "w");
  if (!f)
    return EXIT_FAILURE;
  (void)fwrite(A, sizeof(double), N * N, f);
  (void)fclose(f);
  return EXIT_SUCCESS;
}
