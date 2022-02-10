// generates a list of logical CPUs for `taskset -c` on Intel Xeon Phi 7210
// n is the number of tasks, and 0 <= t < n is the task "ID"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  if (argc != 3) {
    (void)fprintf(stderr, "%s n t\n", *argv);
    return EXIT_FAILURE;
  }
  const int n = atoi(argv[1]);
  if ((n <= 0) || (n > 256))
    return EXIT_FAILURE;
  const int t = atoi(argv[2]);
  if ((t < 0) || (t >= n))
    return EXIT_FAILURE;
  int i = 0, j = 0;
  switch (n) {
  case 1:
    (void)printf("0-255\n");
    break;
  case 2:
    for (i = 2 * t; i < 256; i += 4) {
      if (j++)
        (void)printf(",");
      (void)printf("%d,%d", i, (i + 1));
    }
    (void)printf("\n");
    break;
  case 4:
    for (i = 2 * t; i < 256; i += 8) {
      if (j++)
        (void)printf(",");
      (void)printf("%d,%d", i, (i + 1));
    }
    (void)printf("\n");
    break;
  case 8:
    for (i = 2 * t; i < 256; i += 16) {
      if (j++)
        (void)printf(",");
      (void)printf("%d,%d", i, (i + 1));
    }
    (void)printf("\n");
    break;
  case 16:
    for (i = 2 * t; i < 256; i += 32) {
      if (j++)
        (void)printf(",");
      (void)printf("%d,%d", i, (i + 1));
    }
    (void)printf("\n");
    break;
  case 32:
    for (i = 2 * t; i < 256; i += 64) {
      if (j++)
        (void)printf(",");
      (void)printf("%d,%d", i, (i + 1));
    }
    (void)printf("\n");
    break;
  case 64:
    (void)printf("%d,%d,%d,%d\n", t, (t + 64), (t + 128), (t + 192));
    break;
  case 128:
    (void)printf("%d,%d\n", t, (t + 128));
    break;
  case 256:
    (void)printf("%d\n", t);
    break;
  default:
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
