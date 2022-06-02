// generates a list of logical CPUs for `taskset -c` on 2x Intel Xeon Platinum 8358 CPU
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
  if ((n <= 0) || (n > 128))
    return EXIT_FAILURE;
  const int t = atoi(argv[2]);
  if ((t < 0) || (t >= n))
    return EXIT_FAILURE;
  int i = 0, j = 0;
  switch (n) {
  case 1:
    (void)printf("0-127\n");
    break;
  case 2:
    for (i = 32 * t; j < 2; i += 64) {
      if (j++)
        (void)printf(",");
      (void)printf("%d-%d", i, (i + 31));
    }
    (void)printf("\n");
    break;
  case 4:
    for (i = 16 * t; j < 2; i += 64) {
      if (j++)
        (void)printf(",");
      (void)printf("%d-%d", i, (i + 15));
    }
    (void)printf("\n");
    break;
  case 8:
    for (i = 8 * t; j < 2; i += 64) {
      if (j++)
        (void)printf(",");
      (void)printf("%d-%d", i, (i + 7));
    }
    (void)printf("\n");
    break;
  case 16:
    for (i = 4 * t; j < 2; i += 64) {
      if (j++)
        (void)printf(",");
      (void)printf("%d-%d", i, (i + 3));
    }
    (void)printf("\n");
    break;
  case 32:
    for (i = 2 * t; j < 2; i += 64) {
      if (j++)
        (void)printf(",");
      (void)printf("%d,%d", i, (i + 1));
    }
    (void)printf("\n");
    break;
  case 64:
    for (i = t; j < 2; i += 64) {
      if (j++)
        (void)printf(",");
      (void)printf("%d", i);
    }
    (void)printf("\n");
    break;
  case 128:
    (void)printf("%d\n", t);
    break;
  default:
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
