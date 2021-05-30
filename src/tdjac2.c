#include "rnd.h"
#include "djac2.h"
#include "wnrme.h"

int main(int argc, char *argv[])
{
  alignas(VA) double a11[VDL], a22[VDL], a21[VDL], s[VDL], t[VDL], c[VDL], l1[VDL], l2[VDL];
  fnat n = VDL;
  if (1 == argc) {
    for (fnat i = 0u; i < n; ++i) {
      a11[i] = dfrand(DBL_MAX);
      a22[i] = dfrand(DBL_MAX);
      a21[i] = dfrand(DBL_MAX);
      s[i] = -0.0;
      t[i] = -0.0;
      c[i] = -0.0;
      l1[i] = -0.0;
      l2[i] = -0.0;
    }
  }
  else if (4 == argc) {
    *a11 = atof(argv[1u]);
    *a22 = atof(argv[2u]);
    *a21 = atof(argv[3u]);
    for (fnat i = 1u; i < n; ++i) {
      a11[i] = *a11;
      a22[i] = *a22;
      a21[i] = *a21;
    }
  }
  else {
    (void)fprintf(stderr, "%s [a11 a22 a21]\n", *argv);
    return EXIT_FAILURE;
  }
  (void)printf("djac2=%d\n", djac2_(&n, a11, a22, a21, s, t, c, l1, l2));
  (void)printf("%25s,%25s,%25s,%25s,%25s,%30s,%30s\n", "\"S\"", "\"T\"", "\"C\"", "\"L1\"", "\"L2\"", "\"AE\"", "\"RE\"");
  char a[31] = { '\0' };
  wide ae, re;
  for (fnat i = 0u; i < n; ++i) {
    (void)printf("%25s,", dtoa(a, s[i]));
    (void)printf("%25s,", dtoa(a, t[i]));
    (void)printf("%25s,", dtoa(a, c[i]));
    (void)printf("%25s,", dtoa(a, l1[i]));
    (void)printf("%25s,", dtoa(a, l2[i]));
    re = wrer(a11[i], a22[i], a21[i], s[i], t[i], c[i], l1[i], l2[i], &ae);
    (void)printf("%30s,", xtoa(a, (long double)ae));
    (void)printf("%30s\n", xtoa(a, (long double)re));
  }
  return EXIT_SUCCESS;
}