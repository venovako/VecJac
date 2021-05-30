#include "tzjac2.h"

#include "rnd.h"
#include "zjac2.h"
#include "wnrme.h"

int main(int argc, char *argv[])
{
  alignas(VA) double a11[VDL], a22[VDL], a21r[VDL], a21i[VDL], s[VDL], t[VDL], c[VDL], ca[VDL], sa[VDL], l1[VDL], l2[VDL];
  fnat n = VDL;
  if (1 == argc) {
    for (fnat i = 0u; i < n; ++i) {
      a11[i] = dfrand(DBL_MAX);
      a22[i] = dfrand(DBL_MAX);
      a21r[i] = dfrand(DBL_MAX);
      a21i[i] = dfrand(DBL_MAX);
      s[i] = -0.0;
      t[i] = -0.0;
      c[i] = -0.0;
      ca[i] = -0.0;
      sa[i] = -0.0;
      l1[i] = -0.0;
      l2[i] = -0.0;
    }
  }
  else if (5 == argc) {
    *a11 = atof(argv[1u]);
    *a22 = atof(argv[2u]);
    *a21r = atof(argv[3u]);
    *a21i = atof(argv[4u]);
    for (fnat i = 1u; i < n; ++i) {
      a11[i] = *a11;
      a22[i] = *a22;
      a21r[i] = *a21r;
      a21i[i] = *a21i;
    }
  }
  else {
    (void)fprintf(stderr, "%s [a11 a22 a21r a21i]\n", *argv);
    return EXIT_FAILURE;
  }
  (void)printf("zjac2=%d\n", zjac2_(&n, a11, a22, a21r, a21i, s, t, c, ca, sa, l1, l2));
  (void)printf("%25s,%25s,%25s,%25s,%25s,%25s,%25s,%30s,%30s\n", "\"S\"", "\"T\"", "\"C\"", "\"Ca\"", "\"Sa\"", "\"L1\"", "\"L2\"", "\"AE\"", "\"RE\"");
  char a[31] = { '\0' };
  wide ae, re;
  for (fnat i = 0u; i < n; ++i) {
    (void)printf("%25s,", dtoa(a, s[i]));
    (void)printf("%25s,", dtoa(a, t[i]));
    (void)printf("%25s,", dtoa(a, c[i]));
    (void)printf("%25s,", dtoa(a, ca[i]));
    (void)printf("%25s,", dtoa(a, sa[i]));
    (void)printf("%25s,", dtoa(a, l1[i]));
    (void)printf("%25s,", dtoa(a, l2[i]));
    re = wrec(a11[i], a22[i], a21r[i], a21i[i], s[i], t[i], c[i], ca[i], sa[i], l1[i], l2[i], &ae);
    (void)printf("%30s,", xtoa(a, (long double)ae));
    (void)printf("%30s\n", xtoa(a, (long double)re));
  }
  return EXIT_SUCCESS;
}
