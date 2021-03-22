#include "proba.h"
#include "dnorme.h"

int main(int argc, char *argv[])
{
  static const fnat m = (2u * VDL);
  alignas(VA) double x[m] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0 };
  double e = 0.0, f = 0.0, d = dnorme_(&m, x, &e, &f);
  (void)printf("e=%#.17e f=%#.17e d=%#.17e\n", e, f, d);
  return EXIT_SUCCESS;
}
