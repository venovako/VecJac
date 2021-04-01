#include "proba.h"

#include "dnorme.h"
#include "znorme.h"
#include "dnormx.h"
#include "znormx.h"
#include "ddpscl.h"
#include "zdpscl.h"
#include "zmerge.h"
#include "zsplit.h"
#include "aalloc.h"
#include "dscale.h"
#include "zscale.h"
#include "mtxio.h"

int main(int argc, char *argv[])
{
  static const fnat m = (2u * VDL);
  alignas(VA) double x[m] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0 };
  double e0 = 0.0, f0 = 0.0, e1 = 0.0, f1 = 0.0, d = dnorme_(&m, x, &e0, &f0, &e1, &f1);
  (void)printf("e=(%#.17e,%#.17e) f=(%#.17e,%#.17e) d=%#.17e\n", e0, e1, f0, f1, d);

  /* static const fnat m = VDL, n = VDL, ldA = VDL; */
  /* alignas(VA) double A[VDL][VDL] = { */
  /*   { 0.0, 1.0, 2.0,  3.0,  4.0,  5.0,  6.0,  7.0 }, */
  /*   { 1.0, 2.0, 3.0,  4.0,  5.0,  6.0,  7.0,  8.0 }, */
  /*   { 2.0, 3.0, 4.0,  5.0,  6.0,  7.0,  8.0,  9.0 }, */
  /*   { 3.0, 4.0, 5.0,  6.0,  7.0,  8.0,  9.0, 10.0 }, */
  /*   { 4.0, 5.0, 6.0,  7.0,  8.0,  9.0, 10.0, 11.0 }, */
  /*   { 5.0, 6.0, 7.0,  8.0,  9.0, 10.0, 11.0, 12.0 }, */
  /*   { 6.0, 7.0, 8.0,  9.0, 10.0, 11.0, 12.0, 13.0 }, */
  /*   { 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0 } */
  /* }; */
  /* (void)printf("%#.17e\n", dnormx_(&m, &n, &(A[0u][0u]), &ldA)); */

  /* alignas(VA) double complex A[VDL] = { */
  /*   CMPLX(  1.0, -2.0), */
  /*   CMPLX( -3.0,  4.0), */
  /*   CMPLX(  5.0, -6.0), */
  /*   CMPLX( -7.0,  8.0), */
  /*   CMPLX(  9.0,-10.0), */
  /*   CMPLX(-11.0, 12.0), */
  /*   CMPLX( 13.0,-14.0), */
  /*   CMPLX(-15.0, 16.0) */
  /* }; */
  /* alignas(VA) double Ar[VDL]; */
  /* alignas(VA) double Ai[VDL]; */
  /* const fnat m = VDL, n = 1u, ldA = m, ldAr = m, ldAi = m; */
  /* (void)printf("%d\n", zsplit_(&m, &n, A, &ldA, Ar, &ldAr, Ai, &ldAi)); */
  /* (void)printf("%d\n", zmerge_(&m, &n, Ar, &ldAr, Ai, &ldAi, A, &ldA)); */
  /* for (fnat i = 0u; i < VDL; ++i) */
  /*   (void)printf("%# .17e %# .17e\n", creal(A[i]), cimag(A[i])); */

  return EXIT_SUCCESS;
}
