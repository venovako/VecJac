#include "laev2.h"

#ifndef USE_MKL
extern void LAPACK_C(laev2)(const float complex A[static restrict 1], const float complex B[static restrict 1], const float complex C[static restrict 1], float RT1[static restrict 1], float RT2[static restrict 1], float CS1[static restrict 1], float complex SN1[static restrict 1]);
extern void LAPACK_Z(laev2)(const double complex A[static restrict 1], const double complex B[static restrict 1], const double complex C[static restrict 1], double RT1[static restrict 1], double RT2[static restrict 1], double CS1[static restrict 1], double complex SN1[static restrict 1]);
#endif /* !USE_MKL */

void clevd2(const float a11[static restrict 1], const float a22[static restrict 1], const float a21r[static restrict 1], const float a21i[static restrict 1], float l1[static restrict 1], float l2[static restrict 1], float cs1[static restrict 1], float snr[static restrict 1], float sni[static restrict 1])
{
  float b[2u] = { *a21r, -*a21i };
  float sn1[2u];
  LAPACK_C(laev2)((const float complex*)a11, (const float complex*)b, (const float complex*)a22, l1, l2, cs1, (float complex*)sn1);
  *snr = sn1[0u];
  *sni = sn1[1u];
}

void zlevd2(const double a11[static restrict 1], const double a22[static restrict 1], const double a21r[static restrict 1], const double a21i[static restrict 1], double l1[static restrict 1], double l2[static restrict 1], double cs1[static restrict 1], double snr[static restrict 1], double sni[static restrict 1])
{
  double b[2u] = { *a21r, -*a21i };
  double sn1[2u];
  LAPACK_Z(laev2)((const double complex*)a11, (const double complex*)b, (const double complex*)a22, l1, l2, cs1, (double complex*)sn1);
  *snr = sn1[0u];
  *sni = sn1[1u];
}
