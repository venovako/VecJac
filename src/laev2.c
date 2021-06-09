#include "laev2.h"

void dlevd2(const double a11, const double a22, const double a21, double l1[static restrict 1], double l2[static restrict 1], double cs1[static restrict 1], double sn1[static restrict 1])
{
  LAPACK_D(laev2)(&a11, &a21, &a22, l1, l2, cs1, sn1);
}

void zlevd2(const double a11, const double a22, const double a21r, const double a21i, double l1[static restrict 1], double l2[static restrict 1], double cs1[static restrict 1], double snr[static restrict 1], double sni[static restrict 1])
{
  const double _Complex a = CMPLX(a11, W_ZERO);
  const double _Complex b = CMPLX(a21r, -a21i);
  const double _Complex c = CMPLX(a22, W_ZERO);
  double _Complex sn1;
  LAPACK_Z(laev2)(&a, &b, &c, l1, l2, cs1, &sn1);
  *snr = creal(sn1);
  *sni = cimag(sn1);
}

double dlevd2_pp(const wide cs1, const wide sn1)
{
  return (double)(sn1 / cs1);
}

double zlevd2_pp(const wide cs1, double snr[static restrict 1], double sni[static restrict 1])
{
  const wide sn1 = hypotw(*snr, *sni);
  if (sn1) {
    *snr = (double)(*snr / sn1);
    *sni = (double)(*sni / sn1);
  }
  return (double)(sn1 / cs1);
}
