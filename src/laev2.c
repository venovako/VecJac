#include "laev2.h"

void slevd2(const float a11, const float a22, const float a21, float l1[static restrict 1], float l2[static restrict 1], float cs1[static restrict 1], float sn1[static restrict 1])
{
  LAPACK_S(laev2)(&a11, &a21, &a22, l1, l2, cs1, sn1);
}

void clevd2(const float a11, const float a22, const float a21r, const float a21i, float l1[static restrict 1], float l2[static restrict 1], float cs1[static restrict 1], float snr[static restrict 1], float sni[static restrict 1])
{
  const float _Complex a = CMPLXF(a11, 0.0f);
  const float _Complex b = CMPLXF(a21r, -a21i);
  const float _Complex c = CMPLXF(a22, 0.0f);
  float _Complex sn1
#ifndef NDEBUG
    = CMPLXF(0.0f, 0.0f);
  *l1 = 0.0f;
  *l2 = 0.0f;
  *cs1 = 0.0f
#endif /* !NDEBUG */
    ;
  LAPACK_C(laev2)(&a, &b, &c, l1, l2, cs1, &sn1);
  *snr = crealf(sn1);
  *sni = cimagf(sn1);
}

void dlevd2(const double a11, const double a22, const double a21, double l1[static restrict 1], double l2[static restrict 1], double cs1[static restrict 1], double sn1[static restrict 1])
{
  LAPACK_D(laev2)(&a11, &a21, &a22, l1, l2, cs1, sn1);
}

void zlevd2(const double a11, const double a22, const double a21r, const double a21i, double l1[static restrict 1], double l2[static restrict 1], double cs1[static restrict 1], double snr[static restrict 1], double sni[static restrict 1])
{
  const double _Complex a = CMPLX(a11, 0.0);
  const double _Complex b = CMPLX(a21r, -a21i);
  const double _Complex c = CMPLX(a22, 0.0);
  double _Complex sn1
#ifndef NDEBUG
    = CMPLX(0.0, 0.0);
  *l1 = 0.0;
  *l2 = 0.0;
  *cs1 = 0.0
#endif /* !NDEBUG */
    ;
  LAPACK_Z(laev2)(&a, &b, &c, l1, l2, cs1, &sn1);
  *snr = creal(sn1);
  *sni = cimag(sn1);
}
