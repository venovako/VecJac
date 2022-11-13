#ifndef LAEV2_H
#define LAEV2_H

#include "common.h"

#ifdef USE_INL
// inline version of DLAEV2 from the Fortran sources
// see https://github.com/Reference-LAPACK/lapack for their LICENSE
static inline void _dlaev2(const double A[static restrict 1], const double B[static restrict 1], const double C[static restrict 1], double RT1[static restrict 1], double RT2[static restrict 1], double CS1[static restrict 1], double SN1[static restrict 1])
{
  const double _A = *A;
  const double _B = *B;
  const double _C = *C;

  const double SM = _A + _C;
  const double DF = _A - _C;
  const double ADF = fabs(DF);
  const double TB = _B + _B;
  const double AB = fabs(TB);

  double ACMX, ACMN, RT, CS, CT, TN;
  int SGN1, SGN2;

  if (fabs(_A) > fabs(_C)) {
    ACMX = _A;
    ACMN = _C;
  }
  else {
    ACMX = _C;
    ACMN = _A;
  }

  if (ADF > AB) {
    const double AB_ADF = AB / ADF;
    RT = ADF * sqrt(fma(AB_ADF, AB_ADF, 1.0));
  }
  else if (ADF < AB) {
    const double ADF_AB = ADF / AB;
    RT = AB * sqrt(fma(ADF_AB, ADF_AB, 1.0));
  }
  else
    RT = AB * M_SQRT2;

  if (SM < 0.0) {
    *RT1 = 0.5 * (SM - RT);
    SGN1 = -1;
    *RT2 = (ACMX / *RT1) * ACMN - (_B / *RT1) * _B;
  }
  else if (SM > 0.0) {
    *RT1 = 0.5 * (SM + RT);
    SGN1 = 1;
    *RT2 = (ACMX / *RT1) * ACMN - (_B / *RT1) * _B;
  }
  else {
    *RT1 = 0.5 * RT;
    *RT2 = -*RT1;
    SGN1 = 1;
  }

  if (DF >= 0.0) {
    CS = DF + RT;
    SGN2 = 1;
  }
  else {
    CS = DF - RT;
    SGN2 = -1;
  }

  if (fabs(CS) > AB) {
    CT = -TB / CS;
    *SN1 = 1.0 / sqrt(fma(CT, CT, 1.0));
    *CS1 = *SN1 * CT;
  }
  else {
    if (AB == 0.0) {
      *CS1 = 1.0;
      *SN1 = 0.0;
    }
    else {
      TN = -CS / TB;
      *CS1 = 1.0 / sqrt(fma(TN, TN, 1.0));
      *SN1 = *CS1 * TN;
    }
  }

  if (SGN1 == SGN2) {
    TN = *CS1;
    *CS1 = -*SN1;
    *SN1 = TN;
  }
}

// inline version of SLAEV2 from the Fortran sources
// see https://github.com/Reference-LAPACK/lapack for their LICENSE
static inline void _slaev2(const float A[static restrict 1], const float B[static restrict 1], const float C[static restrict 1], float RT1[static restrict 1], float RT2[static restrict 1], float CS1[static restrict 1], float SN1[static restrict 1])
{
  const float _A = *A;
  const float _B = *B;
  const float _C = *C;

  const float SM = _A + _C;
  const float DF = _A - _C;
  const float ADF = fabsf(DF);
  const float TB = _B + _B;
  const float AB = fabsf(TB);

  float ACMX, ACMN, RT, CS, CT, TN;
  int SGN1, SGN2;

  if (fabsf(_A) > fabsf(_C)) {
    ACMX = _A;
    ACMN = _C;
  }
  else {
    ACMX = _C;
    ACMN = _A;
  }

  if (ADF > AB) {
    const float AB_ADF = AB / ADF;
    RT = ADF * sqrtf(fmaf(AB_ADF, AB_ADF, 1.0f));
  }
  else if (ADF < AB) {
    const float ADF_AB = ADF / AB;
    RT = AB * sqrtf(fmaf(ADF_AB, ADF_AB, 1.0f));
  }
  else
    RT = AB * (float)M_SQRT2;

  if (SM < 0.0f) {
    *RT1 = 0.5f * (SM - RT);
    SGN1 = -1;
    *RT2 = (ACMX / *RT1) * ACMN - (_B / *RT1) * _B;
  }
  else if (SM > 0.0f) {
    *RT1 = 0.5f * (SM + RT);
    SGN1 = 1;
    *RT2 = (ACMX / *RT1) * ACMN - (_B / *RT1) * _B;
  }
  else {
    *RT1 = 0.5f * RT;
    *RT2 = -*RT1;
    SGN1 = 1;
  }

  if (DF >= 0.0f) {
    CS = DF + RT;
    SGN2 = 1;
  }
  else {
    CS = DF - RT;
    SGN2 = -1;
  }

  if (fabsf(CS) > AB) {
    CT = -TB / CS;
    *SN1 = 1.0f / sqrtf(fmaf(CT, CT, 1.0f));
    *CS1 = *SN1 * CT;
  }
  else {
    if (AB == 0.0f) {
      *CS1 = 1.0f;
      *SN1 = 0.0f;
    }
    else {
      TN = -CS / TB;
      *CS1 = 1.0f / sqrtf(fmaf(TN, TN, 1.0f));
      *SN1 = *CS1 * TN;
    }
  }

  if (SGN1 == SGN2) {
    TN = *CS1;
    *CS1 = -*SN1;
    *SN1 = TN;
  }
}

// inline version of ZLAEV2 with B=(2,1) instead of (1,2) from the Fortran sources
// see https://github.com/Reference-LAPACK/lapack for their LICENSE
static inline void _zlaev2(const double A[static restrict 1], const double Br[static restrict 1], const double Bi[static restrict 1], const double C[static restrict 1], double RT1[static restrict 1], double RT2[static restrict 1], double CS1[static restrict 1], double SNr[static restrict 1], double SNi[static restrict 1])
{
  const double _Br = *Br;
  const double _Bi = *Bi;
  const double _B = cr_hypot(_Br, _Bi);
  _dlaev2(A, &_B, C, RT1, RT2, CS1, SNr);
  if (_B == 0.0)
    *SNi = copysign(0.0, *SNr);
  else {
    *SNi = *SNr * (_Bi / _B);
    *SNr *= (_Br / _B);
  }
}

// inline version of CLAEV2 with B=(2,1) instead of (1,2) from the Fortran sources
// see https://github.com/Reference-LAPACK/lapack for their LICENSE
static inline void _claev2(const float A[static restrict 1], const float Br[static restrict 1], const float Bi[static restrict 1], const float C[static restrict 1], float RT1[static restrict 1], float RT2[static restrict 1], float CS1[static restrict 1], float SNr[static restrict 1], float SNi[static restrict 1])
{
  const float _Br = *Br;
  const float _Bi = *Bi;
  const float _B = cr_hypotf(_Br, _Bi);
  _slaev2(A, &_B, C, RT1, RT2, CS1, SNr);
  if (_B == 0.0f)
    *SNi = copysignf(0.0f, *SNr);
  else {
    *SNi = *SNr * (_Bi / _B);
    *SNr *= (_Br / _B);
  }
}
#else /* !USE_INL */
#ifndef USE_MKL
extern void LAPACK_S(laev2)(const float A[static restrict 1], const float B[static restrict 1], const float C[static restrict 1], float RT1[static restrict 1], float RT2[static restrict 1], float CS1[static restrict 1], float SN1[static restrict 1]);
extern void LAPACK_D(laev2)(const double A[static restrict 1], const double B[static restrict 1], const double C[static restrict 1], double RT1[static restrict 1], double RT2[static restrict 1], double CS1[static restrict 1], double SN1[static restrict 1]);
#endif /* !USE_MKL */
extern void _claev2(const float a11[static restrict 1], const float a21r[static restrict 1], const float a21i[static restrict 1], const float a22[static restrict 1], float l1[static restrict 1], float l2[static restrict 1], float cs1[static restrict 1], float snr[static restrict 1], float sni[static restrict 1]);
extern void _zlaev2(const double a11[static restrict 1], const double a21r[static restrict 1], const double a21i[static restrict 1], const double a22[static restrict 1], double l1[static restrict 1], double l2[static restrict 1], double cs1[static restrict 1], double snr[static restrict 1], double sni[static restrict 1]);
#endif /* ?USE_INL */

#endif /* !LAEV2_H */
