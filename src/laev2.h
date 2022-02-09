#ifndef LAEV2_H
#define LAEV2_H

#include "common.h"

#ifndef USE_MKL
extern void LAPACK_S(laev2)(const float A[static restrict 1], const float B[static restrict 1], const float C[static restrict 1], float RT1[static restrict 1], float RT2[static restrict 1], float CS1[static restrict 1], float SN1[static restrict 1]);
extern void LAPACK_D(laev2)(const double A[static restrict 1], const double B[static restrict 1], const double C[static restrict 1], double RT1[static restrict 1], double RT2[static restrict 1], double CS1[static restrict 1], double SN1[static restrict 1]);
#endif /* !USE_MKL */

extern void clevd2(const float a11[static restrict 1], const float a22[static restrict 1], const float a21r[static restrict 1], const float a21i[static restrict 1], float l1[static restrict 1], float l2[static restrict 1], float cs1[static restrict 1], float snr[static restrict 1], float sni[static restrict 1]);
extern void zlevd2(const double a11[static restrict 1], const double a22[static restrict 1], const double a21r[static restrict 1], const double a21i[static restrict 1], double l1[static restrict 1], double l2[static restrict 1], double cs1[static restrict 1], double snr[static restrict 1], double sni[static restrict 1]);

#endif /* !LAEV2_H */
