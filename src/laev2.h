#ifndef LAEV2_H
#define LAEV2_H

#include "common.h"

extern void slevd2(const float a11, const float a22, const float a21, float l1[static restrict 1], float l2[static restrict 1], float cs1[static restrict 1], float sn1[static restrict 1]);
extern void clevd2(const float a11, const float a22, const float a21r, const float a21i, float l1[static restrict 1], float l2[static restrict 1], float cs1[static restrict 1], float snr[static restrict 1], float sni[static restrict 1]);

extern void dlevd2(const double a11, const double a22, const double a21, double l1[static restrict 1], double l2[static restrict 1], double cs1[static restrict 1], double sn1[static restrict 1]);
extern void zlevd2(const double a11, const double a22, const double a21r, const double a21i, double l1[static restrict 1], double l2[static restrict 1], double cs1[static restrict 1], double snr[static restrict 1], double sni[static restrict 1]);

#endif /* !LAEV2_H */
