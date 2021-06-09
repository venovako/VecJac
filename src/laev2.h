#ifndef LAEV2_H
#define LAEV2_H

#include "common.h"

extern void dlevd2(const double a11, const double a22, const double a21, double l1[static restrict 1], double l2[static restrict 1], double cs1[static restrict 1], double sn1[static restrict 1]);
extern void zlevd2(const double a11, const double a22, const double a21r, const double a21i, double l1[static restrict 1], double l2[static restrict 1], double cs1[static restrict 1], double snr[static restrict 1], double sni[static restrict 1]);

extern double dlevd2_pp(const wide cs1, const wide sn1);
extern double zlevd2_pp(const wide cs1, double snr[static restrict 1], double sni[static restrict 1]);

#endif /* !LAEV2_H */
