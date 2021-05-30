#ifndef WNRME_H
#define WNRME_H

#include "common.h"

extern wide wnrmer(const wide a11, const wide a22, const wide a21);
extern wide wnrmec(const wide a11, const wide a22, const wide a21r, const wide a21i);

extern wide wrer(const wide a11, const wide a22, const wide a21, const wide s, const wide t, const wide c, const wide l1, const wide l2, wide ae[static restrict 1], wide an[static restrict 1], wide L1[static restrict 1], wide L2[static restrict 1]);
extern wide wrec(const wide a11, const wide a22, const wide a21r, const wide a21i, const wide s, const wide t, const wide c, const wide ca, const wide sa, const wide l1, const wide l2, wide ae[static restrict 1], wide an[static restrict 1], wide L1[static restrict 1], wide L2[static restrict 1]);

extern wide wlam(const wide L1, const wide L2, const wide l1, const wide l2, wide relmax[static restrict 1], wide relmin[static restrict 1]);

#endif /* !WNRME_H */
