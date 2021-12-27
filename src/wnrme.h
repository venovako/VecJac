#ifndef WNRME_H
#define WNRME_H

#include "common.h"

extern wide wnrmer(const wide a11, const wide a22, const wide a21);
extern wide wnrmec(const wide a11, const wide a22, const wide a21r, const wide a21i);

extern wide wrer(const wide a11, const wide a22, const wide a21, const wide cs, const wide sn, const wide l1, const wide l2, wide ae[static restrict 1], wide an[static restrict 1]);
extern wide wrec(const wide a11, const wide a22, const wide a21r, const wide a21i, const wide cs, const wide ca, const wide sa, const wide l1, const wide l2, wide ae[static restrict 1], wide an[static restrict 1]);

extern wide wrerf(const wide a11, const wide a22, const wide a21, const wide t, const wide c, const wide l1, const wide l2, wide ae[static restrict 1], wide an[static restrict 1]);
extern wide wrecf(const wide a11, const wide a22, const wide a21r, const wide a21i, const wide t, const wide c, const wide ca, const wide sa, const wide l1, const wide l2, wide ae[static restrict 1], wide an[static restrict 1]);

extern wide wlam(wide L1, wide L2, wide l1, wide l2, wide relmax[static restrict 1], wide relmin[static restrict 1]);

#endif /* !WNRME_H */
