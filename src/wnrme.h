#ifndef WNRME_H
#define WNRME_H

#include "common.h"

extern wide wnrmer(const wide a11, const wide a22, const wide a21);
extern wide wnrmec(const wide a11, const wide a22, const wide a21r, const wide a21i);

extern wide wrer(const wide a11, const wide a22, const wide a21, const wide s, const wide t, const wide c, const wide l1, const wide l2, wide ae[static restrict 1]);
extern wide wrec(const wide a11, const wide a22, const wide a21r, const wide a21i, const wide s, const wide t, const wide c, const wide ca, const wide sa, const wide l1, const wide l2, wide ae[static restrict 1]);

#endif /* !WNRME_H */
