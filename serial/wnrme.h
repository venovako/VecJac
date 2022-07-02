#ifndef WNRME_H
#define WNRME_H

#include "common.h"

EXTERN wide wnrmer(const wide a11, const wide a22, const wide a21);
EXTERN wide wnrmec(const wide a11, const wide a22, const wide a21r, const wide a21i);

EXTERN wide worr(const wide cs, const wide sn);
EXTERN wide worc(const wide cs, const wide ca, const wide sa);

EXTERN wide wrer(const wide a11, const wide a22, const wide a21, const wide cs, const wide sn, const wide l1, const wide l2, wide *const ae, wide *const an);
EXTERN wide wrec(const wide a11, const wide a22, const wide a21r, const wide a21i, const wide cs, const wide ca, const wide sa, const wide l1, const wide l2, wide *const ae, wide *const an);

EXTERN wide wlam(wide L1, wide L2, wide l1, wide l2, wide *const relmax, wide *const relmin);

#endif /* !WNRME_H */
