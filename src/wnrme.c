#include "wnrme.h"

// return || A ||_F
wide wnrmer(const wide a11, const wide a22, const wide a21)
{
  return hypotw(hypotw(a11, a22), hypotw(a21, a21));
}

// return || A ||_F
wide wnrmec(const wide a11, const wide a22, const wide a21r, const wide a21i)
{
  return hypotw(hypotw(a11, a22), (W_SQRT2 * hypotw(a21r, a21i)));
}

// absolute error checkers

// return || U L U^T - A ||_F
static wide waer(const wide a11, const wide a22, const wide a21, const wide cs, const wide sn, const wide l1, const wide l2)
{
  const wide c2 = cs * cs;
  const wide t = sn / cs;
  const wide t2 = t * t;
  const wide l = l1 - l2;
  const wide tl = t * l;
  return wnrmer(fmaw(c2, fmaw(l2, t2, l1), -a11), fmaw(c2, fmaw(l1, t2, l2), -a22), fmaw(c2, tl, -a21));
}

// return || U L U^H - A ||_F
static wide waec(const wide a11, const wide a22, const wide a21r, const wide a21i, const wide cs, wide ca, wide sa, const wide l1, const wide l2)
{
  const wide c2 = cs * cs;
  const wide sn = hypotw(ca, sa);
  if (sn != W_ZERO) {
    ca /= sn;
    sa /= sn;
  }
  const wide t = sn / cs;
  const wide t2 = t * t;
  const wide l = l1 - l2;
  const wide tl = t * l;
  return wnrmec(fmaw(c2, fmaw(l2, t2, l1), -a11), fmaw(c2, fmaw(l1, t2, l2), -a22), fmaw(c2, (ca * tl), -a21r), fmaw(c2, (sa * tl), -a21i));
}

// return || U L U^T - A ||_F
static wide waerf(const wide a11, const wide a22, const wide a21, const wide t, const wide c, const wide l1, const wide l2)
{
  const wide l = l1 - l2;
  const wide t2 = t * t;
  const wide c2 = c * c;
  const wide tl = t * l;
  return wnrmer(fmaw(c2, fmaw(l2, t2, l1), -a11), fmaw(c2, fmaw(l1, t2, l2), -a22), fmaw(c2, tl, -a21));
}

// return || U L U^H - A ||_F
static wide waecf(const wide a11, const wide a22, const wide a21r, const wide a21i, const wide t, const wide c, const wide ca, const wide sa, const wide l1, const wide l2)
{
  const wide l = l1 - l2;
  const wide t2 = t * t;
  const wide c2 = c * c;
  const wide tl = t * l;
  return wnrmec(fmaw(c2, fmaw(l2, t2, l1), -a11), fmaw(c2, fmaw(l1, t2, l2), -a22), fmaw(c2, (ca * tl), -a21r), fmaw(c2, (sa * tl), -a21i));
}

// relative error checkers

// return || U L U^T - A ||_F / || A ||_F
wide wrer(const wide a11, const wide a22, const wide a21, const wide cs, const wide sn, const wide l1, const wide l2, wide ae[static restrict 1], wide an[static restrict 1])
{
  *ae = waer(a11, a22, a21, cs, sn, l1, l2);
  *an = wnrmer(a11, a22, a21);
  return fmaxw((*ae / *an), W_ZERO);
}

// return || U L U^H - A ||_F / || A ||_F
wide wrec(const wide a11, const wide a22, const wide a21r, const wide a21i, const wide cs, const wide ca, const wide sa, const wide l1, const wide l2, wide ae[static restrict 1], wide an[static restrict 1])
{
  *ae = waec(a11, a22, a21r, a21i, cs, ca, sa, l1, l2);
  *an = wnrmec(a11, a22, a21r, a21i);
  return fmaxw((*ae / *an), W_ZERO);
}

// return || U L U^T - A ||_F / || A ||_F
wide wrerf(const wide a11, const wide a22, const wide a21, const wide t, const wide c, const wide l1, const wide l2, wide ae[static restrict 1], wide an[static restrict 1])
{
  *ae = waerf(a11, a22, a21, t, c, l1, l2);
  *an = wnrmer(a11, a22, a21);
  return fmaxw((*ae / *an), W_ZERO);
}

// return || U L U^H - A ||_F / || A ||_F
wide wrecf(const wide a11, const wide a22, const wide a21r, const wide a21i, const wide t, const wide c, const wide ca, const wide sa, const wide l1, const wide l2, wide ae[static restrict 1], wide an[static restrict 1])
{
  *ae = waecf(a11, a22, a21r, a21i, t, c, ca, sa, l1, l2);
  *an = wnrmec(a11, a22, a21r, a21i);
  return fmaxw((*ae / *an), W_ZERO);
}

// L1, L2: actual computed eigenvalues; l1, l2: eigenvalues generated for testing
// relmax: relative error on the larger eigenvalue; relmin: relative error on the smaller eigenvalue
// return max{|| L - l ||_F / || l ||_F, 0}
wide wlam(wide L1, wide L2, wide l1, wide l2, wide relmax[static restrict 1], wide relmin[static restrict 1])
{
  L1 = fabsw(L1);
  L2 = fabsw(L2);
  l1 = fabsw(l1);
  l2 = fabsw(l2);

  wide Lmin, Lmax;
  if (L1 <= L2) {
    Lmin = L1;
    Lmax = L2;
  }
  else {
    Lmin = L2;
    Lmax = L1;
  }

  wide lmin, lmax;
  if (l1 <= l2) {
    lmin = l1;
    lmax = l2;
  }
  else {
    lmin = l2;
    lmax = l1;
  }

  const wide
    aelmax = ((Lmax <= lmax) ? (lmax - Lmax) : (Lmax - lmax)),
    aelmin = ((Lmin <= lmin) ? (lmin - Lmin) : (Lmin - lmin));

  *relmax = fmaxw((aelmax / lmax), W_ZERO);
  *relmin = fmaxw((aelmin / lmin), W_ZERO);

  return fmaxw((hypotw(aelmax, aelmin) / hypotw(lmax, lmin)), W_ZERO);
}
