#include "wnrme.h"

// Frobenius norm of A

// return || A ||_F
wide wnrmer(const wide a11, const wide a22, const wide a21)
{
  return hypotw(hypotw(a11, a22), (W_SQRT2 * a21));
}

// return || A ||_F
wide wnrmec(const wide a11, const wide a22, const wide a21r, const wide a21i)
{
  return hypotw(hypotw(a11, a22), (W_SQRT2 * hypotw(a21r, a21i)));
}

// orthogonality checkers

// return cos^2 + sin^2 - 1
wide worr(const wide cs, const wide sn)
{
  return fmaw((sn + W_ONE), (sn - W_ONE), (cs * cs));
}

// return cos^2 + |sin|^2 - 1
wide worc(const wide cs, const wide ca, const wide sa)
{
  const wide _s = hypotw(ca, sa);
  return fmaw((_s + W_ONE), (_s - W_ONE), (cs * cs));
}

// absolute error checkers

// return || U L U^T - A ||_F
static wide waer(const wide a11, const wide a22, const wide a21, const wide cs, const wide sn, const wide l1, const wide l2)
{
  wide d11, d22;
  if (fabsw(cs) >= fabsw(sn)) {
    const wide t = sn / cs;
    const wide t2 = t * t;
    const wide c2 = cs * cs;
    d11 = fmaw(fmaw(l2, t2, l1), c2, -a11);
    d22 = fmaw(fmaw(l1, t2, l2), c2, -a22);
  }
  else {
    const wide g = cs / sn;
    const wide g2 = g * g;
    const wide s2 = sn * sn;
    d11 = fmaw(fmaw(l1, g2, l2), s2, -a11);
    d22 = fmaw(fmaw(l2, g2, l1), s2, -a22);
  }
  const wide d21 = fmaw((l1 - l2), (cs * sn), -a21);
  return wnrmer(d11, d22, d21);
}

// return || U L U^H - A ||_F
static wide waec(const wide a11, const wide a22, const wide a21r, const wide a21i, const wide cs, const wide ca, const wide sa, const wide l1, const wide l2)
{
  const wide _s = hypotw(ca, sa);
  wide d11, d22;
  if (fabsw(cs) >= _s) {
    const wide t = _s / cs;
    const wide t2 = t * t;
    const wide c2 = cs * cs;
    d11 = fmaw(fmaw(l2, t2, l1), c2, -a11);
    d22 = fmaw(fmaw(l1, t2, l2), c2, -a22);
  }
  else {
    const wide g = cs / _s;
    const wide g2 = g * g;
    const wide s2 = _s * _s;
    d11 = fmaw(fmaw(l1, g2, l2), s2, -a11);
    d22 = fmaw(fmaw(l2, g2, l1), s2, -a22);
  }
  const wide l = l1 - l2;
  const wide d21r = fmaw(l, (cs * ca), -a21r);
  const wide d21i = fmaw(l, (cs * sa), -a21i);
  return wnrmec(d11, d22, d21r, d21i);
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
