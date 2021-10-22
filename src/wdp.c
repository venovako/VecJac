#include "wdp.h"
#include "dnorme.h"

double wsq(const fnat n, const double x[static restrict 1])
{
  wide sq = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,x) reduction(+:sq)
#endif /* !_OPENMP */
  for (fnat i = 0u; i < n; ++i) {
    const wide w = x[i];
    sq = fmaw(w, w, sq);
  }
  return (double)sqrtw(sq);
}

double dn2(const fnat n, const double x[static restrict 1])
{
  const fint incx = 1;
  return BLAS_D(nrm2)(&n, x, &incx);
}

double ddp(const fnat n, const double x[static restrict 1])
{
  const fint incx = 1;
  return (double)sqrt(BLAS_D(dot)(&n, x, &incx, x, &incx));
}

double xdp(const fnat n, const double x[static restrict 1])
{
  long double sq = 0.0L;
  for (fnat i = 0u; i < n; ++i) {
    const long double l = x[i];
    sq += (l * l);
  }
  return (double)sqrtl(sq);
}

double dne(const fnat n, const double x[static restrict VDL])
{
  double e0, f0, e1, f1;
  return dnorme_(&n, x, &e0, &f0, &e1, &f1);
}

double dre(const double c, const double e)
{
  const double a = fabs(c - e);
  return ((e == 0.0) ? ((a == 0.0) ? 0.0 : (a / e)) : (a / e));
}
