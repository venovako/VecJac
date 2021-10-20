#include "wdp.h"
#include "dnorme.h"

wide wdsq(const fnat n, const double x[static restrict 1])
{
  wide sq = W_ZERO;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,x) reduction(+:sq)
#endif /* !_OPENMP */
  for (fnat i = 0u; i < n; ++i) {
    const wide w = x[i];
    sq = fmaw(w, w, sq);
  }
  return sqrtw(sq);
}

wide wdn2(const fnat n, const double x[static restrict 1])
{
  const fint incx = 1;
  return BLAS_D(nrm2)(&n, x, &incx);
}

wide wddp(const fnat n, const double x[static restrict 1])
{
  const fint incx = 1;
  return sqrtw(BLAS_D(dot)(&n, x, &incx, x, &incx));
}

wide wdne(const fnat n, const double x[static restrict VDL])
{
  double e0, f0, e1, f1;
  return dnorme_(&n, x, &e0, &f0, &e1, &f1);
}

long double qdnre(const wide c, const wide e)
{
  const wide a = fabsw(c - e);
  return (long double)((e == W_ZERO) ? ((a == W_ZERO) ? W_ZERO : (a / e)) : (a / e));
}
