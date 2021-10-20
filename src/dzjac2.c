#include "dzjac2.h"

fint dzjac2_pp(const fnat n, const double *const restrict s, const double *const restrict c, const double l1[static restrict 1], const double l2[static restrict 1], const unsigned p[static restrict 1], double *const restrict L1, double *const restrict L2)
{
#ifdef _OPENMP
  fint t = 0;
#endif /* _OPENMP */

  if (s && c && L1 && L2) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(n,s,c,p,l1,l2,L1,L2) reduction(max:t)
#endif /* _OPENMP */
    for (fnat i = 0u; i < n; ++i) {
      const double c2 = c[i] * c[i];
      const int _s = -(int)(s[i]);
      const unsigned pi = p[i >> VDLlg];
      if (pi & (1u << (i & VDL_1))) {
        L2[i] = scalbn((c2 * l1[i]), _s);
        L1[i] = scalbn((c2 * l2[i]), _s);
      }
      else {
        L1[i] = scalbn((c2 * l1[i]), _s);
        L2[i] = scalbn((c2 * l2[i]), _s);
      }
#ifdef _OPENMP
      t = imax(t, omp_get_thread_num());
#endif /* _OPENMP */
    }
  }

#ifdef _OPENMP
  return (t + 1);
#else /* !_OPENMP */
  return 0;
#endif /* ?_OPENMP */
}
