#include "zsplit.h"

fint zsplit_(const fnat m[static restrict 1], const fnat n[static restrict 1], const double complex A[static restrict VDL_2], const fnat ldA[static restrict 1], double Ar[static restrict VDL], const fnat ldAr[static restrict 1], double Ai[static restrict VDL], const fnat ldAi[static restrict 1])
{
#ifndef NDEBUG
  if (IS_NOT_VFPENV)
    return -2;
  if (*m & VDL__2)
    return -1;
  if (IS_NOT_ALIGNED(A))
    return -3;
  if (*ldA < *m)
    return -4;
  if (*ldA & VDL__2)
    return -4;
  if (IS_NOT_ALIGNED(Ar))
    return -5;
  if (*ldAr < *m)
    return -6;
  if (*ldAr & VDL_1)
    return -6;
  if (IS_NOT_ALIGNED(Ai))
    return -7;
  if (*ldAi < *m)
    return -8;
  if (*ldAi & VDL_1)
    return -8;
#endif /* !NDEBUG */

  register const VI idx = _mm512_set_epi64(7, 5, 3, 1, 6, 4, 2, 0);

  for (fnat j = 0u; j < *n; ++j) {
    const double complex *const Aj = A + j * (*ldA);
    double *const Arj = Ar + j * (*ldAr);
    double *const Aij = Ai + j * (*ldAi);
    for (fnat i = 0u; i < *m; i += VDL_2) {
      register const VD ri = _mm512_permutexvar_pd(idx, _mm512_load_pd(Aj + i)); VDP(ri);
      _mm256_store_pd((Arj + i), _mm512_extractf64x4_pd(ri, 0x00u));
      _mm256_store_pd((Aij + i), _mm512_extractf64x4_pd(ri, 0x01u));
    }
  }

  return 0;
}
