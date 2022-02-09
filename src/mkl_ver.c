#include "common.h"

/* struct MKLVersion:
  int    MajorVersion;
  int    MinorVersion;
  int    UpdateVersion;
  char * ProductStatus;
  char * Build;
  char * Processor;
  char * Platform;
 */

#ifndef USE_MKL
extern void LAPACK_I(laver)(fint major[static restrict 1], fint minor[static restrict 1], fint patch[static restrict 1]);
#endif /* !USE_MKL */

int main(int argc, char *argv[])
{
#ifdef USE_MKL
  if (argc != 2)
    goto err;
  const char c = (char)toupper(argv[1u][0u]);

  if (c == 'C') {
    const int b = set_cbwr();
    (void)printf("CBWR: %#x; ", (unsigned)b);
    switch (b & ~MKL_CBWR_STRICT) {
    case MKL_CBWR_AVX512_MIC:
      (void)printf("MKL_CBWR_AVX512_MIC");
      break;
    case MKL_CBWR_AVX512:
      (void)printf("MKL_CBWR_AVX512");
      break;
    case MKL_CBWR_AVX512_MIC_E1:
      (void)printf("MKL_CBWR_AVX512_MIC_E1");
      break;
    case MKL_CBWR_AVX512_E1:
      (void)printf("MKL_CBWR_AVX512_E1");
      break;
    default:
      (void)printf("???");
    }
    if (b & MKL_CBWR_STRICT)
      (void)printf(",STRICT");
    (void)printf("\n");
  }
  else if (c == 'V') {
    MKLVersion v;
    mkl_get_version(&v);
    (void)printf("MajorVersion: %d\n", v.MajorVersion);
    (void)printf("MinorVersion: %d\n", v.MinorVersion);
    (void)printf("UpdateVersion: %d\n", v.UpdateVersion);
    (void)printf("ProductStatus: %s\n", v.ProductStatus);
    (void)printf("Build: %s\n", v.Build);
    (void)printf("Processor: %s\n", v.Processor);
    (void)printf("Platform: %s\n", v.Platform);
  }
  else
    goto err;
  return EXIT_SUCCESS;

 err:
  (void)fprintf(stderr, "%s { C | V }\n", *argv);
  return EXIT_FAILURE;
#else /* LAPACK */
  if (argc != 1) {
    (void)fprintf(stderr, "%s\n", *argv);
    return EXIT_FAILURE;
  }
  fint major = 0, minor = 0, patch = 0;
  LAPACK_I(laver)(&major, &minor, &patch);
  (void)printf("LAPACK %lld.%lld.%lld\n", (long long)major, (long long)minor, (long long)patch);
  return EXIT_SUCCESS;
#endif /* ?USE_MKL */
}
