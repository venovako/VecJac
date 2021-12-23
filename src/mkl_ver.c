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

int main(int argc, char *argv[])
{
  if (argc != 2)
    goto err;
  const char c = (char)toupper(argv[1u][0u]);

  if (c == 'C') {
    const int b = set_cbwr();
    switch (b & ~MKL_CBWR_STRICT) {
    case MKL_CBWR_AVX512_MIC:
      (void)printf("CBWR: MKL_CBWR_AVX512_MIC");
      break;
    case MKL_CBWR_AVX512:
      (void)printf("CBWR: MKL_CBWR_AVX512");
      break;
    case MKL_CBWR_AVX512_MIC_E1:
      (void)printf("CBWR: MKL_CBWR_AVX512_MIC_E1");
      break;
    case MKL_CBWR_AVX512_E1:
      (void)printf("CBWR: MKL_CBWR_AVX512_E1");
      break;
    default:
      (void)printf("CBWR: %d", b);
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
  (void)fprintf(stderr, "%s (C|V)\n", *argv);
  return EXIT_FAILURE;
}
