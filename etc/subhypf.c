// icc -std=c18 -fp-model precise -prec-div -prec-sqrt -no-ftz -diag-disable=10441 subhypf.c
#if (defined(__ICC) || defined(__INTEL_COMPILER) || defined(__INTEL_CLANG_COMPILER) || defined(__INTEL_LLVM_COMPILER))
#include <mathimf.h>
#else /* !__ICC */
#include <math.h>
#endif /* ?__ICC */
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
  size_t c = 0u;
  for (float f = FLT_TRUE_MIN; f < FLT_MIN; f = nextafterf(f, FLT_MIN))
    if (hypotf(f, f) == f)
      (void)printf("%zu: %#15.9E (= %#A)\n", ++c, f, f);
  return EXIT_SUCCESS;
}
