// icc -std=c18 -fp-model precise -prec-div -prec-sqrt -no-ftz subhypf.c
#ifdef __ICC
#include <mathimf.h>
#else /* !__ICC */
#include <math.h>
#endif /* ?__ICC */
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
  size_t c = 0u;
  for (float f = FLT_TRUE_MIN; f < FLT_MIN; f = nextafterf(f, FLT_MIN)) {
    if (hypotf(f, f) == f)
      ++c;
  }
  (void)fprintf(stdout, "float = %zu\n", c);
  (void)fflush(stdout);
  return EXIT_SUCCESS;
}
