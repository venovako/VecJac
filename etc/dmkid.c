#if (defined(__ICC) || defined(__INTEL_CLANG_COMPILER) || defined(__INTEL_LLVM_COMPILER))
#include <mathimf.h>
#else /* !__ICC */
#include <math.h>
#endif /* ?__ICC */
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  if (argc != 3) {
    (void)fprintf(stderr, "%s N FN\n", *argv);
    return EXIT_FAILURE;
  }
  const size_t n = (size_t)atoll(argv[1u]);
  if (!n)
    return EXIT_SUCCESS;

  FILE *const f = fopen(argv[2u], "wb");
  if (!f)
    return EXIT_FAILURE;

  for (size_t j = 0u; j < n; ++j)
    for (size_t i = 0u; i < n; ++i) {
      const double d = ((i == j) ? 1.0 : 0.0);
      if (1u != fwrite(&d, sizeof(d), 1u, f))
        return EXIT_FAILURE;
    }

  return (fclose(f) ? EXIT_FAILURE : EXIT_SUCCESS);
}
