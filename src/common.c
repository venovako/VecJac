#include "common.h"

static_assert(CHAR_BIT == 8, "CHAR_BIT != 8");
#ifdef MKL_ILP64
static_assert(sizeof(fint) == 8, "sizeof(fint) != 8");
static_assert(sizeof(fnat) == 8, "sizeof(fnat) != 8");
#else /* LP64 */
static_assert(sizeof(fint) == 4, "sizeof(fint) != 4");
static_assert(sizeof(fnat) == 4, "sizeof(fnat) != 4");
#endif /* ?MKL_ILP64 */
static_assert(sizeof(float) == 4, "sizeof(float) != 4");
static_assert(sizeof(float complex) == 8, "sizeof(float complex) != 8");
static_assert(alignof(float complex) == alignof(float), "alignof(float complex) != alignof(float)");
static_assert(sizeof(double) == 8, "sizeof(double) != 8");
static_assert(sizeof(double complex) == 16, "sizeof(double complex) != 16");
static_assert(alignof(double complex) == alignof(double), "alignof(double complex) != alignof(double)");
static_assert(sizeof(long double) >= 8, "sizeof(long double) < 8");
static_assert(sizeof(long double complex) >= 16, "sizeof(long double complex) < 16");
static_assert(alignof(long double complex) == alignof(long double), "alignof(long double complex) != alignof(long double)");
#ifndef USE_EXTENDED
static_assert(sizeof(__float128) == 16, "sizeof(__float128) != 16");
static_assert(sizeof(__float128 complex) == 32, "sizeof(__float128 complex) != 32");
static_assert(alignof(__float128 complex) == alignof(__float128), "alignof(__float128 complex) != alignof(__float128)");
#endif /* !USE_EXTENDED */

size_t atoz(const char *const s)
{
  char *e = (char*)NULL;
  const size_t z = ((s && *s) ? (size_t)strtoull(s, &e, 0) : (size_t)0u);
  return ((e && *e) ? (size_t)0u : z);
}

char *stoa(char s[static restrict 17], const float x)
{
  int l = sprintf((char*)memset(s, 0, (size_t)17u), "%# -16.9E", (double)x);
  if (l <= 0)
    return (char*)NULL;
  char *d = s + 16;
  char *e = strrchr(s, 'E');
  if (e) {
    for (--d; isblank(*d); --d)
      *d = '\0';
    e += 2;
    l = (int)(strchr(e, '\0') - e);
    if (l >= 2)
      return s;
    d = s + 16;
    e += l;
    for (int i = 0; i < l; ++i)
      *--d = *--e;
    for (--d; isdigit(*d); --d)
      *d = '0';
  }
  else
    for (--d; !*d; --d)
      *d = ' ';
  return s;
}

char *dtoa(char s[static restrict 26], const double x)
{
  int l = sprintf((char*)memset(s, 0, (size_t)26u), "%# -25.17E", x);
  if (l <= 0)
    return (char*)NULL;
  char *d = s + 25;
  char *e = strrchr(s, 'E');
  if (e) {
    for (--d; isblank(*d); --d)
      *d = '\0';
    e += 2;
    l = (int)(strchr(e, '\0') - e);
    if (l >= 3)
      return s;
    d = s + 25;
    e += l;
    for (int i = 0; i < l; ++i)
      *--d = *--e;
    for (--d; isdigit(*d); --d)
      *d = '0';
  }
  else
    for (--d; !*d; --d)
      *d = ' ';
  return s;
}

char *xtoa(char s[static restrict 31], const long double x)
{
  int l = sprintf((char*)memset(s, 0, (size_t)31u), "%# -30.21LE", x);
  if (l <= 0)
    return (char*)NULL;
  char *d = s + 30;
  char *e = strrchr(s, 'E');
  if (e) {
    for (--d; isblank(*d); --d)
      *d = '\0';
    e += 2;
    l = (int)(strchr(e, '\0') - e);
    if (l >= 4)
      return s;
    d = s + 30;
    e += l;
    for (int i = 0; i < l; ++i)
      *--d = *--e;
    for (--d; isdigit(*d); --d)
      *d = '0';
  }
  else
    for (--d; !*d; --d)
      *d = ' ';
  return s;
}

int set_cbwr()
{
#ifdef USE_MKL
  const int b = mkl_cbwr_get_auto_branch();
  return ((mkl_cbwr_set(b) == MKL_CBWR_SUCCESS) ? b : -1);
#else /* !USE_MKL */
  return -2;
#endif /* ?USE_MKL */
}
