#include "common.h"

size_t atoz(const char *const s)
{
  char *e = (char*)NULL;
  const size_t z = ((s && *s) ? (size_t)strtoull(s, &e, 0) : (size_t)0u);
  return ((e && *e) ? (size_t)0u : z);
}

char *stoa(char *const s, const float x)
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

char *dtoa(char *const s, const double x)
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

char *xtoa(char *const s, const long double x)
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
