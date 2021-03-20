#include "vec.h"

int VSprintf(FILE f[static 1], const char *const h, const VS v)
{
  alignas(VA) float s[VSL];

  int ret = (h ? fprintf(f, "\nL: %s\n", h) : 0);
  if (ret < 0) {
    perror("fprintf");
    return -1;
  }

  if (fflush(f)) {
    perror("fflush");
    return -2;
  }

  _mm512_store_ps(s, v);

  char a[17];
  for (unsigned i = 0u; i < VSL; ++i) {
    char *const p = stoa(a, s[i]);
    if (!p) {
      perror("sprintf");
      return -(int)(i + 4u);
    }
    if (20 != fprintf(f, "%X: %s\n", i, p)) {
      perror("fprintf");
      return -(int)(i + 4u);
    }
    ret += 20;
  }

  return (fflush(f) ? -3 : ret);
}

int VDprintf(FILE f[static 1], const char *const h, const VD v)
{
  alignas(VA) double d[VDL];

  int ret = (h ? fprintf(f, "\nL: %s\n", h) : 0);
  if (ret < 0) {
    perror("fprintf");
    return -1;
  }

  if (fflush(f)) {
    perror("fflush");
    return -2;
  }

  _mm512_store_pd(d, v);

  char a[26];
  for (unsigned i = 0u; i < VDL; ++i) {
    char *const p = dtoa(a, d[i]);
    if (!p) {
      perror("sprintf");
      return -(int)(i + 4u);
    }
    if (29 != fprintf(f, "%u: %s\n", i, p)) {
      perror("fprintf");
      return -(int)(i + 4u);
    }
    ret += 29;
  }

  return (fflush(f) ? -3 : ret);
}

int MSprintf(FILE f[static 1], const char *const h, const MS m)
{
  if (fflush(f)) {
    perror("fflush");
    return -1;
  }

  int ret = (h ? fprintf(f, "\n%s: ", h) : 0);
  if (ret < 0) {
    perror("fprintf");
    return -2;
  }

  const unsigned u = MS2U(m);
  for (unsigned i = 0u, o = (1u << VSL_1); i < VSL; ++i) {
    if (1 != fprintf(f, "%c", ((u & o) ? '1' : '0'))) {
      perror("fprintf");
      return -(int)(i + 4u);
    }
    o >>= 1u;
    ++ret;
  }

  if (8 != fprintf(f, " (%04X)\n", u)) {
    perror("fprintf");
    return -(int)(VSL + 4u);
  }
  ret += 8;

  return (fflush(f) ? -3 : ret);
}

int MDprintf(FILE f[static 1], const char *const h, const MD m)
{
  if (fflush(f)) {
    perror("fflush");
    return -1;
  }

  int ret = (h ? fprintf(f, "\n%s: ", h) : 0);
  if (ret < 0) {
    perror("fprintf");
    return -2;
  }

  const unsigned u = MD2U(m);
  for (unsigned i = 0u, o = (1u << VDL_1); i < VDL; ++i) {
    if (1 != fprintf(f, "%c", ((u & o) ? '1' : '0'))) {
      perror("fprintf");
      return -(int)(i + 4u);
    }
    o >>= 1u;
    ++ret;
  }

  if (6 != fprintf(f, " (%02X)\n", u)) {
    perror("fprintf");
    return -(int)(VDL + 4u);
  }
  ret += 6;

  return (fflush(f) ? -3 : ret);
}

size_t VSread(float s[static VSL], const size_t V, FILE f[static 1])
{
  return fread(s, sizeof(float) * VSL, V, f);
}

size_t VDread(double d[static VDL], const size_t V, FILE f[static 1])
{
  return fread(d, sizeof(double) * VDL, V, f);
}

size_t VSwrite(const float s[static VSL], const size_t V, FILE f[static 1])
{
  return fwrite(s, sizeof(float) * VSL, V, f);
}

size_t VDwrite(const double d[static VDL], const size_t V, FILE f[static 1])
{
  return fwrite(d, sizeof(double) * VDL, V, f);
}
