#ifndef tk_vec_name
#error "tk_vec_name missing"
#endif

#ifndef tk_vec_base
#error "tk_vec_base missing"
#endif

#define tk_str(x) #x
#define tk_xstr(x) tk_str(x)
#define tk_strcat2(a, b) a##_##b
#define tk_strcat(a, b) tk_strcat2(a, b)
#define tk_vec_pfx(name) tk_strcat(tk_vec_name, name)

#ifndef tk_vec_limited
static inline void tk_vec_pfx(rasc) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, m0->n, NULL, NULL);
  tk_rvec_t *ranks = tk_rvec_create(L, cols, NULL, NULL);
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    for (size_t c = 0; c < cols; c ++) {
      ranks->a[c].i = (int64_t) c;
      ranks->a[c].d = m0->a[r * cols + c];
    }
    tk_rvec_asc(ranks, 0, ranks->n);
    for (size_t c = 0; c < cols; c ++)
      out->a[r * cols + c] = ranks->a[c].i;
  }
  lua_pop(L, 1);
}
#endif

#ifndef tk_vec_limited
static inline void tk_vec_pfx(rdesc) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, m0->n, NULL, NULL);
  tk_rvec_t *ranks = tk_rvec_create(L, cols, NULL, NULL);
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    for (size_t c = 0; c < cols; c ++) {
      ranks->a[c].i = (int64_t) c;
      ranks->a[c].d = m0->a[r * cols + c];
    }
    tk_rvec_desc(ranks, 0, ranks->n);
    for (size_t c = 0; c < cols; c ++)
      out->a[r * cols + c] = ranks->a[c].i;
  }
  lua_pop(L, 1);
}
#endif

#ifndef tk_vec_limited
static inline void tk_vec_pfx(casc) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  uint64_t rows = m0->n / cols;
  tk_ivec_t *out = tk_ivec_create(L, m0->n, NULL, NULL);
  tk_rvec_t *ranks = tk_rvec_create(L, rows, NULL, NULL);
  for (size_t c = 0; c < cols; c ++) {
    for (uint64_t r = 0; r < rows ; r ++) {
      ranks->a[r].i = (int64_t) r;
      ranks->a[r].d = m0->a[c * rows + r];
    }
    tk_rvec_asc(ranks, 0, ranks->n);
    for (size_t r = 0; r < rows; r ++)
      out->a[c * rows + r] = ranks->a[r].i;
  }
  lua_pop(L, 1);
}
#endif

#ifndef tk_vec_limited
static inline void tk_vec_pfx(cdesc) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  uint64_t rows = m0->n / cols;
  tk_ivec_t *out = tk_ivec_create(L, m0->n, NULL, NULL);
  tk_rvec_t *ranks = tk_rvec_create(L, rows, NULL, NULL);
  for (size_t c = 0; c < cols; c ++) {
    for (uint64_t r = 0; r < rows ; r ++) {
      ranks->a[r].i = (int64_t) r;
      ranks->a[r].d = m0->a[c * rows + r];
    }
    tk_rvec_desc(ranks, 0, ranks->n);
    for (size_t r = 0; r < rows; r ++)
      out->a[c * rows + r] = ranks->a[r].i;
  }
  lua_pop(L, 1);
}
#endif

#ifndef tk_vec_limited
static inline tk_dvec_t *tk_vec_pfx(cmagnitudes) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_dvec_t *out = tk_dvec_create(L, cols, NULL, NULL);
  for (uint64_t c = 0; c < cols; c ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t r = 0; r < m0->n / cols; r ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val * val;
    }
    out->a[c] = sqrt(sum);
  }
  return out;
}
#endif

#ifndef tk_vec_limited
static inline tk_dvec_t *tk_vec_pfx(rmagnitudes) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_dvec_t *out = tk_dvec_create(L, m0->n / cols, NULL, NULL);
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val * val;
    }
    out->a[r] = sqrt(sum);
  }
  return out;
}
#endif

#ifndef tk_vec_limited
static inline tk_ivec_t *tk_vec_pfx(cmaxargs) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, cols, NULL, NULL);
  for (uint64_t c = 0; c < cols; c ++) {
    uint64_t maxr = 0;
    tk_vec_base maxv = m0->a[0 * cols + c];
    for (size_t r = 1; r < m0->n / cols; r ++) {
      if (m0->a[r * cols + c] > maxv) {
        maxr = r;
        maxv = m0->a[r * cols + c];
      }
    }
    out->a[c] = (int64_t) maxr;
  }
  return out;
}
#endif

#ifndef tk_vec_limited
static inline tk_ivec_t *tk_vec_pfx(rmaxargs) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, m0->n / cols, NULL, NULL);
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val;
    }
    out->a[r] = sum;
    uint64_t maxc = 0;
    tk_vec_base maxv = m0->a[r * cols + 0];
    for (uint64_t c = 1; c < cols; c ++) {
      if (m0->a[r * cols + c] > maxv) {
        maxc = c;
        maxv = m0->a[r * cols + c];
      }
    }
    out->a[r] = (int64_t) maxc;
  }
  return out;
}
#endif

#ifndef tk_vec_limited
static inline tk_ivec_t *tk_vec_pfx(cminargs) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, cols, NULL, NULL);
  for (uint64_t c = 0; c < cols; c ++) {
    uint64_t minr = 0;
    tk_vec_base minv = m0->a[0 * cols + c];
    for (size_t r = 1; r < m0->n / cols; r ++) {
      if (m0->a[r * cols + c] < minv) {
        minr = r;
        minv = m0->a[r * cols + c];
      }
    }
    out->a[c] = (int64_t) minr;
  }
  return out;
}
#endif

#ifndef tk_vec_limited
static inline tk_ivec_t *tk_vec_pfx(rminargs) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, m0->n / cols, NULL, NULL);
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val;
    }
    out->a[r] = sum;
    uint64_t minc = 0;
    tk_vec_base minv = m0->a[r * cols + 0];
    for (uint64_t c = 1; c < cols; c ++) {
      if (m0->a[r * cols + c] < minv) {
        minc = c;
        minv = m0->a[r * cols + c];
      }
    }
    out->a[r] = (int64_t) minc;
  }
  return out;
}
#endif

#ifdef tk_vec_lua

#ifndef tk_vec_limited
static inline int tk_vec_pfx(rmagnitudes_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1);
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rmagnitudes)(L, m0, cols);
  return 1;
}
#endif

#ifndef tk_vec_limited
static inline int tk_vec_pfx(cmagnitudes_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1);
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cmagnitudes)(L, m0, cols);
  return 1;
}
#endif

#ifndef tk_vec_limited
static inline int tk_vec_pfx(rminargs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1);
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rminargs)(L, m0, cols);
  return 1;
}
#endif

#ifndef tk_vec_limited
static inline int tk_vec_pfx(cminargs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1);
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cminargs)(L, m0, cols);
  return 1;
}
#endif

#ifndef tk_vec_limited
static inline int tk_vec_pfx(rmaxargs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1);
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rmaxargs)(L, m0, cols);
  return 1;
}
#endif

#ifndef tk_vec_limited
static inline int tk_vec_pfx(cmaxargs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1);
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cmaxargs)(L, m0, cols);
  return 1;
}
#endif

#ifndef tk_vec_limited
static inline int tk_vec_pfx(rasc_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1);
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rasc)(L, m0, cols);
  return 1;
}
#endif

#ifndef tk_vec_limited
static inline int tk_vec_pfx(rdesc_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1);
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rdesc)(L, m0, cols);
  return 1;
}
#endif

#ifndef tk_vec_limited
static inline int tk_vec_pfx(casc_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1);
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(casc)(L, m0, cols);
  return 1;
}
#endif

#ifndef tk_vec_limited
static inline int tk_vec_pfx(cdesc_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1);
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cdesc)(L, m0, cols);
  return 1;
}
#endif

static luaL_Reg tk_vec_pfx(lua_fns)[] =
{
#ifndef tk_vec_limited
  { "cmagnitudes", tk_vec_pfx(cmagnitudes_lua) },
  { "rmagnitudes", tk_vec_pfx(rmagnitudes_lua) },
#endif

#ifndef tk_vec_limited
  { "cmaxargs", tk_vec_pfx(rmaxargs_lua) },
  { "rmaxargs", tk_vec_pfx(cmaxargs_lua) },
#endif

#ifndef tk_vec_limited
  { "cminargs", tk_vec_pfx(rminargs_lua) },
  { "rminargs", tk_vec_pfx(cminargs_lua) },
#endif

#ifndef tk_vec_limited
  { "rasc", tk_vec_pfx(rasc_lua) },
  { "rdesc", tk_vec_pfx(rdesc_lua) },
#endif

#ifndef tk_vec_limited
  { "casc", tk_vec_pfx(casc_lua) },
  { "cdesc", tk_vec_pfx(cdesc_lua) },
#endif

  { NULL, NULL }
};

#endif

#include <santoku/vec.undef.h>
