#if defined(_OPENMP) && !defined(__EMSCRIPTEN__)
#include <omp.h>
#endif
#include <string.h>

#define tk_vec_pfx(name) tk_pp_strcat(tk_vec_name, name)

#ifndef tk_vec_lt
#define tk_vec_lt(a, b) ((a) < (b))
#endif

#ifndef tk_vec_gt
#define tk_vec_gt(a, b) ((a) > (b))
#endif

#ifndef tk_vec_eqx
#define tk_vec_eqx(a, b) (memcmp(&(a), &(b), sizeof(a)) == 0)
#endif

#ifndef tk_vec_eq
#define tk_vec_eq(a, b) tk_vec_eqx(a, b)
#endif

#ifndef tk_vec_err
#define tk_vec_err(L, name, n, ...) \
  tk_lua_verror((L), ((n) + 1), tk_pp_xstr(tk_vec_pfx(name)), __VA_ARGS__)
#endif


static inline void tk_vec_pfx(copy_indexed) (tk_vec_pfx(t) *m0, tk_vec_pfx(t) *m1, tk_ivec_t *indices);
static inline void tk_vec_pfx(scatter_indexed) (tk_vec_pfx(t) *m0, tk_vec_pfx(t) *m1, tk_ivec_t *indices);
#ifndef tk_vec_limited
static inline tk_ivec_t *tk_vec_pfx(rasc) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_ivec_t *tk_vec_pfx(rdesc) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_ivec_t *tk_vec_pfx(casc) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_ivec_t *tk_vec_pfx(cdesc) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_dvec_t *tk_vec_pfx(cmagnitudes) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_dvec_t *tk_vec_pfx(rmagnitudes) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_ivec_t *tk_vec_pfx(cmaxargs) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_ivec_t *tk_vec_pfx(rmaxargs) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_ivec_t *tk_vec_pfx(cminargs) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_ivec_t *tk_vec_pfx(rminargs) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
#endif


static inline int tk_vec_pfx(copy_indexed_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  if (t >= 3 && lua_type(L, 3) == LUA_TUSERDATA) {
    tk_ivec_t *indices = tk_ivec_peekopt(L, 3);
    if (indices != NULL) {
      tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "dest");
      tk_vec_pfx(t) *m1 = tk_vec_pfx(peek)(L, 2, "source");
      bool scatter = t >= 4 && lua_toboolean(L, 4);
      if (scatter)
        tk_vec_pfx(scatter_indexed)(m0, m1, indices);
      else
        tk_vec_pfx(copy_indexed)(m0, m1, indices);
      lua_settop(L, 1);
      return 1;
    }
  }
  return tk_vec_pfx(copy_lua)(L);
}

#ifndef tk_vec_limited

static inline int tk_vec_pfx(rmagnitudes_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rmagnitudes)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(cmagnitudes_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cmagnitudes)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(rminargs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rminargs)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(cminargs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cminargs)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(rmaxargs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rmaxargs)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(cmaxargs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cmaxargs)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(rasc_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rasc)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(rdesc_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rdesc)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(casc_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(casc)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(cdesc_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cdesc)(L, m0, cols);
  return 1;
}

#endif

static inline bool tk_vec_pfx(eq) (tk_vec_pfx(t) *a, tk_vec_pfx(t) *b, uint64_t start, uint64_t end)
{
  if (end > a->n || end > b->n)
    return false;
  for (uint64_t i = start; i < end; i ++)
    if (!tk_vec_eq(a->a[i], b->a[i]))
      return false;
  return true;
}

static inline int tk_vec_pfx(eq_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *a = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(t) *b = tk_vec_pfx(peek)(L, 2, "other");
  uint64_t start = 0, end = a->n;
  bool full = true;
#ifndef tk_vec_limited
  double eps = -1.0;
  if (t == 3 || t == 5)
    eps = tk_lua_checkdouble(L, t, "eps");
#else
  if (t != 2 && t != 4) {
    tk_vec_err(L, eq, 1, "expected 2 or 4 arguments (vec, other or vec, other, start, end)");
    return 0;
  }
#endif
  if (t >= 4) {
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
    full = false;
  }
  if ((full && a->n != b->n) || end > a->n || end > b->n) {
    lua_pushboolean(L, false);
    return 1;
  }
  bool r;
#ifndef tk_vec_limited
  if (eps >= 0.0) {
    r = true;
    for (uint64_t i = start; i < end; i ++)
      if (fabs((double) a->a[i] - (double) b->a[i]) > eps) {
        r = false;
        break;
      }
  } else {
    r = tk_vec_pfx(eq)(a, b, start, end);
  }
#else
  r = tk_vec_pfx(eq)(a, b, start, end);
#endif
  lua_pushboolean(L, r);
  return 1;
}

#if !defined(tk_vec_limited) && defined(tk_vec_peekbase)

static inline int tk_vec_pfx(where_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *v = tk_vec_pfx(peek)(L, 1, "vector");
  const char *cmp = luaL_checkstring(L, 2);
  tk_vec_base x = tk_vec_peekbase(L, 3);
  uint64_t start = 0, end = v->n;
  if (t >= 5) {
    start = tk_lua_checkunsigned(L, 4, "start");
    end = tk_lua_checkunsigned(L, 5, "end");
  }
  if (end > v->n)
    end = v->n;
  int c;
  if (strcmp(cmp, "lt") == 0) c = 0;
  else if (strcmp(cmp, "le") == 0) c = 1;
  else if (strcmp(cmp, "eq") == 0) c = 2;
  else if (strcmp(cmp, "ne") == 0) c = 3;
  else if (strcmp(cmp, "ge") == 0) c = 4;
  else if (strcmp(cmp, "gt") == 0) c = 5;
  else {
    tk_vec_err(L, where, 1, "cmp must be one of lt, le, eq, ne, ge, gt");
    return 0;
  }
  tk_ivec_t *out = tk_ivec_create(L, 0);
  for (uint64_t i = start; i < end; i ++) {
    tk_vec_base a = v->a[i];
    bool m;
    switch (c) {
      case 0: m = tk_vec_lt(a, x); break;
      case 1: m = !tk_vec_gt(a, x); break;
      case 2: m = tk_vec_eq(a, x); break;
      case 3: m = !tk_vec_eq(a, x); break;
      case 4: m = !tk_vec_lt(a, x); break;
      default: m = tk_vec_gt(a, x); break;
    }
    if (m)
      tk_ivec_push(out, (int64_t) i);
  }
  return 1;
}

static inline void tk_vec_pfx(fill_segments) (tk_vec_pfx(t) *v, tk_ivec_t *offsets, tk_vec_pfx(t) *values)
{
  uint64_t nd = offsets->n > 0 ? offsets->n - 1 : 0;
  for (uint64_t d = 0; d < nd && d < values->n; d ++) {
    int64_t s = offsets->a[d], e = offsets->a[d + 1];
    if (s < 0)
      s = 0;
    if (e > (int64_t) v->n)
      e = (int64_t) v->n;
    for (int64_t i = s; i < e; i ++)
      v->a[i] = values->a[d];
  }
}

static inline int tk_vec_pfx(fill_segments_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *v = tk_vec_pfx(peek)(L, 1, "vector");
  tk_ivec_t *offsets = tk_ivec_peek(L, 2, "offsets");
  tk_vec_pfx(t) *values = tk_vec_pfx(peek)(L, 3, "values");
  tk_vec_pfx(fill_segments)(v, offsets, values);
  lua_settop(L, 1);
  return 1;
}

#endif

static inline void tk_vec_pfx(persist) (lua_State *L, tk_vec_pfx(t) *v, FILE *fh)
{
  uint64_t n64 = (uint64_t) v->n;
  tk_lua_fwrite(L, (char *) &n64, sizeof(uint64_t), 1, fh);
  tk_lua_fwrite(L, (char *) v->a, sizeof(tk_vec_base) * v->n, 1, fh);
}

static inline int tk_vec_pfx(persist_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  FILE *fh = tk_lua_fopen(L, luaL_checkstring(L, 2), "w");
  tk_vec_pfx(persist)(L, m0, fh);
  tk_lua_fclose(L, fh);
  return 0;
}

static luaL_Reg tk_vec_pfx(lua_mt_ext_fns)[] =
{
  { "copy", tk_vec_pfx(copy_indexed_lua) },
  { "persist", tk_vec_pfx(persist_lua) },
  { "eq", tk_vec_pfx(eq_lua) },
#if !defined(tk_vec_limited) && defined(tk_vec_peekbase)
  { "where", tk_vec_pfx(where_lua) },
  { "fill_segments", tk_vec_pfx(fill_segments_lua) },
#endif
#ifndef tk_vec_limited
#endif
  { NULL, NULL }
};

static inline void tk_vec_pfx(suppress_unused_lua_mt_ext_fns) (void)
  { (void) tk_vec_pfx(lua_mt_ext_fns); }



#include <santoku/parallel/tpl.h>
#include <santoku/vec/ext/tpl_para.h>


#define TK_GENERATE_SINGLE
#include <santoku/parallel/tpl.h>
#include <santoku/vec/ext/tpl_para.h>
#undef TK_GENERATE_SINGLE

#include <santoku/vec/undef.h>
