#if defined(_OPENMP) && !defined(__EMSCRIPTEN__)
#include <omp.h>
#endif

#define tk_vec_pfx(name) tk_pp_strcat(tk_vec_name, name)


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
  bool tostr = lua_type(L, 2) == LUA_TNIL;
  FILE *fh;
  if (tostr)
    fh = tk_lua_tmpfile(L);
  else
    fh = tk_lua_fopen(L, luaL_checkstring(L, 2), "w");
  tk_vec_pfx(persist)(L, m0, fh);
  if (!tostr) {
    tk_lua_fclose(L, fh);
    return 0;
  } else {
    size_t len;
    char *data = tk_lua_fslurp(L, fh, &len);
    if (data) {
      tk_cvec_t *cv = tk_cvec_create(L, len);
      memcpy(cv->a, data, len);
      free(data);
      return 1;
    } else {
      tk_lua_fclose(L, fh);
      return 0;
    }
  }
}

static luaL_Reg tk_vec_pfx(lua_mt_ext_fns)[] =
{
  { "copy", tk_vec_pfx(copy_indexed_lua) },
  { "persist", tk_vec_pfx(persist_lua) },
#ifndef tk_vec_limited
  { "cmagnitudes", tk_vec_pfx(cmagnitudes_lua) },
  { "rmagnitudes", tk_vec_pfx(rmagnitudes_lua) },
  { "cmaxargs", tk_vec_pfx(cmaxargs_lua) },
  { "rmaxargs", tk_vec_pfx(rmaxargs_lua) },
  { "cminargs", tk_vec_pfx(cminargs_lua) },
  { "rminargs", tk_vec_pfx(rminargs_lua) },
  { "rasc", tk_vec_pfx(rasc_lua) },
  { "rdesc", tk_vec_pfx(rdesc_lua) },
  { "casc", tk_vec_pfx(casc_lua) },
  { "cdesc", tk_vec_pfx(cdesc_lua) },
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
