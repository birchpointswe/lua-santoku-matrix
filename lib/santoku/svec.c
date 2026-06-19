#include <santoku/iuset.h>
#include <santoku/svec.h>
#include <santoku/ivec.h>
#include <santoku/cvec.h>
#include <string.h>

static inline int tk_svec_set_jaccard_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_svec_t *a = tk_svec_peek(L, 1, "a");
  tk_svec_t *b = tk_svec_peek(L, 2, "b");
  tk_dvec_t *weights = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_dvec_peek(L, 3, "weights") : NULL;
  double inter_w, sum_a, sum_b;
  tk_svec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  double sim = tk_svec_set_jaccard(inter_w, sum_a, sum_b);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_svec_set_overlap_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_svec_t *a = tk_svec_peek(L, 1, "a");
  tk_svec_t *b = tk_svec_peek(L, 2, "b");
  tk_dvec_t *weights = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_dvec_peek(L, 3, "weights") : NULL;
  double inter_w, sum_a, sum_b;
  tk_svec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  double sim = tk_svec_set_overlap(inter_w, sum_a, sum_b);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_svec_set_dice_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_svec_t *a = tk_svec_peek(L, 1, "a");
  tk_svec_t *b = tk_svec_peek(L, 2, "b");
  tk_dvec_t *weights = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_dvec_peek(L, 3, "weights") : NULL;
  double inter_w, sum_a, sum_b;
  tk_svec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  double sim = tk_svec_set_dice(inter_w, sum_a, sum_b);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_svec_set_tversky_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_svec_t *a = tk_svec_peek(L, 1, "a");
  tk_svec_t *b = tk_svec_peek(L, 2, "b");

  tk_dvec_t *weights = NULL;
  double alpha = 0.5;
  double beta = 0.5;

  if (nargs >= 3) {
    if (lua_isnumber(L, 3)) {
      alpha = luaL_checknumber(L, 3);
      beta = (nargs >= 4) ? luaL_checknumber(L, 4) : 0.5;
    } else if (!lua_isnil(L, 3)) {
      weights = tk_dvec_peek(L, 3, "weights");
      alpha = (nargs >= 4) ? luaL_checknumber(L, 4) : 0.5;
      beta = (nargs >= 5) ? luaL_checknumber(L, 5) : 0.5;
    }
  }

  double inter_w, sum_a, sum_b;
  tk_svec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  double sim = tk_svec_set_tversky(inter_w, sum_a, sum_b, alpha, beta);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_svec_set_find_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_svec_t *vec = tk_svec_peek(L, 1, "vector");
  int32_t value = (int32_t)luaL_checkinteger(L, 2);
  int32_t start = (nargs >= 3) ? (int32_t)luaL_checkinteger(L, 3) : 0;
  int32_t end = (nargs >= 4) ? (int32_t)luaL_checkinteger(L, 4) : (int32_t)vec->n;
  if (start < 0) start = 0;
  if (end > (int32_t)vec->n) end = (int32_t)vec->n;
  if (start > end) start = end;

  int32_t result = tk_svec_set_find(vec->a, start, end, value);
  lua_pushinteger(L, result);
  return 1;
}

static inline int tk_svec_set_insert_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_svec_t *vec = tk_svec_peek(L, 1, "vector");
  int32_t pos = (int32_t)luaL_checkinteger(L, 2);
  int32_t value = (int32_t)luaL_checkinteger(L, 3);
  tk_svec_set_insert(vec, pos, value);
  return 0;
}

static inline int tk_svec_set_intersect_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_svec_t *a = tk_svec_peek(L, 1, "a");
  tk_svec_t *b = tk_svec_peek(L, 2, "b");
  tk_svec_t *out = NULL;
  if (!lua_isnil(L, 3))
    out = tk_svec_peek(L, 3, "output");
  tk_svec_set_intersect(L, a, b, out);
  return out == NULL ? 1 : 0;
}

static inline int tk_svec_set_union_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_svec_t *a = tk_svec_peek(L, 1, "a");
  tk_svec_t *b = tk_svec_peek(L, 2, "b");
  tk_svec_t *out = NULL;
  if (!lua_isnil(L, 3))
    out = tk_svec_peek(L, 3, "output");
  tk_svec_set_union(L, a, b, out);
  return out == NULL ? 1 : 0;
}

static inline int tk_svec_lookup_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_svec_t *indices = tk_svec_peek(L, 1, "indices");
  tk_svec_t *source = tk_svec_peek(L, 2, "source");
  tk_svec_lookup(indices, source);
  return 0;
}

static inline int tk_svec_to_ivec_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_svec_t *v = tk_svec_peek(L, 1, "svec");
  tk_svec_to_ivec(L, v);
  return 1;
}

static inline int tk_svec_to_dvec_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_svec_t *v = tk_svec_peek(L, 1, "svec");
  tk_svec_to_dvec(L, v);
  return 1;
}

static inline int tk_svec_index_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_svec_t *v = tk_svec_peek(L, 1, "vector");
  tk_iumap_t *m = tk_iumap_create(L, v->n);
  if (!m)
    return luaL_error(L, "index: allocation failed");
  int kha;
  uint32_t khi;
  for (int64_t i = 0; i < (int64_t) v->n; i ++) {
    khi = tk_iumap_put(m, (int64_t) v->a[i], &kha);
    if (kha < 0)
      return luaL_error(L, "index: allocation failed");
    tk_iumap_setval(m, khi, i);
  }
  return 1;
}

static luaL_Reg tk_svec_lua_mt_ext2_fns[] =
{
  { "set_jaccard", tk_svec_set_jaccard_lua },
  { "set_overlap", tk_svec_set_overlap_lua },
  { "set_dice", tk_svec_set_dice_lua },
  { "set_tversky", tk_svec_set_tversky_lua },
  { "set_find", tk_svec_set_find_lua },
  { "set_insert", tk_svec_set_insert_lua },
  { "set_intersect", tk_svec_set_intersect_lua },
  { "set_union", tk_svec_set_union_lua },
  { "lookup", tk_svec_lookup_lua },
  { "index", tk_svec_index_lua },
  { "to_ivec", tk_svec_to_ivec_lua },
  { "to_dvec", tk_svec_to_dvec_lua },
  { NULL, NULL }
};

static luaL_Reg tk_svec_lua_ext_fns[] =
{
  { NULL, NULL }
};

int luaopen_santoku_svec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_svec_lua_fns);
  luaL_register(L, NULL, tk_svec_lua_ext_fns);
  tk_svec_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_svec_lua_mt_fns);
  luaL_register(L, NULL, tk_svec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_svec_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
