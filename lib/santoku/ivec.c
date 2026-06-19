#include <santoku/iuset.h>
#include <santoku/ivec.h>
#include <santoku/svec.h>
#include <santoku/cvec.h>
#include <string.h>

static inline int tk_ivec_set_jaccard_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_dvec_t *weights = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_dvec_peek(L, 3, "weights") : NULL;
  double inter_w, sum_a, sum_b;
  tk_ivec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  double sim = tk_ivec_set_jaccard(inter_w, sum_a, sum_b);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_ivec_set_overlap_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_dvec_t *weights = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_dvec_peek(L, 3, "weights") : NULL;
  double inter_w, sum_a, sum_b;
  tk_ivec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  double sim = tk_ivec_set_overlap(inter_w, sum_a, sum_b);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_ivec_set_dice_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_dvec_t *weights = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_dvec_peek(L, 3, "weights") : NULL;
  double inter_w, sum_a, sum_b;
  tk_ivec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  double sim = tk_ivec_set_dice(inter_w, sum_a, sum_b);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_ivec_set_tversky_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");

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
  tk_ivec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  double sim = tk_ivec_set_tversky(inter_w, sum_a, sum_b, alpha, beta);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_ivec_set_find_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *vec = tk_ivec_peek(L, 1, "vector");
  int64_t value = luaL_checkinteger(L, 2);
  int64_t start = (nargs >= 3) ? luaL_checkinteger(L, 3) : 0;
  int64_t end = (nargs >= 4) ? luaL_checkinteger(L, 4) : (int64_t)vec->n;
  if (start < 0) start = 0;
  if (end > (int64_t)vec->n) end = (int64_t)vec->n;
  if (start > end) start = end;

  int64_t result = tk_ivec_set_find(vec->a, start, end, value);
  lua_pushinteger(L, result);
  return 1;
}

static inline int tk_ivec_set_insert_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_ivec_t *vec = tk_ivec_peek(L, 1, "vector");
  int64_t pos = luaL_checkinteger(L, 2);
  int64_t value = luaL_checkinteger(L, 3);
  tk_ivec_set_insert(vec, pos, value);
  return 0;
}

static inline int tk_ivec_set_intersect_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_ivec_t *out = NULL;
  if (!lua_isnil(L, 3))
    out = tk_ivec_peek(L, 3, "output");
  tk_ivec_set_intersect(L, a, b, out);
  return out == NULL ? 1 : 0;
}

static inline int tk_ivec_set_union_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_ivec_t *out = NULL;
  if (!lua_isnil(L, 3))
    out = tk_ivec_peek(L, 3, "output");
  tk_ivec_set_union(L, a, b, out);
  return out == NULL ? 1 : 0;
}

static inline int tk_ivec_lookup_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_ivec_t *indices = tk_ivec_peek(L, 1, "indices");
  tk_ivec_t *source = tk_ivec_peek(L, 2, "source");
  tk_ivec_lookup(indices, source);
  return 0;
}

static inline int tk_ivec_index_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_ivec_t *v = tk_ivec_peek(L, 1, "vector");
  tk_iumap_from_ivec(L, v);
  return 1;
}

static inline int tk_ivec_bincount_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_ivec_t *v = tk_ivec_peek(L, 1, "vector");
  uint64_t n_bins = tk_lua_checkunsigned(L, 2, "n_bins");
  tk_dvec_t *w = t >= 3 && !lua_isnil(L, 3) ? tk_dvec_peek(L, 3, "weights") : NULL;
  if (w != NULL && w->n != v->n)
    return tk_lua_verror(L, 2, "bincount", "weights length differs from vector length");
  if (w != NULL) {
    tk_dvec_t *out = tk_dvec_create(L, n_bins);
    for (uint64_t b = 0; b < n_bins; b ++)
      out->a[b] = 0.0;
    for (uint64_t i = 0; i < v->n; i ++) {
      int64_t b = v->a[i];
      if (b < 0 || (uint64_t) b >= n_bins)
        return tk_lua_verror(L, 2, "bincount", "value out of range");
      out->a[b] += w->a[i];
    }
  } else {
    tk_ivec_t *out = tk_ivec_create(L, n_bins);
    for (uint64_t b = 0; b < n_bins; b ++)
      out->a[b] = 0;
    for (uint64_t i = 0; i < v->n; i ++) {
      int64_t b = v->a[i];
      if (b < 0 || (uint64_t) b >= n_bins)
        return tk_lua_verror(L, 2, "bincount", "value out of range");
      out->a[b] ++;
    }
  }
  return 1;
}

static inline int tk_ivec_to_dvec_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_ivec_t *v = tk_ivec_peek(L, 1, "ivec");
  tk_ivec_to_dvec(L, v);
  return 1;
}

static inline int tk_ivec_to_svec_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_ivec_t *v = tk_ivec_peek(L, 1, "ivec");
  tk_ivec_to_svec(L, v);
  return 1;
}

static luaL_Reg tk_ivec_lua_mt_ext2_fns[] =
{
  { "set_jaccard", tk_ivec_set_jaccard_lua },
  { "set_overlap", tk_ivec_set_overlap_lua },
  { "set_dice", tk_ivec_set_dice_lua },
  { "set_tversky", tk_ivec_set_tversky_lua },
  { "set_find", tk_ivec_set_find_lua },
  { "set_insert", tk_ivec_set_insert_lua },
  { "set_intersect", tk_ivec_set_intersect_lua },
  { "set_union", tk_ivec_set_union_lua },
  { "lookup", tk_ivec_lookup_lua },
  { "index", tk_ivec_index_lua },
  { "bincount", tk_ivec_bincount_lua },
  { "to_dvec", tk_ivec_to_dvec_lua },
  { "to_svec", tk_ivec_to_svec_lua },
  { NULL, NULL }
};

static luaL_Reg tk_ivec_lua_ext_fns[] =
{
  { NULL, NULL }
};

int luaopen_santoku_ivec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_ivec_lua_fns);
  luaL_register(L, NULL, tk_ivec_lua_ext_fns);
  tk_ivec_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_ivec_lua_mt_fns);
  luaL_register(L, NULL, tk_ivec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_ivec_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
