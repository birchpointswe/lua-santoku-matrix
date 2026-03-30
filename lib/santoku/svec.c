#include <santoku/iuset.h>
#include <santoku/svec.h>
#include <santoku/ivec.h>
#include <santoku/cvec.h>
#include <string.h>

static inline int tk_svec_copy_pkeys_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_svec_t *m0 = tk_svec_peek(L, 1, "dest");
  tk_pvec_t *m1 = tk_pvec_peek(L, 2, "source");
  int64_t start, end, dest;
  if (t == 2) {
    start = 0;
    end = (int64_t) m1->n;
    dest = (int64_t) m0->n;
  } else if (t == 3) {
    start = 0;
    end = (int64_t) m1->n;
    dest = tk_lua_checkinteger(L, 3, "dest");
  } else {
    start = tk_lua_checkinteger(L, 3, "start");
    end = tk_lua_checkinteger(L, 4, "end");
    dest = tk_lua_checkinteger(L, 5, "dest");
  }
  tk_svec_copy_pkeys(m0, m1, start, end, dest);
  return 0;
}

static inline int tk_svec_copy_rkeys_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_svec_t *m0 = tk_svec_peek(L, 1, "dest");
  tk_rvec_t *m1 = tk_rvec_peek(L, 2, "source");
  int64_t start, end, dest;
  if (t == 2) {
    start = 0;
    end = (int64_t) m1->n;
    dest = (int64_t) m0->n;
  } else if (t == 3) {
    start = 0;
    end = (int64_t) m1->n;
    dest = tk_lua_checkinteger(L, 3, "dest");
  } else {
    start = tk_lua_checkinteger(L, 3, "start");
    end = tk_lua_checkinteger(L, 4, "end");
    dest = tk_lua_checkinteger(L, 5, "dest");
  }
  tk_svec_copy_rkeys(m0, m1, start, end, dest);
  return 0;
}

static inline int tk_svec_copy_pvalues_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_svec_t *m0 = tk_svec_peek(L, 1, "dest");
  tk_pvec_t *m1 = tk_pvec_peek(L, 2, "source");
  int64_t start, end, dest;
  if (t == 2) {
    start = 0;
    end = (int64_t) m1->n;
    dest = (int64_t) m0->n;
  } else if (t == 3) {
    start = 0;
    end = (int64_t) m1->n;
    dest = tk_lua_checkinteger(L, 3, "dest");
  } else {
    start = tk_lua_checkinteger(L, 3, "start");
    end = tk_lua_checkinteger(L, 4, "end");
    dest = tk_lua_checkinteger(L, 5, "dest");
  }
  tk_svec_copy_pvalues(m0, m1, start, end, dest);
  return 0;
}

static inline int tk_svec_copy_rvalues_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_svec_t *m0 = tk_svec_peek(L, 1, "dest");
  tk_rvec_t *m1 = tk_rvec_peek(L, 2, "source");
  int64_t start, end, dest;
  if (t == 2) {
    start = 0;
    end = (int64_t) m1->n;
    dest = (int64_t) m0->n;
  } else if (t == 3) {
    start = 0;
    end = (int64_t) m1->n;
    dest = tk_lua_checkinteger(L, 3, "dest");
  } else {
    start = tk_lua_checkinteger(L, 3, "start");
    end = tk_lua_checkinteger(L, 4, "end");
    dest = tk_lua_checkinteger(L, 5, "dest");
  }
  tk_svec_copy_rvalues(m0, m1, start, end, dest);
  return 0;
}

static inline int tk_svec_set_stats_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_svec_t *a = tk_svec_peek(L, 1, "a");
  tk_svec_t *b = tk_svec_peek(L, 2, "b");
  tk_dvec_t *weights = lua_isnil(L, 3) ? NULL : tk_dvec_peek(L, 3, "weights");
  double inter_w, sum_a, sum_b;
  tk_svec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  lua_pushnumber(L, inter_w);
  lua_pushnumber(L, sum_a);
  lua_pushnumber(L, sum_b);
  return 3;
}

static inline int tk_svec_set_similarity_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_svec_t *a = tk_svec_peek(L, 1, "a");
  tk_svec_t *b = tk_svec_peek(L, 2, "b");
  tk_dvec_t *weights = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_dvec_peek(L, 3, "weights") : NULL;
  tk_svec_sim_type_t type = TK_SVEC_JACCARD;
  if (nargs >= 4 && lua_isstring(L, 4)) {
    const char *type_str = lua_tostring(L, 4);
    if (strcmp(type_str, "jaccard") == 0) {
      type = TK_SVEC_JACCARD;
    } else if (strcmp(type_str, "overlap") == 0) {
      type = TK_SVEC_OVERLAP;
    } else if (strcmp(type_str, "dice") == 0) {
      type = TK_SVEC_DICE;
    } else if (strcmp(type_str, "tversky") == 0) {
      type = TK_SVEC_TVERSKY;
    }
  }
  double tversky_alpha = (nargs >= 5) ? luaL_checknumber(L, 5) : 0.5;
  double tversky_beta = (nargs >= 6) ? luaL_checknumber(L, 6) : 0.5;
  double sim = tk_svec_set_similarity(a->a, a->n, b->a, b->n, weights, type, tversky_alpha, tversky_beta);
  lua_pushnumber(L, sim);
  return 1;
}

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

static inline int tk_svec_set_weights_by_rank_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_svec_t *features = tk_svec_peek(L, 1, "features");
  tk_dvec_t *weights = (nargs >= 2 && !lua_isnil(L, 2)) ? tk_dvec_peek(L, 2, "weights") : NULL;
  tk_svec_t *ranks = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_svec_peek(L, 3, "ranks") : NULL;
  uint64_t n_ranks = (nargs >= 4) ? tk_lua_checkunsigned(L, 4, "n_ranks") : 1;
  tk_dvec_t *weights_by_rank = tk_dvec_create(L, n_ranks);
  tk_svec_set_weights_by_rank(features->a, features->n, weights, ranks, n_ranks, weights_by_rank->a);
  weights_by_rank->n = n_ranks;
  return 1;
}

static inline int tk_svec_set_similarity_by_rank_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_dvec_t *wacc = tk_dvec_peek(L, 1, "wacc");
  int64_t vsid = tk_lua_checkinteger(L, 2, "vsid");
  tk_dvec_t *q_weights_by_rank = tk_dvec_peek(L, 3, "q_weights_by_rank");
  tk_dvec_t *e_weights_by_rank = tk_dvec_peek(L, 4, "e_weights_by_rank");
  uint64_t n_ranks = tk_lua_checkunsigned(L, 5, "n_ranks");
  int64_t rank_decay_window = (nargs >= 6) ? tk_lua_checkinteger(L, 6, "rank_decay_window") : -1;
  double rank_decay_sigma = (nargs >= 7) ? luaL_checknumber(L, 7) : 10.0;
  double rank_decay_floor = (nargs >= 8) ? luaL_checknumber(L, 8) : 0.01;
  tk_svec_sim_type_t type = TK_SVEC_JACCARD;
  if (nargs >= 9 && lua_isstring(L, 9)) {
    const char *type_str = lua_tostring(L, 9);
    if (strcmp(type_str, "jaccard") == 0) {
      type = TK_SVEC_JACCARD;
    } else if (strcmp(type_str, "overlap") == 0) {
      type = TK_SVEC_OVERLAP;
    } else if (strcmp(type_str, "dice") == 0) {
      type = TK_SVEC_DICE;
    } else if (strcmp(type_str, "tversky") == 0) {
      type = TK_SVEC_TVERSKY;
    }
  }
  double tversky_alpha = (nargs >= 10) ? luaL_checknumber(L, 10) : 0.5;
  double tversky_beta = (nargs >= 11) ? luaL_checknumber(L, 11) : 0.5;
  double sim = tk_svec_set_similarity_by_rank(
    wacc, vsid, q_weights_by_rank->a,
    e_weights_by_rank->a, n_ranks,
    rank_decay_window,
    rank_decay_sigma,
    rank_decay_floor, type,
    tversky_alpha,
    tversky_beta);
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

static inline int tk_svec_scores_elbow_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_svec_t *scores = tk_svec_peek(L, 1, "svec");
  const char *method = luaL_checkstring(L, 2);
  int32_t alpha = (int32_t)luaL_optnumber(L, 3, 0.0);
  int32_t val;
  size_t idx;
  if (strcmp(method, "lmethod") == 0) {
    idx = tk_svec_scores_lmethod(scores, &val);
  } else if (strcmp(method, "max_gap") == 0) {
    idx = tk_svec_scores_max_gap(scores, &val);
  } else if (strcmp(method, "max_curvature") == 0) {
    idx = tk_svec_scores_max_curvature(scores, &val);
  } else if (strcmp(method, "kneedle") == 0) {
    double sensitivity = (alpha > 0) ? (double)alpha : 1.0;
    idx = tk_svec_scores_kneedle(scores, sensitivity, &val);
  } else if (strcmp(method, "plateau") == 0) {
    int32_t tolerance = (alpha > 0) ? alpha : 0;
    idx = tk_svec_scores_plateau(scores, tolerance, &val);
  } else if (strcmp(method, "otsu") == 0) {
    idx = tk_svec_scores_otsu(scores, &val);
  } else if (strcmp(method, "first_gap") == 0) {
    double ratio = (alpha > 0) ? (double)alpha : 3.0;
    idx = tk_svec_scores_first_gap(scores, ratio, &val);
  } else {
    return luaL_error(L, "unknown elbow method: %s", method);
  }
  lua_pushinteger(L, (lua_Integer)val);
  lua_pushinteger(L, (lua_Integer)(idx + 1));
  return 2;
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

static luaL_Reg tk_svec_lua_mt_ext2_fns[] =
{
  { "copy_pkeys", tk_svec_copy_pkeys_lua },
  { "copy_rkeys", tk_svec_copy_rkeys_lua },
  { "copy_pvalues", tk_svec_copy_pvalues_lua },
  { "copy_rvalues", tk_svec_copy_rvalues_lua },
  { "set_stats", tk_svec_set_stats_lua },
  { "set_similarity", tk_svec_set_similarity_lua },
  { "set_jaccard", tk_svec_set_jaccard_lua },
  { "set_overlap", tk_svec_set_overlap_lua },
  { "set_dice", tk_svec_set_dice_lua },
  { "set_tversky", tk_svec_set_tversky_lua },
  { "set_weights_by_rank", tk_svec_set_weights_by_rank_lua },
  { "set_similarity_by_rank", tk_svec_set_similarity_by_rank_lua },
  { "set_find", tk_svec_set_find_lua },
  { "set_insert", tk_svec_set_insert_lua },
  { "set_intersect", tk_svec_set_intersect_lua },
  { "set_union", tk_svec_set_union_lua },
  { "lookup", tk_svec_lookup_lua },
  { "scores_elbow", tk_svec_scores_elbow_lua },
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
