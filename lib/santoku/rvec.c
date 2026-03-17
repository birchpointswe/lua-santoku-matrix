#define TK_RVEC_INIT
#include <santoku/iuset.h>
#include <santoku/rvec.h>
#include <santoku/ivec/ext.h>

static inline int tk_rvec_split_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  tk_ivec_t *K = tk_ivec_peek(L, 2, "keys");
  tk_dvec_t *V = tk_dvec_peek(L, 3, "values");
  bool append = tk_lua_optboolean(L, 4, "append", true);
  if (tk_rvec_split(P, K, V, append) != 0)
    tk_lua_verror(L, 2, "split", "allocation error");
  return 0;
}

static inline int tk_rvec_get_lua (lua_State *L) {
  lua_settop(L, 2);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  uint64_t i = tk_lua_checkunsigned(L, 2, "idx");
  if (i > P->n)
    return 0;
  lua_pushinteger(L, P->a[i].i);
  lua_pushnumber(L, P->a[i].d);
  return 2;
}

static inline int tk_rvec_set_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  uint64_t i = tk_lua_checkunsigned(L, 2, "idx");
  int64_t a = tk_lua_checkinteger(L, 3, "i");
  double d = tk_lua_checkdouble(L, 4, "d");
  tk_rvec_ensure(P, i + 1);
  P->a[i] = (tk_rank_t) { a, d };
  return 0;
}

static inline int tk_rvec_push_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  int64_t a = tk_lua_checkinteger(L, 2, "i");
  double d = tk_lua_checkdouble(L, 3, "d");
  tk_rvec_push(P, (tk_rank_t) { a, d });
  return 0;
}

static inline int tk_rvec_hmax_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  int64_t a = tk_lua_checkinteger(L, 2, "i");
  double d = tk_lua_checkdouble(L, 3, "d");
  uint64_t m = tk_lua_checkunsigned(L, 4, "max");
  tk_rvec_hmax(P, m, (tk_rank_t) { a, d });
  return 0;
}

static inline int tk_rvec_hmax_init_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  tk_rvec_hmax_init(P);
  return 0;
}

static inline int tk_rvec_hmax_push_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  int64_t a = tk_lua_checkinteger(L, 2, "i");
  double d = tk_lua_checkdouble(L, 3, "d");
  tk_rvec_hmax_push(P, (tk_rank_t) { a, d });
  return 0;
}

static inline int tk_rvec_hmax_pop_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  if (!P->n)
    return 0;
  tk_rank_t r = tk_rvec_hmax_pop(P);
  lua_pushinteger(L, r.i);
  lua_pushnumber(L, r.d);
  return 2;
}

static inline int tk_rvec_hmin_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  int64_t a = tk_lua_checkinteger(L, 2, "i");
  double d = tk_lua_checkdouble(L, 3, "d");
  uint64_t m = tk_lua_checkunsigned(L, 4, "max");
  tk_rvec_hmin(P, m, (tk_rank_t) { a, d });
  return 0;
}

static inline int tk_rvec_hmin_init_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  tk_rvec_hmin_init(P);
  return 0;
}

static inline int tk_rvec_hmin_push_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  int64_t a = tk_lua_checkinteger(L, 2, "i");
  double d = tk_lua_checkdouble(L, 3, "d");
  tk_rvec_hmin_push(P, (tk_rank_t) { a, d });
  return 0;
}

static inline int tk_rvec_hmin_pop_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  if (!P->n)
    return 0;
  tk_rank_t r = tk_rvec_hmin_pop(P);
  lua_pushinteger(L, r.i);
  lua_pushnumber(L, r.d);
  return 2;
}

static inline int tk_rvec_keys_lua (lua_State *L) {
  lua_settop(L, 2);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  tk_ivec_t *out = lua_isnil(L, 2) ? NULL : tk_ivec_peek(L, 2, "out");
  tk_rvec_keys(L, P, out);
  if (out)
    return 0;
  return 1;
}

static inline int tk_rvec_values_lua (lua_State *L) {
  lua_settop(L, 2);
  tk_rvec_t *P = tk_rvec_peek(L, 1, "rvec");
  tk_dvec_t *out = lua_isnil(L, 2) ? NULL : tk_dvec_peek(L, 2, "out");
  tk_rvec_values(L, P, out);
  if (out)
    return 0;
  return 1;
}

static inline int tk_rvec_scores_elbow_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_rvec_t *scores = tk_rvec_peek(L, 1, "rvec");
  const char *method = luaL_checkstring(L, 2);
  double alpha = luaL_optnumber(L, 3, 0.0);
  double val;
  size_t idx;
  if (strcmp(method, "lmethod") == 0) {
    idx = tk_rvec_scores_lmethod(scores, &val);
  } else if (strcmp(method, "max_gap") == 0) {
    idx = tk_rvec_scores_max_gap(scores, &val);
  } else if (strcmp(method, "max_curvature") == 0) {
    idx = tk_rvec_scores_max_curvature(scores, &val);
  } else if (strcmp(method, "kneedle") == 0) {
    double sensitivity = (alpha > 0.0) ? alpha : 1.0;
    idx = tk_rvec_scores_kneedle(scores, sensitivity, &val);
  } else if (strcmp(method, "plateau") == 0) {
    double tolerance = (alpha > 0.0) ? alpha : 1e-3;
    idx = tk_rvec_scores_plateau(scores, tolerance, &val);
  } else if (strcmp(method, "otsu") == 0) {
    idx = tk_rvec_scores_otsu(scores, &val);
  } else if (strcmp(method, "first_gap") == 0) {
    double ratio = (alpha > 0.0) ? alpha : 3.0;
    idx = tk_rvec_scores_first_gap(scores, ratio, &val);
  } else {
    return luaL_error(L, "unknown elbow method: %s", method);
  }
  lua_pushnumber(L, (lua_Number)val);
  lua_pushinteger(L, (lua_Integer)(idx + 1));
  return 2;
}

static inline int tk_rvec_rankings_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "scores");
  uint64_t n_visible = tk_lua_checkunsigned(L, 2, "n_visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 3, "n_hidden");
  tk_rvec_rankings(L, scores, n_visible, n_hidden);
  return 1;
}

static luaL_Reg tk_rvec_lua_ext_fns[] =
{
  { "rankings", tk_rvec_rankings_lua },
  { NULL, NULL }
};

static luaL_Reg tk_rvec_lua_mt_ext2_fns[] =
{
  { "get", tk_rvec_get_lua },
  { "set", tk_rvec_set_lua },
  { "push", tk_rvec_push_lua },
  { "keys", tk_rvec_keys_lua },
  { "split", tk_rvec_split_lua },
  { "values", tk_rvec_values_lua },
  { "each", tk_rvec_each_lua },
  { "ieach", tk_rvec_ieach_lua },
  { "hmax", tk_rvec_hmax_lua },
  { "hmax_init", tk_rvec_hmax_init_lua },
  { "hmax_pop", tk_rvec_hmax_pop_lua },
  { "hmax_push", tk_rvec_hmax_push_lua },
  { "hmin", tk_rvec_hmin_lua },
  { "hmin_init", tk_rvec_hmin_init_lua },
  { "hmin_pop", tk_rvec_hmin_pop_lua },
  { "hmin_push", tk_rvec_hmin_push_lua },
  { "scores_elbow", tk_rvec_scores_elbow_lua },
  { NULL, NULL }
};

void tk_rvec_init_mt (lua_State *L)
{
  luaL_register(L, NULL, tk_rvec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_rvec_lua_mt_ext2_fns);
}

int luaopen_santoku_rvec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_rvec_lua_fns);
  luaL_register(L, NULL, tk_rvec_lua_ext_fns);
  tk_rvec_create(L, 0, 0, 0);
  lua_pop(L, 1);
  return 1;
}
