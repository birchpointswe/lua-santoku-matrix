#define TK_PVEC_INIT
#include <santoku/iuset.h>
#include <santoku/pvec.h>

static inline int tk_pvec_get_lua (lua_State *L) {
  lua_settop(L, 2);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  uint64_t i = tk_lua_checkunsigned(L, 2, "idx");
  if (i > P->n)
    return 0;
  lua_pushinteger(L, P->a[i].i);
  lua_pushinteger(L, P->a[i].p);
  return 2;
}

static inline int tk_pvec_set_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  uint64_t i = tk_lua_checkunsigned(L, 2, "idx");
  int64_t a = tk_lua_checkinteger(L, 3, "i");
  int64_t b = tk_lua_checkinteger(L, 4, "p");
  tk_pvec_ensure(P, i + 1);
  P->a[i] = (tk_pair_t) { a, b };
  return 0;
}

static inline int tk_pvec_push_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  int64_t a = tk_lua_checkinteger(L, 2, "i");
  int64_t b = tk_lua_checkinteger(L, 3, "p");
  tk_pvec_push(P, (tk_pair_t) { a, b });
  return 0;
}

static inline int tk_pvec_hmax_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  int64_t a = tk_lua_checkinteger(L, 2, "i");
  int64_t b = tk_lua_checkinteger(L, 3, "p");
  uint64_t m = tk_lua_checkunsigned(L, 4, "max");
  tk_pvec_hmax(P, m, (tk_pair_t) { a, b });
  return 0;
}

static inline int tk_pvec_hmax_init_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  tk_pvec_hmax_init(P);
  return 0;
}

static inline int tk_pvec_hmax_push_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  int64_t a = tk_lua_checkinteger(L, 2, "i");
  int64_t b = tk_lua_checkinteger(L, 3, "p");
  tk_pvec_hmax_push(P, (tk_pair_t) { a, b });
  return 0;
}

static inline int tk_pvec_hmax_pop_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  if (!P->n)
    return 0;
  tk_pair_t r = tk_pvec_hmax_pop(P);
  lua_pushinteger(L, r.i);
  lua_pushinteger(L, r.p);
  return 2;
}

static inline int tk_pvec_hmin_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  int64_t a = tk_lua_checkinteger(L, 2, "i");
  int64_t b = tk_lua_checkinteger(L, 3, "p");
  uint64_t m = tk_lua_checkunsigned(L, 4, "max");
  tk_pvec_hmin(P, m, (tk_pair_t) { a, b });
  return 0;
}

static inline int tk_pvec_hmin_init_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  tk_pvec_hmin_init(P);
  return 0;
}

static inline int tk_pvec_hmin_push_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  int64_t a = tk_lua_checkinteger(L, 2, "i");
  int64_t b = tk_lua_checkinteger(L, 3, "p");
  tk_pvec_hmin_push(P, (tk_pair_t) { a, b });
  return 0;
}

static inline int tk_pvec_hmin_pop_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  if (!P->n)
    return 0;
  tk_pair_t r = tk_pvec_hmin_pop(P);
  lua_pushinteger(L, r.i);
  lua_pushinteger(L, r.p);
  return 2;
}

static inline int tk_pvec_keys_lua (lua_State *L) {
  lua_settop(L, 2);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  tk_ivec_t *out = lua_isnil(L, 2) ? NULL : tk_ivec_peek(L, 2, "out");
  tk_pvec_keys(L, P, out);
  if (out)
    return 0;
  return 1;
}

static inline int tk_pvec_values_lua (lua_State *L) {
  lua_settop(L, 2);
  tk_pvec_t *P = tk_pvec_peek(L, 1, "pvec");
  tk_ivec_t *out = lua_isnil(L, 2) ? NULL : tk_ivec_peek(L, 2, "out");
  tk_pvec_values(L, P, out);
  if (out)
    return 0;
  return 1;
}

static inline int tk_pvec_scores_elbow_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_pvec_t *scores = tk_pvec_peek(L, 1, "pvec");
  const char *method = luaL_checkstring(L, 2);
  int64_t alpha = (int64_t)luaL_optnumber(L, 3, 0.0);
  int64_t val;
  size_t idx;
  if (strcmp(method, "lmethod") == 0) {
    idx = tk_pvec_scores_lmethod(scores, &val);
  } else if (strcmp(method, "max_gap") == 0) {
    idx = tk_pvec_scores_max_gap(scores, &val);
  } else if (strcmp(method, "max_curvature") == 0) {
    idx = tk_pvec_scores_max_curvature(scores, &val);
  } else if (strcmp(method, "kneedle") == 0) {
    double sensitivity = (alpha > 0) ? (double)alpha : 1.0;
    idx = tk_pvec_scores_kneedle(scores, sensitivity, &val);
  } else if (strcmp(method, "plateau") == 0) {
    int64_t tolerance = (alpha > 0) ? alpha : 0;
    idx = tk_pvec_scores_plateau(scores, tolerance, &val);
  } else if (strcmp(method, "otsu") == 0) {
    idx = tk_pvec_scores_otsu(scores, &val);
  } else if (strcmp(method, "first_gap") == 0) {
    double ratio = (alpha > 0) ? (double)alpha : 3.0;
    idx = tk_pvec_scores_first_gap(scores, ratio, &val);
  } else {
    return luaL_error(L, "unknown elbow method: %s", method);
  }
  lua_pushinteger(L, (lua_Integer)val);
  lua_pushinteger(L, (lua_Integer)(idx + 1));
  return 2;
}

static luaL_Reg tk_pvec_lua_mt_ext2_fns[] =
{
  { "get", tk_pvec_get_lua },
  { "set", tk_pvec_set_lua },
  { "push", tk_pvec_push_lua },
  { "keys", tk_pvec_keys_lua },
  { "values", tk_pvec_values_lua },
  { "each", tk_pvec_each_lua },
  { "ieach", tk_pvec_ieach_lua },
  { "hmax", tk_pvec_hmax_lua },
  { "hmax_init", tk_pvec_hmax_init_lua },
  { "hmax_pop", tk_pvec_hmax_pop_lua },
  { "hmax_push", tk_pvec_hmax_push_lua },
  { "hmin", tk_pvec_hmin_lua },
  { "hmin_init", tk_pvec_hmin_init_lua },
  { "hmin_pop", tk_pvec_hmin_pop_lua },
  { "hmin_push", tk_pvec_hmin_push_lua },
  { "scores_elbow", tk_pvec_scores_elbow_lua },
  { NULL, NULL }
};

void tk_pvec_init_mt (lua_State *L)
{
  luaL_register(L, NULL, tk_pvec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_pvec_lua_mt_ext2_fns);
}

int luaopen_santoku_pvec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_pvec_lua_fns);
  tk_pvec_create(L, 0, 0, 0);
  lua_pop(L, 1);
  return 1;
}
