#include <santoku/iuset.h>
#include <santoku/dvec.h>

static inline int tk_dvec_multiply_bits_lua (lua_State *L)
{
  lua_settop(L, 5);
  tk_dvec_t *P = tk_dvec_peek(L, 1, "dvec");
  tk_ivec_t *raw_features = tk_ivec_peek(L, 2, "bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  tk_dvec_multiply_bits(L, P, raw_features, n_samples, n_features, n_hidden);
  return 1;
}

static inline int tk_dvec_center_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_dvec_t *M = tk_dvec_peek(L, 1, "dvec");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_dims = tk_lua_checkunsigned(L, 3, "n_dims");
  tk_dvec_center(M->a, n_samples, n_dims);
  lua_pushvalue(L, 1);
  return 1;
}

static inline int tk_dvec_rnorml2_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_dvec_t *M = tk_dvec_peek(L, 1, "dvec");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_dims = tk_lua_checkunsigned(L, 3, "n_dims");
  tk_dvec_rnorml2(M->a, n_samples, n_dims);
  lua_pushvalue(L, 1);
  return 1;
}

static inline int tk_dvec_scores_kaiser_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_kaiser(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_scores_max_curvature_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_max_curvature(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_scores_lmethod_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_lmethod(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_scores_max_gap_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_max_gap(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_scores_max_drop_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_max_drop(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_scores_max_acceleration_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_max_acceleration(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static inline int tk_dvec_scores_tolerance_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double tolerance = luaL_optnumber(L, 2, 1e-3);
  size_t start, end;
  tk_dvec_scores_tolerance(scores->a, scores->n, tolerance, &start, &end);
  double start_val = (start < scores->n) ? scores->a[start] : 0.0;
  double end_val = (end < scores->n) ? scores->a[end] : 0.0;
  lua_pushnumber(L, (lua_Number)start_val);
  lua_pushnumber(L, (lua_Number)start);
  lua_pushnumber(L, (lua_Number)end_val);
  lua_pushnumber(L, (lua_Number)end);
  return 4;
}

static inline int tk_dvec_scores_plateau_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double tolerance = luaL_optnumber(L, 2, 1e-3);
  double val;
  size_t idx = tk_dvec_scores_plateau(scores->a, scores->n, tolerance, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)idx);
  return 2;
}

static luaL_Reg tk_dvec_lua_mt_ext2_fns[] =
{
  { "multiply_bits", tk_dvec_multiply_bits_lua },
  { "center", tk_dvec_center_lua },
  { "rnorml2", tk_dvec_rnorml2_lua },
  { "scores_kaiser", tk_dvec_scores_kaiser_lua },
  { "scores_max_curvature", tk_dvec_scores_max_curvature_lua },
  { "scores_lmethod", tk_dvec_scores_lmethod_lua },
  { "scores_max_gap", tk_dvec_scores_max_gap_lua },
  { "scores_max_drop", tk_dvec_scores_max_drop_lua },
  { "scores_max_acceleration", tk_dvec_scores_max_acceleration_lua },
  { "scores_tolerance", tk_dvec_scores_tolerance_lua },
  { "scores_plateau", tk_dvec_scores_plateau_lua },
  { NULL, NULL }
};

int luaopen_santoku_dvec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_dvec_lua_fns);
  tk_dvec_create(L, 0, 0, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_dvec_lua_mt_fns);
  luaL_register(L, NULL, tk_dvec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_dvec_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
