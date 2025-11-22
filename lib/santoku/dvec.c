#include <santoku/iuset.h>
#include <santoku/ivec.h>
#include <santoku/dvec.h>

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
  lua_pushnumber(L, (lua_Number)(idx + 1));
  return 2;
}

static inline int tk_dvec_scores_max_curvature_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_max_curvature(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)(idx + 1));
  return 2;
}

static inline int tk_dvec_scores_lmethod_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_lmethod(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)(idx + 1));
  return 2;
}

static inline int tk_dvec_scores_max_gap_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_max_gap(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)(idx + 1));
  return 2;
}

static inline int tk_dvec_scores_max_drop_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_max_drop(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)(idx + 1));
  return 2;
}

static inline int tk_dvec_scores_max_acceleration_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *scores = tk_dvec_peek(L, 1, "dvec");
  double val;
  size_t idx = tk_dvec_scores_max_acceleration(scores->a, scores->n, &val);
  lua_pushnumber(L, (lua_Number)val);
  lua_pushnumber(L, (lua_Number)(idx + 1));
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
  lua_pushnumber(L, (lua_Number)(start + 1));
  lua_pushnumber(L, (lua_Number)end_val);
  lua_pushnumber(L, (lua_Number)(end + 1));
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
  lua_pushnumber(L, (lua_Number)(idx + 1));
  return 2;
}

static inline int tk_dvec_multiply_blas_lua (lua_State *L)
{
  lua_settop(L, 6);
  tk_dvec_t *a = tk_dvec_peek(L, 1, "a");
  tk_dvec_t *b = tk_dvec_peek(L, 2, "b");
  tk_dvec_t *c = tk_dvec_peek(L, 3, "c");
  uint64_t k = tk_lua_checkunsigned(L, 4, "k");
  bool transpose_a = lua_toboolean(L, 5);
  bool transpose_b = lua_toboolean(L, 6);
  tk_dvec_multiply_override(a, b, c, k, transpose_a, transpose_b);
  lua_pushvalue(L, 3);
  return 1;
}

static inline int tk_dvec_rsums_blas_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *m = tk_dvec_peek(L, 1, "dvec");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_dvec_rsums_override(L, m, cols);
  return 1;
}

static inline int tk_dvec_csums_blas_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *m = tk_dvec_peek(L, 1, "dvec");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_dvec_csums_override(L, m, cols);
  return 1;
}

static inline int tk_dvec_rmags_blas_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *m = tk_dvec_peek(L, 1, "dvec");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_dvec_rmags_override(L, m, cols);
  return 1;
}

static inline int tk_dvec_cmags_blas_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *m = tk_dvec_peek(L, 1, "dvec");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_dvec_cmags_override(L, m, cols);
  return 1;
}

static inline int tk_dvec_mtx_select_lua (lua_State *L)
{
  int n_args = lua_gettop(L);
  tk_dvec_t *src_matrix = tk_dvec_peek(L, 1, "src_matrix");
  tk_ivec_t *selected_features = lua_isnil(L, 2) ? NULL : tk_ivec_peek(L, 2, "selected_features");
  tk_ivec_t *sample_ids = NULL;
  uint64_t n_features;
  tk_dvec_t *dest = NULL;
  uint64_t dest_sample = 0;
  uint64_t dest_stride = 0;

  if (n_args == 4) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
  } else if (n_args == 5) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
    dest = tk_dvec_peek(L, 5, "dest");
  } else if (n_args == 6) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
    dest = tk_dvec_peek(L, 5, "dest");
    dest_sample = tk_lua_checkunsigned(L, 6, "dest_sample");
  } else if (n_args == 7) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
    dest = tk_dvec_peek(L, 5, "dest");
    dest_sample = tk_lua_checkunsigned(L, 6, "dest_sample");
    dest_stride = tk_lua_checkunsigned(L, 7, "dest_stride");
  } else {
    n_features = tk_lua_checkunsigned(L, 3, "n_features");
  }

  tk_dvec_t *result = tk_dvec_mtx_select(src_matrix, selected_features, sample_ids, n_features, dest, dest_sample, dest_stride);
  if (result == NULL)
    return luaL_error(L, "mtx_select failed");
  return 0;
}

static inline int tk_dvec_mtx_extend_lua (lua_State *L)
{
  int nargs = lua_gettop(L);

  if (nargs == 4) {
    tk_dvec_t *base = tk_dvec_peek(L, 1, "base_matrix");
    tk_dvec_t *ext = tk_dvec_peek(L, 2, "ext_matrix");
    uint64_t n_base_features = tk_lua_checkunsigned(L, 3, "n_base_features");
    uint64_t n_ext_features = tk_lua_checkunsigned(L, 4, "n_ext_features");
    tk_dvec_t *result = tk_dvec_mtx_extend(base, ext, n_base_features, n_ext_features);
    if (!result)
      return luaL_error(L, "mtx_extend: allocation failed");
  } else if (nargs == 6 || nargs == 7) {
    tk_dvec_t *base = tk_dvec_peek(L, 1, "base_matrix");
    tk_dvec_t *ext = tk_dvec_peek(L, 2, "ext_matrix");
    tk_ivec_t *aids = tk_ivec_peek(L, 3, "aids");
    tk_ivec_t *bids = tk_ivec_peek(L, 4, "bids");
    uint64_t n_base_features = tk_lua_checkunsigned(L, 5, "n_base_features");
    uint64_t n_ext_features = tk_lua_checkunsigned(L, 6, "n_ext_features");
    bool project = false;
    if (nargs == 7)
      project = lua_toboolean(L, 7);
    if (tk_dvec_mtx_extend_mapped(base, ext, aids, bids, n_base_features, n_ext_features, project) != 0)
      return luaL_error(L, "mtx_extend_mapped: allocation failed");
  } else {
    return luaL_error(L, "mtx_extend expects 4, 6, or 7 arguments, got %d", nargs);
  }
  return 0;
}

static inline int tk_dvec_mtx_top_variance_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_dvec_t *matrix = tk_dvec_peek(L, 1, "matrix");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t top_k = lua_isnil(L, 4) ? n_features : tk_lua_checkunsigned(L, 4, "top_k");
  tk_dvec_mtx_top_variance(L, matrix, n_samples, n_features, top_k);
  return 2;
}

static inline int tk_dvec_mtx_top_skewness_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_dvec_t *matrix = tk_dvec_peek(L, 1, "matrix");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t top_k = lua_isnil(L, 4) ? n_features : tk_lua_checkunsigned(L, 4, "top_k");
  tk_dvec_mtx_top_skewness(L, matrix, n_samples, n_features, top_k);
  return 2;
}

static inline int tk_dvec_mtx_top_esber_lua (lua_State *L)
{
  lua_settop(L, 5);
  tk_dvec_t *matrix = tk_dvec_peek(L, 1, "matrix");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t top_k = lua_isnil(L, 4) ? n_features : tk_lua_checkunsigned(L, 4, "top_k");
  uint64_t n_bins = lua_isnil(L, 5) ? 0 : tk_lua_checkunsigned(L, 5, "n_bins");
  tk_dvec_mtx_top_esber(L, matrix, n_samples, n_features, top_k, n_bins);
  return 2;
}

static inline int tk_dvec_mtx_top_bimodality_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_dvec_t *matrix = tk_dvec_peek(L, 1, "matrix");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t top_k = lua_isnil(L, 4) ? n_features : tk_lua_checkunsigned(L, 4, "top_k");
  tk_dvec_mtx_top_bimodality(L, matrix, n_samples, n_features, top_k);
  return 2;
}

static inline int tk_dvec_mtx_top_dip_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_dvec_t *matrix = tk_dvec_peek(L, 1, "matrix");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t top_k = lua_isnil(L, 4) ? n_features : tk_lua_checkunsigned(L, 4, "top_k");
  tk_dvec_mtx_top_dip(L, matrix, n_samples, n_features, top_k);
  return 2;
}

static luaL_Reg tk_dvec_lua_mt_ext2_fns[] =
{
  { "center", tk_dvec_center_lua },
  { "rnorml2", tk_dvec_rnorml2_lua },
  { "multiply", tk_dvec_multiply_blas_lua },
  { "rsums", tk_dvec_rsums_blas_lua },
  { "csums", tk_dvec_csums_blas_lua },
  { "rmags", tk_dvec_rmags_blas_lua },
  { "cmags", tk_dvec_cmags_blas_lua },
  { "scores_kaiser", tk_dvec_scores_kaiser_lua },
  { "scores_max_curvature", tk_dvec_scores_max_curvature_lua },
  { "scores_lmethod", tk_dvec_scores_lmethod_lua },
  { "scores_max_gap", tk_dvec_scores_max_gap_lua },
  { "scores_max_drop", tk_dvec_scores_max_drop_lua },
  { "scores_max_acceleration", tk_dvec_scores_max_acceleration_lua },
  { "scores_tolerance", tk_dvec_scores_tolerance_lua },
  { "scores_plateau", tk_dvec_scores_plateau_lua },
  { "mtx_select", tk_dvec_mtx_select_lua },
  { "mtx_extend", tk_dvec_mtx_extend_lua },
  { "mtx_top_variance", tk_dvec_mtx_top_variance_lua },
  { "mtx_top_skewness", tk_dvec_mtx_top_skewness_lua },
  { "mtx_top_esber", tk_dvec_mtx_top_esber_lua },
  { "mtx_top_bimodality", tk_dvec_mtx_top_bimodality_lua },
  { "mtx_top_dip", tk_dvec_mtx_top_dip_lua },
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
