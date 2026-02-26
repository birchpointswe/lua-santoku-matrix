#include <santoku/iuset.h>
#include <santoku/ivec.h>
#include <santoku/cvec.h>
#include <string.h>

static inline int tk_ivec_bits_top_chi2_lua (lua_State *L)
{
  lua_settop(L, 8);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_visible : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  char *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_isnil(L, 2)) {
  } else if (tk_lua_testuserdata(L, 2, "tk_cvec_t")) {
    tk_cvec_t *cvec = tk_cvec_peek(L, 2, "codes");
    codes = cvec->a;
  } else {
    labels = tk_ivec_peek(L, 2, "labels");
  }
  tk_ivec_bits_top_chi2(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_ivec_bits_top_mi_lua (lua_State *L)
{
  lua_settop(L, 8);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_visible : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  char *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_isnil(L, 2)) {
  } else if (tk_lua_testuserdata(L, 2, "tk_cvec_t")) {
    tk_cvec_t *cvec = tk_cvec_peek(L, 2, "codes");
    codes = cvec->a;
  } else {
    labels = tk_ivec_peek(L, 2, "labels");
  }
  tk_ivec_bits_top_mi(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_ivec_bits_top_bns_lua (lua_State *L)
{
  lua_settop(L, 8);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_visible : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  char *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_isnil(L, 2)) {
  } else if (tk_lua_testuserdata(L, 2, "tk_cvec_t")) {
    tk_cvec_t *cvec = tk_cvec_peek(L, 2, "codes");
    codes = cvec->a;
  } else {
    labels = tk_ivec_peek(L, 2, "labels");
  }
  tk_ivec_bits_top_bns(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_ivec_bits_top_entropy_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 3, "hidden");
  uint64_t top_k = lua_isnil(L, 4) ? n_hidden : tk_lua_checkunsigned(L, 4, "top_k");
  tk_ivec_bits_top_entropy(L, set_bits, n_samples, n_hidden, top_k);
  return 2;
}

static inline int tk_ivec_bits_top_idf_lua (lua_State *L)
{
  lua_settop(L, 6);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 3, "features");
  uint64_t top_k = lua_isnil(L, 4) ? n_visible : tk_lua_checkunsigned(L, 4, "top_k");
  double min_df = tk_lua_optnumber(L, 5, "min_df", 0.0);
  double max_df = tk_lua_optnumber(L, 6, "max_df", 1.0);
  if (!tk_ivec_bits_top_idf(L, set_bits, n_samples, n_visible, min_df, max_df, top_k))
    return tk_lua_verror(L, 2, "bits_top_idf", "allocation failed");
  return 2;
}

static inline int tk_ivec_bits_top_df_lua (lua_State *L)
{
  lua_settop(L, 6);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 3, "features");
  uint64_t top_k = lua_isnil(L, 4) ? n_visible : tk_lua_checkunsigned(L, 4, "top_k");
  double min_df = tk_lua_optnumber(L, 5, "min_df", 0.0);
  double max_df = tk_lua_optnumber(L, 6, "max_df", 1.0);
  if (!tk_ivec_bits_top_df(L, set_bits, n_samples, n_visible, min_df, max_df, top_k))
    return tk_lua_verror(L, 2, "bits_top_df", "allocation failed");
  return 2;
}

static inline int tk_ivec_bits_top_reg_f_lua (lua_State *L)
{
  lua_settop(L, 8);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  tk_dvec_t *targets = tk_dvec_peek(L, 2, "targets");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  tk_ivec_bits_top_reg_f(L, set_bits, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_ivec_bits_top_reg_pearson_lua (lua_State *L)
{
  lua_settop(L, 8);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  tk_dvec_t *targets = tk_dvec_peek(L, 2, "targets");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  tk_ivec_bits_top_reg_pearson(L, set_bits, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_ivec_bits_top_reg_r2_lua (lua_State *L)
{
  lua_settop(L, 8);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  tk_dvec_t *targets = tk_dvec_peek(L, 2, "targets");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  tk_ivec_bits_top_reg_r2(L, set_bits, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_ivec_bits_top_reg_cohen_lua (lua_State *L)
{
  lua_settop(L, 8);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  tk_dvec_t *targets = tk_dvec_peek(L, 2, "targets");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  tk_ivec_bits_top_reg_cohen(L, set_bits, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_ivec_bits_bipartite_lua (lua_State *L)
{
  int n_args = lua_gettop(L);
  tk_ivec_t *src = tk_ivec_peek(L, 1, "src_bits");
  uint64_t n_docs = tk_lua_checkunsigned(L, 2, "n_docs");
  uint64_t n_labels = tk_lua_checkunsigned(L, 3, "n_labels");
  const char *mode = luaL_checkstring(L, 4);
  tk_ivec_t *tokens = (n_args >= 5 && !lua_isnil(L, 5)) ? tk_ivec_peek(L, 5, "tokens") : NULL;
  uint64_t n_tokens = (n_args >= 6) ? tk_lua_checkunsigned(L, 6, "n_tokens") : 0;
  uint64_t n_features = 0;
  tk_ivec_t *out = tk_ivec_bits_bipartite(L, src, n_docs, n_labels, mode, tokens, n_tokens, &n_features);
  if (!out)
    return luaL_error(L, "bits_bipartite: allocation failed");
  lua_pushinteger(L, (lua_Integer)n_features);
  return 2;
}

static inline int tk_ivec_bits_select_lua (lua_State *L)
{
  int n_args = lua_gettop(L);
  tk_ivec_t *src_bits = tk_ivec_peek(L, 1, "src_bits");
  tk_ivec_t *selected_features = lua_isnil(L, 2) ? NULL : tk_ivec_peek(L, 2, "selected_features");
  tk_ivec_t *sample_ids = NULL;
  uint64_t n_features;
  tk_ivec_t *dest = NULL;
  uint64_t dest_sample = 0;
  uint64_t dest_stride = 0;

  if (n_args == 4) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
  } else if (n_args == 5) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
    dest = tk_ivec_peek(L, 5, "dest");
  } else if (n_args == 6) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
    dest = tk_ivec_peek(L, 5, "dest");
    dest_sample = tk_lua_checkunsigned(L, 6, "dest_sample");
  } else if (n_args == 7) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
    dest = tk_ivec_peek(L, 5, "dest");
    dest_sample = tk_lua_checkunsigned(L, 6, "dest_sample");
    dest_stride = tk_lua_checkunsigned(L, 7, "dest_stride");
  } else {
    n_features = tk_lua_checkunsigned(L, 3, "n_features");
  }

  tk_ivec_t *result = tk_ivec_bits_select(src_bits, selected_features, sample_ids, n_features, dest, dest_sample, dest_stride);
  if (result == NULL)
    return luaL_error(L, "bits_select failed");
  lua_settop(L, 1);
  return 1;
}

static inline int tk_ivec_bits_to_cvec_lua (lua_State *L)
{
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  if (tk_lua_testuserdata(L, 4, "tk_ivec_t")) {
    uint64_t n_visible = tk_lua_checkunsigned(L, 3, "n_visible");
    tk_ivec_t *offsets = tk_ivec_peek(L, 4, "offsets");
    tk_ivec_t *features = tk_ivec_peek(L, 5, "features");
    bool flip_interleave = lua_toboolean(L, 6);
    uint64_t max_k = tk_ivec_bits_to_cvec_grouped(L, set_bits, n_samples, n_visible, offsets, features, flip_interleave);
    lua_pushinteger(L, (lua_Integer)max_k);
    return 2;
  } else {
    uint64_t n_features = tk_lua_checkunsigned(L, 3, "features");
    bool flip_interleave = lua_toboolean(L, 4);
    tk_ivec_bits_to_cvec(L, set_bits, n_samples, n_features, flip_interleave);
    return 1;
  }
}

static inline int tk_ivec_bits_from_cvec_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_cvec_t *cvec = tk_cvec_peek(L, 1, "bitmap");
  const char *bm = cvec->a;
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "features");
  tk_ivec_bits_from_cvec(L, bm, n_samples, n_features);
  return 1;
}

static inline int tk_ivec_bits_extend_cvec_helper (
  lua_State *L,
  tk_ivec_t *base,
  tk_cvec_t *ext_cvec,
  uint64_t n_base_features,
  uint64_t n_ext_features
) {
  tk_ivec_t *ext = tk_cvec_bits_to_ivec(L, ext_cvec, n_ext_features);
  if (!ext)
    return -1;
  tk_ivec_t *result = tk_ivec_bits_extend(base, ext, n_base_features, n_ext_features);
  tk_ivec_destroy(ext);
  lua_remove(L, -1);
  if (!result)
    return -1;
  return 0;
}

static inline int tk_ivec_bits_extend_lua (lua_State *L)
{
  int nargs = lua_gettop(L);

  if (nargs == 4) {
    tk_ivec_t *base = tk_ivec_peek(L, 1, "base_bits");
    uint64_t n_base_features = tk_lua_checkunsigned(L, 3, "features");
    uint64_t n_ext_features = tk_lua_checkunsigned(L, 4, "extended");
    tk_ivec_t *ext_ivec = tk_ivec_peekopt(L, 2);
    if (ext_ivec) {
      tk_ivec_t *result = tk_ivec_bits_extend(base, ext_ivec, n_base_features, n_ext_features);
      if (!result)
        return luaL_error(L, "bits_extend: allocation failed");
    } else {
      tk_cvec_t *ext_cvec = tk_cvec_peek(L, 2, "ext_bits");
      if (tk_ivec_bits_extend_cvec_helper(L, base, ext_cvec, n_base_features, n_ext_features) != 0)
        return luaL_error(L, "bits_extend: allocation failed");
    }
  } else if (nargs == 6 || nargs == 7) {
    tk_ivec_t *base = tk_ivec_peek(L, 1, "base_bits");
    tk_ivec_t *aids = tk_ivec_peek(L, 3, "aids");
    tk_ivec_t *bids = tk_ivec_peek(L, 4, "bids");
    uint64_t n_base_features = tk_lua_checkunsigned(L, 5, "features");
    uint64_t n_ext_features = tk_lua_checkunsigned(L, 6, "extended");

    bool project = false;
    if (nargs == 7)
      project = lua_toboolean(L, 7);
    tk_ivec_t *ext_ivec = tk_ivec_peekopt(L, 2);
    if (ext_ivec) {
      if (tk_ivec_bits_extend_mapped(base, ext_ivec, aids, bids, n_base_features, n_ext_features, project) != 0)
        return luaL_error(L, "bits_extend_mapped: allocation failed");
    } else {
      tk_cvec_t *ext_cvec = tk_cvec_peek(L, 2, "ext_bits");
      tk_ivec_t *ext = tk_cvec_bits_to_ivec(L, ext_cvec, n_ext_features);
      if (!ext)
        return luaL_error(L, "bits_to_ivec: allocation failed");
      if (tk_ivec_bits_extend_mapped(base, ext, aids, bids, n_base_features, n_ext_features, project) != 0) {
        tk_ivec_destroy(ext);
        lua_remove(L, -1);
        return luaL_error(L, "bits_extend_mapped: allocation failed");
      }
      tk_ivec_destroy(ext);
      lua_remove(L, -1);
    }
  } else {
    return luaL_error(L, "bits_extend expects 4, 6, or 7 arguments, got %d", nargs);
  }
  lua_settop(L, 1);
  return 1;
}

static inline int tk_ivec_copy_pkeys_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_ivec_t *m0 = tk_ivec_peek(L, 1, "dest");
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
  tk_ivec_copy_pkeys(m0, m1, start, end, dest);
  return 0;
}

static inline int tk_ivec_copy_rkeys_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_ivec_t *m0 = tk_ivec_peek(L, 1, "dest");
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
  tk_ivec_copy_rkeys(m0, m1, start, end, dest);
  return 0;
}

static inline int tk_ivec_copy_pvalues_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_ivec_t *m0 = tk_ivec_peek(L, 1, "dest");
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
  tk_ivec_copy_pvalues(m0, m1, start, end, dest);
  return 0;
}

static inline int tk_ivec_copy_rvalues_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_ivec_t *m0 = tk_ivec_peek(L, 1, "dest");
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
  tk_ivec_copy_rvalues(m0, m1, start, end, dest);
  return 0;
}

static inline int tk_ivec_set_stats_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_dvec_t *weights = lua_isnil(L, 3) ? NULL : tk_dvec_peek(L, 3, "weights");
  double inter_w, sum_a, sum_b;
  tk_ivec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  lua_pushnumber(L, inter_w);
  lua_pushnumber(L, sum_a);
  lua_pushnumber(L, sum_b);
  return 3;
}

static inline int tk_ivec_set_similarity_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_dvec_t *weights = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_dvec_peek(L, 3, "weights") : NULL;
  tk_ivec_sim_type_t type = TK_IVEC_JACCARD;
  if (nargs >= 4 && lua_isstring(L, 4)) {
    const char *type_str = lua_tostring(L, 4);
    if (strcmp(type_str, "jaccard") == 0) {
      type = TK_IVEC_JACCARD;
    } else if (strcmp(type_str, "overlap") == 0) {
      type = TK_IVEC_OVERLAP;
    } else if (strcmp(type_str, "dice") == 0) {
      type = TK_IVEC_DICE;
    } else if (strcmp(type_str, "tversky") == 0) {
      type = TK_IVEC_TVERSKY;
    }
  }
  double tversky_alpha = (nargs >= 5) ? luaL_checknumber(L, 5) : 0.5;
  double tversky_beta = (nargs >= 6) ? luaL_checknumber(L, 6) : 0.5;
  double sim = tk_ivec_set_similarity(a->a, a->n, b->a, b->n, weights, type, tversky_alpha, tversky_beta);
  lua_pushnumber(L, sim);
  return 1;
}

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

static inline int tk_ivec_set_weights_by_rank_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *features = tk_ivec_peek(L, 1, "features");
  tk_dvec_t *weights = (nargs >= 2 && !lua_isnil(L, 2)) ? tk_dvec_peek(L, 2, "weights") : NULL;
  tk_ivec_t *ranks = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_ivec_peek(L, 3, "ranks") : NULL;
  uint64_t n_ranks = (nargs >= 4) ? tk_lua_checkunsigned(L, 4, "n_ranks") : 1;
  tk_dvec_t *weights_by_rank = tk_dvec_create(L, n_ranks, 0, 0);
  tk_ivec_set_weights_by_rank(features->a, features->n, weights, ranks, n_ranks, weights_by_rank->a);
  weights_by_rank->n = n_ranks;
  return 1;
}

static inline int tk_ivec_set_similarity_by_rank_lua (lua_State *L)
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
  tk_ivec_sim_type_t type = TK_IVEC_JACCARD;
  if (nargs >= 9 && lua_isstring(L, 9)) {
    const char *type_str = lua_tostring(L, 9);
    if (strcmp(type_str, "jaccard") == 0) {
      type = TK_IVEC_JACCARD;
    } else if (strcmp(type_str, "overlap") == 0) {
      type = TK_IVEC_OVERLAP;
    } else if (strcmp(type_str, "dice") == 0) {
      type = TK_IVEC_DICE;
    } else if (strcmp(type_str, "tversky") == 0) {
      type = TK_IVEC_TVERSKY;
    }
  }
  double tversky_alpha = (nargs >= 10) ? luaL_checknumber(L, 10) : 0.5;
  double tversky_beta = (nargs >= 11) ? luaL_checknumber(L, 11) : 0.5;
  double sim = tk_ivec_set_similarity_by_rank(
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

static inline int tk_ivec_scores_elbow_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_ivec_t *scores = tk_ivec_peek(L, 1, "ivec");
  const char *method = luaL_checkstring(L, 2);
  int64_t alpha = (int64_t)luaL_optnumber(L, 3, 0.0);
  int64_t val;
  size_t idx;
  if (strcmp(method, "lmethod") == 0) {
    idx = tk_ivec_scores_lmethod(scores, &val);
  } else if (strcmp(method, "max_gap") == 0) {
    idx = tk_ivec_scores_max_gap(scores, &val);
  } else if (strcmp(method, "max_curvature") == 0) {
    idx = tk_ivec_scores_max_curvature(scores, &val);
  } else if (strcmp(method, "kneedle") == 0) {
    double sensitivity = (alpha > 0) ? (double)alpha : 1.0;
    idx = tk_ivec_scores_kneedle(scores, sensitivity, &val);
  } else if (strcmp(method, "plateau") == 0) {
    int64_t tolerance = (alpha > 0) ? alpha : 0;
    idx = tk_ivec_scores_plateau(scores, tolerance, &val);
  } else if (strcmp(method, "otsu") == 0) {
    idx = tk_ivec_scores_otsu(scores, &val);
  } else if (strcmp(method, "first_gap") == 0) {
    double ratio = (alpha > 0) ? (double)alpha : 3.0;
    idx = tk_ivec_scores_first_gap(scores, ratio, &val);
  } else {
    return luaL_error(L, "unknown elbow method: %s", method);
  }
  lua_pushinteger(L, (lua_Integer)val);
  lua_pushinteger(L, (lua_Integer)(idx + 1));
  return 2;
}

static inline int tk_ivec_bits_to_csr_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_labels = tk_lua_checkunsigned(L, 3, "n_labels");
  if (n_labels == 0)
    return luaL_error(L, "n_labels must be > 0");
  tk_ivec_t *offsets = tk_ivec_create(L, n_samples + 1, NULL, NULL);
  memset(offsets->a, 0, (n_samples + 1) * sizeof(int64_t));
  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t v = set_bits->a[i];
    if (v < 0) continue;
    uint64_t s = (uint64_t)v / n_labels;
    if (s < n_samples)
      offsets->a[s + 1]++;
  }
  for (uint64_t i = 1; i <= n_samples; i++)
    offsets->a[i] += offsets->a[i - 1];
  uint64_t total = (uint64_t)offsets->a[n_samples];
  tk_ivec_t *labels = tk_ivec_create(L, total, NULL, NULL);
  tk_ivec_t *counts = tk_ivec_create(L, n_samples, NULL, NULL);
  memset(counts->a, 0, n_samples * sizeof(int64_t));
  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t v = set_bits->a[i];
    if (v < 0) continue;
    uint64_t s = (uint64_t)v / n_labels;
    uint64_t label = (uint64_t)v % n_labels;
    if (s < n_samples) {
      uint64_t pos = (uint64_t)offsets->a[s] + (uint64_t)counts->a[s];
      labels->a[pos] = (int64_t)label;
      counts->a[s]++;
    }
  }
  lua_pop(L, 1);
  return 2;
}

static inline int tk_ivec_bits_transpose_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_ivec_t *bits = tk_ivec_peek(L, 1, "bits");
  uint64_t n_rows = tk_lua_checkunsigned(L, 2, "n_rows");
  uint64_t n_cols = tk_lua_checkunsigned(L, 3, "n_cols");
  if (n_cols == 0)
    return luaL_error(L, "n_cols must be > 0");
  tk_ivec_t *out = tk_ivec_create(L, bits->n, NULL, NULL);
  for (uint64_t i = 0; i < bits->n; i++) {
    int64_t v = bits->a[i];
    if (v < 0) {
      out->a[i] = v;
      continue;
    }
    uint64_t row = (uint64_t)v / n_cols;
    uint64_t col = (uint64_t)v % n_cols;
    out->a[i] = (int64_t)(col * n_rows + row);
  }
  tk_ivec_asc(out, 0, out->n);
  return 1;
}

static inline int tk_ivec_to_dvec_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_ivec_t *v = tk_ivec_peek(L, 1, "ivec");
  tk_ivec_to_dvec(L, v);
  return 1;
}

static luaL_Reg tk_ivec_lua_mt_ext2_fns[] =
{
  { "bits_to_csr", tk_ivec_bits_to_csr_lua },
  { "bits_transpose", tk_ivec_bits_transpose_lua },
  { "copy_pkeys", tk_ivec_copy_pkeys_lua },
  { "copy_rkeys", tk_ivec_copy_rkeys_lua },
  { "copy_pvalues", tk_ivec_copy_pvalues_lua },
  { "copy_rvalues", tk_ivec_copy_rvalues_lua },
  { "bits_top_chi2", tk_ivec_bits_top_chi2_lua },
  { "bits_top_mi", tk_ivec_bits_top_mi_lua },
  { "bits_top_bns", tk_ivec_bits_top_bns_lua },
  { "bits_top_entropy", tk_ivec_bits_top_entropy_lua },
  { "bits_top_idf", tk_ivec_bits_top_idf_lua },
  { "bits_top_df", tk_ivec_bits_top_df_lua },
  { "bits_top_reg_f", tk_ivec_bits_top_reg_f_lua },
  { "bits_top_reg_pearson", tk_ivec_bits_top_reg_pearson_lua },
  { "bits_top_reg_r2", tk_ivec_bits_top_reg_r2_lua },
  { "bits_top_reg_cohen", tk_ivec_bits_top_reg_cohen_lua },
  { "bits_select", tk_ivec_bits_select_lua },
  { "bits_bipartite", tk_ivec_bits_bipartite_lua },
  { "bits_to_cvec", tk_ivec_bits_to_cvec_lua },
  { "bits_extend", tk_ivec_bits_extend_lua },
  { "set_stats", tk_ivec_set_stats_lua },
  { "set_similarity", tk_ivec_set_similarity_lua },
  { "set_jaccard", tk_ivec_set_jaccard_lua },
  { "set_overlap", tk_ivec_set_overlap_lua },
  { "set_dice", tk_ivec_set_dice_lua },
  { "set_tversky", tk_ivec_set_tversky_lua },
  { "set_weights_by_rank", tk_ivec_set_weights_by_rank_lua },
  { "set_similarity_by_rank", tk_ivec_set_similarity_by_rank_lua },
  { "set_find", tk_ivec_set_find_lua },
  { "set_insert", tk_ivec_set_insert_lua },
  { "set_intersect", tk_ivec_set_intersect_lua },
  { "set_union", tk_ivec_set_union_lua },
  { "lookup", tk_ivec_lookup_lua },
  { "index", tk_ivec_index_lua },
  { "scores_elbow", tk_ivec_scores_elbow_lua },
  { "to_dvec", tk_ivec_to_dvec_lua },
  { NULL, NULL }
};

static inline int tk_ivec_from_rvec_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_rvec_t *R = tk_rvec_peek(L, 1, "rvec");
  tk_ivec_from_rvec(L, R);
  return 1;
}

static luaL_Reg tk_ivec_lua_ext_fns[] =
{
  { "bits_from_cvec", tk_ivec_bits_from_cvec_lua },
  { "from_rvec", tk_ivec_from_rvec_lua },
  { NULL, NULL }
};

int luaopen_santoku_ivec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_ivec_lua_fns);
  luaL_register(L, NULL, tk_ivec_lua_ext_fns);
  tk_ivec_create(L, 0, 0, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_ivec_lua_mt_fns);
  luaL_register(L, NULL, tk_ivec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_ivec_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
