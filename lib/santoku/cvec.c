#include <santoku/iuset.h>
#include <santoku/cvec.h>
#include <santoku/ivec.h>
#include <stdbool.h>

static inline int tk_cvec_bits_flip_interleave_lua (lua_State *L) {
  lua_settop(L, 2);
  tk_cvec_t *v = tk_cvec_peek(L, 1, "cvec");
  uint64_t n_features = tk_lua_checkunsigned(L, 2, "n_features");
  if (!tk_cvec_bits_flip_interleave(v, n_features))
    return tk_lua_verror(L, 2, "bits_flip_interleave", "allocation failed");
  return 0;
}

static inline int tk_cvec_bits_to_ivec_lua (lua_State *L) {
  lua_settop(L, 2);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  uint64_t n_features = tk_lua_checkunsigned(L, 2, "n_features");
  tk_cvec_bits_to_ivec(L, bitmap, n_features);
  return 1;
}

static inline int tk_cvec_bits_from_ivec_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  tk_cvec_bits_from_ivec(L, set_bits, n_samples, n_features);
  return 1;
}

static inline tk_cvec_t *tk_cvec_bits_extend_ivec_helper (
  lua_State *L,
  tk_cvec_t *base,
  tk_ivec_t *ext_ivec,
  uint64_t n_base_features,
  uint64_t n_ext_features,
  uint64_t n_samples
) {
  tk_cvec_t *ext = tk_cvec_bits_from_ivec(L, ext_ivec, n_samples, n_ext_features);
  tk_cvec_t *result = tk_cvec_bits_extend(base, ext, n_base_features, n_ext_features);
  tk_cvec_destroy(ext);
  lua_remove(L, -1);
  return result;
}

static inline int tk_cvec_bits_extend_lua (lua_State *L) {
  int nargs = lua_gettop(L);

  if (nargs == 4) {
    tk_cvec_t *base = tk_cvec_peek(L, 1, "base");
    uint64_t n_base_features = tk_lua_checkunsigned(L, 3, "n_base_features");
    uint64_t n_ext_features = tk_lua_checkunsigned(L, 4, "n_ext_features");
    tk_cvec_t *ext_cvec = tk_cvec_peekopt(L, 2);
    if (ext_cvec) {
      if (!tk_cvec_bits_extend(base, ext_cvec, n_base_features, n_ext_features))
        return tk_lua_verror(L, 2, "bits_extend", "allocation failed");
    } else {
      tk_ivec_t *ext_ivec = tk_ivec_peek(L, 2, "ext");
      uint64_t n_samples = base->n / TK_CVEC_BITS_BYTES(n_base_features);
      if (!tk_cvec_bits_extend_ivec_helper(L, base, ext_ivec, n_base_features, n_ext_features, n_samples))
        return tk_lua_verror(L, 2, "bits_extend", "allocation failed");
    }
  } else if (nargs == 6 || nargs == 7) {
    tk_cvec_t *base = tk_cvec_peek(L, 1, "base");
    tk_ivec_t *aids = tk_ivec_peek(L, 3, "aids");
    tk_ivec_t *bids = tk_ivec_peek(L, 4, "bids");
    uint64_t n_base_features = tk_lua_checkunsigned(L, 5, "n_base_features");
    uint64_t n_ext_features = tk_lua_checkunsigned(L, 6, "n_ext_features");

    bool project = false;
    if (nargs == 7) {
      project = lua_toboolean(L, 7);
    }

    tk_cvec_t *ext_cvec = tk_cvec_peekopt(L, 2);
    if (ext_cvec) {
      if (tk_cvec_bits_extend_mapped(base, ext_cvec, aids, bids, n_base_features, n_ext_features, project) != 0)
        return tk_lua_verror(L, 2, "bits_extend_mapped", "allocation failed");
    } else {
      tk_ivec_t *ext_ivec = tk_ivec_peek(L, 2, "ext");
      uint64_t n_samples = ext_ivec->n > 0 ?
        ((uint64_t)ext_ivec->a[ext_ivec->n - 1] / n_ext_features + 1) : 0;
      tk_cvec_t *ext = tk_cvec_bits_from_ivec(L, ext_ivec, n_samples, n_ext_features);
      if (!ext)
        return tk_lua_verror(L, 2, "bits_from_ivec", "allocation failed");
      if (tk_cvec_bits_extend_mapped(base, ext, aids, bids, n_base_features, n_ext_features, project) != 0) {
        tk_cvec_destroy(ext);
        lua_remove(L, -1);
        return tk_lua_verror(L, 2, "bits_extend_mapped", "allocation failed");
      }
      tk_cvec_destroy(ext);
      lua_remove(L, -1);
    }
  } else {
    return luaL_error(L, "bits_extend expects 4, 6, or 7 arguments, got %d", nargs);
  }
  return 0;
}

static inline int tk_cvec_bits_select_lua (lua_State *L) {
  int n_args = lua_gettop(L);
  tk_cvec_t *src_bitmap = tk_cvec_peek(L, 1, "src_bitmap");
  tk_ivec_t *selected_features = lua_isnil(L, 2) ? NULL : tk_ivec_peek(L, 2, "selected_features");
  tk_ivec_t *sample_ids = NULL;
  uint64_t n_features;
  tk_cvec_t *dest = NULL;
  uint64_t dest_sample = 0;
  uint64_t dest_stride = 0;

  if (n_args == 4) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
  } else if (n_args == 5) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
    dest = tk_cvec_peek(L, 5, "dest");
  } else if (n_args == 6) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
    dest = tk_cvec_peek(L, 5, "dest");
    dest_sample = tk_lua_checkunsigned(L, 6, "dest_sample");
  } else if (n_args == 7) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
    dest = tk_cvec_peek(L, 5, "dest");
    dest_sample = tk_lua_checkunsigned(L, 6, "dest_sample");
    dest_stride = tk_lua_checkunsigned(L, 7, "dest_stride");
  } else {
    n_features = tk_lua_checkunsigned(L, 3, "n_features");
  }

  if (!tk_cvec_bits_select(src_bitmap, selected_features, sample_ids, n_features, dest, dest_sample, dest_stride))
    return tk_lua_verror(L, 2, "bits_select", "allocation failed");
  return 0;
}

static inline int tk_cvec_bits_top_chi2_lua (lua_State *L) {
  lua_settop(L, 9);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  tk_cvec_t *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_isnil(L, 2)) {
  } else if (tk_lua_testuserdata(L, 2, "tk_cvec_t")) {
    codes = tk_cvec_peek(L, 2, "codes");
  } else {
    labels = tk_ivec_peek(L, 2, "labels");
  }
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(tk_lua_optstring(L, 8, "pool", "max"));
  tk_cvec_bits_top_chi2(L, bitmap, codes, labels, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_cvec_bits_top_mi_lua (lua_State *L) {
  lua_settop(L, 9);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  tk_cvec_t *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_isnil(L, 2)) {
  } else if (tk_lua_testuserdata(L, 2, "tk_cvec_t")) {
    codes = tk_cvec_peek(L, 2, "codes");
  } else {
    labels = tk_ivec_peek(L, 2, "labels");
  }
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(tk_lua_optstring(L, 8, "pool", "max"));
  tk_cvec_bits_top_mi(L, bitmap, codes, labels, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_cvec_bits_top_bns_lua (lua_State *L) {
  lua_settop(L, 9);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  tk_cvec_t *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_isnil(L, 2)) {
  } else if (tk_lua_testuserdata(L, 2, "tk_cvec_t")) {
    codes = tk_cvec_peek(L, 2, "codes");
  } else {
    labels = tk_ivec_peek(L, 2, "labels");
  }
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(tk_lua_optstring(L, 8, "pool", "max"));
  tk_cvec_bits_top_bns(L, bitmap, codes, labels, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_cvec_bits_top_entropy_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_cvec_t *codes = tk_cvec_peek(L, 1, "codes");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 3, "n_hidden");
  uint64_t top_k = lua_isnil(L, 4) ? n_hidden : tk_lua_checkunsigned(L, 4, "top_k");
  tk_cvec_bits_top_entropy(L, codes, n_samples, n_hidden, top_k);
  return 2;
}

static inline int tk_cvec_bits_top_idf_lua (lua_State *L) {
  lua_settop(L, 6);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t top_k = lua_isnil(L, 4) ? n_features : tk_lua_checkunsigned(L, 4, "top_k");
  double min_df = tk_lua_optnumber(L, 5, "min_df", 0.0);
  double max_df = tk_lua_optnumber(L, 6, "max_df", 1.0);
  if (!tk_cvec_bits_top_idf(L, bitmap, n_samples, n_features, min_df, max_df, top_k))
    return tk_lua_verror(L, 2, "bits_top_idf", "allocation failed");
  return 2;
}

static inline int tk_cvec_bits_top_df_lua (lua_State *L) {
  lua_settop(L, 6);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t top_k = lua_isnil(L, 4) ? n_features : tk_lua_checkunsigned(L, 4, "top_k");
  double min_df = tk_lua_optnumber(L, 5, "min_df", 0.0);
  double max_df = tk_lua_optnumber(L, 6, "max_df", 1.0);
  if (!tk_cvec_bits_top_df(L, bitmap, n_samples, n_features, min_df, max_df, top_k))
    return tk_lua_verror(L, 2, "bits_top_df", "allocation failed");
  return 2;
}

static inline int tk_cvec_bits_top_reg_pearson_lua (lua_State *L) {
  lua_settop(L, 8);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  tk_dvec_t *targets = tk_dvec_peek(L, 2, "targets");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  tk_cvec_bits_top_reg_pearson(L, bitmap, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_cvec_bits_top_reg_r2_lua (lua_State *L) {
  lua_settop(L, 8);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  tk_dvec_t *targets = tk_dvec_peek(L, 2, "targets");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  tk_cvec_bits_top_reg_r2(L, bitmap, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_cvec_bits_top_reg_cohen_lua (lua_State *L) {
  lua_settop(L, 8);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  tk_dvec_t *targets = tk_dvec_peek(L, 2, "targets");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  tk_cvec_bits_top_reg_cohen(L, bitmap, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_cvec_bits_top_reg_auc_lua (lua_State *L) {
  lua_settop(L, 8);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  tk_dvec_t *targets = tk_dvec_peek(L, 2, "targets");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  tk_cvec_bits_top_reg_auc(L, bitmap, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static inline int tk_cvec_bits_popcount_lua (lua_State *L) {
  lua_settop(L, 2);
  tk_cvec_t *vec = tk_cvec_peek(L, 1, "cvec");
  uint64_t n_bits = tk_lua_checkunsigned(L, 2, "n_bits");
  uint64_t count = tk_cvec_bits_popcount((const uint8_t *)vec->a, n_bits);
  lua_pushinteger(L, (lua_Integer)count);
  return 1;
}

static inline int tk_cvec_bits_hamming_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_cvec_t *a = tk_cvec_peek(L, 1, "a");
  tk_cvec_t *b = tk_cvec_peek(L, 2, "b");
  uint64_t n_bits = tk_lua_checkunsigned(L, 3, "n_bits");
  uint64_t distance = tk_cvec_bits_hamming((const uint8_t *)a->a, (const uint8_t *)b->a, n_bits);
  lua_pushinteger(L, (lua_Integer)distance);
  return 1;
}

static inline int tk_cvec_bits_hamming_mask_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_cvec_t *a = tk_cvec_peek(L, 1, "a");
  tk_cvec_t *b = tk_cvec_peek(L, 2, "b");
  tk_cvec_t *mask = tk_cvec_peek(L, 3, "mask");
  uint64_t n_bits = tk_lua_checkunsigned(L, 4, "n_bits");
  uint64_t distance = tk_cvec_bits_hamming_mask((const uint8_t *)a->a, (const uint8_t *)b->a, (const uint8_t *)mask->a, n_bits);
  lua_pushinteger(L, (lua_Integer)distance);
  return 1;
}

static inline int tk_cvec_bits_and_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_cvec_t *out = tk_cvec_peek(L, 1, "out");
  tk_cvec_t *a = tk_cvec_peek(L, 2, "a");
  tk_cvec_t *b = tk_cvec_peek(L, 3, "b");
  uint64_t n_bits = tk_lua_checkunsigned(L, 4, "n_bits");
  tk_cvec_bits_and((uint8_t *)out->a, (const uint8_t *)a->a, (const uint8_t *)b->a, n_bits);
  return 0;
}

static inline int tk_cvec_bits_or_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_cvec_t *out = tk_cvec_peek(L, 1, "out");
  tk_cvec_t *a = tk_cvec_peek(L, 2, "a");
  tk_cvec_t *b = tk_cvec_peek(L, 3, "b");
  uint64_t n_bits = tk_lua_checkunsigned(L, 4, "n_bits");
  tk_cvec_bits_or((uint8_t *)out->a, (const uint8_t *)a->a, (const uint8_t *)b->a, n_bits);
  return 0;
}

static inline int tk_cvec_bits_xor_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_cvec_t *out = tk_cvec_peek(L, 1, "out");
  tk_cvec_t *a = tk_cvec_peek(L, 2, "a");
  tk_cvec_t *b = tk_cvec_peek(L, 3, "b");
  uint64_t n_bits = tk_lua_checkunsigned(L, 4, "n_bits");
  tk_cvec_bits_xor((uint8_t *)out->a, (const uint8_t *)a->a, (const uint8_t *)b->a, n_bits);
  return 0;
}

static inline int tk_cvec_bits_to_ascii_lua (lua_State *L) {
  int nargs = lua_gettop(L);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");

  uint64_t start_bit = 0;
  uint64_t end_bit = bitmap->n * 8;


  if (nargs >= 3) {
    start_bit = tk_lua_checkunsigned(L, 2, "start_bit");
    end_bit = tk_lua_checkunsigned(L, 3, "end_bit");
  }


  if (start_bit > end_bit) {
    luaL_error(L, "start_bit must be <= end_bit");
  }
  if (end_bit > bitmap->n * 8) {
    luaL_error(L, "end_bit must be <= bitmap size in bits");
  }

  tk_cvec_bits_to_ascii(L, bitmap, start_bit, end_bit);
  return 1;
}

static inline int tk_cvec_bits_transpose_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  uint64_t n_rows = tk_lua_checkunsigned(L, 2, "n_rows");
  uint64_t n_cols = tk_lua_checkunsigned(L, 3, "n_cols");
  if (n_cols == 0)
    return luaL_error(L, "n_cols must be > 0");
  uint64_t in_bytes_per_row = TK_CVEC_BITS_BYTES(n_cols);
  uint64_t out_bytes_per_row = TK_CVEC_BITS_BYTES(n_rows);
  uint64_t out_size = n_cols * out_bytes_per_row;
  tk_cvec_t *out = tk_cvec_create(L, out_size, NULL, NULL);
  memset(out->a, 0, out_size);
  for (uint64_t row = 0; row < n_rows; row++) {
    for (uint64_t col = 0; col < n_cols; col++) {
      uint64_t in_byte = row * in_bytes_per_row + col / CHAR_BIT;
      uint64_t in_bit = col % CHAR_BIT;
      if ((uint8_t)bitmap->a[in_byte] & (1 << in_bit)) {
        uint64_t out_byte = col * out_bytes_per_row + row / CHAR_BIT;
        uint64_t out_bit = row % CHAR_BIT;
        out->a[out_byte] |= (1 << out_bit);
      }
    }
  }
  return 1;
}

static luaL_Reg tk_cvec_lua_mt_ext2_fns[] =
{
  { "bits_flip_interleave", tk_cvec_bits_flip_interleave_lua },
  { "bits_to_ivec", tk_cvec_bits_to_ivec_lua },
  { "bits_select", tk_cvec_bits_select_lua },
  { "bits_extend", tk_cvec_bits_extend_lua },
  { "bits_top_chi2", tk_cvec_bits_top_chi2_lua },
  { "bits_top_mi", tk_cvec_bits_top_mi_lua },
  { "bits_top_entropy", tk_cvec_bits_top_entropy_lua },
  { "bits_top_idf", tk_cvec_bits_top_idf_lua },
  { "bits_top_df", tk_cvec_bits_top_df_lua },
  { "bits_top_reg_pearson", tk_cvec_bits_top_reg_pearson_lua },
  { "bits_top_reg_r2", tk_cvec_bits_top_reg_r2_lua },
  { "bits_top_reg_cohen", tk_cvec_bits_top_reg_cohen_lua },
  { "bits_top_reg_auc", tk_cvec_bits_top_reg_auc_lua },
  { "bits_popcount", tk_cvec_bits_popcount_lua },
  { "bits_hamming", tk_cvec_bits_hamming_lua },
  { "bits_hamming_mask", tk_cvec_bits_hamming_mask_lua },
  { "bits_and", tk_cvec_bits_and_lua },
  { "bits_or", tk_cvec_bits_or_lua },
  { "bits_xor", tk_cvec_bits_xor_lua },
  { "bits_to_ascii", tk_cvec_bits_to_ascii_lua },
  { "bits_transpose", tk_cvec_bits_transpose_lua },
  { "bits_top_bns", tk_cvec_bits_top_bns_lua },
  { NULL, NULL }
};

static inline int tk_cvec_bits_bytes_lua (lua_State *L) {
  lua_settop(L, 1);
  uint64_t n_bits = tk_lua_checkunsigned(L, 1, "n_bits");
  lua_pushinteger(L, (lua_Integer) TK_CVEC_BITS_BYTES(n_bits));
  return 1;
}

static inline int tk_cvec_bits_byte_lua (lua_State *L) {
  lua_settop(L, 1);
  uint64_t n_bits = tk_lua_checkunsigned(L, 1, "n_bits");
  lua_pushinteger(L, (lua_Integer) TK_CVEC_BITS_BYTE(n_bits));
  return 1;
}

static inline int tk_cvec_bits_bit_lua (lua_State *L) {
  lua_settop(L, 1);
  uint64_t n_bits = tk_lua_checkunsigned(L, 1, "n_bits");
  lua_pushinteger(L, (lua_Integer) TK_CVEC_BITS_BIT(n_bits));
  return 1;
}

static luaL_Reg tk_cvec_lua_ext_fns[] =
{
  { "bits_from_ivec", tk_cvec_bits_from_ivec_lua },
  { "bits_bytes", tk_cvec_bits_bytes_lua },
  { "bits_byte", tk_cvec_bits_byte_lua },
  { "bits_bit", tk_cvec_bits_bit_lua },
  { NULL, NULL }
};

int luaopen_santoku_cvec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_cvec_lua_fns);
  luaL_register(L, NULL, tk_cvec_lua_ext_fns);
  tk_cvec_create(L, 0, 0, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_cvec_lua_mt_fns);
  luaL_register(L, NULL, tk_cvec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_cvec_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
