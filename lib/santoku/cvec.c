#include <santoku/cvec.h>
#include <santoku/ivec.h>
#include <santoku/iuset.h>

static inline int tk_cvec_bits_flip_interleave_lua (lua_State *L) {
  lua_settop(L, 2);
  tk_cvec_t *v = tk_cvec_peek(L, 1, "cvec");
  uint64_t n_features = tk_lua_checkunsigned(L, 2, "n_features");
  tk_cvec_bits_flip_interleave(v, n_features);
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

static inline int tk_cvec_bits_rearrange_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  tk_ivec_t *ids = tk_ivec_peek(L, 2, "ids");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  tk_cvec_bits_rearrange(bitmap, ids, n_features);
  return 0;
}

static inline int tk_cvec_bits_extend_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_cvec_t *base = tk_cvec_peek(L, 1, "base");
  tk_cvec_t *ext = tk_cvec_peek(L, 2, "ext");
  uint64_t n_base_features = tk_lua_checkunsigned(L, 3, "n_base_features");
  uint64_t n_ext_features = tk_lua_checkunsigned(L, 4, "n_ext_features");
  tk_cvec_bits_extend(base, ext, n_base_features, n_ext_features);
  return 0;
}

static inline int tk_cvec_bits_filter_lua (lua_State *L) {
  int n_args = lua_gettop(L);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  tk_ivec_t *selected_features = lua_isnil(L, 2) ? NULL : tk_ivec_peek(L, 2, "selected_features");
  tk_ivec_t *sample_ids = NULL;
  uint64_t n_features;

  if (n_args == 4) {

    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 4, "n_features");
  } else {

    n_features = tk_lua_checkunsigned(L, 3, "n_features");
  }

  tk_cvec_bits_filter(bitmap, selected_features, sample_ids, n_features);
  return 0;
}

static inline int tk_cvec_bits_copy_lua (lua_State *L) {
  int n_args = lua_gettop(L);
  tk_cvec_t *dest = tk_cvec_peek(L, 1, "dest");
  tk_cvec_t *src_bitmap = tk_cvec_peek(L, 2, "src_bitmap");
  tk_ivec_t *selected_features = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "selected_features");
  tk_ivec_t *sample_ids = NULL;
  uint64_t n_features;

  if (n_args == 5) {

    sample_ids = lua_isnil(L, 4) ? NULL : tk_ivec_peek(L, 4, "sample_ids");
    n_features = tk_lua_checkunsigned(L, 5, "n_features");
  } else {

    n_features = tk_lua_checkunsigned(L, 4, "n_features");
  }

  tk_cvec_bits_copy(dest, src_bitmap, selected_features, sample_ids, n_features);
  return 0;
}

static inline int tk_cvec_bits_score_chi2_lua (lua_State *L) {
  lua_settop(L, 6);
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
  unsigned int n_threads = tk_threads_getn(L, 6, "n_threads", NULL);
  tk_cvec_bits_score_chi2(L, bitmap, codes, labels, n_samples, n_features, n_hidden, n_threads);
  return 1;
}

static inline int tk_cvec_bits_score_mi_lua (lua_State *L) {
  lua_settop(L, 6);
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
  unsigned int n_threads = tk_threads_getn(L, 6, "n_threads", NULL);
  tk_cvec_bits_score_mi(L, bitmap, codes, labels, n_samples, n_features, n_hidden, n_threads);
  return 1;
}

static inline int tk_cvec_bits_score_entropy_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_cvec_t *codes = tk_cvec_peek(L, 1, "codes");
  unsigned int n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  unsigned int n_hidden = tk_lua_checkunsigned(L, 3, "n_hidden");
  unsigned int n_threads = tk_threads_getn(L, 4, "n_threads", NULL);
  tk_cvec_bits_score_entropy(L, codes, n_samples, n_hidden, n_threads);
  return 1;
}

static inline int tk_cvec_bits_top_chi2_lua (lua_State *L) {
  lua_settop(L, 7);
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
  uint64_t top_k = tk_lua_checkunsigned(L, 6, "top_k");
  unsigned int n_threads = tk_threads_getn(L, 7, "n_threads", NULL);
  tk_cvec_bits_top_chi2(L, bitmap, codes, labels, n_samples, n_features, n_hidden, top_k, n_threads);
  return 2;
}

static inline int tk_cvec_bits_top_mi_lua (lua_State *L) {
  lua_settop(L, 7);
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
  uint64_t top_k = tk_lua_checkunsigned(L, 6, "top_k");
  unsigned int n_threads = tk_threads_getn(L, 7, "n_threads", NULL);
  tk_cvec_bits_top_mi(L, bitmap, codes, labels, n_samples, n_features, n_hidden, top_k, n_threads);
  return 2;
}

static inline int tk_cvec_bits_top_entropy_lua (lua_State *L) {
  lua_settop(L, 5);
  tk_cvec_t *codes = tk_cvec_peek(L, 1, "codes");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 3, "n_hidden");
  uint64_t top_k = tk_lua_checkunsigned(L, 4, "top_k");
  unsigned int n_threads = tk_threads_getn(L, 5, "n_threads", NULL);
  tk_cvec_bits_top_entropy(L, codes, n_samples, n_hidden, top_k, n_threads);
  return 1;
}

static luaL_Reg tk_cvec_lua_mt_ext2_fns[] =
{
  { "bits_flip_interleave", tk_cvec_bits_flip_interleave_lua },
  { "bits_to_ivec", tk_cvec_bits_to_ivec_lua },
  { "bits_rearrange", tk_cvec_bits_rearrange_lua },
  { "bits_extend", tk_cvec_bits_extend_lua },
  { "bits_filter", tk_cvec_bits_filter_lua },
  { "bits_copy", tk_cvec_bits_copy_lua },
  { "bits_score_chi2", tk_cvec_bits_score_chi2_lua },
  { "bits_score_mi", tk_cvec_bits_score_mi_lua },
  { "bits_score_entropy", tk_cvec_bits_score_entropy_lua },
  { "bits_top_chi2", tk_cvec_bits_top_chi2_lua },
  { "bits_top_mi", tk_cvec_bits_top_mi_lua },
  { "bits_top_entropy", tk_cvec_bits_top_entropy_lua },
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
