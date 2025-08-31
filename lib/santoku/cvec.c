#include <santoku/cvec.h>
#include <santoku/cvec/ext.h>
#include <santoku/ivec.h>
#include <santoku/threads.h>

static inline int tk_cvec_bits_flip_interleave_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_cvec_t *v = tk_cvec_peek(L, 1, "cvec");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  tk_cvec_bits_flip_interleave(v, n_samples, n_features);
  return 0;
}

static inline int tk_cvec_bits_to_ivec_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  tk_cvec_bits_to_ivec(L, bitmap, n_samples, n_features);
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
  lua_settop(L, 5);
  tk_cvec_t *base = tk_cvec_peek(L, 1, "base");
  tk_cvec_t *ext = tk_cvec_peek(L, 2, "ext");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_base_features = tk_lua_checkunsigned(L, 4, "n_base_features");
  uint64_t n_ext_features = tk_lua_checkunsigned(L, 5, "n_ext_features");
  tk_cvec_bits_extend(base, ext, n_samples, n_base_features, n_ext_features);
  return 0;
}

static inline int tk_cvec_bits_filter_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  tk_ivec_t *selected_features = tk_ivec_peek(L, 2, "selected_features");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  tk_cvec_bits_filter(bitmap, selected_features, n_samples, n_features);
  return 0;
}

static luaL_Reg tk_cvec_lua_mt_ext2_fns[] =
{
  { "bits_flip_interleave", tk_cvec_bits_flip_interleave_lua },
  { "bits_to_ivec", tk_cvec_bits_to_ivec_lua },
  { "bits_rearrange", tk_cvec_bits_rearrange_lua },
  { "bits_extend", tk_cvec_bits_extend_lua },
  { "bits_filter", tk_cvec_bits_filter_lua },
  { NULL, NULL }
};

static luaL_Reg tk_cvec_lua_ext_fns[] =
{
  { "bits_from_ivec", tk_cvec_bits_from_ivec_lua },
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
