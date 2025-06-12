#define _GNU_SOURCE

#include <santoku/dvec.h>
#include <santoku/threads.h>

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

static luaL_Reg tk_dvec_lua_mt_ext2_fns[] =
{
  { "multiply_bits", tk_dvec_multiply_bits_lua },
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
