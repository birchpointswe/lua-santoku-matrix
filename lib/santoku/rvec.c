#include <santoku/rvec.h>

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

static inline int tk_rvec_keys_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_rvec_keys(L, tk_rvec_peek(L, 1, "rvec"));
  return 1;
}

static inline int tk_rvec_values_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_rvec_values(L, tk_rvec_peek(L, 1, "rvec"));
  return 1;
}

static luaL_Reg tk_rvec_lua_mt_ext2_fns[] =
{
  { "get", tk_rvec_get_lua },
  { "set", tk_rvec_set_lua },
  { "keys", tk_rvec_keys_lua },
  { "values", tk_rvec_values_lua },
  { "each", tk_rvec_each0_lua },
  { "ieach", tk_rvec_ieach0_lua },
  { NULL, NULL }
};

int luaopen_santoku_rvec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_rvec_lua_fns);
  tk_rvec_create(L, 0, 0, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_rvec_lua_mt_fns);
  luaL_register(L, NULL, tk_rvec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_rvec_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
