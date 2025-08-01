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

static inline int tk_pvec_keys_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_pvec_keys(L, tk_pvec_peek(L, 1, "pvec"));
  return 1;
}

static inline int tk_pvec_values_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_pvec_values(L, tk_pvec_peek(L, 1, "pvec"));
  return 1;
}

static luaL_Reg tk_pvec_lua_mt_ext2_fns[] =
{
  { "get", tk_pvec_get_lua },
  { "set", tk_pvec_set_lua },
  { "push", tk_pvec_push_lua },
  { "keys", tk_pvec_keys_lua },
  { "values", tk_pvec_values_lua },
  { "each", tk_pvec_each0_lua },
  { "ieach", tk_pvec_ieach0_lua },
  { NULL, NULL }
};

int luaopen_santoku_pvec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_pvec_lua_fns);
  tk_pvec_create(L, 0, 0, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_pvec_lua_mt_fns);
  luaL_register(L, NULL, tk_pvec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_pvec_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
