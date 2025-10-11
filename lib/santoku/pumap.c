#include <santoku/iuset.h>
#include <santoku/pumap.h>
#include <santoku/cvec.h>

static inline int tk_pumap_setval_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_pumap_t *h = tk_pumap_peek(L, 1, "pumap");
  uint32_t khi = tk_lua_checkunsigned(L, 2, "i");
  int64_t i = tk_lua_checkinteger(L, 3, "i");
  int64_t p = tk_lua_checkinteger(L, 4, "p");
  tk_pumap_setval(h, khi, (tk_pair_t) { i, p });
  return 0;
}

static inline int tk_pumap_val_lua (lua_State *L) {
  lua_settop(L, 2);
  tk_pumap_t *h = tk_pumap_peek(L, 1, "pumap");
  uint32_t khi = tk_lua_checkunsigned(L, 2, "i");
  tk_pair_t p = tk_pumap_val(h, khi);
  lua_pushinteger(L, p.i);
  lua_pushinteger(L, p.p);
  return 2;
}

static luaL_Reg tk_pumap_lua_mt_ext2_fns[] =
{
  { "setval", tk_pumap_setval_lua },
  { "val", tk_pumap_val_lua },
  { NULL, NULL }
};

int luaopen_santoku_pumap (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_pumap_lua_fns);
  tk_pumap_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_pumap_lua_mt_fns);
  luaL_register(L, NULL, tk_pumap_lua_ext_fns);
  luaL_register(L, NULL, tk_pumap_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
