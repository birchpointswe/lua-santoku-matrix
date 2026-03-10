#include <santoku/dmmap.h>

#if !defined(__EMSCRIPTEN__)

int luaopen_santoku_dmmap (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_dmmap_lua_fns);
  tk_dmmap_create(L, "/dev/null", 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_dmmap_lua_mt_fns);
  luaL_register(L, NULL, tk_dmmap_lua_ext_fns);
  lua_pop(L, 2);
  return 1;
}

#endif
