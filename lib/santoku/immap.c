#include <santoku/immap.h>

#if !defined(__EMSCRIPTEN__)

int luaopen_santoku_immap (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_immap_lua_fns);
  tk_immap_create(L, "/dev/null", 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_immap_lua_mt_fns);
  luaL_register(L, NULL, tk_immap_lua_ext_fns);
  lua_pop(L, 2);
  return 1;
}

#endif
