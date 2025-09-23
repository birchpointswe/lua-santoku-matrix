#include <santoku/iumap.h>
#include <santoku/cvec.h>

int luaopen_santoku_iumap (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_iumap_lua_fns);
  tk_iumap_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_iumap_lua_mt_fns);
  luaL_register(L, NULL, tk_iumap_lua_ext_fns);
  lua_pop(L, 2);
  return 1;
}
