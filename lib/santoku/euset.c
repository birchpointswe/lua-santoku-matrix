#include <santoku/euset.h>
#include <santoku/cvec.h>

int luaopen_santoku_euset (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_euset_lua_fns);
  tk_euset_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_euset_lua_mt_fns);
  luaL_register(L, NULL, tk_euset_lua_ext_fns);
  lua_pop(L, 2);
  return 1;
}
