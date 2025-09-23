#include <santoku/iuset.h>
#include <santoku/cvec.h>

int luaopen_santoku_iuset (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_iuset_lua_fns);
  tk_iuset_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_iuset_lua_mt_fns);
  luaL_register(L, NULL, tk_iuset_lua_ext_fns);
  lua_pop(L, 2);
  return 1;
}
