#include <santoku/cvec.h>
#include <santoku/threads.h>

int luaopen_santoku_cvec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_cvec_lua_fns);
  tk_cvec_create(L, 0, 0, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_cvec_lua_mt_fns);
  luaL_register(L, NULL, tk_cvec_lua_mt_ext_fns);
  lua_pop(L, 2);
  return 1;
}
