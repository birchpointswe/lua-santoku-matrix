#include <santoku/iuset.h>
#include <santoku/cvec.h>
#include <santoku/ivec.h>
#include <stdbool.h>

















static luaL_Reg tk_cvec_lua_mt_ext2_fns[] =
{
  { NULL, NULL }
};




static luaL_Reg tk_cvec_lua_ext_fns[] =
{
  { NULL, NULL }
};

int luaopen_santoku_cvec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_cvec_lua_fns);
  luaL_register(L, NULL, tk_cvec_lua_ext_fns);
  tk_cvec_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_cvec_lua_mt_fns);
  luaL_register(L, NULL, tk_cvec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_cvec_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
