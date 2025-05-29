#include <santoku/lua/utils.h>
#include <santoku/klib.h>

#include <math.h>
#include <float.h>
#include <errno.h>

#define tk_vec_name tk_dvec
#define tk_vec_base double
#define tk_vec_pushbase(...) lua_pushnumber(__VA_ARGS__)
#define tk_vec_peekbase(...) luaL_checknumber(__VA_ARGS__)
#define tk_vec_lua
#include <santoku/vec.template.h>

int luaopen_santoku_dvec_base (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_dvec_lua_fns);
  return 1;
}
