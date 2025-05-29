#include <santoku/ivec/base.h>
#include <santoku/rvec/base.h>
#include <santoku/dvec/base.h>

#define tk_vec_name tk_dvec
#define tk_vec_base double
#define tk_vec_pushbase(...) lua_pushnumber(__VA_ARGS__)
#define tk_vec_peekbase(...) luaL_checknumber(__VA_ARGS__)
#define tk_vec_lua
#include <santoku/vec.ext.template.h>

int luaopen_santoku_dvec_ext (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_dvec_lua_fns);
  return 1;
}
