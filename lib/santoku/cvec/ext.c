#include <santoku/ivec/base.h>
#include <santoku/rvec/base.h>
#include <santoku/dvec/base.h>

#define tk_vec_name tk_cvec
#define tk_vec_base double
#define tk_vec_luaopen luaopen_santoku_cvec_ext
#define tk_vec_limited
#define tk_vec_lua
#include <santoku/vec.ext.template.h>

int luaopen_santoku_cvec_ext (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_cvec_lua_fns);
  return 1;
}
