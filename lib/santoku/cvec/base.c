#include <santoku/lua/utils.h>
#include <santoku/klib.h>

#include <math.h>
#include <float.h>
#include <errno.h>

#define tk_vec_name tk_cvec
#define tk_vec_base char
#define tk_vec_limited
#define tk_vec_lua
#include <santoku/vec.template.h>

int luaopen_santoku_cvec_base (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_cvec_lua_fns);
  return 1;
}
