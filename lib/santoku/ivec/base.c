#include <santoku/lua/utils.h>
#include <santoku/klib.h>

#include <math.h>
#include <float.h>
#include <errno.h>

#define tk_vec_name tk_ivec
#define tk_vec_base int64_t
#define tk_vec_pushbase(...) lua_pushinteger(__VA_ARGS__)
#define tk_vec_peekbase(...) luaL_checkinteger(__VA_ARGS__)
#define tk_vec_lua
#include <santoku/vec.template.h>

int luaopen_santoku_ivec_base (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_ivec_lua_fns);
  return 1;
}
