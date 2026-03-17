#ifndef TK_IVEC_BASE_H
#define TK_IVEC_BASE_H

#include <santoku/lua/utils.h>
#include <santoku/klib.h>

#include <math.h>
#include <float.h>
#include <errno.h>

#define tk_vec_name tk_ivec
#define tk_vec_base int64_t
#define tk_vec_pushbase(...) lua_pushinteger(__VA_ARGS__)
#define tk_vec_peekbase(...) luaL_checkinteger(__VA_ARGS__)
#define tk_vec_abs(...) llabs(__VA_ARGS__)
#include <santoku/vec/tpl.h>

#endif
