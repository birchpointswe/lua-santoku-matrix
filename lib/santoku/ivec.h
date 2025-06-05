#ifndef TK_IVEC_H
#define TK_IVEC_H

#include <santoku/vec/base.h>

#define tk_vec_name tk_ivec
#define tk_vec_base int64_t
#define tk_vec_pushbase(...) lua_pushinteger(__VA_ARGS__)
#define tk_vec_peekbase(...) luaL_checkinteger(__VA_ARGS__)
#define tk_vec_abs(...) fabs(__VA_ARGS__)
#include <santoku/vec/ext/tpl.h>

#include <santoku/ivec/ext.h>

#endif
