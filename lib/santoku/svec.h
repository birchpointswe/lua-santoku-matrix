#ifndef TK_SVEC_H
#define TK_SVEC_H

#include <santoku/vec/base.h>

#define tk_vec_name tk_svec
#define tk_vec_base int32_t
#define tk_vec_pushbase(...) lua_pushinteger(__VA_ARGS__)
#define tk_vec_peekbase(...) ((int32_t)luaL_checkinteger(__VA_ARGS__))
#define tk_vec_abs(...) abs(__VA_ARGS__)
#include <santoku/vec/ext/tpl.h>

#include <santoku/svec/ext.h>

#endif
