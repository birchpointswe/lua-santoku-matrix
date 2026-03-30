#ifndef TK_SVEC_BASE_H
#define TK_SVEC_BASE_H

#include <santoku/lua/utils.h>
#include <santoku/klib.h>

#include <math.h>
#include <float.h>
#include <errno.h>
#include <stdlib.h>

#define tk_vec_name tk_svec
#define tk_vec_base int32_t
#define tk_vec_pushbase(...) lua_pushinteger(__VA_ARGS__)
#define tk_vec_peekbase(...) ((int32_t)luaL_checkinteger(__VA_ARGS__))
#define tk_vec_abs(...) abs(__VA_ARGS__)
#define tk_vec_module "santoku.svec"
#include <santoku/vec/tpl.h>

#endif
