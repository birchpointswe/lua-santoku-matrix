#ifndef TK_FVEC_BASE_H
#define TK_FVEC_BASE_H

#include <santoku/lua/utils.h>
#include <santoku/klib.h>

#include <math.h>
#include <float.h>
#include <errno.h>

#define tk_vec_name tk_fvec
#define tk_vec_base float
#define tk_vec_pushbase(...) lua_pushnumber(__VA_ARGS__)
#define tk_vec_peekbase(...) ((float)luaL_checknumber(__VA_ARGS__))
#define tk_vec_abs(...) fabsf(__VA_ARGS__)
#include <santoku/vec/tpl.h>

#endif
