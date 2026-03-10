#ifndef TK_FVEC_H
#define TK_FVEC_H

#include <santoku/vec/base.h>

#define tk_vec_name tk_fvec
#define tk_vec_base float
#define tk_vec_pushbase(...) lua_pushnumber(__VA_ARGS__)
#define tk_vec_peekbase(...) ((float)luaL_checknumber(__VA_ARGS__))
#define tk_vec_abs(...) fabsf(__VA_ARGS__)
#include <santoku/vec/ext/tpl.h>

#include <santoku/fvec/ext.h>

#endif
