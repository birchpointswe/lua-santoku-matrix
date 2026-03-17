#ifndef TK_DVEC_BASE_H
#define TK_DVEC_BASE_H

#include <santoku/lua/utils.h>
#include <santoku/klib.h>

#include <math.h>
#include <float.h>
#include <errno.h>

#define tk_vec_name tk_dvec
#define tk_vec_base double
#define tk_vec_pushbase(...) lua_pushnumber(__VA_ARGS__)
#define tk_vec_peekbase(...) luaL_checknumber(__VA_ARGS__)
#define tk_vec_abs(...) fabs(__VA_ARGS__)
#ifdef TK_DVEC_INIT
extern void tk_dvec_init_mt(lua_State *L);
#define tk_vec_init tk_dvec_init_mt
#endif
#include <santoku/vec/tpl.h>

#endif
