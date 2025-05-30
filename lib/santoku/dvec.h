#ifndef TK_DVEC_H
#define TK_DVEC_H

#include <santoku/ivec/base.h>
#include <santoku/rvec/base.h>
#include <santoku/dvec/base.h>

#define tk_vec_name tk_dvec
#define tk_vec_base double
#define tk_vec_pushbase(...) lua_pushnumber(__VA_ARGS__)
#define tk_vec_peekbase(...) luaL_checknumber(__VA_ARGS__)
#define tk_vec_abs(...) fabs(__VA_ARGS__)
#include <santoku/vec.ext.template.h>

#endif
