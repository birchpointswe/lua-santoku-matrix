#ifndef TK_CVEC_BASE_H
#define TK_CVEC_BASE_H

#include <santoku/lua/utils.h>
#include <santoku/klib.h>

#include <math.h>
#include <float.h>
#include <errno.h>

#define tk_vec_name tk_cvec
#define tk_vec_base char
#define tk_vec_limited
#ifdef TK_CVEC_INIT
extern void tk_cvec_init_mt(lua_State *L);
#define tk_vec_init tk_cvec_init_mt
#endif
#include <santoku/vec/tpl.h>

#endif
