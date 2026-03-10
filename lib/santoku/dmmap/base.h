#ifndef TK_DMMAP_BASE_H
#define TK_DMMAP_BASE_H

#include <santoku/dvec/base.h>

#define tk_mmap_name tk_dmmap
#define tk_mmap_base double
#define tk_mmap_pushbase(...) lua_pushnumber(__VA_ARGS__)
#define tk_mmap_peekbase(...) luaL_checknumber(__VA_ARGS__)
#define tk_mmap_vec_name tk_dvec
#include <santoku/mmap/tpl.h>

#endif
