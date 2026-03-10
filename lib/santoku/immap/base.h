#ifndef TK_IMMAP_BASE_H
#define TK_IMMAP_BASE_H

#include <santoku/ivec/base.h>

#define tk_mmap_name tk_immap
#define tk_mmap_base int64_t
#define tk_mmap_pushbase(...) lua_pushinteger(__VA_ARGS__)
#define tk_mmap_peekbase(...) luaL_checkinteger(__VA_ARGS__)
#define tk_mmap_vec_name tk_ivec
#include <santoku/mmap/tpl.h>

#endif
