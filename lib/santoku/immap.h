#ifndef TK_IMMAP_H
#define TK_IMMAP_H

#include <santoku/immap/base.h>

#define tk_mmap_name tk_immap
#define tk_mmap_base int64_t
#define tk_mmap_pushbase(...) lua_pushinteger(__VA_ARGS__)
#define tk_mmap_peekbase(...) luaL_checkinteger(__VA_ARGS__)
#define tk_mmap_vec_name tk_ivec
#include <santoku/mmap/ext/tpl.h>

#endif
