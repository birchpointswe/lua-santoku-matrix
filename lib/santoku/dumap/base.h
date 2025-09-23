#ifndef TK_DUMAP_BASE_H
#define TK_DUMAP_BASE_H

#include <santoku/klib.h>
#include <santoku/dvec/base.h>

#define tk_umap_name tk_dumap
#define tk_umap_key int64_t
#define tk_umap_value double
#define tk_umap_peekkey(...) tk_lua_checkinteger(__VA_ARGS__)
#define tk_umap_peekvalue(...) tk_lua_checknumber(__VA_ARGS__)
#define tk_umap_pushkey(...) lua_pushinteger(__VA_ARGS__)
#define tk_umap_pushvalue(...) lua_pushnumber(__VA_ARGS__)
#define tk_umap_eq(a, b) ((a) == (b))
#define tk_umap_hash(a) (kh_int64_hash_func((uint64_t) a))
#include <santoku/umap/tpl.h>

#endif
