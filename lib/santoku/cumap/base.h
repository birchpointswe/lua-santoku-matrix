#ifndef TK_CUMAP_BASE_H
#define TK_CUMAP_BASE_H

#include <santoku/klib.h>
#include <santoku/cvec/base.h>

#define tk_umap_name tk_cumap
#define tk_umap_key int64_t
#define tk_umap_value const char *
#define tk_umap_peekkey(...) tk_lua_checkinteger(__VA_ARGS__)
#define tk_umap_pushkey(...) lua_pushinteger(__VA_ARGS__)
#define tk_umap_peekvalue(...) tk_lua_checkstring(__VA_ARGS__)
#define tk_umap_pushvalue(...) lua_pushstring(__VA_ARGS__)
#define tk_umap_eq(a, b) ((a) == (b))
#define tk_umap_hash(a) (kh_int64_hash_func((uint64_t) a))
#include <santoku/umap/tpl.h>

#endif
