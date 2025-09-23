#ifndef TK_CUMAP_H
#define TK_CUMAP_H

#include <santoku/cvec/base.h>
#include <santoku/cumap/base.h>

#define tk_umap_name tk_cumap
#define tk_umap_key int64_t
#define tk_umap_value const char *
#define tk_umap_peekkey(...) tk_lua_checkinteger(__VA_ARGS__)
#define tk_umap_pushkey(...) lua_pushinteger(__VA_ARGS__)
#define tk_umap_peekvalue(...) tk_lua_checkstring(__VA_ARGS__)
#define tk_umap_pushvalue(...) lua_pushstring(__VA_ARGS__)
#define tk_umap_eq(a, b) ((a) == (b))
#define tk_umap_hash(a) (kh_int64_hash_func((uint64_t) a))
#define tk_umap_no_persist
#include <santoku/umap/ext/tpl.h>

#include <santoku/cumap/ext.h>

#endif
