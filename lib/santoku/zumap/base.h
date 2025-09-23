#ifndef TK_ZUMAP_BASE_H
#define TK_ZUMAP_BASE_H

#include <santoku/klib.h>
#include <santoku/ivec/base.h>

#define tk_umap_name tk_zumap
#define tk_umap_key const char *
#define tk_umap_value int64_t
#define tk_umap_peekkey(...) tk_lua_checkstring(__VA_ARGS__)
#define tk_umap_peekvalue(...) tk_lua_checkinteger(__VA_ARGS__)
#define tk_umap_pushkey(...) lua_pushstring(__VA_ARGS__)
#define tk_umap_pushvalue(...) lua_pushinteger(__VA_ARGS__)
#define tk_umap_eq(a, b) (strcmp((a), (b)) == 0)
#define tk_umap_hash(a) (kh_str_hash_func(a))
#include <santoku/umap/tpl.h>

#endif
