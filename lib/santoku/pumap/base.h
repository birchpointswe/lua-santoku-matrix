#ifndef TK_PUMAP_BASE_H
#define TK_PUMAP_BASE_H

#include <santoku/klib.h>
#include <santoku/pvec/base.h>

#define tk_umap_name tk_pumap
#define tk_umap_key int64_t
#define tk_umap_value tk_pair_t
#define tk_umap_peekkey(...) tk_lua_checkinteger(__VA_ARGS__)
#define tk_umap_pushkey(...) lua_pushinteger(__VA_ARGS__)
#define tk_umap_eq(a, b) ((a) == (b))
#define tk_umap_hash(a) (kh_int64_hash_func((uint64_t) a))
#include <santoku/umap/tpl.h>

#endif