#ifndef TK_DUSET_BASE_H
#define TK_DUSET_BASE_H

#include <santoku/klib.h>
#include <santoku/dvec/base.h>

#define tk_umap_name tk_duset
#define tk_umap_key double
#define tk_umap_peekkey(...) tk_lua_checknumber(__VA_ARGS__)
#define tk_umap_pushkey(...) lua_pushnumber(__VA_ARGS__)
#define tk_umap_eq(a, b) ((a) == (b))
#define tk_umap_hash(a) (tk_hash_double(a))
#include <santoku/umap/tpl.h>

#endif