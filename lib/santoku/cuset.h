// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2024 Birch Point SWE
#ifndef TK_CUSET_H
#define TK_CUSET_H

#include <santoku/cvec/base.h>
#include <santoku/cuset/base.h>

#define tk_umap_name tk_cuset
#define tk_umap_key const char *
#define tk_umap_peekkey(...) tk_lua_checkstring(__VA_ARGS__)
#define tk_umap_pushkey(...) lua_pushstring(__VA_ARGS__)
#define tk_umap_eq(a, b) (strcmp((a), (b)) == 0)
#define tk_umap_hash(a) (kh_str_hash_func(a))
#define tk_umap_no_persist
#include <santoku/umap/ext/tpl.h>

#endif
