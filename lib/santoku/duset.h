// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2024 Birch Point SWE
#ifndef TK_DUSET_H
#define TK_DUSET_H

#include <santoku/cvec/base.h>
#include <santoku/duset/base.h>

#define tk_umap_name tk_duset
#define tk_umap_key double
#define tk_umap_peekkey(...) tk_lua_checknumber(__VA_ARGS__)
#define tk_umap_pushkey(...) lua_pushnumber(__VA_ARGS__)
#define tk_umap_eq(a, b) ((a) == (b))
#define tk_umap_hash(a) (tk_hash_double(a))
#include <santoku/umap/ext/tpl.h>

#include <santoku/dvec.h>

#endif
