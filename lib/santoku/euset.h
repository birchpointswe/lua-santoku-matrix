#ifndef TK_EUSET_H
#define TK_EUSET_H

#include <santoku/cvec/base.h>
#include <santoku/euset/base.h>

#define tk_umap_name tk_euset
#define tk_umap_key tk_edge_t
#define tk_umap_eq(a, b) tk_edge_eq(a, b)
#define tk_umap_hash(a) tk_edge_hash(a)
#include <santoku/umap/ext/tpl.h>

#include <santoku/euset/ext.h>

#endif
