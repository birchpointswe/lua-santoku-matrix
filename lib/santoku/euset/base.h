#ifndef TK_EUSET_BASE_H
#define TK_EUSET_BASE_H

#include <santoku/cvec/base.h>
#include <santoku/evec/base.h>

#define tk_edge_eq(a, b) ((a).u == (b).u && ((a).v == (b).v))
#define tk_edge_hash(a) (tk_hash_128((uint64_t) kh_int64_hash_func((a).u), (uint64_t) kh_int64_hash_func((a).v)))

#define tk_umap_name tk_euset
#define tk_umap_key tk_edge_t
#define tk_umap_eq(a, b) tk_edge_eq(a, b)
#define tk_umap_hash(a) tk_edge_hash(a)
#include <santoku/umap/tpl.h>

#endif
