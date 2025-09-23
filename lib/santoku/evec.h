#ifndef TK_EVEC_H
#define TK_EVEC_H

#include <santoku/vec/base.h>
#include <santoku/evec/base.h>

#define tk_vec_name tk_evec
#define tk_vec_base tk_edge_t
#define tk_vec_lt(a, b) ((a).w < (b).w)
#define tk_vec_gt(a, b) ((a).w > (b).w)
#define tk_vec_ltx(a,b)  ((a).u < (b).u || ((a).u == (b).u && ((a).v < (b).v || ((a).v == (b).v && (a).w < (b).w))))
#define tk_vec_gtx(a,b)  ((a).u > (b).u || ((a).u == (b).u && ((a).v > (b).v || ((a).v == (b).v && (a).w > (b).w))))
#define tk_vec_eqx(a, b) ((a).u == (b).u && (a).v == (b).v && (a).w == (b).w)
#define tk_vec_limited
#include <santoku/vec/ext/tpl.h>

#include <santoku/evec/ext.h>

#endif
