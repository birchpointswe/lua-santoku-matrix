#ifndef TK_PVEC_H
#define TK_PVEC_H

#include <santoku/vec/base.h>

#define tk_vec_name tk_pvec
#define tk_vec_base tk_pair_t
#define tk_vec_lt(a, b) ((a).p < (b).p)
#define tk_vec_gt(a, b) ((a).p > (b).p)
#define tk_vec_ltx(a,b)  ((a).i < (b).i || ((a).i == (b).i && (a).p < (b).p))
#define tk_vec_gtx(a,b)  ((a).i > (b).i || ((a).i == (b).i && (a).p > (b).p))
#define tk_vec_eqx(a, b) ((a).i == (b).i && (a).p == (b).p)
#define tk_vec_limited
#include <santoku/vec/ext/tpl.h>

#include <santoku/pvec/ext.h>

#endif
