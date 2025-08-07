#ifndef TK_RVEC_H
#define TK_RVEC_H

#include <santoku/vec/base.h>

#define tk_vec_name tk_rvec
#define tk_vec_base tk_rank_t
#define tk_vec_lt(a, b) ((a).d < (b).d)
#define tk_vec_gt(a, b) ((a).d > (b).d)
#define tk_vec_ltx(a,b)  ((a).i < (b).i || ((a).i == (b).i && (a).d < (b).d))
#define tk_vec_gtx(a,b)  ((a).i > (b).i || ((a).i == (b).i && (a).d > (b).d))
#define tk_vec_eqx(a, b) ((a).i == (b).i && (a).d == (b).d)
#define tk_vec_limited
#include <santoku/vec/ext/tpl.h>

#include <santoku/rvec/ext.h>

#endif
