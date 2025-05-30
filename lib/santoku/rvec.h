#ifndef TK_RVEC_H
#define TK_RVEC_H

#include <santoku/ivec/base.h>
#include <santoku/rvec/base.h>
#include <santoku/dvec/base.h>

#define tk_vec_name tk_rvec
#define tk_vec_base tk_rank_t
#define tk_vec_lt(a, b) ((a).d < (b).d)
#define tk_vec_gt(a, b) ((a).d > (b).d)
#define tk_vec_limited
#include <santoku/vec.ext.template.h>

#endif
