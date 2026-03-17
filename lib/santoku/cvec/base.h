#ifndef TK_CVEC_BASE_H
#define TK_CVEC_BASE_H

#include <santoku/lua/utils.h>
#include <santoku/klib.h>

#include <math.h>
#include <float.h>
#include <errno.h>

#define tk_vec_name tk_cvec
#define tk_vec_base char
#define tk_vec_limited
#define tk_vec_module "santoku.cvec"
#include <santoku/vec/tpl.h>

#endif
