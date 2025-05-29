#ifndef TK_CVEC_H
#define TK_CVEC_H

#include <santoku/ivec/base.h>
#include <santoku/rvec/base.h>
#include <santoku/dvec/base.h>

#include <santoku/cvec/ext.h>

#define tk_vec_name tk_cvec
#define tk_vec_base char
#define tk_vec_limited
#include <santoku/vec.ext.template.h>

#endif

