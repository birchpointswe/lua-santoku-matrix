#ifndef TK_RVEC_BASE_H
#define TK_RVEC_BASE_H

typedef struct { int64_t i; double d; } tk_rank_t;
#define tk_rank(a, b) ((tk_rank_t) { a, b })

#define tk_vec_name tk_rvec
#define tk_vec_base tk_rank_t
#define tk_vec_lt(a, b) ((a).d < (b).d)
#define tk_vec_gt(a, b) ((a).d > (b).d)
#define tk_vec_limited
#include <santoku/vec/tpl.h>

static inline void tk_rvec_hasc (tk_rvec_t *v, tk_rank_t r)
{
  if (v->n < v->m) {
    v->a[v->n ++] = r;
    if (v->n == v->m)
      ks_heapmake(tk_rvec_asc, v->m, v->a);
  } else if (r.d < v->a[0].d) {
    v->a[0] = r;
    ks_heapadjust(tk_rvec_asc, 0, v->m, v->a);
  }
}

static inline void tk_rvec_hdesc (tk_rvec_t *v, tk_rank_t r)
{
  if (v->n < v->m) {
    v->a[v->n ++] = r;
    if (v->n == v->m)
      ks_heapmake(tk_rvec_desc, v->m, v->a);
  } else if (r.d > v->a[0].d) {
    v->a[0] = r;
    ks_heapadjust(tk_rvec_desc, 0, v->m, v->a);
  }
}

#endif
