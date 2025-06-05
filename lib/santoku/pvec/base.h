#ifndef TK_PVEC_BASE_H
#define TK_PVEC_BASE_H

typedef struct { int64_t i; int64_t p; } tk_pair_t;
#define tk_pair(a, b) ((tk_pair_t) { a, b })

#define tk_vec_name tk_pvec
#define tk_vec_base tk_pair_t
#define tk_vec_lt(a, b) ((a).p < (b).p)
#define tk_vec_gt(a, b) ((a).p > (b).p)
#define tk_vec_limited
#include <santoku/vec/tpl.h>

static inline void tk_pvec_hasc (tk_pvec_t *v, tk_pair_t r)
{
  if (v->n < v->m) {
    v->a[v->n ++] = r;
    if (v->n == v->m)
      ks_heapmake(tk_pvec_asc, v->m, v->a);
  } else if (r.p > v->a[0].p) {
    v->a[0] = r;
    ks_heapadjust(tk_pvec_asc, 0, v->m, v->a);
  }
}

static inline void tk_pvec_hdesc (tk_pvec_t *v, tk_pair_t r)
{
  if (v->n < v->m) {
    v->a[v->n ++] = r;
    if (v->n == v->m)
      ks_heapmake(tk_pvec_desc, v->m, v->a);
  } else if (r.p < v->a[0].p) {
    v->a[0] = r;
    ks_heapadjust(tk_pvec_desc, 0, v->m, v->a);
  }
}

#endif

