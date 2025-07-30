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

static inline void tk_rvec_hmax (tk_rvec_t *v, tk_rank_t r)
{
  if (v->n < v->m) {
    v->a[v->n ++] = r;
    if (v->n == v->m)
      ks_heapmake(tk_rvec_desc, v->m, v->a);
  } else if (r.d < v->a[0].d) {
    v->a[0] = r;
    ks_heapadjust(tk_rvec_desc, 0, v->m, v->a);
  }
}

static inline void tk_rvec_hmin (tk_rvec_t *v, tk_rank_t r)
{
  if (v->n < v->m) {
    v->a[v->n ++] = r;
    if (v->n == v->m)
      ks_heapmake(tk_rvec_asc, v->m, v->a);
  } else if (r.d > v->a[0].d) {
    v->a[0] = r;
    ks_heapadjust(tk_rvec_asc, 0, v->m, v->a);
  }
}

static inline void tk_rvec_hmax_init (tk_rvec_t *v) {
  ks_heapmake(tk_rvec_desc, v->n, v->a);
}

static inline void tk_rvec_hmin_init (tk_rvec_t *v) {
  ks_heapmake(tk_rvec_asc, v->n, v->a);
}

static inline tk_rank_t tk_rvec_hmax_pop (tk_rvec_t *v) {
  tk_rank_t top = v->a[0];
  v->a[0] = v->a[--v->n];
  ks_heapadjust(tk_rvec_desc, 0, v->n, v->a);
  return top;
}

static inline tk_rank_t tk_rvec_hmin_pop (tk_rvec_t *v) {
  tk_rank_t top = v->a[0];
  v->a[0] = v->a[--v->n];
  ks_heapadjust(tk_rvec_asc, 0, v->n, v->a);
  return top;
}

static inline void tk_rvec_hmax_push (tk_rvec_t *v, tk_rank_t r) {
  size_t i = v->n;
  tk_rvec_push(v, r);
  while (i && tk_rvec_lt(v->a[i], v->a[(i - 1) >> 1])) {
    size_t p = (i - 1) >> 1;
    tk_rank_t tmp = v->a[i];
    v->a[i] = v->a[p];
    v->a[p] = tmp;
    i = p;
  }
}

static inline void tk_rvec_hmin_push (tk_rvec_t *v, tk_rank_t r) {
  size_t i = v->n;
  tk_rvec_push(v, r);
  while (i && tk_rvec_gt(v->a[i], v->a[(i - 1) >> 1])) {
    size_t p = (i - 1) >> 1;
    tk_rank_t tmp = v->a[i];
    v->a[i] = v->a[p];
    v->a[p] = tmp;
    i = p;
  }
}

#endif
