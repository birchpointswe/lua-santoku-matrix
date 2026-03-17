#ifndef TK_EVEC_BASE_H
#define TK_EVEC_BASE_H

typedef struct { int64_t u, v; double w; } tk_edge_t;
#define tk_edge(u, v, w) (((u) < (v)) ? ((tk_edge_t) { (u), (v), (w) }) : ((tk_edge_t) { (v), (u), (w) }))

#define tk_vec_name tk_evec
#define tk_vec_base tk_edge_t
#define tk_vec_lt(a, b) ((a).w < (b).w)
#define tk_vec_gt(a, b) ((a).w > (b).w)
#define tk_vec_ltx(a,b)  ((a).u < (b).u || ((a).u == (b).u && ((a).v < (b).v || ((a).v == (b).v && (a).w < (b).w))))
#define tk_vec_gtx(a,b)  ((a).u > (b).u || ((a).u == (b).u && ((a).v > (b).v || ((a).v == (b).v && (a).w > (b).w))))
#define tk_vec_eqx(a, b) ((a).u == (b).u && (a).v == (b).v && (a).w == (b).w)
#define tk_vec_limited
#include <santoku/vec/tpl.h>

static inline void tk_evec_hmax (tk_evec_t *v, size_t k, tk_edge_t r)
{
  if (v->n < k) {
    tk_evec_push(v, r);
    if (v->n == k)
      ks_heapmake(tk_evec_asc, k, v->a);
  } else if (r.w < v->a[0].w) {
    v->a[0] = r;
    ks_heapadjust(tk_evec_asc, 0, k, v->a);
  }
}

static inline void tk_evec_hmin (tk_evec_t *v, size_t k, tk_edge_t r)
{
  if (v->n < k) {
    tk_evec_push(v, r);
    if (v->n == k)
      ks_heapmake(tk_evec_desc, k, v->a);
  } else if (r.w > v->a[0].w) {
    v->a[0] = r;
    ks_heapadjust(tk_evec_desc, 0, k, v->a);
  }
}

static inline void tk_evec_hmax_init (tk_evec_t *v) {
  ks_heapmake(tk_evec_desc, v->n, v->a);
}

static inline void tk_evec_hmin_init (tk_evec_t *v) {
  ks_heapmake(tk_evec_asc, v->n, v->a);
}

static inline tk_edge_t tk_evec_hmax_pop (tk_evec_t *v) {
  tk_edge_t top = v->a[0];
  v->a[0] = v->a[--v->n];
  ks_heapadjust(tk_evec_desc, 0, v->n, v->a);
  return top;
}

static inline tk_edge_t tk_evec_hmin_pop (tk_evec_t *v) {
  tk_edge_t top = v->a[0];
  v->a[0] = v->a[--v->n];
  ks_heapadjust(tk_evec_asc, 0, v->n, v->a);
  return top;
}

static inline void tk_evec_hmax_push (tk_evec_t *v, tk_edge_t r) {
  size_t i = v->n;
  tk_evec_push(v, r);
  while (i && tk_evec_lt(v->a[i], v->a[(i - 1) >> 1])) {
    size_t p = (i - 1) >> 1;
    tk_edge_t tmp = v->a[i];
    v->a[i] = v->a[p];
    v->a[p] = tmp;
    i = p;
  }
}

static inline void tk_evec_hmin_push (tk_evec_t *v, tk_edge_t r) {
  size_t i = v->n;
  tk_evec_push(v, r);
  while (i && tk_evec_gt(v->a[i], v->a[(i - 1) >> 1])) {
    size_t p = (i - 1) >> 1;
    tk_edge_t tmp = v->a[i];
    v->a[i] = v->a[p];
    v->a[p] = tmp;
    i = p;
  }
}

#endif
