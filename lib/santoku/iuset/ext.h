#ifndef TK_IUSET_EXT_H
#define TK_IUSET_EXT_H

#include <santoku/ivec/base.h>
#include <santoku/iuset/base.h>

static inline void tk_iuset_union (tk_iuset_t *a, tk_iuset_t *b)
{
  int kha;
  int64_t x;
  tk_umap_foreach_keys(b, x, ({
    tk_iuset_put(a, x, &kha);
  }));
}

static inline double tk_iuset_jaccard (tk_iuset_t *a, tk_iuset_t *b)
{
  uint64_t intersection = tk_iuset_size(a);
  uint64_t union_count = 0;
  union_count = tk_iuset_size(a) + tk_iuset_size(b) - intersection;
  if (union_count == 0)
    return 0.0;
  return (double) intersection / (double) union_count;
}

static inline void tk_iuset_intersect (tk_iuset_t *a, tk_iuset_t *b)
{
  uint32_t i;
  int64_t k;
  tk_umap_foreach_iters(a, i, ({
    k = tk_iuset_key(a, i);
    if (!tk_iuset_contains(b, k))
      tk_iuset_del(a, i);
  }))
}

static inline void tk_iuset_subtract (tk_iuset_t *a, tk_iuset_t *b)
{
  uint32_t i;
  int64_t x;
  tk_umap_foreach_keys(b, x, ({
    i = tk_iuset_get(a, x);
    if (tk_iuset_exist(a, i))
      tk_iuset_del(a, i);
  }))
}

static inline tk_iuset_t *tk_iuset_from_ivec (lua_State *L, tk_ivec_t *v)
{
  int kha;
  tk_iuset_t *s = tk_iuset_create(L, v->n);
  if (!s)
    return NULL;
  for (uint64_t i = 0; i < v->n; i ++) {
    tk_iuset_put(s, v->a[i], &kha);
    if (kha < 0) {
      tk_iuset_destroy(s);
      return NULL;
    }
  }
  return s;
}

static inline tk_ivec_t *tk_iuset_keys (lua_State *L, tk_iuset_t *S)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_iuset_size(S));
  int64_t k;
  out->n = 0;
  tk_umap_foreach_keys(S, k, ({
    out->a[out->n ++] = k;
  }));
  return out;
}

static inline tk_ivec_t *tk_ivec_from_iuset (lua_State *L, tk_iuset_t *s)
{
  tk_ivec_t *v = tk_ivec_create(L, tk_iuset_size(s));
  int64_t x;
  v->n = 0;
  tk_umap_foreach_keys(s, x, ({
    v->a[v->n ++] = x;
  }))
  return v;
}

#endif
