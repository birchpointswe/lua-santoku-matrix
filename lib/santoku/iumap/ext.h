#ifndef TK_IUMAP_EXT_H
#define TK_IUMAP_EXT_H

#include <santoku/cvec/base.h>

static inline tk_iumap_t *tk_iumap_from_ivec (lua_State *L, tk_ivec_t *V)
{
  int kha;
  uint32_t khi;
  tk_iumap_t *M = tk_iumap_create(L, V->n);
  if (!M)
    return NULL;
  for (int64_t i = 0; i < (int64_t) V->n; i ++) {
    khi = tk_iumap_put(M, V->a[i], &kha);
    if (kha < 0) {
      tk_iumap_destroy(M);
      return NULL;
    }
    tk_iumap_setval(M, khi, i);
  }
  return M;
}

static inline tk_ivec_t *tk_iumap_values (lua_State *L, tk_iumap_t *M)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_iumap_size(M));
  int64_t v;
  out->n = 0;
  tk_umap_foreach_values(M, v, ({
    out->a[out->n ++] = v;
  }));
  return out;
}

static inline tk_ivec_t *tk_iumap_keys (lua_State *L, tk_iumap_t *M)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_iumap_size(M));
  int64_t k;
  out->n = 0;
  tk_umap_foreach_keys(M, k, ({
    out->a[out->n ++] = k;
  }));
  return out;
}

static inline void tk_iumap_inc (tk_iumap_t *map, int64_t key)
{
  int absent;
  uint32_t k = tk_iumap_put(map, key, &absent);
  if (absent) {
    tk_iumap_setval(map, k, 1);
  } else {
    tk_iumap_setval(map, k, tk_iumap_val(map, k) + 1);
  }
}

static inline void tk_iumap_add (tk_iumap_t *map, int64_t key, int64_t val)
{
  int absent;
  uint32_t k = tk_iumap_put(map, key, &absent);
  if (absent) {
    tk_iumap_setval(map, k, val);
  } else {
    tk_iumap_setval(map, k, tk_iumap_val(map, k) + val);
  }
}

static inline int64_t tk_iumap_get_or (tk_iumap_t *map, int64_t key, int64_t default_val)
{
  uint32_t k = tk_iumap_get(map, key);
  if (k == tk_iumap_end(map))
    return default_val;
  return tk_iumap_val(map, k);
}

static inline void tk_iumap_union (tk_iumap_t *a, tk_iumap_t *b)
{
  int absent;
  uint32_t k;
  int64_t key, value;
  tk_umap_foreach(b, key, value, ({
    k = tk_iumap_put(a, key, &absent);
    if (absent) {
      tk_iumap_setval(a, k, value);
    }

  }))
}

static inline void tk_iumap_intersect (tk_iumap_t *a, tk_iumap_t *b)
{
  uint32_t i;
  int64_t key;
  tk_umap_foreach_iters(a, i, ({
    key = tk_iumap_key(a, i);
    if (!tk_iumap_contains(b, key))
      tk_iumap_del(a, i);
  }))
}

static inline void tk_iumap_subtract (tk_iumap_t *a, tk_iumap_t *b)
{
  uint32_t i;
  int64_t key;
  tk_umap_foreach_keys(b, key, ({
    i = tk_iumap_get(a, key);
    if (tk_iumap_exist(a, i))
      tk_iumap_del(a, i);
  }))
}


#endif
