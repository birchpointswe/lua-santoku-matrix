#ifndef TK_DUMAP_EXT_H
#define TK_DUMAP_EXT_H

#include <santoku/ivec.h>
#include <santoku/dvec.h>

static inline tk_dvec_t *tk_dumap_values (lua_State *L, tk_dumap_t *M)
{
  tk_dvec_t *out = tk_dvec_create(L, tk_dumap_size(M), 0, 0);
  double v;
  out->n = 0;
  tk_umap_foreach_values(M, v, ({
    out->a[out->n ++] = v;
  }));
  return out;
}

static inline tk_ivec_t *tk_dumap_keys (lua_State *L, tk_dumap_t *M)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_dumap_size(M), 0, 0);
  int64_t k;
  out->n = 0;
  tk_umap_foreach_keys(M, k, ({
    out->a[out->n ++] = k;
  }));
  return out;
}

#endif