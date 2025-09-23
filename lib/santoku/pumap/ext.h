#ifndef TK_PUMAP_EXT_H
#define TK_PUMAP_EXT_H

#include <santoku/ivec.h>
#include <santoku/pvec.h>

static inline tk_pvec_t *tk_pumap_values (lua_State *L, tk_pumap_t *M)
{
  tk_pvec_t *out = tk_pvec_create(L, tk_pumap_size(M), 0, 0);
  tk_pair_t v;
  out->n = 0;
  tk_umap_foreach_values(M, v, ({
    out->a[out->n ++] = v;
  }));
  return out;
}

static inline tk_ivec_t *tk_pumap_keys (lua_State *L, tk_pumap_t *M)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_pumap_size(M), 0, 0);
  int64_t k;
  out->n = 0;
  tk_umap_foreach_keys(M, k, ({
    out->a[out->n ++] = k;
  }));
  return out;
}

#endif