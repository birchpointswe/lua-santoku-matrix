#ifndef TK_EUSET_EXT_H
#define TK_EUSET_EXT_H

#include <santoku/evec/base.h>

static inline tk_evec_t *tk_euset_keys (lua_State *L, tk_euset_t *M)
{
  tk_evec_t *out = tk_evec_create(L, tk_euset_size(M), 0, 0);
  tk_edge_t k;
  out->n = 0;
  tk_umap_foreach_keys(M, k, ({
    out->a[out->n ++] = k;
  }));
  return out;
}


#endif
