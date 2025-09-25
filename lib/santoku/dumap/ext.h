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

static inline double tk_dumap_correlation (
  tk_dumap_t *nbrs_a,
  tk_dumap_t *nbrs_b
) {
  uint32_t i;
  uint32_t n = 0;
  int64_t a, ra, rb;
  double sum_squared_diff = 0.0;
  tk_umap_foreach(nbrs_a, a, ra, ({
    i = tk_dumap_get(nbrs_b, a);
    if (i == tk_dumap_end(nbrs_b))
      continue;
    n ++;
    rb = tk_dumap_val(nbrs_b, i);
    double diff = ra - rb;
    sum_squared_diff += diff * diff;
  }))
  if (n > 1) {
    double correlation = 1.0 - (6.0 * sum_squared_diff) / (n * (n * n - 1));
    return correlation;
  }
  return 1.0;
}

#endif
