// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2024 Birch Point SWE
#ifndef TK_IUMAP_EXT_H
#define TK_IUMAP_EXT_H

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

static inline int64_t tk_iumap_get_or (tk_iumap_t *map, int64_t key, int64_t default_val)
{
  uint32_t k = tk_iumap_get(map, key);
  if (k == tk_iumap_end(map))
    return default_val;
  return tk_iumap_val(map, k);
}

#endif
