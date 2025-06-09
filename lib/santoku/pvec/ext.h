#ifndef TK_PVEC_EXT_H
#define TK_PVEC_EXT_H

static inline tk_pvec_t *tk_pvec_from_ivec (
  lua_State *L,
  tk_ivec_t *I
) {
  tk_pvec_t *P = tk_pvec_create(L, I->n, 0, 0);
  for (int64_t i = 0; i < (int64_t) I->n; i ++)
    P->a[i] = tk_pair(i, I->a[i]);
  return P;
}

static inline tk_ivec_t *tk_pvec_keys (
  lua_State *L,
  tk_pvec_t *P
) {
  tk_ivec_t *I = tk_ivec_create(L, P->n, 0, 0);
  for (uint64_t i = 0; i < P->n; i ++)
    I->a[i] = P->a[i].i;
  return I;
}

static inline tk_ivec_t *tk_pvec_values (
  lua_State *L,
  tk_pvec_t *P
) {
  tk_ivec_t *I = tk_ivec_create(L, P->n, 0, 0);
  for (uint64_t i = 0; i < P->n; i ++)
    I->a[i] = P->a[i].p;
  return I;
}

#endif

