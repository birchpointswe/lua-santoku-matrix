#ifndef TK_PVEC_EXT_H
#define TK_PVEC_EXT_H

#include <santoku/dumap.h>

static inline tk_pvec_t *tk_pvec_from_ivec (
  lua_State *L,
  tk_ivec_t *ivec
) {
  tk_pvec_t *P = tk_pvec_create(L, ivec->n, 0, 0);
  for (int64_t i = 0; i < (int64_t) ivec->n; i ++)
    P->a[i] = tk_pair(i, ivec->a[i]);
  return P;
}

static inline tk_ivec_t *tk_pvec_keys (
  lua_State *L,
  tk_pvec_t *P,
  tk_ivec_t *out
) {
  tk_ivec_t *result = out ? out : tk_ivec_create(L, P->n, 0, 0);
  if (out)
    tk_ivec_ensure(result, P->n);
  for (uint64_t i = 0; i < P->n; i ++)
    result->a[i] = P->a[i].i;
  if (out)
    result->n = P->n;
  return result;
}

static inline tk_ivec_t *tk_pvec_values (
  lua_State *L,
  tk_pvec_t *P,
  tk_ivec_t *out
) {
  tk_ivec_t *result = out ? out : tk_ivec_create(L, P->n, 0, 0);
  if (out)
    tk_ivec_ensure(result, P->n);
  for (uint64_t i = 0; i < P->n; i ++)
    result->a[i] = P->a[i].p;
  if (out)
    result->n = P->n;
  return result;
}

static inline int tk_pvec_each_lua_iter (lua_State *L)
{
  lua_settop(L, 0);
  tk_pvec_t *m0 = tk_pvec_peek(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_pair_t v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, v.i);
  lua_pushinteger(L, v.p);
  return 2;
}

static inline int tk_pvec_ieach_lua_iter (lua_State *L)
{
  lua_settop(L, 0);
  tk_pvec_t *m0 = tk_pvec_peek(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_pair_t v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, (int64_t) n);
  lua_pushinteger(L, v.i);
  lua_pushinteger(L, v.p);
  return 3;
}

static inline int tk_pvec_each_lua (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_pvec_each_lua_iter, 2);
  return 1;
}

static inline int tk_pvec_ieach_lua (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_pvec_ieach_lua_iter, 2);
  return 1;
}

static inline void tk_pvec_ranks (
  tk_pvec_t *v,
  tk_dumap_t *r
) {
  tk_dumap_clear(r);
  int kha;
  uint32_t khi;
  for (size_t i = 0; i < v->n; i ++) {
    double rank = (double) i;
    size_t count = 1;
    while (i + 1 < v->n && v->a[i].p == v->a[i + 1].p) {
      count ++;
      i ++;
    }
    double average_rank = (rank + (rank + count - 1)) / 2.0;
    for (size_t j = 0; j < count; j ++) {
      khi = tk_dumap_put(r, v->a[i - j].i, &kha);
      tk_dumap_setval(r, khi, average_rank);
    }
  }
}

#endif
