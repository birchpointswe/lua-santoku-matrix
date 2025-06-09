#ifndef TK_IUMAP_H
#define TK_IUMAP_H

#include <santoku/klib.h>
#include <santoku/ivec.h>

KHASH_INIT(tk_iumap, int64_t, int64_t, 1, kh_int64_hash_func, kh_int64_hash_equal)
typedef khash_t(tk_iumap) tk_iumap_t;

#define tk_iumap_put(...) kh_put(tk_iumap, __VA_ARGS__)
#define tk_iumap_get(...) kh_get(tk_iumap, __VA_ARGS__)
#define tk_iumap_del(...) kh_del(tk_iumap, __VA_ARGS__)
#define tk_iumap_exist(...) kh_exist(__VA_ARGS__)
#define tk_iumap_key(...) kh_key(__VA_ARGS__)
#define tk_iumap_value(...) kh_value(__VA_ARGS__)
#define tk_iumap_begin(...) kh_begin(__VA_ARGS__)
#define tk_iumap_end(...) kh_end(__VA_ARGS__)
#define tk_iumap_size(...) kh_size(__VA_ARGS__)
#define tk_iumap_resize(...) kh_resize(tk_iumap, __VA_ARGS__)
#define tk_iumap_clear(...) kh_clear(tk_iumap, __VA_ARGS__)
#define tk_iumap_destroy(...) kh_destroy(tk_iumap, __VA_ARGS__)
#define tk_iumap_create() kh_init(tk_iumap)
#define tk_iumap_foreach(...) kh_foreach(__VA_ARGS__)

static inline tk_iumap_t *tk_iumap_from_ivec (tk_ivec_t *V)
{
  int kha;
  khint_t khi;
  tk_iumap_t *M = tk_iumap_create();
  for (int64_t i = 0; i < (int64_t) V->n; i ++) {
    khi = tk_iumap_put(M, V->a[i], &kha);
    tk_iumap_value(M, khi) = i;
  }
  return M;
}

static inline tk_ivec_t *tk_iumap_values (lua_State *L, tk_iumap_t *M)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_iumap_size(M), 0, 0);
  int64_t k, v;
  out->n = 0;
  tk_iumap_foreach(M, k, v, ({
    out->a[out->n ++] = v;
  }));
  return out;
}

static inline tk_ivec_t *tk_iumap_keys (lua_State *L, tk_iumap_t *M)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_iumap_size(M), 0, 0);
  int64_t k, v;
  out->n = 0;
  tk_iumap_foreach(M, k, v, ({
    out->a[out->n ++] = k;
  }));
  return out;
}








#endif
