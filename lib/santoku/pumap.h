#ifndef TK_PUMAP_H
#define TK_PUMAP_H

#include <santoku/klib.h>
#include <santoku/ivec.h>
#include <santoku/pvec.h>

KHASH_INIT(tk_pumap, int64_t, tk_pair_t, 1, kh_int64_hash_func, kh_int64_hash_equal)
typedef khash_t(tk_pumap) tk_pumap_t;

#define tk_pumap_put(...) kh_put(tk_pumap, __VA_ARGS__)
#define tk_pumap_get(...) kh_get(tk_pumap, __VA_ARGS__)
#define tk_pumap_del(...) kh_del(tk_pumap, __VA_ARGS__)
#define tk_pumap_exist(...) kh_exist(__VA_ARGS__)
#define tk_pumap_key(...) kh_key(__VA_ARGS__)
#define tk_pumap_value(...) kh_value(__VA_ARGS__)
#define tk_pumap_begin(...) kh_begin(__VA_ARGS__)
#define tk_pumap_end(...) kh_end(__VA_ARGS__)
#define tk_pumap_size(...) kh_size(__VA_ARGS__)
#define tk_pumap_resize(...) kh_resize(tk_pumap, __VA_ARGS__)
#define tk_pumap_clear(...) kh_clear(tk_pumap, __VA_ARGS__)
#define tk_pumap_destroy(...) kh_destroy(tk_pumap, __VA_ARGS__)
#define tk_pumap_create() kh_init(tk_pumap)
#define tk_pumap_foreach(...) kh_foreach(__VA_ARGS__)

#define tk_pumap_foreach_keys(h, kvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(kvar) = kh_key(h,__i); \
		code;	\
	} }

#define tk_pumap_foreach_values(h, vvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(vvar) = kh_value(h,__i); \
		code;	\
	} }

static inline tk_pvec_t *tk_pumap_values (lua_State *L, tk_pumap_t *M)
{
  tk_pvec_t *out = tk_pvec_create(L, tk_pumap_size(M), 0, 0);
  tk_pair_t v;
  out->n = 0;
  tk_pumap_foreach_values(M, v, ({
    out->a[out->n ++] = v;
  }));
  return out;
}

static inline tk_ivec_t *tk_pumap_keys (lua_State *L, tk_pumap_t *M)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_pumap_size(M), 0, 0);
  int64_t k;
  out->n = 0;
  tk_pumap_foreach_keys(M, k, ({
    out->a[out->n ++] = k;
  }));
  return out;
}

#endif
