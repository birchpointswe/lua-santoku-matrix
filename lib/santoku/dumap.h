#ifndef TK_DUMAP_H
#define TK_DUMAP_H

#include <santoku/klib.h>
#include <santoku/ivec.h>
#include <santoku/dvec.h>

#define tk_dumap_hash(k) (kh_int64_hash_func((uint64_t)k))
KHASH_INIT(tk_dumap, int64_t, double, 1, tk_dumap_hash, kh_int64_hash_equal)
typedef khash_t(tk_dumap) tk_dumap_t;

#define tk_dumap_put(...) kh_put(tk_dumap, __VA_ARGS__)
#define tk_dumap_get(...) kh_get(tk_dumap, __VA_ARGS__)
#define tk_dumap_del(...) kh_del(tk_dumap, __VA_ARGS__)
#define tk_dumap_exist(...) kh_exist(__VA_ARGS__)
#define tk_dumap_key(...) kh_key(__VA_ARGS__)
#define tk_dumap_value(...) kh_value(__VA_ARGS__)
#define tk_dumap_begin(...) kh_begin(__VA_ARGS__)
#define tk_dumap_end(...) kh_end(__VA_ARGS__)
#define tk_dumap_size(...) kh_size(__VA_ARGS__)
#define tk_dumap_resize(...) kh_resize(tk_dumap, __VA_ARGS__)
#define tk_dumap_clear(...) kh_clear(tk_dumap, __VA_ARGS__)
#define tk_dumap_destroy(...) kh_destroy(tk_dumap, __VA_ARGS__)
#define tk_dumap_create() kh_init(tk_dumap)
#define tk_dumap_foreach(...) kh_foreach(__VA_ARGS__)

#define tk_dumap_foreach_keys(h, kvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(kvar) = kh_key(h,__i); \
		code;	\
	} }

#define tk_dumap_foreach_values(h, vvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(vvar) = kh_value(h,__i); \
		code;	\
	} }

static inline tk_dvec_t *tk_dumap_values (lua_State *L, tk_dumap_t *M)
{
  tk_dvec_t *out = tk_dvec_create(L, tk_dumap_size(M), 0, 0);
  double v;
  out->n = 0;
  tk_dumap_foreach_values(M, v, ({
    out->a[out->n ++] = v;
  }));
  return out;
}

static inline tk_ivec_t *tk_dumap_keys (lua_State *L, tk_dumap_t *M)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_dumap_size(M), 0, 0);
  int64_t k;
  out->n = 0;
  tk_dumap_foreach_keys(M, k, ({
    out->a[out->n ++] = k;
  }));
  return out;
}

#endif
