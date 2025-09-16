#ifndef TK_IUMAP_H
#define TK_IUMAP_H

#include <santoku/klib.h>
#include <santoku/ivec/base.h>

#define tk_iumap_hash(k) (kh_int64_hash_func((uint64_t)k))
KHASH_INIT(tk_iumap, int64_t, int64_t, 1, tk_iumap_hash, kh_int64_hash_equal)
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

#define tk_iumap_foreach_keys(h, kvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(kvar) = kh_key(h,__i); \
		code;	\
	} }

#define tk_iumap_foreach_values(h, vvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(vvar) = kh_value(h,__i); \
		code;	\
	} }

static inline tk_ivec_t *tk_iumap_values (lua_State *L, tk_iumap_t *M)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_iumap_size(M), 0, 0);
  int64_t v;
  out->n = 0;
  tk_iumap_foreach_values(M, v, ({
    out->a[out->n ++] = v;
  }));
  return out;
}

static inline tk_ivec_t *tk_iumap_keys (lua_State *L, tk_iumap_t *M)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_iumap_size(M), 0, 0);
  int64_t k;
  out->n = 0;
  tk_iumap_foreach_keys(M, k, ({
    out->a[out->n ++] = k;
  }));
  return out;
}

static inline void tk_iumap_inc (tk_iumap_t *map, int64_t key)
{
  int absent;
  khint_t k = tk_iumap_put(map, key, &absent);
  if (absent) {
    tk_iumap_value(map, k) = 1;
  } else {
    tk_iumap_value(map, k)++;
  }
}

static inline void tk_iumap_add (tk_iumap_t *map, int64_t key, int64_t val)
{
  int absent;
  khint_t k = tk_iumap_put(map, key, &absent);
  if (absent) {
    tk_iumap_value(map, k) = val;
  } else {
    tk_iumap_value(map, k) += val;
  }
}

static inline int64_t tk_iumap_get_or (tk_iumap_t *map, int64_t key, int64_t default_val)
{
  khint_t k = tk_iumap_get(map, key);
  if (k == tk_iumap_end(map)) {
    return default_val;
  }
  return tk_iumap_value(map, k);
}

static inline int tk_ivec_bits_rearrange (
  tk_ivec_t *m0,
  tk_ivec_t *ids,
  uint64_t n_features
) {
  tk_ivec_asc(m0, 0, m0->n);
  tk_iumap_t *remap = tk_iumap_create();
  int kha;
  khint_t khi;
  for (int64_t i = 0; i < (int64_t) ids->n; i ++) {
    khi = tk_iumap_put(remap, ids->a[i], &kha);
    tk_iumap_value(remap, khi) = i;
  }
  size_t write = 0;
  for (int64_t i = 0; i < (int64_t) m0->n; i ++) {
    int64_t b = m0->a[i];
    int64_t s0 = b / (int64_t) n_features;
    khi = tk_iumap_get(remap, s0);
    if (khi == tk_iumap_end(remap))
      continue;
    int64_t s1 = tk_iumap_value(remap, khi);
    int64_t f = b % (int64_t) n_features;
    m0->a[write ++] = s1 * (int64_t) n_features + f;
  }
  m0->n = write;
  tk_iumap_destroy(remap);
  tk_ivec_asc(m0, 0, m0->n);
  return 0;
}

#endif
