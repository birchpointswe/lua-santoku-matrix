#ifndef TK_IUSET_H
#define TK_IUSET_H

#include <santoku/klib.h>
#include <santoku/ivec.h>

KHASH_INIT(tk_iuset, int64_t, char, 0, kh_int64_hash_func, kh_int64_hash_equal)
typedef khash_t(tk_iuset) tk_iuset_t;

#define tk_iuset_put(...) kh_put(tk_iuset, __VA_ARGS__)
#define tk_iuset_get(...) kh_get(tk_iuset, __VA_ARGS__)
#define tk_iuset_del(...) kh_del(tk_iuset, __VA_ARGS__)
#define tk_iuset_exist(...) kh_exist(__VA_ARGS__)
#define tk_iuset_key(...) kh_key(__VA_ARGS__)
#define tk_iuset_value(...) kh_value(__VA_ARGS__)
#define tk_iuset_begin(...) kh_begin(__VA_ARGS__)
#define tk_iuset_end(...) kh_end(__VA_ARGS__)
#define tk_iuset_size(...) kh_size(__VA_ARGS__)
#define tk_iuset_resize(...) kh_resize(tk_iuset, __VA_ARGS__)
#define tk_iuset_clear(...) kh_clear(tk_iuset, __VA_ARGS__)
#define tk_iuset_destroy(...) kh_destroy(tk_iuset, __VA_ARGS__)
#define tk_iuset_create() kh_init(tk_iuset)
#define tk_iuset_contains(h, v) (tk_iuset_get(h, v) != tk_iuset_end(h))

#define tk_iuset_foreach(h, kvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(kvar) = kh_key(h,__i); \
		code;	\
	} }

static inline tk_iuset_t *tk_iuset_from_ivec (tk_ivec_t *v)
{
  int kha;
  tk_iuset_t *s = tk_iuset_create();
  for (uint64_t i = 0; i < v->n; i ++)
    tk_iuset_put(s, v->a[i], &kha);
  return s;
}








#endif
