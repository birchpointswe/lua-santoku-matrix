#ifndef TK_DUSET_H
#define TK_DUSET_H

#include <santoku/lua/utils.h>
#include <santoku/klib.h>
#include <santoku/dvec.h>

#define tk_duset_double_equal(a, b) ((a) == (b))
KHASH_INIT(tk_duset, double, char, 0, tk_lua_hash_double, tk_duset_double_equal)
typedef khash_t(tk_duset) tk_duset_t;

#define tk_duset_put(...) kh_put(tk_duset, __VA_ARGS__)
#define tk_duset_get(...) kh_get(tk_duset, __VA_ARGS__)
#define tk_duset_del(...) kh_del(tk_duset, __VA_ARGS__)
#define tk_duset_exist(...) kh_exist(__VA_ARGS__)
#define tk_duset_key(...) kh_key(__VA_ARGS__)
#define tk_duset_value(...) kh_value(__VA_ARGS__)
#define tk_duset_begin(...) kh_begin(__VA_ARGS__)
#define tk_duset_end(...) kh_end(__VA_ARGS__)
#define tk_duset_size(...) kh_size(__VA_ARGS__)
#define tk_duset_resize(...) kh_resize(tk_duset, __VA_ARGS__)
#define tk_duset_clear(...) kh_clear(tk_duset, __VA_ARGS__)
#define tk_duset_destroy(...) kh_destroy(tk_duset, __VA_ARGS__)
#define tk_duset_create() kh_init(tk_duset)

#define tk_duset_foreach(h, kvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(kvar) = kh_key(h,__i); \
		code;	\
	} }

static void tk_duset_dump (tk_duset_t *s, tk_dvec_t *v)
{
  double x;
  tk_duset_foreach(s, x, ({
    tk_dvec_push(v, x);
  }))
}








#endif
