#ifndef TK_CUSET_H
#define TK_CUSET_H

#include <santoku/lua/utils.h>
#include <santoku/klib.h>

KHASH_SET_INIT_STR(tk_cuset)
typedef khash_t(tk_cuset) tk_cuset_t;

#define tk_cuset_put(...) kh_put(tk_cuset, __VA_ARGS__)
#define tk_cuset_get(...) kh_get(tk_cuset, __VA_ARGS__)
#define tk_cuset_del(...) kh_del(tk_cuset, __VA_ARGS__)
#define tk_cuset_exist(...) kh_exist(__VA_ARGS__)
#define tk_cuset_key(...) kh_key(__VA_ARGS__)
#define tk_cuset_value(...) kh_value(__VA_ARGS__)
#define tk_cuset_begin(...) kh_begin(__VA_ARGS__)
#define tk_cuset_end(...) kh_end(__VA_ARGS__)
#define tk_cuset_size(...) kh_size(__VA_ARGS__)
#define tk_cuset_resize(...) kh_resize(tk_cuset, __VA_ARGS__)
#define tk_cuset_clear(...) kh_clear(tk_cuset, __VA_ARGS__)
#define tk_cuset_destroy(...) kh_destroy(tk_cuset, __VA_ARGS__)
#define tk_cuset_create() kh_init(tk_cuset)

#define tk_cuset_foreach(h, kvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(kvar) = kh_key(h,__i); \
		code;	\
	} }

#endif
