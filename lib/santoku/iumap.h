#ifndef TK_IUMAP_H
#define TK_IUMAP_H

#include <santoku/klib.h>

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








#endif
