#ifndef TK_ZUMAP_H
#define TK_ZUMAP_H

#include <santoku/klib.h>

KHASH_INIT(tk_zumap, const char *, int64_t, 1, kh_str_hash_func, kh_str_hash_equal)
typedef khash_t(tk_zumap) tk_zumap_t;

#define tk_zumap_put(...) kh_put(tk_zumap, __VA_ARGS__)
#define tk_zumap_get(...) kh_get(tk_zumap, __VA_ARGS__)
#define tk_zumap_del(...) kh_del(tk_zumap, __VA_ARGS__)
#define tk_zumap_exist(...) kh_exist(__VA_ARGS__)
#define tk_zumap_key(...) kh_key(__VA_ARGS__)
#define tk_zumap_value(...) kh_value(__VA_ARGS__)
#define tk_zumap_begin(...) kh_begin(__VA_ARGS__)
#define tk_zumap_end(...) kh_end(__VA_ARGS__)
#define tk_zumap_size(...) kh_size(__VA_ARGS__)
#define tk_zumap_resize(...) kh_resize(tk_zumap, __VA_ARGS__)
#define tk_zumap_clear(...) kh_clear(tk_zumap, __VA_ARGS__)
#define tk_zumap_destroy(...) kh_destroy(tk_zumap, __VA_ARGS__)
#define tk_zumap_create() kh_init(tk_zumap)
#define tk_zumap_foreach(...) kh_foreach(__VA_ARGS__)

#endif
