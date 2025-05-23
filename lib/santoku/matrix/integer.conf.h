#ifndef TK_MATRIX_INTEGER_CONF_H
#define TK_MATRIX_INTEGER_CONF_H

#define tk_base_t int64_t
#define tk_sort(...) ks_introsort(int64_t, __VA_ARGS__)
#define tk_matrix_pushbase(...) lua_pushinteger(__VA_ARGS__)
#define tk_matrix_peekbase(...) luaL_checkinteger(__VA_ARGS__)
#define TK_MT "santoku_matrix_integer"
#define TK_OPEN luaopen_santoku_matrix_integer

#endif
