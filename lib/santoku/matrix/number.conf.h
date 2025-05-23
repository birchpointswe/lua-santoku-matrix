#ifndef TK_MATRIX_NUMBER_CONF_H
#define TK_MATRIX_NUMBER_CONF_H

#define tk_base_t double
#define tk_sort(...) ks_introsort(double, __VA_ARGS__)
#define tk_matrix_pushbase(...) lua_pushnumber(__VA_ARGS__)
#define tk_matrix_peekbase(...) luaL_checknumber(__VA_ARGS__)
#define TK_MT "santoku_matrix_number"
#define TK_OPEN luaopen_santoku_matrix_number

#endif
