#ifndef TK_CVEC_H
#define TK_CVEC_H

#include <santoku/vec/base.h>

#define tk_vec_name tk_cvec
#define tk_vec_base char
#define tk_vec_limited
#define tk_vec_pushbase(L, x) lua_pushinteger(L, (unsigned char)(x))
#define tk_vec_peekbase(L, idx) tk_cvec_checkbyte(L, idx)

static inline unsigned char tk_cvec_checkbyte(lua_State *L, int idx) {
  lua_Integer val = luaL_checkinteger(L, idx);
  if (val < 0 || val > 255) {
    luaL_error(L, "byte value must be between 0 and 255, got %lld", (long long)val);
  }
  return (unsigned char)val;
}

#include <santoku/vec/ext/tpl.h>

#include <santoku/cvec/ext.h>

#endif
