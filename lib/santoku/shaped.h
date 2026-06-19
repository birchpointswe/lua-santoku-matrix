#ifndef TK_SHAPED_H
#define TK_SHAPED_H

#include <santoku/lua/utils.h>
#include <string.h>

typedef enum {
  TK_TAG_NONE = 0,
  TK_TAG_U8,
  TK_TAG_I32,
  TK_TAG_I64,
  TK_TAG_F32,
  TK_TAG_F64,
  TK_TAG_BITS
} tk_tag_t;

static inline const char *tk_tag_name (tk_tag_t t)
{
  switch (t) {
    case TK_TAG_U8: return "u8";
    case TK_TAG_I32: return "i32";
    case TK_TAG_I64: return "i64";
    case TK_TAG_F32: return "f32";
    case TK_TAG_F64: return "f64";
    case TK_TAG_BITS: return "bits";
    default: return "none";
  }
}

static inline tk_tag_t tk_tag_from_string (const char *s)
{
  if (s == NULL) return TK_TAG_NONE;
  if (strcmp(s, "u8") == 0) return TK_TAG_U8;
  if (strcmp(s, "i32") == 0) return TK_TAG_I32;
  if (strcmp(s, "i64") == 0) return TK_TAG_I64;
  if (strcmp(s, "f32") == 0) return TK_TAG_F32;
  if (strcmp(s, "f64") == 0) return TK_TAG_F64;
  if (strcmp(s, "bits") == 0) return TK_TAG_BITS;
  return TK_TAG_NONE;
}

static inline size_t tk_tag_size (tk_tag_t t)
{
  switch (t) {
    case TK_TAG_U8: return 1;
    case TK_TAG_I32: return 4;
    case TK_TAG_I64: return 8;
    case TK_TAG_F32: return 4;
    case TK_TAG_F64: return 8;
    case TK_TAG_BITS: return 1;
    default: return 0;
  }
}



static inline void tk_lua_require_mod (lua_State *L, const char *mod)
{
  lua_getglobal(L, "require");
  lua_pushstring(L, mod);
  lua_call(L, 1, 0);
}



static inline void tk_lua_extend_mt (lua_State *L, const char *mt, luaL_Reg *fns)
{
  luaL_getmetatable(L, mt);
  lua_getfield(L, -1, "__index");
  luaL_register(L, NULL, fns);
  lua_pop(L, 2);
}



static inline void tk_eph_init (lua_State *L, int i)
{
  lua_pushvalue(L, i);
  lua_newtable(L);
  lua_setfenv(L, -2);
  lua_pop(L, 1);
}

static inline void tk_eph_anchor (lua_State *L, int ip, int ic, void *ptr)
{
  lua_pushvalue(L, ic);
  lua_pushvalue(L, ip);
  lua_getfenv(L, -1);
  lua_pushlightuserdata(L, ptr);
  lua_pushvalue(L, -4);
  lua_settable(L, -3);
  lua_pop(L, 3);
}

static inline void tk_eph_release (lua_State *L, int ip, void *ptr)
{
  lua_pushvalue(L, ip);
  lua_getfenv(L, -1);
  lua_pushlightuserdata(L, ptr);
  lua_pushnil(L);
  lua_settable(L, -3);
  lua_pop(L, 2);
}

static inline void tk_eph_get (lua_State *L, int ip, void *ptr)
{
  lua_pushvalue(L, ip);
  lua_getfenv(L, -1);
  lua_pushlightuserdata(L, ptr);
  lua_gettable(L, -2);
  lua_replace(L, -3);
  lua_pop(L, 1);
}

#endif
