#include <santoku/iuset.h>
#include <santoku/ivec.h>
#include <santoku/fvec.h>
#include <santoku/dvec.h>












static inline int tk_dvec_round_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_dvec_t *v = tk_dvec_peek(L, 1, "dvec");
  uint64_t start = (t >= 2) ? tk_lua_checkunsigned(L, 2, "start") : 0;
  uint64_t end = (t >= 3) ? tk_lua_checkunsigned(L, 3, "end") : v->n;
  tk_dvec_round(v, start, end);
  lua_pushvalue(L, 1);
  return 1;
}

static inline int tk_dvec_trunc_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_dvec_t *v = tk_dvec_peek(L, 1, "dvec");
  uint64_t start = (t >= 2) ? tk_lua_checkunsigned(L, 2, "start") : 0;
  uint64_t end = (t >= 3) ? tk_lua_checkunsigned(L, 3, "end") : v->n;
  tk_dvec_trunc(v, start, end);
  lua_pushvalue(L, 1);
  return 1;
}

static inline int tk_dvec_floor_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_dvec_t *v = tk_dvec_peek(L, 1, "dvec");
  uint64_t start = (t >= 2) ? tk_lua_checkunsigned(L, 2, "start") : 0;
  uint64_t end = (t >= 3) ? tk_lua_checkunsigned(L, 3, "end") : v->n;
  tk_dvec_floor(v, start, end);
  lua_pushvalue(L, 1);
  return 1;
}

static inline int tk_dvec_ceil_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_dvec_t *v = tk_dvec_peek(L, 1, "dvec");
  uint64_t start = (t >= 2) ? tk_lua_checkunsigned(L, 2, "start") : 0;
  uint64_t end = (t >= 3) ? tk_lua_checkunsigned(L, 3, "end") : v->n;
  tk_dvec_ceil(v, start, end);
  lua_pushvalue(L, 1);
  return 1;
}

static inline int tk_dvec_to_ivec_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_dvec_t *v = tk_dvec_peek(L, 1, "dvec");
  tk_dvec_to_ivec(L, v);
  return 1;
}

static inline int tk_dvec_to_fvec_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_dvec_t *v = tk_dvec_peek(L, 1, "dvec");
  tk_fvec_t *out = lua_isnil(L, 2) ? NULL : tk_fvec_peek(L, 2, "out");
  tk_dvec_to_fvec(L, v, out);
  return out == NULL ? 1 : 0;
}






static luaL_Reg tk_dvec_lua_mt_ext2_fns[] =
{
  { "round", tk_dvec_round_lua },
  { "trunc", tk_dvec_trunc_lua },
  { "floor", tk_dvec_floor_lua },
  { "ceil", tk_dvec_ceil_lua },
  { "to_ivec", tk_dvec_to_ivec_lua },
  { "to_fvec", tk_dvec_to_fvec_lua },
  { NULL, NULL }
};

int luaopen_santoku_dvec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_dvec_lua_fns);
  tk_dvec_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_dvec_lua_mt_fns);
  luaL_register(L, NULL, tk_dvec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_dvec_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
