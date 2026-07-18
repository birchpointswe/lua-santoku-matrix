#include <santoku/iuset.h>
#include <santoku/ivec.h>
#include <santoku/fvec.h>
#include <santoku/dvec.h>
#include <math.h>












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
















static inline int tk_dvec_group_gauge_lua (lua_State *L)
{
  lua_settop(L, 6);
  tk_dvec_t *pc = tk_dvec_peek(L, 1, "colsumsq");
  tk_ivec_t *go = tk_ivec_peek(L, 2, "group_offsets");
  luaL_checktype(L, 3, LUA_TTABLE);
  double n = luaL_checknumber(L, 4);
  tk_fvec_t *w = NULL;
  double floorv = 1e-6;
  int exps_idx = 0;
  if (!lua_isnil(L, 5)) {
    luaL_checktype(L, 5, LUA_TTABLE);
    lua_getfield(L, 5, "w");
    if (lua_isnil(L, -1)) lua_pop(L, 1);
    else w = tk_fvec_peek(L, -1, "w");
    lua_getfield(L, 5, "floor");
    if (!lua_isnil(L, -1)) floorv = luaL_checknumber(L, -1);
    lua_pop(L, 1);
    lua_getfield(L, 5, "exps");
    if (lua_isnil(L, -1)) lua_pop(L, 1);
    else exps_idx = lua_gettop(L);
  }
  tk_fvec_t *out_arg = lua_isnil(L, 6) ? NULL : tk_fvec_peek(L, 6, "out");
  uint64_t G = go->n > 0 ? go->n - 1 : 0;
  uint64_t total = G > 0 ? (uint64_t) go->a[G] : 0;
  if (total > pc->n)
    return luaL_error(L, "group_gauge: offsets exceed colsumsq length");
  if (w && w->n < total)
    return luaL_error(L, "group_gauge: w shorter than columns");
  if (out_arg && out_arg->m < total)
    return luaL_error(L, "group_gauge: out buffer shorter than columns");
  tk_fvec_t *out;
  if (out_arg) {
    out = out_arg;
    lua_pushvalue(L, 6);
  } else {
    out = tk_fvec_create(L, total);
  }
  out->n = total;
  for (uint64_t g = 0; g < G; g ++) {
    int64_t lo = go->a[g], hi = go->a[g + 1];
    lua_rawgeti(L, 3, (int) g + 1);
    double sg = lua_type(L, -1) == LUA_TNUMBER ? lua_tonumber(L, -1) : 1.0;
    lua_pop(L, 1);
    double eg = 1.0;
    if (exps_idx) {
      lua_rawgeti(L, exps_idx, (int) g + 1);
      if (lua_type(L, -1) == LUA_TNUMBER) eg = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }
    double wssq = 0.0;
    if (w) {
      double logsum = 0.0;
      for (int64_t c = lo; c < hi; c ++) {
        double wv = (double) w->a[c];
        if (wv < floorv) wv = floorv;
        logsum += log(wv);
      }
      double geo = hi > lo ? exp(logsum / (double) (hi - lo)) : 1.0;
      for (int64_t c = lo; c < hi; c ++) {
        double wv = (double) w->a[c];
        if (wv < floorv) wv = floorv;
        double cs = pow(wv / geo, eg);
        out->a[c] = (float) cs;
        wssq += cs * cs * pc->a[c];
      }
    } else {
      for (int64_t c = lo; c < hi; c ++) {
        out->a[c] = 1.0f;
        wssq += pc->a[c];
      }
    }
    double mult = wssq > 0.0 ? sg * sqrt(n / wssq) : 0.0;
    for (int64_t c = lo; c < hi; c ++)
      out->a[c] = (float) ((double) out->a[c] * mult);
  }
  return 1;
}

static luaL_Reg tk_dvec_lua_mt_ext2_fns[] =
{
  { "group_gauge", tk_dvec_group_gauge_lua },
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
