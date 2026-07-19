#include <santoku/iuset.h>
#include <santoku/ivec.h>
#include <santoku/dvec.h>
#include <santoku/fvec.h>
#include <math.h>












static inline int tk_fvec_round_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_fvec_t *v = tk_fvec_peek(L, 1, "fvec");
  uint64_t start = (t >= 2) ? tk_lua_checkunsigned(L, 2, "start") : 0;
  uint64_t end = (t >= 3) ? tk_lua_checkunsigned(L, 3, "end") : v->n;
  tk_fvec_round(v, start, end);
  lua_pushvalue(L, 1);
  return 1;
}

static inline int tk_fvec_trunc_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_fvec_t *v = tk_fvec_peek(L, 1, "fvec");
  uint64_t start = (t >= 2) ? tk_lua_checkunsigned(L, 2, "start") : 0;
  uint64_t end = (t >= 3) ? tk_lua_checkunsigned(L, 3, "end") : v->n;
  tk_fvec_trunc(v, start, end);
  lua_pushvalue(L, 1);
  return 1;
}

static inline int tk_fvec_floor_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_fvec_t *v = tk_fvec_peek(L, 1, "fvec");
  uint64_t start = (t >= 2) ? tk_lua_checkunsigned(L, 2, "start") : 0;
  uint64_t end = (t >= 3) ? tk_lua_checkunsigned(L, 3, "end") : v->n;
  tk_fvec_floor(v, start, end);
  lua_pushvalue(L, 1);
  return 1;
}

static inline int tk_fvec_ceil_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_fvec_t *v = tk_fvec_peek(L, 1, "fvec");
  uint64_t start = (t >= 2) ? tk_lua_checkunsigned(L, 2, "start") : 0;
  uint64_t end = (t >= 3) ? tk_lua_checkunsigned(L, 3, "end") : v->n;
  tk_fvec_ceil(v, start, end);
  lua_pushvalue(L, 1);
  return 1;
}






static inline int tk_fvec_colscale_lua (lua_State *L)
{
  lua_settop(L, 5);
  tk_fvec_t *w = tk_fvec_peek(L, 1, "fvec");
  tk_dvec_t *pc = tk_dvec_peek(L, 2, "colsumsq");
  double e = luaL_checknumber(L, 3);
  double floorv = luaL_optnumber(L, 4, 1e-6);
  uint64_t k = w->n;
  if (k > pc->n) k = pc->n;


  tk_fvec_t *cs;
  if (lua_isnil(L, 5)) {
    cs = tk_fvec_create(L, k);
  } else {
    cs = tk_fvec_peek(L, 5, "out");
    tk_fvec_ensure(cs, k);
    lua_pushvalue(L, 5);
  }
  cs->n = k;
  double logsum = 0.0;
  #pragma omp parallel for reduction(+:logsum) schedule(static)
  for (uint64_t r = 0; r < k; r ++) {
    double wv = (double) w->a[r];
    if (wv < floorv) wv = floorv;
    logsum += log(wv);
  }
  double g = k > 0 ? exp(logsum / (double) k) : 1.0;
  double wssq = 0.0;
  #pragma omp parallel for reduction(+:wssq) schedule(static)
  for (uint64_t r = 0; r < k; r ++) {
    double wv = (double) w->a[r];
    if (wv < floorv) wv = floorv;
    double m = pow(wv / g, e);
    cs->a[r] = (float) m;
    wssq += m * m * pc->a[r];
  }
  lua_pushnumber(L, wssq);
  return 2;
}

static inline int tk_fvec_to_ivec_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_fvec_t *v = tk_fvec_peek(L, 1, "fvec");
  tk_fvec_to_ivec(L, v);
  return 1;
}

static inline int tk_fvec_to_dvec_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_fvec_t *v = tk_fvec_peek(L, 1, "fvec");
  tk_dvec_t *out = lua_isnil(L, 2) ? NULL : tk_dvec_peek(L, 2, "out");
  tk_fvec_to_dvec(L, v, out);
  return out == NULL ? 1 : 0;
}






static luaL_Reg tk_fvec_lua_mt_ext2_fns[] =
{
  { "round", tk_fvec_round_lua },
  { "trunc", tk_fvec_trunc_lua },
  { "floor", tk_fvec_floor_lua },
  { "ceil", tk_fvec_ceil_lua },
  { "to_ivec", tk_fvec_to_ivec_lua },
  { "to_dvec", tk_fvec_to_dvec_lua },
  { "colscale", tk_fvec_colscale_lua },
  { NULL, NULL }
};

int luaopen_santoku_fvec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_fvec_lua_fns);
  tk_fvec_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_fvec_lua_mt_fns);
  luaL_register(L, NULL, tk_fvec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_fvec_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
