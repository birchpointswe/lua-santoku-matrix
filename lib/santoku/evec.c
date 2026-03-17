#define TK_EVEC_INIT
#include <santoku/iuset.h>
#include <santoku/evec.h>

static inline int tk_evec_get_lua (lua_State *L) {
  lua_settop(L, 2);
  tk_evec_t *P = tk_evec_peek(L, 1, "evec");
  uint64_t i = tk_lua_checkunsigned(L, 2, "idx");
  if (i > P->n)
    return 0;
  lua_pushinteger(L, P->a[i].u);
  lua_pushinteger(L, P->a[i].v);
  lua_pushnumber(L, P->a[i].w);
  return 3;
}

static inline int tk_evec_set_lua (lua_State *L) {
  lua_settop(L, 5);
  tk_evec_t *P = tk_evec_peek(L, 1, "evec");
  uint64_t i = tk_lua_checkunsigned(L, 2, "idx");
  int64_t u = tk_lua_checkinteger(L, 3, "u");
  int64_t v = tk_lua_checkinteger(L, 4, "v");
  double w = tk_lua_checkdouble(L, 5, "w");
  tk_evec_ensure(P, i + 1);
  P->a[i] = tk_edge(u, v, w);
  return 0;
}

static inline int tk_evec_push_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_evec_t *P = tk_evec_peek(L, 1, "evec");
  int64_t u = tk_lua_checkinteger(L, 2, "u");
  int64_t v = tk_lua_checkinteger(L, 3, "v");
  double w = tk_lua_checkdouble(L, 4, "w");
  tk_evec_push(P, tk_edge(u, v, w));
  return 0;
}

static inline int tk_evec_hmax_lua (lua_State *L) {
  lua_settop(L, 5);
  tk_evec_t *P = tk_evec_peek(L, 1, "evec");
  int64_t u = tk_lua_checkinteger(L, 2, "u");
  int64_t v = tk_lua_checkinteger(L, 3, "v");
  double w = tk_lua_checkdouble(L, 4, "w");
  uint64_t m = tk_lua_checkunsigned(L, 5, "max");
  tk_evec_hmax(P, m, tk_edge(u, v, w));
  return 0;
}

static inline int tk_evec_hmax_init_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_evec_t *P = tk_evec_peek(L, 1, "evec");
  tk_evec_hmax_init(P);
  return 0;
}

static inline int tk_evec_hmax_push_lua (lua_State *L) {
  lua_settop(L, 4);
  tk_evec_t *P = tk_evec_peek(L, 1, "evec");
  int64_t u = tk_lua_checkinteger(L, 2, "u");
  int64_t v = tk_lua_checkinteger(L, 3, "v");
  double w = tk_lua_checkdouble(L, 4, "w");
  tk_evec_hmax_push(P, tk_edge(u, v, w));
  return 0;
}

static inline int tk_evec_hmax_pop_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_evec_t *P = tk_evec_peek(L, 1, "evec");
  if (!P->n)
    return 0;
  tk_edge_t r = tk_evec_hmax_pop(P);
  lua_pushinteger(L, r.u);
  lua_pushinteger(L, r.v);
  lua_pushnumber(L, r.w);
  return 3;
}

static inline int tk_evec_hmin_lua (lua_State *L) {
  lua_settop(L, 5);
  tk_evec_t *P = tk_evec_peek(L, 1, "evec");
  int64_t u = tk_lua_checkinteger(L, 2, "u");
  int64_t v = tk_lua_checkinteger(L, 3, "v");
  double w = tk_lua_checkdouble(L, 4, "w");
  uint64_t m = tk_lua_checkunsigned(L, 5, "max");
  tk_evec_hmin(P, m, tk_edge(u, v, w));
  return 0;
}

static inline int tk_evec_hmin_init_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_evec_t *P = tk_evec_peek(L, 1, "evec");
  tk_evec_hmin_init(P);
  return 0;
}

static inline int tk_evec_hmin_push_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_evec_t *P = tk_evec_peek(L, 1, "evec");
  int64_t u = tk_lua_checkinteger(L, 2, "u");
  int64_t v = tk_lua_checkinteger(L, 3, "v");
  double w = tk_lua_checkdouble(L, 4, "w");
  tk_evec_hmin_push(P, tk_edge(u, v, w));
  return 0;
}

static inline int tk_evec_hmin_pop_lua (lua_State *L) {
  lua_settop(L, 1);
  tk_evec_t *P = tk_evec_peek(L, 1, "evec");
  if (!P->n)
    return 0;
  tk_edge_t r = tk_evec_hmin_pop(P);
  lua_pushinteger(L, r.u);
  lua_pushinteger(L, r.v);
  lua_pushnumber(L, r.w);
  return 3;
}

static luaL_Reg tk_evec_lua_mt_ext2_fns[] =
{
  { "get", tk_evec_get_lua },
  { "set", tk_evec_set_lua },
  { "push", tk_evec_push_lua },
  { "each", tk_evec_each_lua },
  { "ieach", tk_evec_ieach_lua },
  { "hmax", tk_evec_hmax_lua },
  { "hmax_init", tk_evec_hmax_init_lua },
  { "hmax_pop", tk_evec_hmax_pop_lua },
  { "hmax_push", tk_evec_hmax_push_lua },
  { "hmin", tk_evec_hmin_lua },
  { "hmin_init", tk_evec_hmin_init_lua },
  { "hmin_pop", tk_evec_hmin_pop_lua },
  { "hmin_push", tk_evec_hmin_push_lua },
  { NULL, NULL }
};

void tk_evec_init_mt (lua_State *L)
{
  luaL_register(L, NULL, tk_evec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_evec_lua_mt_ext2_fns);
}

int luaopen_santoku_evec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_evec_lua_fns);
  tk_evec_create(L, 0, 0, 0);
  lua_pop(L, 1);
  return 1;
}
