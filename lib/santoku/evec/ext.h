#ifndef TK_EVEC_EXT_H
#define TK_EVEC_EXT_H

static inline int tk_evec_each_lua_iter (lua_State *L)
{
  lua_settop(L, 0);
  tk_evec_t *m0 = tk_evec_peek(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_edge_t v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, v.u);
  lua_pushinteger(L, v.v);
  lua_pushnumber(L, v.w);
  return 3;
}

static inline int tk_evec_ieach_lua_iter (lua_State *L)
{
  lua_settop(L, 0);
  tk_evec_t *m0 = tk_evec_peek(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_edge_t v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, (int64_t) n);
  lua_pushinteger(L, v.u);
  lua_pushinteger(L, v.v);
  lua_pushnumber(L, v.w);
  return 4;
}

static inline int tk_evec_each_lua (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_evec_each_lua_iter, 2);
  return 1;
}

static inline int tk_evec_ieach_lua (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_evec_ieach_lua_iter, 2);
  return 1;
}

#endif
