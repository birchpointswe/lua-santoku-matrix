#if !defined(__EMSCRIPTEN__)

#include <santoku/lua/utils.h>

#ifndef tk_mmap_name
#error "tk_mmap_name not defined"
#endif

#ifndef tk_mmap_base
#error "tk_mmap_base not defined"
#endif

#define tk_mmap_pfx(name) tk_pp_strcat(tk_mmap_name, name)
#define tk_mmap_mt tk_pp_xstr(tk_mmap_pfx(t))
#define tk_mmap_err(L, name, n, ...) \
  tk_lua_verror((L), ((n) + 1), tk_pp_xstr(tk_mmap_pfx(name)), __VA_ARGS__)

static inline void tk_mmap_pfx(persist) (lua_State *L, tk_mmap_pfx(t) *m, FILE *fh)
{
  tk_lua_fwrite(L, (char *) m->a, sizeof(tk_mmap_base) * m->n, 1, fh);
}

static inline int tk_mmap_pfx(persist_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_mmap_pfx(t) *m = tk_mmap_pfx(peek)(L, 1, "mmap");
  const char *path = luaL_checkstring(L, 2);
  FILE *fh = tk_lua_fopen(L, path, "w");
  tk_mmap_pfx(persist)(L, m, fh);
  tk_lua_fclose(L, fh);
  return 0;
}

static luaL_Reg tk_mmap_pfx(lua_ext_fns)[] =
{
  { "persist", tk_mmap_pfx(persist_lua) },
  { NULL, NULL }
};

static inline void tk_mmap_pfx(suppress_unused_lua_ext_fns) (void)
  { (void) tk_mmap_pfx(lua_ext_fns); }

#undef tk_mmap_pfx
#undef tk_mmap_mt
#undef tk_mmap_err

#endif
