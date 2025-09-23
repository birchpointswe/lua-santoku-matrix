#define tk_umap_pfx(name) tk_pp_strcat(tk_umap_name, name)

#ifndef tk_umap_no_persist
static inline void tk_umap_pfx(persist) (lua_State *L, tk_umap_pfx(t) *h, FILE *fh)
{
  uint32_t n = tk_umap_pfx(size)(h);
  tk_lua_fwrite(L, (char *) &n, sizeof(uint32_t), 1, fh);
#ifdef tk_umap_value
  tk_umap_key k;
  tk_umap_value v;
  tk_umap_foreach(h, k, v, ({
    tk_lua_fwrite(L, (char *) &k, sizeof(tk_umap_key), 1, fh);
    tk_lua_fwrite(L, (char *) &v, sizeof(tk_umap_value), 1, fh);
  }))
#else
  tk_umap_key k;
  tk_umap_foreach_keys(h, k, ({
    tk_lua_fwrite(L, (char *) &k, sizeof(tk_umap_key), 1, fh);
  }))
#endif
}
#endif

#ifndef tk_umap_no_persist
static inline int tk_umap_pfx(persist_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  bool tostr = lua_type(L, 2) == LUA_TNIL;
  FILE *fh;
  if (tostr)
    fh = tk_lua_tmpfile(L);
  else
    fh = tk_lua_fopen(L, luaL_checkstring(L, 2), "w");
  tk_umap_pfx(persist)(L, h, fh);
  if (!tostr) {
    tk_lua_fclose(L, fh);
    return 0;
  } else {
    size_t len;
    char *data = tk_lua_fslurp(L, fh, &len);
    if (data) {
      tk_cvec_create(L, len, data, 0);
      return 1;
    } else {
      tk_lua_fclose(L, fh);
      return 0;
    }
  }
}
#endif

static luaL_Reg tk_umap_pfx(lua_ext_fns)[] =
{
#ifndef tk_umap_no_persist
  { "persist", tk_umap_pfx(persist_lua) },
#endif
  { NULL, NULL }
};

static inline void tk_umap_pfx(suppress_unused_lua_ext_fns) (void)
  { (void) tk_umap_pfx(lua_ext_fns); }

#include <santoku/umap/undef.h>
