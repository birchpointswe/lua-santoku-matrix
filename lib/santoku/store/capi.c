#include <santoku/lua/utils.h>
#include <santoku/fvec.h>
#include <santoku/dvec.h>
#include <santoku/ivec.h>
#include <santoku/svec.h>
#include <santoku/cvec.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#if !defined(__EMSCRIPTEN__)
#include <sys/mman.h>
#include <sys/statvfs.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#endif

#define TK_STORE_MT "tk_store_t"
#define TK_STORE_ALIGN 64

typedef enum { TK_STORE_F32, TK_STORE_F64, TK_STORE_I64, TK_STORE_I32, TK_STORE_U8 } tk_store_kind_t;

typedef struct {
  bool disk;
  int state;
  char *base;
  size_t len;
  bool mapped;
  uint64_t n_pieces, cap;
  uint64_t *ns;
  uint8_t *kinds;
} tk_store_t;

static inline tk_store_t *tk_store_peek (lua_State *L, int i)
{
  return (tk_store_t *) luaL_checkudata(L, i, TK_STORE_MT);
}

static inline size_t tk_store_esz (uint8_t kind)
{
  switch (kind) {
    case TK_STORE_F32: return sizeof(float);
    case TK_STORE_F64: return sizeof(double);
    case TK_STORE_I64: return sizeof(int64_t);
    case TK_STORE_I32: return sizeof(int32_t);
    default: return 1;
  }
}

static inline void tk_store_release (tk_store_t *s)
{
  if (s->state != 1)
    return;
#if !defined(__EMSCRIPTEN__)
  if (s->mapped) {
    if (s->base) munmap(s->base, s->len);
  } else {
    free(s->base);
  }
#else
  free(s->base);
#endif
  s->base = NULL;
  s->len = 0;
  s->state = 2;
}

static int tk_store_gc_lua (lua_State *L)
{
  tk_store_t *s = tk_store_peek(L, 1);
  tk_store_release(s);
  free(s->ns);
  free(s->kinds);
  s->ns = NULL;
  s->kinds = NULL;
  return 0;
}

static inline void tk_store_view_set (lua_State *L, int vi, uint8_t kind, char *a, uint64_t n)
{
  switch (kind) {
    case TK_STORE_F32: { tk_fvec_t *v = tk_fvec_peek(L, vi, "view"); v->a = (float *) a; v->n = v->m = n; break; }
    case TK_STORE_F64: { tk_dvec_t *v = tk_dvec_peek(L, vi, "view"); v->a = (double *) a; v->n = v->m = n; break; }
    case TK_STORE_I64: { tk_ivec_t *v = tk_ivec_peek(L, vi, "view"); v->a = (int64_t *) a; v->n = v->m = n; break; }
    case TK_STORE_I32: { tk_svec_t *v = tk_svec_peek(L, vi, "view"); v->a = (int32_t *) a; v->n = v->m = n; break; }
    default: { tk_cvec_t *v = tk_cvec_peek(L, vi, "view"); v->a = a; v->n = v->m = n; break; }
  }
}

static inline void tk_store_view_mark (lua_State *L, int vi, uint8_t kind)
{
  switch (kind) {
    case TK_STORE_F32: { tk_fvec_t *v = tk_fvec_peek(L, vi, "view"); free(v->a); v->a = NULL; v->n = v->m = 0; v->lua_managed = 3; break; }
    case TK_STORE_F64: { tk_dvec_t *v = tk_dvec_peek(L, vi, "view"); free(v->a); v->a = NULL; v->n = v->m = 0; v->lua_managed = 3; break; }
    case TK_STORE_I64: { tk_ivec_t *v = tk_ivec_peek(L, vi, "view"); free(v->a); v->a = NULL; v->n = v->m = 0; v->lua_managed = 3; break; }
    case TK_STORE_I32: { tk_svec_t *v = tk_svec_peek(L, vi, "view"); free(v->a); v->a = NULL; v->n = v->m = 0; v->lua_managed = 3; break; }
    default: { tk_cvec_t *v = tk_cvec_peek(L, vi, "view"); free(v->a); v->a = NULL; v->n = v->m = 0; v->lua_managed = 3; break; }
  }
}

static int tk_store_declare (lua_State *L, uint8_t kind)
{
  lua_settop(L, 2);
  tk_store_t *s = tk_store_peek(L, 1);
  if (s->state != 0)
    return tk_lua_verror(L, 2, "store", "pieces must be declared before open");
  uint64_t n = tk_lua_checkunsigned(L, 2, "n");
  if (s->n_pieces == s->cap) {
    uint64_t cap = s->cap ? s->cap * 2 : 8;
    uint64_t *ns = (uint64_t *) realloc(s->ns, cap * sizeof(uint64_t));
    if (!ns) return tk_lua_verror(L, 2, "store", "out of memory");
    s->ns = ns;
    uint8_t *kinds = (uint8_t *) realloc(s->kinds, cap);
    if (!kinds) return tk_lua_verror(L, 2, "store", "out of memory");
    s->kinds = kinds;
    s->cap = cap;
  }
  switch (kind) {
    case TK_STORE_F32: tk_fvec_create(L, 0); break;
    case TK_STORE_F64: tk_dvec_create(L, 0); break;
    case TK_STORE_I64: tk_ivec_create(L, 0); break;
    case TK_STORE_I32: tk_svec_create(L, 0); break;
    default: tk_cvec_create(L, 0); break;
  }
  int vi = lua_gettop(L);
  tk_store_view_mark(L, vi, kind);
  lua_newtable(L);
  lua_pushvalue(L, 1);
  lua_setfield(L, -2, "store");
  lua_setfenv(L, vi);
  lua_getfenv(L, 1);
  lua_pushvalue(L, vi);
  lua_rawseti(L, -2, (int) s->n_pieces + 1);
  lua_pop(L, 1);
  s->ns[s->n_pieces] = n;
  s->kinds[s->n_pieces] = kind;
  s->n_pieces ++;
  lua_settop(L, vi);
  return 1;
}

static int tk_store_fvec_lua (lua_State *L) { return tk_store_declare(L, TK_STORE_F32); }
static int tk_store_dvec_lua (lua_State *L) { return tk_store_declare(L, TK_STORE_F64); }
static int tk_store_ivec_lua (lua_State *L) { return tk_store_declare(L, TK_STORE_I64); }
static int tk_store_svec_lua (lua_State *L) { return tk_store_declare(L, TK_STORE_I32); }
static int tk_store_cvec_lua (lua_State *L) { return tk_store_declare(L, TK_STORE_U8); }

static int tk_store_open_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_store_t *s = tk_store_peek(L, 1);
  if (s->state != 0)
    return tk_lua_verror(L, 2, "store", "open called twice or after close");
  size_t total = 0;
  for (uint64_t i = 0; i < s->n_pieces; i ++) {
    total = (total + TK_STORE_ALIGN - 1) & ~((size_t) TK_STORE_ALIGN - 1);
    total += s->ns[i] * tk_store_esz(s->kinds[i]);
  }
  char *base = NULL;
  bool mapped = false;
#if !defined(__EMSCRIPTEN__)
  if (s->disk && total > 0) {
    const char *dir = getenv("TMPDIR");
    if (!dir || !*dir) dir = "/var/tmp";
    struct statvfs sv;
    if (statvfs(dir, &sv) != 0)
      return tk_lua_verror(L, 2, "store", strerror(errno));
    if ((uint64_t) sv.f_bavail * (uint64_t) sv.f_frsize < total)
      return tk_lua_verror(L, 2, "store", "not enough free disk space to spill");
    char path[4096];
    snprintf(path, sizeof path, "%s/tkstore.XXXXXX", dir);
    int fd = mkstemp(path);
    if (fd < 0)
      return tk_lua_verror(L, 2, "store", strerror(errno));
    unlink(path);
    if (ftruncate(fd, (off_t) total) != 0) {
      close(fd);
      return tk_lua_verror(L, 2, "store", strerror(errno));
    }
    void *p = mmap(NULL, total, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);
    if (p == MAP_FAILED)
      return tk_lua_verror(L, 2, "store", strerror(errno));
    base = (char *) p;
    mapped = true;
  }
#endif
  if (!mapped && total > 0) {
    size_t kb = total >> 10;
    if (kb > (size_t) lua_gc(L, LUA_GCCOUNT, 0))
      lua_gc(L, LUA_GCCOLLECT, 0);
    else if (kb > 0)
      lua_gc(L, LUA_GCSTEP, (int) kb);
    base = (char *) calloc(1, total);
    if (!base)
      return tk_lua_verror(L, 2, "store", "out of memory");
  }
  s->base = base;
  s->len = total;
  s->mapped = mapped;
  s->state = 1;
  lua_getfenv(L, 1);
  int views = lua_gettop(L);
  size_t off = 0;
  for (uint64_t i = 0; i < s->n_pieces; i ++) {
    off = (off + TK_STORE_ALIGN - 1) & ~((size_t) TK_STORE_ALIGN - 1);
    lua_rawgeti(L, views, (int) i + 1);
    tk_store_view_set(L, lua_gettop(L), s->kinds[i], base + off, s->ns[i]);
    lua_pop(L, 1);
    off += s->ns[i] * tk_store_esz(s->kinds[i]);
  }
  lua_settop(L, 1);
  return 1;
}

static int tk_store_close_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_store_t *s = tk_store_peek(L, 1);
  if (s->state == 1) {
    lua_getfenv(L, 1);
    int views = lua_gettop(L);
    for (uint64_t i = 0; i < s->n_pieces; i ++) {
      lua_rawgeti(L, views, (int) i + 1);
      tk_store_view_set(L, lua_gettop(L), s->kinds[i], NULL, 0);
      lua_pop(L, 1);
    }
  }
  tk_store_release(s);
  s->state = 2;
  return 0;
}

static int tk_store_bytes_lua (lua_State *L)
{
  tk_store_t *s = tk_store_peek(L, 1);
  lua_pushinteger(L, (lua_Integer) s->len);
  return 1;
}

static int tk_store_on_disk_lua (lua_State *L)
{
  tk_store_t *s = tk_store_peek(L, 1);
  lua_pushboolean(L, s->mapped);
  return 1;
}

static luaL_Reg tk_store_mt_fns[] = {
  { "fvec", tk_store_fvec_lua },
  { "dvec", tk_store_dvec_lua },
  { "ivec", tk_store_ivec_lua },
  { "svec", tk_store_svec_lua },
  { "cvec", tk_store_cvec_lua },
  { "open", tk_store_open_lua },
  { "close", tk_store_close_lua },
  { "bytes", tk_store_bytes_lua },
  { "on_disk", tk_store_on_disk_lua },
  { NULL, NULL }
};

static int tk_store_create_lua (lua_State *L)
{
  lua_settop(L, 1);
  bool disk = lua_isnil(L, 1) ? false : tk_lua_foptboolean(L, 1, "store", "disk", false);
  tk_store_t *s = (tk_store_t *) lua_newuserdata(L, sizeof(tk_store_t));
  memset(s, 0, sizeof(tk_store_t));
  s->disk = disk;
  if (luaL_newmetatable(L, TK_STORE_MT)) {
    lua_newtable(L);
    luaL_register(L, NULL, tk_store_mt_fns);
    lua_setfield(L, -2, "__index");
    lua_pushcfunction(L, tk_store_gc_lua);
    lua_setfield(L, -2, "__gc");
  }
  lua_setmetatable(L, -2);
  lua_newtable(L);
  lua_setfenv(L, -2);
  return 1;
}

static luaL_Reg tk_store_fns[] = {
  { "create", tk_store_create_lua },
  { NULL, NULL }
};

int luaopen_santoku_store_capi (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_store_fns);
  return 1;
}
