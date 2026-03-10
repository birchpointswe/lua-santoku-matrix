#if !defined(__EMSCRIPTEN__)

#include <santoku/lua/utils.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#ifndef tk_mmap_name
#error "tk_mmap_name not defined"
#endif

#ifndef tk_mmap_base
#error "tk_mmap_base not defined"
#endif

#ifndef tk_mmap_vec_name
#error "tk_mmap_vec_name not defined"
#endif

#define tk_mmap_pfx(name) tk_pp_strcat(tk_mmap_name, name)
#define tk_mmap_mt tk_pp_xstr(tk_mmap_pfx(t))
#define tk_mmap_err(L, name, n, ...) \
  tk_lua_verror((L), ((n) + 1), tk_pp_xstr(tk_mmap_pfx(name)), __VA_ARGS__)

#define tk_mmap_vec_pfx(name) tk_pp_strcat(tk_mmap_vec_name, name)
#define tk_mmap_vec_t tk_mmap_vec_pfx(t)

typedef struct {
  size_t n;
  tk_mmap_base *a;
  size_t map_len;
  bool lua_managed;
} tk_mmap_pfx(t);

static inline tk_mmap_pfx(t) *tk_mmap_pfx(create) (lua_State *, const char *, size_t);

static inline tk_mmap_pfx(t) *tk_mmap_pfx(peek) (lua_State *L, int i, const char *name)
{
  return (tk_mmap_pfx(t) *) tk_lua_checkuserdata(L, i, tk_mmap_mt, name);
}

static inline tk_mmap_pfx(t) *tk_mmap_pfx(peekopt) (lua_State *L, int i)
{
  return (tk_mmap_pfx(t) *) tk_lua_testuserdata(L, i, tk_mmap_mt);
}

static inline void tk_mmap_pfx(destroy) (tk_mmap_pfx(t) *m)
{
  if (m->a && m->a != MAP_FAILED) {
    munmap(m->a, m->map_len);
    m->a = NULL;
    m->map_len = 0;
    m->n = 0;
  }
}

static inline int tk_mmap_pfx(msync) (tk_mmap_pfx(t) *m)
{
  if (m->a && m->map_len > 0)
    return msync(m->a, m->map_len, MS_SYNC);
  return 0;
}

static inline size_t tk_mmap_pfx(size) (tk_mmap_pfx(t) *m)
{
  return m->n;
}

static inline tk_mmap_base tk_mmap_pfx(get) (tk_mmap_pfx(t) *m, size_t i)
{
  return m->a[i];
}

static inline void tk_mmap_pfx(set) (tk_mmap_pfx(t) *m, size_t i, tk_mmap_base v)
{
  m->a[i] = v;
}

static inline tk_mmap_vec_t tk_mmap_pfx(vec) (tk_mmap_pfx(t) *m)
{
  return (tk_mmap_vec_t) { .a = m->a, .n = m->n, .m = m->n, .lua_managed = 0 };
}

static inline int tk_mmap_pfx(destroy_lua) (lua_State *L)
{
  tk_mmap_pfx(destroy)(tk_mmap_pfx(peek)(L, 1, "mmap"));
  return 0;
}

static inline int tk_mmap_pfx(gc_lua) (lua_State *L)
{
  tk_mmap_pfx(t) *m = tk_mmap_pfx(peek)(L, 1, "mmap");
  if (m->lua_managed)
    tk_mmap_pfx(destroy)(m);
  return 0;
}

static inline int tk_mmap_pfx(msync_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_mmap_pfx(t) *m = tk_mmap_pfx(peek)(L, 1, "mmap");
  if (tk_mmap_pfx(msync)(m) != 0)
    return tk_mmap_err(L, msync, 1, strerror(errno));
  return 0;
}

static inline int tk_mmap_pfx(size_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_mmap_pfx(t) *m = tk_mmap_pfx(peek)(L, 1, "mmap");
  lua_pushinteger(L, (int64_t) tk_mmap_pfx(size)(m));
  return 1;
}

static inline int tk_mmap_pfx(vec_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_mmap_pfx(t) *m = tk_mmap_pfx(peek)(L, 1, "mmap");
  tk_mmap_vec_t *v = tk_mmap_vec_pfx(create)(L, 0, NULL, NULL);
  v->a = m->a;
  v->n = m->n;
  v->m = m->n;
  v->lua_managed = 0;
  return 1;
}

#ifdef tk_mmap_pushbase
static inline int tk_mmap_pfx(get_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_mmap_pfx(t) *m = tk_mmap_pfx(peek)(L, 1, "mmap");
  uint64_t i = tk_lua_checkunsigned(L, 2, "idx");
  tk_mmap_pushbase(L, tk_mmap_pfx(get)(m, i));
  return 1;
}
#endif

#ifdef tk_mmap_peekbase
static inline int tk_mmap_pfx(set_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_mmap_pfx(t) *m = tk_mmap_pfx(peek)(L, 1, "mmap");
  uint64_t i = tk_lua_checkunsigned(L, 2, "idx");
  tk_mmap_base v = tk_mmap_peekbase(L, 3);
  tk_mmap_pfx(set)(m, i, v);
  return 0;
}
#endif

static luaL_Reg tk_mmap_pfx(lua_mt_fns)[] =
{
  { "destroy", tk_mmap_pfx(destroy_lua) },
  { "msync", tk_mmap_pfx(msync_lua) },
  { "size", tk_mmap_pfx(size_lua) },
  { "vec", tk_mmap_pfx(vec_lua) },
#ifdef tk_mmap_pushbase
  { "get", tk_mmap_pfx(get_lua) },
#endif
#ifdef tk_mmap_peekbase
  { "set", tk_mmap_pfx(set_lua) },
#endif
  { NULL, NULL }
};

static inline void tk_mmap_pfx(suppress_unused_lua_mt_fns) (void)
  { (void) tk_mmap_pfx(lua_mt_fns); }

static inline tk_mmap_pfx(t) *tk_mmap_pfx(create) (lua_State *L, const char *path, size_t n)
{
  size_t map_len = n * sizeof(tk_mmap_base);
  int fd = open(path, O_RDWR | O_CREAT | O_TRUNC, 0666);
  if (fd < 0)
    return (tk_mmap_pfx(t) *)(intptr_t) tk_mmap_err(L, create, 1, strerror(errno));
  if (map_len > 0) {
    if (ftruncate(fd, (off_t) map_len) != 0) {
      close(fd);
      return (tk_mmap_pfx(t) *)(intptr_t) tk_mmap_err(L, create, 1, strerror(errno));
    }
  }
  tk_mmap_pfx(t) *m = tk_lua_newuserdata(L, tk_mmap_pfx(t), tk_mmap_mt, tk_mmap_pfx(lua_mt_fns), tk_mmap_pfx(gc_lua));
  m->lua_managed = true;
  if (map_len > 0) {
    m->a = (tk_mmap_base *) mmap(NULL, map_len, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (m->a == MAP_FAILED) {
      close(fd);
      return (tk_mmap_pfx(t) *)(intptr_t) tk_mmap_err(L, create, 1, strerror(errno));
    }
    memset(m->a, 0, map_len);
  } else {
    m->a = NULL;
  }
  close(fd);
  m->n = n;
  m->map_len = map_len;
  return m;
}

static inline tk_mmap_pfx(t) *tk_mmap_pfx(open) (lua_State *L, const char *path)
{
  int fd = open(path, O_RDWR);
  if (fd < 0)
    return (tk_mmap_pfx(t) *)(intptr_t) tk_mmap_err(L, open, 1, strerror(errno));
  struct stat st;
  if (fstat(fd, &st) != 0) {
    close(fd);
    return (tk_mmap_pfx(t) *)(intptr_t) tk_mmap_err(L, open, 1, strerror(errno));
  }
  size_t map_len = (size_t) st.st_size;
  size_t n = map_len / sizeof(tk_mmap_base);
  tk_mmap_pfx(t) *m = tk_lua_newuserdata(L, tk_mmap_pfx(t), tk_mmap_mt, tk_mmap_pfx(lua_mt_fns), tk_mmap_pfx(gc_lua));
  m->lua_managed = true;
  if (map_len > 0) {
    m->a = (tk_mmap_base *) mmap(NULL, map_len, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (m->a == MAP_FAILED) {
      close(fd);
      return (tk_mmap_pfx(t) *)(intptr_t) tk_mmap_err(L, open, 1, strerror(errno));
    }
  } else {
    m->a = NULL;
  }
  close(fd);
  m->n = n;
  m->map_len = map_len;
  return m;
}

static inline int tk_mmap_pfx(create_lua) (lua_State *L)
{
  lua_settop(L, 2);
  const char *path = luaL_checkstring(L, 1);
  uint64_t n = tk_lua_checkunsigned(L, 2, "size");
  tk_mmap_pfx(create)(L, path, n);
  return 1;
}

static inline int tk_mmap_pfx(open_lua) (lua_State *L)
{
  lua_settop(L, 1);
  const char *path = luaL_checkstring(L, 1);
  tk_mmap_pfx(open)(L, path);
  return 1;
}

static luaL_Reg tk_mmap_pfx(lua_fns)[] =
{
  { "create", tk_mmap_pfx(create_lua) },
  { "open", tk_mmap_pfx(open_lua) },
  { NULL, NULL }
};

static inline void tk_mmap_pfx(suppress_unused_lua_fns) (void)
  { (void) tk_mmap_pfx(lua_fns); }

#include <santoku/mmap/undef.h>

#endif
