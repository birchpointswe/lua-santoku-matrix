#ifndef TK_UMAP_H
#define TK_UMAP_H

#include <santoku/lua/utils.h>
#include <santoku/klib.h>
#include <assert.h>

#define tk_umap_foreach(...) kh_foreach(__VA_ARGS__)

#define tk_umap_foreach_keys(h, kvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(kvar) = kh_key(h,__i); \
		code;	\
	} }

#define tk_umap_foreach_values(h, vvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(vvar) = kh_value(h,__i); \
		code;	\
	} }

#define tk_umap_foreach_iters(h, ivar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); ++__i) { \
		if (!kh_exist(h,__i)) continue; \
		(ivar) = __i; \
		code;	\
	} }

#endif

#ifndef tk_umap_name
#error "tk_umap_name not defined"
#endif

#define tk_umap_pfx(name) tk_pp_strcat(tk_umap_name, name)
#define tk_umap_err(L, name, n, ...) \
  tk_lua_verror((L), ((n) + 1), tk_pp_xstr(tk_umap_pfx(name)), __VA_ARGS__)

#ifndef tk_umap_key
#error "tk_umap_key not defined"
#endif

#ifndef tk_umap_hash
#error "tk_umap_hash not defined"
#endif

#ifndef tk_umap_eq
#error "tk_umap_eq not defined"
#endif

#define tk_umap_mt tk_pp_xstr(tk_umap_pfx(t))

#ifdef tk_umap_value
KHASH_INIT(tk_umap_name, tk_umap_key, tk_umap_value, 1, tk_umap_hash, tk_umap_eq)
#else
KHASH_INIT(tk_umap_name, tk_umap_key, char, 0, tk_umap_hash, tk_umap_eq)
#endif

#define tk_umap_t(name) khash_t(name)
#define tk_umap_init(name, h, m) kh_init(name, h, m)
#define tk_umap_destroy(name, h) kh_destroy(name, h)
#define tk_umap_resize(name, h, m) kh_resize(name, h, m)
#define tk_umap_size(h) kh_size(h)
#define tk_umap_put(name, h, k, absent) kh_put(name, h, k, absent)
#define tk_umap_get(name, h, k) kh_get(name, h, k)
#define tk_umap_del(name, h, i) kh_del(name, h, i)
#define tk_umap_begin(name, h) kh_begin(h)
#define tk_umap_end(name, h) kh_end(h)
#define tk_umap_clear(name, h) kh_clear(name, h)

typedef tk_umap_t(tk_umap_name) tk_umap_pfx(t);

static inline tk_umap_pfx(t) *tk_umap_pfx(create) (lua_State *, uint32_t);

static inline tk_umap_pfx(t) *tk_umap_pfx(peek) (lua_State *L, int i, const char *name)
{
  return (tk_umap_pfx(t) *) tk_lua_checkuserdata(L, i, tk_umap_mt, name);
}

static inline tk_umap_pfx(t) *tk_umap_pfx(peekopt) (lua_State *L, int i)
{
  return (tk_umap_pfx(t) *) tk_lua_testuserdata(L, i, tk_umap_mt);
}

static inline void tk_umap_pfx(resize) (tk_umap_pfx(t) *h, uint32_t m) {
  tk_umap_resize(tk_umap_name, h, m);
}

static inline void tk_umap_pfx(shrink) (tk_umap_pfx(t) *h) {
  tk_umap_pfx(resize)(h, tk_umap_size(h));
}

static inline uint32_t tk_umap_pfx(put) (tk_umap_pfx(t) *h, tk_umap_key k, int *absent)
{
  return (uint32_t) tk_umap_put(tk_umap_name, h, k, absent);
}

static inline uint32_t tk_umap_pfx(get) (tk_umap_pfx(t) *h, tk_umap_key k)
{
  return (uint32_t) tk_umap_get(tk_umap_name, h, k);
}

static inline void tk_umap_pfx(del) (tk_umap_pfx(t) *h, uint32_t i)
{
  tk_umap_del(tk_umap_name, h, i);
}

static inline tk_umap_key tk_umap_pfx(key) (tk_umap_pfx(t) *h, uint32_t i)
{
  return kh_key(h, i);
}

#ifdef tk_umap_value
static inline tk_umap_value tk_umap_pfx(val) (tk_umap_pfx(t) *h, uint32_t i)
{
  return kh_value(h, i);
}
#endif

static inline void tk_umap_pfx(setkey) (tk_umap_pfx(t) *h, uint32_t i, tk_umap_key k)
{
  kh_key(h, i) = k;
}

#ifdef tk_umap_value
static inline void tk_umap_pfx(setval) (tk_umap_pfx(t) *h, uint32_t i, tk_umap_value v)
{
  kh_value(h, i) = v;
}
#endif

static inline uint32_t tk_umap_pfx(begin) (tk_umap_pfx(t) *h)
{
  (void) h;
  return tk_umap_begin(tk_umap_name, h);
}

static inline uint32_t tk_umap_pfx(end) (tk_umap_pfx(t) *h)
{
  (void) h;
  return tk_umap_end(tk_umap_name, h);
}

static inline bool tk_umap_pfx(exist) (tk_umap_pfx(t) *h, uint32_t k)
{
  return kh_exist(h, k);
}

static inline bool tk_umap_pfx(contains) (tk_umap_pfx(t) *h, tk_umap_key k)
{
  uint32_t i = tk_umap_pfx(get)(h, k);
  return i != tk_umap_pfx(end)(h);
}

static inline void tk_umap_pfx(destroy) (tk_umap_pfx(t) *h)
{
#if defined(tk_umap_destroy_key) && defined(tk_umap_key)
  for (uint32_t i = tk_umap_pfx(begin)(h); i < tk_umap_pfx(end)(h); i ++)
    if (tk_umap_pfx(exist)(h, i))
      tk_umap_destroy_key(tk_umap_key(h, i));
#endif
#if defined(tk_umap_destroy_value) && defined(tk_umap_value)
  for (uint32_t i = tk_umap_pfx(begin)(h); i < tk_umap_pfx(end)(h); i ++)
    if (tk_umap_pfx(exist)(h, i))
      tk_umap_destroy_value(tk_umap_val(h, i));
#endif
  bool lua_managed = h->lua_managed;
  tk_umap_destroy(tk_umap_name, h);
  memset(h, 0, sizeof(tk_umap_pfx(t)));
  h->lua_managed = lua_managed;
}

static inline uint32_t tk_umap_pfx(size) (tk_umap_pfx(t) *h)
{
  return kh_size(h);
}

static inline void tk_umap_pfx(clear) (tk_umap_pfx(t) *h)
{
  return tk_umap_clear(tk_umap_name, h);
}

#if defined(tk_umap_pushkey) && defined(tk_umap_pushvalue)
static inline int tk_umap_pfx(each_lua_iter) (lua_State *L)
{
  lua_settop(L, 0);
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, lua_upvalueindex(1), "umap");
  uint32_t i = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  uint32_t end = tk_umap_pfx(end)(h);


  while (i < end && !tk_umap_pfx(exist)(h, i)) {
    i++;
  }

  if (i >= end)
    return 0;


  lua_pushinteger(L, i + 1);
  lua_replace(L, lua_upvalueindex(2));


  tk_umap_pushkey(L, tk_umap_pfx(key)(h, i));
  tk_umap_pushvalue(L, tk_umap_pfx(val)(h, i));
  return 2;
}

static inline int tk_umap_pfx(each_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  lua_pushvalue(L, 1);
  lua_pushinteger(L, tk_umap_pfx(begin)(h));
  lua_pushcclosure(L, tk_umap_pfx(each_lua_iter), 2);
  return 1;
}
#endif

#ifdef tk_umap_pushkey
static inline int tk_umap_pfx(keach_lua_iter) (lua_State *L)
{
  lua_settop(L, 0);
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, lua_upvalueindex(1), "umap");
  uint32_t i = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  uint32_t end = tk_umap_pfx(end)(h);


  while (i < end && !tk_umap_pfx(exist)(h, i)) {
    i++;
  }

  if (i >= end)
    return 0;


  lua_pushinteger(L, i + 1);
  lua_replace(L, lua_upvalueindex(2));


  tk_umap_pushkey(L, tk_umap_pfx(key)(h, i));
  return 1;
}

static inline int tk_umap_pfx(keach_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  lua_pushvalue(L, 1);
  lua_pushinteger(L, tk_umap_pfx(begin)(h));
  lua_pushcclosure(L, tk_umap_pfx(keach_lua_iter), 2);
  return 1;
}
#endif

static inline int tk_umap_pfx(ieach_lua_iter) (lua_State *L)
{
  lua_settop(L, 0);
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, lua_upvalueindex(1), "umap");
  uint32_t i = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  uint32_t end = tk_umap_pfx(end)(h);


  while (i < end && !tk_umap_pfx(exist)(h, i)) {
    i++;
  }

  if (i >= end)
    return 0;


  lua_pushinteger(L, i + 1);
  lua_replace(L, lua_upvalueindex(2));


  lua_pushinteger(L, i);
  return 1;
}

static inline int tk_umap_pfx(ieach_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  lua_pushvalue(L, 1);
  lua_pushinteger(L, tk_umap_pfx(begin)(h));
  lua_pushcclosure(L, tk_umap_pfx(ieach_lua_iter), 2);
  return 1;
}

#ifdef tk_umap_pushvalue
static inline int tk_umap_pfx(veach_lua_iter) (lua_State *L)
{
  lua_settop(L, 0);
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, lua_upvalueindex(1), "umap");
  uint32_t i = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  uint32_t end = tk_umap_pfx(end)(h);


  while (i < end && !tk_umap_pfx(exist)(h, i)) {
    i++;
  }

  if (i >= end)
    return 0;


  lua_pushinteger(L, i + 1);
  lua_replace(L, lua_upvalueindex(2));


  tk_umap_pushvalue(L, tk_umap_pfx(val)(h, i));
  return 1;
}

static inline int tk_umap_pfx(veach_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  lua_pushvalue(L, 1);
  lua_pushinteger(L, tk_umap_pfx(begin)(h));
  lua_pushcclosure(L, tk_umap_pfx(veach_lua_iter), 2);
  return 1;
}
#endif

static inline int tk_umap_pfx(clear_lua) (lua_State *L)
{
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  tk_umap_pfx(clear)(h);
  return 0;
}

static inline int tk_umap_pfx(del_lua) (lua_State *L)
{
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  uint32_t i = tk_lua_checkunsigned(L, 2, "i");
  tk_umap_pfx(del)(h, i);
  return 0;
}

static inline int tk_umap_pfx(exist_lua) (lua_State *L)
{
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  uint32_t i = tk_lua_checkunsigned(L, 2, "i");
  lua_pushboolean(L, tk_umap_pfx(exist)(h, i));
  return 1;
}

#ifdef tk_umap_peekkey
static inline int tk_umap_pfx(get_lua) (lua_State *L)
{
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  tk_umap_key k = tk_umap_peekkey(L, 2, "key");
  lua_pushinteger(L, tk_umap_pfx(get)(h, k));
  return 1;
}
#endif

#ifdef tk_umap_peekkey
static inline int tk_umap_pfx(put_lua) (lua_State *L)
{
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  tk_umap_key k = tk_umap_peekkey(L, 2, "key");
  int absent;
  uint32_t iter = tk_umap_pfx(put)(h, k, &absent);
  lua_pushinteger(L, iter);
  return 1;
}
#endif

#ifdef tk_umap_peekkey
static inline int tk_umap_pfx(contains_lua) (lua_State *L)
{
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  tk_umap_key k = tk_umap_peekkey(L, 2, "key");
  lua_pushboolean(L, tk_umap_pfx(contains)(h, k));
  return 1;
}
#endif

#ifdef tk_umap_peekkey
static inline int tk_umap_pfx(setkey_lua) (lua_State *L)
{
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  uint32_t i = tk_lua_checkunsigned(L, 2, "i");
  tk_umap_key k = tk_umap_peekkey(L, 3, "key");
  tk_umap_pfx(setkey)(h, i, k);
  return 0;
}
#endif

#ifdef tk_umap_peekvalue
static inline int tk_umap_pfx(setval_lua) (lua_State *L)
{
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  uint32_t i = tk_lua_checkunsigned(L, 2, "i");
  tk_umap_value v = tk_umap_peekvalue(L, 3, "value");
  tk_umap_pfx(setval)(h, i, v);
  return 0;
}
#endif

#ifdef tk_umap_pushkey
static inline int tk_umap_pfx(key_lua) (lua_State *L)
{
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  uint32_t i = tk_lua_checkunsigned(L, 2, "i");
  tk_umap_key k = tk_umap_pfx(key)(h, i);
  tk_umap_pushkey(L, k);
  return 1;
}
#endif

#ifdef tk_umap_pushvalue
static inline int tk_umap_pfx(val_lua) (lua_State *L)
{
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  uint32_t i = tk_lua_checkunsigned(L, 2, "i");
  tk_umap_value v = tk_umap_pfx(val)(h, i);
  tk_umap_pushvalue(L, v);
  return 1;
}
#endif

static inline int tk_umap_pfx(destroy_lua) (lua_State *L)
{
  tk_umap_pfx(destroy)(tk_umap_pfx(peek)(L, 1, "umap"));
  return 0;
}

static inline int tk_umap_pfx(gc_lua) (lua_State *L)
{
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  if (h->lua_managed) {
    tk_umap_pfx(destroy)(h);
  }
  return 0;
}

static inline tk_umap_pfx(t) *tk_umap_pfx(load) (lua_State *L, FILE *fh)
{
  uint32_t n;
  tk_lua_fread(L, &n, sizeof(uint32_t), 1, fh);
  tk_umap_pfx(t) *h = tk_umap_pfx(create)(L, n);
#ifdef tk_umap_value
  uint32_t khi;
  tk_umap_key k;
  tk_umap_value v;
  int a;
  for (uint32_t i = 0; i < n; i ++) {
    tk_lua_fread(L, &k, sizeof(tk_umap_key), 1, fh);
    tk_lua_fread(L, &v, sizeof(tk_umap_value), 1, fh);
    khi = tk_umap_pfx(put)(h, k, &a);
    assert(a);
    tk_umap_pfx(setval)(h, khi, v);
  }
#else
  tk_umap_key k;
  int a;
  for (uint32_t i = 0; i < n; i ++) {
    tk_lua_fread(L, &k, sizeof(tk_umap_key), 1, fh);
    tk_umap_pfx(put)(h, k, &a);
    assert(a);
  }
#endif
  return h;
}

static inline int tk_umap_pfx(load_lua) (lua_State *L)
{
  lua_settop(L, 2);
  size_t len;
  const char *data = luaL_checklstring(L, 1, &len);
  bool isstr = lua_type(L, 2) == LUA_TBOOLEAN && lua_toboolean(L, 2);
  FILE *fh = isstr ? tk_lua_fmemopen(L, (char *) data, len, "r") : tk_lua_fopen(L, data, "r");
  tk_umap_pfx(load)(L, fh);
  return 1;
}

static inline int tk_umap_pfx(size_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  lua_pushinteger(L, (int64_t) tk_umap_pfx(size)(h));
  return 1;
}

static inline int tk_umap_pfx(resize_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  uint64_t m = tk_lua_checkunsigned(L, 2, "size");
  tk_umap_pfx(resize)(h, m);
  return 0;
}

static inline int tk_umap_pfx(shrink_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_umap_pfx(t) *h = tk_umap_pfx(peek)(L, 1, "umap");
  tk_umap_pfx(shrink)(h);
  return 0;
}


static luaL_Reg tk_umap_pfx(lua_mt_fns)[] =
{
  { "clear", tk_umap_pfx(clear_lua) },
  { "del", tk_umap_pfx(del_lua) },
  { "destroy", tk_umap_pfx(destroy_lua) },
  { "exist", tk_umap_pfx(exist_lua) },
  { "resize", tk_umap_pfx(resize_lua) },
  { "shrink", tk_umap_pfx(shrink_lua) },
  { "size", tk_umap_pfx(size_lua) },

#ifdef tk_umap_peekkey
  { "get", tk_umap_pfx(get_lua) },
  { "put", tk_umap_pfx(put_lua) },
  { "contains", tk_umap_pfx(contains_lua) },
  { "setkey", tk_umap_pfx(setkey_lua) },
#endif
#ifdef tk_umap_peekvalue
  { "setval", tk_umap_pfx(setval_lua) },
#endif
#ifdef tk_umap_pushkey
  { "key", tk_umap_pfx(key_lua) },
  { "keach", tk_umap_pfx(keach_lua) },
#endif
#ifdef tk_umap_pushvalue
  { "val", tk_umap_pfx(val_lua) },
  { "veach", tk_umap_pfx(veach_lua) },
#endif
#if defined(tk_umap_pushkey) && defined(tk_umap_pushvalue)
  { "each", tk_umap_pfx(each_lua) },
#endif
  { "ieach", tk_umap_pfx(ieach_lua) },

  { NULL, NULL }
};

static inline tk_umap_pfx(t) *tk_umap_pfx(create) (lua_State *L, uint32_t n)
{
  tk_umap_pfx(t) *h = L == NULL
    ? malloc(sizeof(tk_umap_pfx(t)))
    : tk_lua_newuserdata(L, tk_umap_pfx(t), tk_umap_mt, tk_umap_pfx(lua_mt_fns), tk_umap_pfx(gc_lua));
  if (h == NULL)
    return h;
  bool lua_managed = L != NULL;
  tk_umap_init(tk_umap_name, h, lua_managed);
  if (n)
    tk_umap_pfx(resize)(h, n);
  return h;
}

static inline tk_umap_pfx(t) *tk_umap_pfx(register) (lua_State *L, tk_umap_pfx(t) *h0)
{
  tk_umap_pfx(t) *h1 = tk_lua_newuserdata(L, tk_umap_pfx(t), tk_umap_mt, tk_umap_pfx(lua_mt_fns), tk_umap_pfx(gc_lua));
  memcpy(h1, h0, sizeof(tk_umap_pfx(t)));
  tk_umap_pfx(destroy)(h0);
  return h1;
}

static inline int tk_umap_pfx(create_lua) (lua_State *L)
{
  lua_settop(L, 1);
  uint32_t n = tk_lua_optunsigned(L, 1, "size", 0);
  tk_umap_pfx(create)(L, n);
  return 1;
}

static inline void tk_umap_pfx(suppress_unused_lua_mt_fns) (void)
  { (void) tk_umap_pfx(lua_mt_fns); }

static luaL_Reg tk_umap_pfx(lua_fns)[] =
{
  { "create", tk_umap_pfx(create_lua) },
  { "load", tk_umap_pfx(load_lua) },
  { NULL, NULL }
};

static inline void tk_umap_pfx(suppress_unused_lua_fns) (void)
  { (void) tk_umap_pfx(lua_fns); }

#include <santoku/umap/undef.h>
