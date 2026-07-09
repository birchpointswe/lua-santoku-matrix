#include <santoku/lua/utils.h>
#include <santoku/klib.h>
#if defined(_OPENMP) && !defined(__EMSCRIPTEN__)
#include <omp.h>
#endif
#if !defined(__EMSCRIPTEN__)
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#endif

#ifndef tk_vec_lt
#define tk_vec_lt(a, b) ((a) < (b))
#endif

#ifndef tk_vec_gt
#define tk_vec_gt(a, b) ((a) > (b))
#endif

#ifndef tk_vec_eqx
#define tk_vec_eqx(a, b) (memcmp(&(a), &(b), sizeof(a)) == 0)
#endif

#ifndef tk_vec_eq
#define tk_vec_eq(a, b) tk_vec_eqx(a, b)
#endif

#ifndef tk_vec_mt
#define tk_vec_mt tk_pp_xstr(tk_vec_pfx(t))
#endif

#define tk_vec_pfx(name) tk_pp_strcat(tk_vec_name, name)
#define tk_vec_err(L, name, n, ...) \
  tk_lua_verror((L), ((n) + 1), tk_pp_xstr(tk_vec_pfx(name)), __VA_ARGS__)

typedef kvec_t(tk_vec_base) tk_vec_pfx(t);
#define tk_vec_ksort(...) KSORT_INIT(__VA_ARGS__)
tk_vec_ksort(tk_vec_pfx(asc), tk_vec_base, tk_vec_lt)
tk_vec_ksort(tk_vec_pfx(desc), tk_vec_base, tk_vec_gt)

#define tk_vec_shuffle(...) ks_shuffle(__VA_ARGS__)
#define tk_vec_introsort(...) ks_introsort(__VA_ARGS__)
#define tk_vec_ksmall(...) ks_ksmall(__VA_ARGS__)


static inline void tk_vec_pfx(transpose) (tk_vec_pfx(t) *m0, tk_vec_pfx(t) *m1, uint64_t cols);
#ifndef tk_vec_limited
static inline tk_vec_pfx(t) *tk_vec_pfx(csums) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_vec_pfx(t) *tk_vec_pfx(rsums) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_vec_pfx(t) *tk_vec_pfx(cmaxs) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_vec_pfx(t) *tk_vec_pfx(rmaxs) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_vec_pfx(t) *tk_vec_pfx(cmins) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline tk_vec_pfx(t) *tk_vec_pfx(rmins) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t cols);
static inline void tk_vec_pfx(scale) (tk_vec_pfx(t) *m0, tk_vec_base scale, uint64_t start, uint64_t end);
static inline void tk_vec_pfx(add) (tk_vec_pfx(t) *m0, tk_vec_base add, uint64_t start, uint64_t end);
static inline void tk_vec_pfx(add_scaled) (tk_vec_pfx(t) *m0, tk_vec_base add, uint64_t start, uint64_t end);
static inline void tk_vec_pfx(scalev) (tk_vec_pfx(t) *m0, tk_vec_pfx(t) *m1, uint64_t start, uint64_t end);
static inline void tk_vec_pfx(addv) (tk_vec_pfx(t) *m0, tk_vec_pfx(t) *m1, uint64_t start, uint64_t end);
static inline void tk_vec_pfx(abs) (tk_vec_pfx(t) *m0, uint64_t start, uint64_t end);
static inline double tk_vec_pfx(dot) (tk_vec_pfx(t) *a, tk_vec_pfx(t) *b);
static inline void tk_vec_pfx(multiply) (tk_vec_pfx(t) *a, tk_vec_pfx(t) *b, tk_vec_pfx(t) *c, uint64_t k, bool transpose_a, bool transpose_b);
static inline void tk_vec_pfx(pow) (tk_vec_pfx(t) *v, double exponent, uint64_t start, uint64_t end);
static inline void tk_vec_pfx(log) (tk_vec_pfx(t) *v, uint64_t start, uint64_t end);
static inline void tk_vec_pfx(exp) (tk_vec_pfx(t) *v, uint64_t start, uint64_t end);
static inline void tk_vec_pfx(fill_indices) (tk_vec_pfx(t) *v);
#endif
static inline void tk_vec_pfx(fill) (tk_vec_pfx(t) *v, tk_vec_base x, uint64_t start, uint64_t end);

static inline tk_vec_pfx(t) tk_vec_name (tk_vec_base *a, uint64_t n)
{
  return (tk_vec_pfx(t)) { .a = a, .n = n, .m = n };
}

static inline tk_vec_pfx(t) *tk_vec_pfx(peek) (lua_State *L, int i, const char *name)
{
  return (tk_vec_pfx(t) *) tk_lua_checkuserdata(L, i, tk_vec_mt, name);
}

static inline tk_vec_pfx(t) *tk_vec_pfx(peekopt) (lua_State *L, int i)
{
  return (tk_vec_pfx(t) *) tk_lua_testuserdata(L, i, tk_vec_mt);
}

static inline void tk_vec_pfx(destroy) (tk_vec_pfx(t) *r)
{
#ifdef tk_vec_destroy_item
  for (uint64_t i = 0; i < r->n; i ++)
    tk_vec_destroy_item(r->a[i]);
#endif
  bool lua_managed = r->lua_managed;
#if !defined(__EMSCRIPTEN__)
  if (r->lua_managed == 2) {
    if (r->a) munmap(r->a, r->m * sizeof(tk_vec_base));
    r->a = NULL;
    r->n = r->m = 0;
    r->lua_managed = -1;
    return;
  }
#endif
  kv_destroy(*r);
  if (!lua_managed)
    free(r);
  else {
    r->a = NULL;
    r->n = r->m = 0;
  }
}

static inline bool tk_vec_pfx(lt) (tk_vec_base a, tk_vec_base b)
{
  return tk_vec_lt(a, b);
}

static inline bool tk_vec_pfx(gt) (tk_vec_base a, tk_vec_base b)
{
  return tk_vec_gt(a, b);
}

static inline tk_vec_pfx(t) *tk_vec_pfx(create) (lua_State *, size_t);

static inline int tk_vec_pfx(resize) (
  tk_vec_pfx(t) *m0,
  size_t m,
  bool setn
) {
  if (m0->lua_managed == 2)
    return -1;
  tk_vec_pfx(t) v = *m0;
  if (m == 0) {
    free(v.a);
    v.a = NULL;
    v.m = v.n = 0;
    *m0 = v;
    return 0;
  } else {
    int rc = kv_resize(tk_vec_base, v, m);
    if (rc != 0)
      return -1;
    v.m = m;
    if (setn)
      v.n = m;
    *m0 = v;
    return 0;
  }
}

static inline int tk_vec_pfx(setn) (
  tk_vec_pfx(t) *m0,
  size_t n
) {
  if (n > m0->m)
    return tk_vec_pfx(resize)(m0, n, true);
  else {
    m0->n = n;
    return 0;
  }
}

static inline int tk_vec_pfx(ensure) (
  tk_vec_pfx(t) *m0,
  size_t m
) {
  if (m0->m >= m)
    return 0;
  return tk_vec_pfx(resize)(m0, m, false);
}

static inline int tk_vec_pfx(reverse) (
  tk_vec_pfx(t) *m0,
  int64_t start,
  int64_t end
) {
  if (start < 0 || start >= end || start >= (int64_t) m0->n)
    return 0;
  if (end >= (int64_t) m0->n)
    end = (int64_t) m0->n;
  tk_vec_base b;
  for (int64_t i = start; i < start + (end - start) / 2; i ++) {
    b = m0->a[i];
    m0->a[i] = m0->a[start + end - 1 - i];
    m0->a[start + end - 1 - i] = b;
  }
  return 0;
}

static inline int tk_vec_pfx(copy) (
  tk_vec_pfx(t) *m0,
  tk_vec_pfx(t) *m1,
  int64_t start,
  int64_t end,
  int64_t dest
) {
  if (start < 0 || start >= end || start >= (int64_t) m1->n)
    return 0;
  if (end >= (int64_t) m1->n)
    end = (int64_t) m1->n;
  uint64_t m = (uint64_t) dest + (uint64_t) (end - start);
  if (tk_vec_pfx(ensure)(m0, m) != 0)
    return -1;
  if (m0 == m1)
    memmove(m0->a + dest, m1->a + start, sizeof(tk_vec_base) * (uint64_t) (end - start));
  else
    memcpy(m0->a + dest, m1->a + start, sizeof(tk_vec_base) * (uint64_t) (end - start));
  if (m0->n < m)
    m0->n = m;
  return 0;
}


static inline uint64_t tk_vec_pfx(capacity) (tk_vec_pfx(t) *m0)
{
  return m0->m;
}

static inline uint64_t tk_vec_pfx(size) (tk_vec_pfx(t) *m0)
{
  return m0->n;
}

static inline void tk_vec_pfx(shrink) (
  tk_vec_pfx(t) *v
) {
  tk_vec_pfx(resize)(v, v->n, true);
}

static inline void tk_vec_pfx(clear) (tk_vec_pfx(t) *v) {
  v->n = 0;
}
static inline void tk_vec_pfx(clear_range) (tk_vec_pfx(t) *v, int64_t start, int64_t end) {
  if (start < 0) start = 0;
  if (end < 0 || end > (int64_t) v->n) end = (int64_t) v->n;
  if (start < end)
    memset(v->a + start, 0, (size_t)(end - start) * sizeof(tk_vec_base));
}
static inline void tk_vec_pfx(zero) (tk_vec_pfx(t) *v) {
  memset(v->a, 0, v->m * sizeof(tk_vec_base));
}
static inline tk_vec_base tk_vec_pfx(get) (tk_vec_pfx(t) *v, uint64_t i) {
  return v->a[i];
}
static inline int64_t tk_vec_pfx(find) (tk_vec_pfx(t) *v, tk_vec_base x) {
  for (uint64_t i = 0; i < v->n; i++) {
    if (tk_vec_eq(v->a[i], x))
      return (int64_t) i;
  }
  return -1;
}
static inline int tk_vec_pfx(push) (tk_vec_pfx(t) *v, tk_vec_base x) {
  if (v->lua_managed == 2)
    return -1;
  tk_vec_pfx(t) v0 = *v;
  int rc = kv_push(tk_vec_base, v0, x);
  *v = v0;
  return rc;
}

static inline int tk_vec_pfx(insert) (
  tk_vec_pfx(t) *v,
  uint64_t idx,
  tk_vec_base x
) {
  if (idx > v->n) idx = v->n;
  int rc = tk_vec_pfx(ensure)(v, v->n + 1);
  if (rc != 0)
    return -1;
  if (idx < v->n)
    memmove(v->a + idx + 1, v->a + idx, (v->n - idx) * sizeof(tk_vec_base));
  v->a[idx] = x;
  v->n++;
  return 0;
}

static inline int tk_vec_pfx(set) (tk_vec_pfx(t) *v, uint64_t i, tk_vec_base x) {
  int rc = tk_vec_pfx(ensure)(v, i + 1);
  if (rc != 0)
    return -1;
  v->a[i] = x;
  if (i + 1 > v->n)
    v->n = i + 1;
  return 0;
}

#ifdef tk_vec_pushbase

static inline void tk_vec_pfx(table) (lua_State *L, tk_vec_pfx(t) *m0, uint64_t start, uint64_t end) {
  lua_newtable(L);
  for (uint64_t i = start; i < end; i ++) {
    lua_pushinteger(L, (int64_t) i + 1);
    tk_vec_pushbase(L, m0->a[i]);
    lua_settable(L, -3);
  }
}

static inline void tk_vec_pfx(ctable) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols,
  uint64_t start,
  uint64_t end
) {
  bool done = false;
  lua_newtable(L);
  for (uint64_t c = start; !done && c < end; c ++) {
    lua_pushinteger(L, (int64_t) c + 1);
    lua_newtable(L);
    for (uint64_t r = 0; r < m0->n / cols; r ++) {
      if (r * cols + c >= m0->n)
        continue;
      lua_pushinteger(L, (int64_t) r + 1);
      tk_vec_pushbase(L, m0->a[r * cols + c]);
      lua_settable(L, -3);
    }
    lua_settable(L, -3);
  }
}

static inline void tk_vec_pfx(rtable) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols,
  uint64_t start,
  uint64_t end
) {
  bool done = false;
  lua_newtable(L);
  for (uint64_t r = start; !done && r < end; r ++) {
    lua_pushinteger(L, (int64_t) r + 1);
    lua_newtable(L);
    for (uint64_t c = 0; c < cols; c ++) {
      if (r * cols + c >= m0->n) {
        done = true;
        break;
      }
      lua_pushinteger(L, (int64_t) c + 1);
      tk_vec_pushbase(L, m0->a[r * cols + c]);
      lua_settable(L, -3);
    }
    lua_settable(L, -3);
  }
}

static inline int tk_vec_pfx(each_lua_iter) (lua_State *L)
{
  lua_settop(L, 0);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_vec_base v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  tk_vec_pushbase(L, v);
  return 1;
}

static inline int tk_vec_pfx(ieach_lua_iter) (lua_State *L)
{
  lua_settop(L, 0);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_vec_base v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, (int64_t) n);
  tk_vec_pushbase(L, v);
  return 2;
}

static inline int tk_vec_pfx(each_lua) (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_vec_pfx(each_lua_iter), 2);
  return 1;
}

static inline int tk_vec_pfx(ieach_lua) (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_vec_pfx(ieach_lua_iter), 2);
  return 1;
}

#endif

static inline void tk_vec_pfx(shuffle) (tk_vec_pfx(t) *v, uint64_t s, uint64_t e) {
  if (s >= e || s >= v->n) return;
  if (e > v->n) e = v->n;
  uint64_t n = e - s;
  if (n <= 1) return;
  tk_vec_base *a = v->a + s;
  for (uint64_t i = n; i > 1; --i) {
    uint64_t j = tk_fast_index((unsigned int)i);
    tk_vec_base tmp = a[j];
    a[j] = a[i-1];
    a[i-1] = tmp;
  }
}

static inline void tk_vec_pfx(asc) (tk_vec_pfx(t) *v, uint64_t s, uint64_t e) {
  tk_vec_introsort(tk_vec_pfx(asc), e - s, v->a + s);
}

static inline void tk_vec_pfx(desc) (tk_vec_pfx(t) *v, uint64_t s, uint64_t e) {
  tk_vec_introsort(tk_vec_pfx(desc), e - s, v->a + s);
}

static inline uint64_t tk_vec_pfx(uasc) (tk_vec_pfx(t) *v, uint64_t s, uint64_t e) {
  tk_vec_introsort(tk_vec_pfx(asc), e - s, v->a + s);
  if (e - s == 0)
    return e;
  uint64_t write = s + 1;
  for (uint64_t read = s + 1; read < e; read ++)
    if (!tk_vec_eq(v->a[read - 1], v->a[read]))
      v->a[write ++] = v->a[read];
  return write;
}

static inline uint64_t tk_vec_pfx(udesc) (tk_vec_pfx(t) *v, uint64_t s, uint64_t e) {
  tk_vec_introsort(tk_vec_pfx(desc), e - s, v->a + s);
  if (e - s == 0)
    return e;
  uint64_t write = s + 1;
  for (uint64_t read = s + 1; read < e; read ++)
    if (!tk_vec_eq(v->a[read - 1], v->a[read]))
      v->a[write ++] = v->a[read];
  return write;
}

static inline void tk_vec_pfx(kasc) (tk_vec_pfx(t) *v, size_t k, uint64_t s, uint64_t e) {
  tk_vec_ksmall(tk_vec_pfx(asc), e - s, v->a + s, k);
}

static inline void tk_vec_pfx(kdesc) (tk_vec_pfx(t) *v, size_t k, uint64_t s, uint64_t e) {
  tk_vec_ksmall(tk_vec_pfx(desc), e - s, v->a + s, k);
}

#ifndef tk_vec_limited

static inline double tk_vec_pfx(magnitude) (
  tk_vec_pfx(t) *m0
) {
  tk_vec_base sum = 0.0;
  for (size_t i = 0; i < m0->n; i ++) {
    tk_vec_base val = m0->a[i];
    sum += val * val;
  }
  return sqrt(sum);
}

static inline tk_vec_base tk_vec_pfx(sum) (
  tk_vec_pfx(t) *m0
) {
  tk_vec_base sum = 0.0;
  for (size_t i = 0; i < m0->n; i ++) {
    tk_vec_base val = m0->a[i];
    sum += val;
  }
  return sum;
}

static inline tk_vec_base tk_vec_pfx(max) (
  tk_vec_pfx(t) *m0,
  uint64_t *colp
) {
  if (!m0->n)
    return 0;
  uint64_t maxcol = 0;
  tk_vec_base maxval = m0->a[0];
  for (size_t i = 1; i < m0->n; i ++) {
    if (m0->a[i] > maxval) {
      maxcol = i;
      maxval = m0->a[i];
    }
  }
  *colp = maxcol;
  return maxval;
}

static inline tk_vec_base tk_vec_pfx(min) (
  tk_vec_pfx(t) *m0,
  uint64_t *colp
) {
  if (!m0->n)
    return 0;
  uint64_t mincol = 0;
  tk_vec_base minval = m0->a[0];
  for (size_t i = 1; i < m0->n; i ++) {
    if (m0->a[i] < minval) {
      mincol = i;
      minval = m0->a[i];
    }
  }
  *colp = mincol;
  return minval;
}

#endif

static inline int tk_vec_pfx(destroy_lua) (lua_State *L)
{
  tk_vec_pfx(destroy)(tk_vec_pfx(peek)(L, 1, "vector"));
  return 0;
}

static inline int tk_vec_pfx(gc_lua) (lua_State *L)
{
  tk_vec_pfx(t) *v = tk_vec_pfx(peek)(L, 1, "vector");
  if (v->lua_managed) {
    tk_vec_pfx(destroy)(v);
  }
  return 0;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(load) (lua_State *L, FILE *fh)
{
  uint64_t n64;
  tk_lua_fread(L, &n64, sizeof(uint64_t), 1, fh);
  size_t n = (size_t) n64;
  tk_vec_pfx(t) *v = tk_vec_pfx(create)(L, n);
  tk_lua_fread(L, v->a, sizeof(tk_vec_base), n, fh);
  return v;
}

static inline int tk_vec_pfx(load_lua) (lua_State *L)
{
  lua_settop(L, 1);
  const char *data = luaL_checkstring(L, 1);
  FILE *fh = tk_lua_fopen(L, data, "r");
  tk_vec_pfx(load)(L, fh);
  tk_lua_fclose(L, fh);
  return 1;
}

static inline int tk_vec_pfx(from_raw_lua) (lua_State *L)
{
  lua_settop(L, 1);
  size_t len;
  const char *data = luaL_checklstring(L, 1, &len);
  size_t n = len / sizeof(tk_vec_base);
  if (n * sizeof(tk_vec_base) != len)
    return luaL_error(L, "from_raw: data length not a multiple of element size");
  tk_vec_pfx(t) *v = tk_vec_pfx(create)(L, n);
  memcpy(v->a, data, len);
  v->n = n;
  return 1;
}

static inline int tk_vec_pfx(create_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peekopt)(L, 1);
  uint64_t m = tk_lua_optunsigned(L, 2, "size", 0);
  if (m0 != NULL) {
    tk_vec_pfx(t) *m1 = tk_vec_pfx(create)(L, m > 0 ? m : m0->n);
    tk_vec_pfx(copy)(m1, m0, 0, (int64_t) m0->n, 0);
    return 1;
#ifndef tk_vec_limited
  } else if (lua_type(L, 1) == LUA_TTABLE) {
    tk_vec_pfx(t) *m0 = tk_vec_pfx(create)(L, m);
    uint64_t l = lua_objlen(L, 1);
    tk_vec_pfx(ensure)(m0, l);
    m0->n = l;
    for (uint64_t i = 0; i < l; i ++) {
      lua_pushinteger(L, (int64_t) i + 1);
      lua_gettable(L, 1);
      m0->a[i] = tk_vec_peekbase(L, -1);
      lua_pop(L, 1);
    }
    return 1;
#endif
  } else if (lua_type(L, 1) == LUA_TNUMBER || lua_type(L, 1) == LUA_TNIL) {
    uint64_t m = tk_lua_optunsigned(L, 1, "size", 0);
    tk_vec_pfx(create)(L, m);
    return 1;
  } else {
    tk_vec_err(L, create, 1, "expected either nil, a table, or a number");
    return 0;
  }
}

#ifndef tk_vec_limited
static inline int tk_vec_pfx(dot_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "a");
  tk_vec_pfx(t) *m1 = tk_vec_pfx(peek)(L, 2, "b");
  lua_pushnumber(L, tk_vec_pfx(dot)(m0, m1));
  return 1;
}

static inline int tk_vec_pfx(multiply_lua) (lua_State *L)
{
  lua_settop(L, 6);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "a");
  tk_vec_pfx(t) *m1 = tk_vec_pfx(peek)(L, 2, "b");
  tk_vec_pfx(t) *m2 = tk_vec_pfx(peek)(L, 3, "c");
  uint64_t k = tk_lua_checkunsigned(L, 4, "middle_dimension");
  bool transpose_a = lua_toboolean(L, 5);
  bool transpose_b = lua_toboolean(L, 6);
  tk_vec_pfx(multiply)(m0, m1, m2, k, transpose_a, transpose_b);
  return 0;
}
#endif

static inline int tk_vec_pfx(reverse_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "source");
  int64_t start, end;
  if (t == 1) {
    start = 0;
    end = (int64_t) m0->n;
  } else {
    start = tk_lua_checkinteger(L, 2, "start");
    end = tk_lua_checkinteger(L, 3, "end");
  }
  tk_vec_pfx(reverse)(m0, start, end);
  return 0;
}

static inline int tk_vec_pfx(copy_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "dest");
  tk_vec_pfx(t) *m1 = tk_vec_pfx(peek)(L, 2, "source");
  int64_t start, end, dest;
  if (t == 2) {
    start = 0;
    end = (int64_t) m1->n;
    dest = (int64_t) m0->n;
  } else if (t == 3) {
    start = 0;
    end = (int64_t) m1->n;
    dest = tk_lua_checkinteger(L, 3, "dest");
  } else {
    start = tk_lua_checkinteger(L, 3, "start");
    end = tk_lua_checkinteger(L, 4, "end");
    dest = tk_lua_checkinteger(L, 5, "dest");
  }
  tk_vec_pfx(copy)(m0, m1, start, end, dest);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(capacity_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  lua_pushinteger(L, (int64_t) tk_vec_pfx(capacity)(m0));
  return 1;
}

static inline int tk_vec_pfx(setn_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t n = tk_lua_checkunsigned(L, 2, "n");
  tk_vec_pfx(setn)(m0, n);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(size_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  lua_pushinteger(L, (int64_t) tk_vec_pfx(size)(m0));
  return 1;
}

static inline int tk_vec_pfx(resize_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t m = tk_lua_checkunsigned(L, 2, "size");
  tk_vec_pfx(resize)(m0, m, true);
  return 0;
}

static inline int tk_vec_pfx(ensure_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t m = tk_lua_checkunsigned(L, 2, "size");
  tk_vec_pfx(ensure)(m0, m);
  return 0;
}

static inline int tk_vec_pfx(clear_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  if (t == 1) {
    tk_vec_pfx(clear)(m0);
  } else {
    int64_t start = tk_lua_optinteger(L, 2, "start", 0);
    int64_t end = tk_lua_optinteger(L, 3, "end", (int64_t) m0->n);
    tk_vec_pfx(clear_range)(m0, start, end);
  }
  return 0;
}

static inline int tk_vec_pfx(zero_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(zero)(m0);
  return 1;
}

static inline int tk_vec_pfx(shrink_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(shrink)(m0);
  return 0;
}

static inline int tk_vec_pfx(transpose_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "dest");
  tk_vec_pfx(t) *m1 = tk_vec_pfx(peek)(L, 2, "source");
  uint64_t cols = tk_lua_checkunsigned(L, 3, "cols");
  tk_vec_pfx(transpose)(m0, m1, cols);
  return 0;
}

static inline int tk_vec_pfx(shuffle_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start = tk_lua_optunsigned(L, 2, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 3, "end", m0->n);
  tk_vec_pfx(shuffle)(m0, start, end);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(asc_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start = tk_lua_optunsigned(L, 2, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 3, "end", m0->n);
  tk_vec_pfx(asc)(m0, start, end);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(desc_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start = tk_lua_optunsigned(L, 2, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 3, "end", m0->n);
  tk_vec_pfx(desc)(m0, start, end);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(uasc_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start = tk_lua_optunsigned(L, 2, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 3, "end", m0->n);
  uint64_t end0 = tk_vec_pfx(uasc)(m0, start, end);
  if (start == 0 && end == m0->n)
    m0->n = end0;
  lua_pushinteger(L, (int64_t) end0);
  return 1;
}

static inline int tk_vec_pfx(udesc_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start = tk_lua_optunsigned(L, 2, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 3, "end", m0->n);
  uint64_t end0 = tk_vec_pfx(udesc)(m0, start, end);
  if (start == 0 && end == m0->n)
    m0->n = end0;
  lua_pushinteger(L, (int64_t) end0);
  return 1;
}

static inline int tk_vec_pfx(kasc_lua) (lua_State *L)
{
  lua_settop(L, 4);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t k = tk_lua_checkunsigned(L, 2, "k");
  uint64_t start = tk_lua_optunsigned(L, 3, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 4, "end", m0->n);
  tk_vec_pfx(kasc)(m0, k, start, end);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(kdesc_lua) (lua_State *L)
{
  lua_settop(L, 4);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t k = tk_lua_checkunsigned(L, 2, "k");
  uint64_t start = tk_lua_optunsigned(L, 3, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 4, "end", m0->n);
  tk_vec_pfx(kdesc)(m0, k, start, end);
  lua_settop(L, 1);
  return 1;
}

#ifdef tk_vec_pushbase
static inline int tk_vec_pfx(get_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t i = tk_lua_checkunsigned(L, 2, "idx");
  tk_vec_base x = tk_vec_pfx(get)(m0, i);
  tk_vec_pushbase(L, x);
  return 1;
}
#endif

#ifdef tk_vec_peekbase
static inline int tk_vec_pfx(find_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base x = tk_vec_peekbase(L, 2);
  int64_t idx = tk_vec_pfx(find)(m0, x);
  if (idx < 0) {
    lua_pushnil(L);
  } else {
    lua_pushinteger(L, idx);
  }
  return 1;
}
#endif

#ifdef tk_vec_peekbase
static inline int tk_vec_pfx(set_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t i = tk_lua_checkunsigned(L, 2, "idx");
  tk_vec_base x = tk_vec_peekbase(L, 3);
  tk_vec_pfx(set)(m0, i, x);
  return 0;
}
#endif

#ifdef tk_vec_pushbase

static inline int tk_vec_pfx(table_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start, end;
  if (t == 1) {
    start = 0;
    end =  m0->n;
  } else if (t == 3) {
    start = tk_lua_checkunsigned(L, 2, "start");
    end = tk_lua_checkunsigned(L, 3, "end");
  } else {
    tk_vec_err(L, table, 1, "expected either 1 or 3 arguments (vec and optionally start/end indices)");
    return 0;
  }
  tk_vec_pfx(table)(L, m0, start, end);
  return 1;
}

static inline int tk_vec_pfx(ctable_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  uint64_t start, end;
  if (t == 2) {
    start = 0;
    end = cols;
  } else if (t == 4) {
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, ctable, 1, "expected either 2 or 4 arguments (vec, cols and optionally start/end rows)");
    return 0;
  }
  tk_vec_pfx(ctable)(L, m0, cols, start, end);
  return 1;
}

static inline int tk_vec_pfx(rtable_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  uint64_t start, end;
  if (t == 2) {
    start = 0;
    end = (m0->n + cols - 1) / cols;
  } else if (t == 4) {
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, rtable, 1, "expected either 2 or 4 arguments (vec, cols and optionally start/end cols)");
    return 0;
  }
  tk_vec_pfx(rtable)(L, m0, cols, start, end);
  return 1;
}

#endif

#ifndef tk_vec_limited

static inline int tk_vec_pfx(push_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base x = tk_vec_peekbase(L, 2);
  tk_vec_pfx(push)(m0, x);
  return 0;
}

static inline int tk_vec_pfx(insert_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t idx = tk_lua_checkunsigned(L, 2, "idx");
  tk_vec_base x = tk_vec_peekbase(L, 3);
  tk_vec_pfx(insert)(m0, idx, x);
  return 0;
}

static inline int tk_vec_pfx(magnitude_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  lua_pushnumber(L, tk_vec_pfx(magnitude)(m0));
  return 1;
}

static inline int tk_vec_pfx(scale_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base scale;
  uint64_t start, end;
  if (t == 2) {
    scale = tk_vec_peekbase(L, 2);
    start = 0;
    end = m0->n;
  } else if (t == 4) {
    scale = tk_vec_peekbase(L, 2);
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, scale, 1, "expected either 2 or 4 arguments (vec, scale, or vec, scale, start, end)");
    return 0;
  }
  tk_vec_pfx(scale)(m0, scale, start, end);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(add_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base add;
  uint64_t start, end;
  if (t == 2) {
    add = tk_vec_peekbase(L, 2);
    start = 0;
    end = m0->n;
  } else if (t == 4) {
    add = tk_vec_peekbase(L, 2);
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, add, 1, "expected either 2 or 4 arguments (vec, add or vec, add, start, end)");
    return 0;
  }
  tk_vec_pfx(add)(m0, add, start, end);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(add_scaled_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base add;
  uint64_t start, end;
  if (t == 2) {
    add = tk_vec_peekbase(L, 2);
    start = 0;
    end = m0->n;
  } else if (t == 4) {
    add = tk_vec_peekbase(L, 2);
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, add_scaled, 1, "expected 2 or 4 arguments (vec, x or vec, x, start, end); adds x*i at each index i");
    return 0;
  }
  tk_vec_pfx(add_scaled)(m0, add, start, end);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(scalev_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(t) *m1 = tk_vec_pfx(peek)(L, 2, "scale");
  uint64_t start, end;
  if (t == 2) {
    start = 0;
    end = m0->n;
  } else if (t == 4) {
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, scalev, 1, "expected either 2 or 4 arguments (vec, scale, or vec, scale, start, end)");
    return 0;
  }
  tk_vec_pfx(scalev)(m0, m1, start, end);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(addv_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(t) *m1 = tk_vec_pfx(peek)(L, 2, "vector");
  uint64_t start, end;
  if (t == 2) {
    start = 0;
    end = m0->n;
  } else if (t == 4) {
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, addv, 1, "expected either 2 or 4 arguments (vec, add or vec, add, start, end)");
    return 0;
  }
  tk_vec_pfx(addv)(m0, m1, start, end);
  lua_settop(L, 1);
  return 1;
}

#ifdef tk_vec_abs
static inline int tk_vec_pfx(abs_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start, end;
  if (t == 1) {
    start = 0;
    end = m0->n;
  } else if (t == 2) {
    start = tk_lua_checkunsigned(L, 2, "start");
    end = tk_lua_checkunsigned(L, 3, "end");
  } else {
    tk_vec_err(L, abs, 1, "expected either 1 or 3 arguments (vec or vec, start, end)");
    return 0;
  }
  tk_vec_pfx(abs)(m0, start, end);
  lua_settop(L, 1);
  return 1;
}
#endif

static inline int tk_vec_pfx(exp_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start, end;
  if (t == 1) {
    start = 0;
    end = m0->n;
  } else if (t == 2) {
    start = tk_lua_checkunsigned(L, 2, "start");
    end = tk_lua_checkunsigned(L, 3, "end");
  } else {
    tk_vec_err(L, abs, 1, "expected either 1 or 3 arguments (vec or vec, start, end)");
    return 0;
  }
  tk_vec_pfx(exp)(m0, start, end);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(log_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start, end;
  if (t == 1) {
    start = 0;
    end = m0->n;
  } else if (t == 2) {
    start = tk_lua_checkunsigned(L, 2, "start");
    end = tk_lua_checkunsigned(L, 3, "end");
  } else {
    tk_vec_err(L, abs, 1, "expected either 1 or 3 arguments (vec or vec, start, end)");
    return 0;
  }
  tk_vec_pfx(log)(m0, start, end);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(pow_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  double e;
  uint64_t start, end;
  if (t == 2) {
    e = tk_lua_checkdouble(L, 2, "e");
    start = 0;
    end = m0->n;
  } else if (t == 4) {
    e = tk_lua_checkdouble(L, 2, "e");
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, pow, 1, "expected either 2 or 4 arguments (vec, e or vec, e, start, end)");
    return 0;
  }
  tk_vec_pfx(pow)(m0, e, start, end);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(min_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  if (!m0->n)
    return 0;
  uint64_t minc;
  tk_vec_base minv = tk_vec_pfx(min)(m0, &minc);
  tk_vec_pushbase(L, minv);
  lua_pushinteger(L, (int64_t) minc);
  return 2;
}

static inline int tk_vec_pfx(rmins_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rmins)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(cmins_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cmins)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(max_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  if (!m0->n)
    return 0;
  uint64_t maxc;
  tk_vec_base maxv = tk_vec_pfx(max)(m0, &maxc);
  tk_vec_pushbase(L, maxv);
  lua_pushinteger(L, (int64_t) maxc);
  return 2;
}


static inline int tk_vec_pfx(cumsum_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, m0->n + 1);
  out->n = m0->n + 1;
  tk_vec_base acc = 0;
  out->a[0] = 0;
  for (uint64_t i = 0; i < m0->n; i ++) {
    acc += m0->a[i];
    out->a[i + 1] = acc;
  }
  return 1;
}

static inline int tk_vec_pfx(rmaxs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rmaxs)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(cmaxs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cmaxs)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(sum_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pushbase(L, tk_vec_pfx(sum)(m0));
  return 1;
}

static inline int tk_vec_pfx(csums_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(csums)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(rsums_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rsums)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(fill_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *v = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base x = tk_vec_peekbase(L, 2);
  uint64_t start, end;
  if (t == 2) {
    start = 0;
    end = v->n;
  } else if (t == 4) {
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, fill, 1, "expected either 2 or 4 arguments (vec, value or vec, value, start, end)");
    return 0;
  }
  tk_vec_pfx(fill)(v, x, start, end);
  lua_settop(L, 1);
  return 1;
}

static inline int tk_vec_pfx(fill_indices_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *v = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(fill_indices)(v);
  return 1;
}

#endif

static inline int tk_vec_pfx(raw_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  size_t len = sizeof(tk_vec_base) * m0->n;
  lua_pushlstring(L, (char *) m0->a, len);
  return 1;
}

#if !defined(__EMSCRIPTEN__)
static inline int tk_vec_pfx(mmap_sync_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *v = tk_vec_pfx(peek)(L, 1, "vector");
  if (v->lua_managed == 2 && v->a) {
    if (msync(v->a, v->m * sizeof(tk_vec_base), MS_SYNC) != 0)
      return tk_vec_err(L, mmap_sync, 1, strerror(errno));
  }
  return 0;
}
#endif

static luaL_Reg tk_vec_pfx(lua_mt_fns)[] =
{
  { "copy", tk_vec_pfx(copy_lua) },
  { "reverse", tk_vec_pfx(reverse_lua) },
  { "destroy", tk_vec_pfx(destroy_lua) },
  { "size", tk_vec_pfx(size_lua) },
  { "capacity", tk_vec_pfx(capacity_lua) },
  { "resize", tk_vec_pfx(resize_lua) },
  { "setn", tk_vec_pfx(setn_lua) },
  { "ensure", tk_vec_pfx(ensure_lua) },
  { "shrink", tk_vec_pfx(shrink_lua) },
  { "clear", tk_vec_pfx(clear_lua) },
  { "zero", tk_vec_pfx(zero_lua) },
  { "raw", tk_vec_pfx(raw_lua) },
#if !defined(__EMSCRIPTEN__)
  { "mmap_sync", tk_vec_pfx(mmap_sync_lua) },
#endif

#ifdef tk_vec_peekbase
  { "get", tk_vec_pfx(get_lua) },
  { "find", tk_vec_pfx(find_lua) },
#endif

#ifdef tk_vec_pushbase
  { "set", tk_vec_pfx(set_lua) },
  { "table", tk_vec_pfx(table_lua) },
  { "ctable", tk_vec_pfx(ctable_lua) },
  { "rtable", tk_vec_pfx(rtable_lua) },
  { "each", tk_vec_pfx(each_lua) },
  { "ieach", tk_vec_pfx(ieach_lua) },
#endif

  { "shuffle", tk_vec_pfx(shuffle_lua) },
  { "asc", tk_vec_pfx(asc_lua) },
  { "desc", tk_vec_pfx(desc_lua) },
  { "uasc", tk_vec_pfx(uasc_lua) },
  { "udesc", tk_vec_pfx(udesc_lua) },
  { "kasc", tk_vec_pfx(kasc_lua) },
  { "kdesc", tk_vec_pfx(kdesc_lua) },

#ifndef tk_vec_limited
  { "push", tk_vec_pfx(push_lua) },
  { "insert", tk_vec_pfx(insert_lua) },
  { "add", tk_vec_pfx(add_lua) },
  { "add_scaled", tk_vec_pfx(add_scaled_lua) },
  { "scale", tk_vec_pfx(scale_lua) },
  { "addv", tk_vec_pfx(addv_lua) },
  { "scalev", tk_vec_pfx(scalev_lua) },
  { "exp", tk_vec_pfx(exp_lua) },
  { "log", tk_vec_pfx(log_lua) },
  { "pow", tk_vec_pfx(pow_lua) },

#ifdef tk_vec_abs
  { "abs", tk_vec_pfx(abs_lua) },
#endif
  { "dot", tk_vec_pfx(dot_lua) },
  { "magnitude", tk_vec_pfx(magnitude_lua) },
  { "sum", tk_vec_pfx(sum_lua) },
  { "max", tk_vec_pfx(max_lua) },
  { "min", tk_vec_pfx(min_lua) },
  { "cumsum", tk_vec_pfx(cumsum_lua) },
  { "fill", tk_vec_pfx(fill_lua) },
  { "fill_indices", tk_vec_pfx(fill_indices_lua) },
#endif

  { NULL, NULL }
};

static inline void tk_vec_pfx(suppress_unused_lua_mt_fns) (void)
  { (void) tk_vec_pfx(lua_mt_fns); }

#ifdef tk_vec_module
static inline void tk_vec_pfx(ensure_init) (lua_State *L)
{
  lua_getmetatable(L, -1);
  lua_getfield(L, -1, "__ext_init");
  if (lua_isnil(L, -1)) {
    lua_pop(L, 1);
    lua_pushboolean(L, 1);
    lua_setfield(L, -2, "__ext_init");
    lua_pop(L, 1);
    lua_getfield(L, LUA_REGISTRYINDEX, "_LOADED");
    lua_getfield(L, -1, tk_vec_module);
    if (lua_isnil(L, -1)) {
      lua_pop(L, 2);
      lua_getglobal(L, "require");
      lua_pushstring(L, tk_vec_module);
      lua_pcall(L, 1, 1, 0);
      lua_pop(L, 1);
    } else {
      lua_pop(L, 2);
    }
  } else {
    lua_pop(L, 2);
  }
}
#endif

static inline tk_vec_pfx(t) *tk_vec_pfx(create) (lua_State *L, size_t n)
{
  tk_vec_pfx(t) *v = L == NULL
    ? malloc(sizeof(tk_vec_pfx(t)))
    : tk_lua_newuserdata(L, tk_vec_pfx(t), tk_vec_mt, tk_vec_pfx(lua_mt_fns), tk_vec_pfx(gc_lua));
  if (v == NULL)
    return v;
#ifdef tk_vec_module
  if (L != NULL)
    tk_vec_pfx(ensure_init)(L);
#endif
  tk_vec_pfx(t) v0 = *v;
  bool lua_managed = L != NULL;
  kv_init(v0, lua_managed);
  if (kv_resize(tk_vec_base, v0, n) != 0) {
    if (L)
      tk_error(L, "vec_create resize", ENOMEM);
    else if (!lua_managed)
      free(v);
    return NULL;
  }
  v0.n = n;
  *v = v0;
  return v;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(register) (lua_State *L, tk_vec_pfx(t) *v)
{
  tk_vec_pfx(t) *x = tk_lua_newuserdata(L, tk_vec_pfx(t), tk_vec_mt, tk_vec_pfx(lua_mt_fns), tk_vec_pfx(gc_lua));
  *x = *v;
  free(v);
#ifdef tk_vec_module
  tk_vec_pfx(ensure_init)(L);
#endif
  return x;
}

#if !defined(__EMSCRIPTEN__)

static inline tk_vec_pfx(t) *tk_vec_pfx(mmap_create) (lua_State *L, const char *path, size_t n)
{
  size_t map_len = n * sizeof(tk_vec_base);
  int fd = open(path, O_RDWR | O_CREAT | O_TRUNC, 0666);
  if (fd < 0)
    return (tk_vec_pfx(t) *)(intptr_t) tk_vec_err(L, mmap_create, 1, strerror(errno));
  if (map_len > 0 && ftruncate(fd, (off_t) map_len) != 0) {
    close(fd);
    return (tk_vec_pfx(t) *)(intptr_t) tk_vec_err(L, mmap_create, 1, strerror(errno));
  }
  tk_vec_pfx(t) *v = tk_lua_newuserdata(L, tk_vec_pfx(t), tk_vec_mt, tk_vec_pfx(lua_mt_fns), tk_vec_pfx(gc_lua));
#ifdef tk_vec_module
  tk_vec_pfx(ensure_init)(L);
#endif
  v->lua_managed = 2;
  if (map_len > 0) {
    v->a = (tk_vec_base *) mmap(NULL, map_len, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (v->a == MAP_FAILED) {
      close(fd);
      return (tk_vec_pfx(t) *)(intptr_t) tk_vec_err(L, mmap_create, 1, strerror(errno));
    }
  } else {
    v->a = NULL;
  }
  close(fd);
  v->n = n;
  v->m = n;
  return v;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(mmap_open) (lua_State *L, const char *path)
{
  int fd = open(path, O_RDWR);
  if (fd < 0)
    return (tk_vec_pfx(t) *)(intptr_t) tk_vec_err(L, mmap_open, 1, strerror(errno));
  struct stat st;
  if (fstat(fd, &st) != 0) {
    close(fd);
    return (tk_vec_pfx(t) *)(intptr_t) tk_vec_err(L, mmap_open, 1, strerror(errno));
  }
  size_t map_len = (size_t) st.st_size;
  size_t n = map_len / sizeof(tk_vec_base);
  tk_vec_pfx(t) *v = tk_lua_newuserdata(L, tk_vec_pfx(t), tk_vec_mt, tk_vec_pfx(lua_mt_fns), tk_vec_pfx(gc_lua));
#ifdef tk_vec_module
  tk_vec_pfx(ensure_init)(L);
#endif
  v->lua_managed = 2;
  if (map_len > 0) {
    v->a = (tk_vec_base *) mmap(NULL, map_len, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (v->a == MAP_FAILED) {
      close(fd);
      return (tk_vec_pfx(t) *)(intptr_t) tk_vec_err(L, mmap_open, 1, strerror(errno));
    }
  } else {
    v->a = NULL;
  }
  close(fd);
  v->n = n;
  v->m = n;
  return v;
}

static inline int tk_vec_pfx(mmap_create_lua) (lua_State *L)
{
  lua_settop(L, 2);
  const char *path = luaL_checkstring(L, 1);
  uint64_t n = (uint64_t)luaL_checkinteger(L, 2);
  tk_vec_pfx(mmap_create)(L, path, n);
  return 1;
}

static inline int tk_vec_pfx(mmap_open_lua) (lua_State *L)
{
  lua_settop(L, 1);
  const char *path = luaL_checkstring(L, 1);
  tk_vec_pfx(mmap_open)(L, path);
  return 1;
}

#endif

static luaL_Reg tk_vec_pfx(lua_fns)[] =
{
  { "create", tk_vec_pfx(create_lua) },
  { "load", tk_vec_pfx(load_lua) },
  { "from_raw", tk_vec_pfx(from_raw_lua) },
#if !defined(__EMSCRIPTEN__)
  { "mmap_create", tk_vec_pfx(mmap_create_lua) },
  { "mmap_open", tk_vec_pfx(mmap_open_lua) },
#endif
  { NULL, NULL }
};

static inline void tk_vec_pfx(suppress_unused_lua_fns) (void)
  { (void) tk_vec_pfx(lua_fns); }


#define TK_GENERATE_SINGLE
#include <santoku/parallel/tpl.h>
#include <santoku/vec/tpl_para.h>
#undef TK_GENERATE_SINGLE

#include <santoku/parallel/tpl.h>
#include <santoku/vec/tpl_para.h>

#include <santoku/vec/undef.h>
