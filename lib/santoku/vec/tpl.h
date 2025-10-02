#include <santoku/lua/utils.h>
#include <santoku/klib.h>

#ifndef tk_vec_lt
#define tk_vec_lt(a, b) ((a) < (b))
#endif

#ifndef tk_vec_gt
#define tk_vec_gt(a, b) ((a) > (b))
#endif

#ifndef tk_vec_ltx
#define tk_vec_ltx(a, b) (memcmp(&(a), &(b), sizeof(a)) < 0)
#endif

#ifndef tk_vec_gtx
#define tk_vec_gtx(a, b) (memcmp(&(a), &(b), sizeof(a)) > 0)
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
tk_vec_ksort(tk_vec_pfx(xasc), tk_vec_base, tk_vec_ltx)
tk_vec_ksort(tk_vec_pfx(xdesc), tk_vec_base, tk_vec_gtx)

#define tk_vec_shuffle(...) ks_shuffle(__VA_ARGS__)
#define tk_vec_introsort(...) ks_introsort(__VA_ARGS__)
#define tk_vec_ksmall(...) ks_ksmall(__VA_ARGS__)

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
  kv_destroy(*r);
  if (!lua_managed)
    free(r);
  else {
    memset(r, 0, sizeof(tk_vec_pfx(t)));
    r->lua_managed = lua_managed;
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

static inline tk_vec_pfx(t) *tk_vec_pfx(create) (lua_State *, size_t, tk_vec_base *, tk_vec_pfx(t) *);

static inline tk_vec_pfx(t) *tk_vec_pfx(resize) (
  tk_vec_pfx(t) *m0,
  size_t m,
  bool setn
) {
  tk_vec_pfx(t) v = *m0;
  if (m == 0) {
    free(v.a);
    v.a = NULL;
    v.m = v.n = 0;
    *m0 = v;
    return m0;
  } else {
    kv_resize(tk_vec_base, v, m);
    if (m0->a == NULL)
      return NULL;
    v.m = m;
    if (setn)
      v.n = m;
    *m0 = v;
    return m0;
  }
}

static inline void tk_vec_pfx(setn) (
  tk_vec_pfx(t) *m0,
  size_t n
) {
  if (n > m0->m)
    tk_vec_pfx(resize)(m0, n, true);
  else
    m0->n = n;
}

static inline void tk_vec_pfx(ensure) (
  tk_vec_pfx(t) *m0,
  size_t m
) {
  if (m0->m >= m)
    return;
  tk_vec_pfx(resize)(m0, m, false);
}

static inline void tk_vec_pfx(copy) (
  tk_vec_pfx(t) *m0,
  tk_vec_pfx(t) *m1,
  int64_t start,
  int64_t end,
  int64_t dest
) {
  if (start < 0 || start >= end || start >= (int64_t) m1->n)
    return;
  if (end >= (int64_t) m1->n)
    end = (int64_t) m1->n;
  uint64_t m = (uint64_t) dest + (uint64_t) (end - start);
  tk_vec_pfx(ensure)(m0, m);
  if (m0 == m1)
    memmove(m0->a + dest, m1->a + start, sizeof(tk_vec_base) * (uint64_t) (end - start));
  else
    memcpy(m0->a + dest, m1->a + start, sizeof(tk_vec_base) * (uint64_t) (end - start));
  if (m0->n < m)
    m0->n = m;
}

static inline void tk_vec_pfx(transpose) (
  tk_vec_pfx(t) *m0,
  tk_vec_pfx(t) *m1,
  uint64_t cols
) {
  tk_vec_pfx(ensure)(0, m1->n);
  uint64_t rows = m1->n / cols;
  for (uint64_t r = 0; r < rows; r ++)
    for (uint64_t c = 0; c < cols; c ++)
      m0->a[c * rows + r] = m1->a[r * cols + c];
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
static inline void tk_vec_pfx(push) (tk_vec_pfx(t) *v, tk_vec_base x) {
  tk_vec_pfx(t) v0 = *v;
  kv_push(tk_vec_base, v0, x);
  *v = v0;
}

static inline void tk_vec_pfx(insert) (
  tk_vec_pfx(t) *v,
  uint64_t idx,
  tk_vec_base x
) {
  if (idx > v->n) idx = v->n;
  tk_vec_pfx(ensure)(v, v->n + 1);
  if (idx < v->n)
    memmove(v->a + idx + 1, v->a + idx, (v->n - idx) * sizeof(tk_vec_base));
  v->a[idx] = x;
  v->n++;
}

static inline void tk_vec_pfx(set) (tk_vec_pfx(t) *v, uint64_t i, tk_vec_base x) {
  tk_vec_pfx(ensure)(v, i + 1);
  v->a[i] = x;
  if (i + 1 > v->n)
    v->n = i + 1;
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

static inline void tk_vec_pfx(shuffle) (tk_vec_pfx(t) *v) {
  tk_vec_shuffle(tk_vec_pfx(asc), v->n, v->a);
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

static inline uint64_t tk_vec_pfx(xasc) (tk_vec_pfx(t) *v, uint64_t s, uint64_t e) {
  tk_vec_introsort(tk_vec_pfx(xasc), e - s, v->a + s);
  if (e - s == 0)
    return e;
  uint64_t write = s + 1;
  for (uint64_t read = s + 1; read < e; read ++)
    if (!tk_vec_eqx(v->a[read - 1], v->a[read]))
      v->a[write ++] = v->a[read];
  return write;
}

static inline uint64_t tk_vec_pfx(xdesc) (tk_vec_pfx(t) *v, uint64_t s, uint64_t e) {
  tk_vec_introsort(tk_vec_pfx(xdesc), e - s, v->a + s);
  if (e - s == 0)
    return e;
  uint64_t write = s + 1;
  for (uint64_t read = s + 1; read < e; read ++)
    if (!tk_vec_eqx(v->a[read - 1], v->a[read]))
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

static inline tk_vec_pfx(t) *tk_vec_pfx(csums) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, cols, NULL, NULL);
  for (uint64_t c = 0; c < cols; c ++) {
    tk_vec_base sum = 0;
    for (uint64_t r = 0; r < m0->n / cols; r ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val;
    }
    out->a[c] = sum;
  }
  return out;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(rsums) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, m0->n / cols, NULL, NULL);
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val;
    }
    out->a[r] = sum;
  }
  return out;
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

static inline tk_vec_pfx(t) *tk_vec_pfx(cmaxs) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, cols, NULL, NULL);
  for (uint64_t c = 0; c < cols; c ++) {
    tk_vec_base maxv = m0->a[0 * cols + c];
    for (size_t r = 1; r < m0->n / cols; r ++) {
      if (m0->a[r * cols + c] > maxv) {
        maxv = m0->a[r * cols + c];
      }
    }
    out->a[c] = maxv;
  }
  return out;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(rmaxs) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, m0->n / cols, NULL, NULL);
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val;
    }
    out->a[r] = sum;
    tk_vec_base maxv = m0->a[r * cols + 0];
    for (uint64_t c = 1; c < cols; c ++) {
      if (m0->a[r * cols + c] > maxv) {
        maxv = m0->a[r * cols + c];
      }
    }
    out->a[r] = maxv;
  }
  return out;
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

static inline tk_vec_pfx(t) *tk_vec_pfx(cmins) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, cols, NULL, NULL);
  for (uint64_t c = 0; c < cols; c ++) {
    tk_vec_base minv = m0->a[0 * cols + c];
    for (size_t r = 1; r < m0->n / cols; r ++) {
      if (m0->a[r * cols + c] < minv) {
        minv = m0->a[r * cols + c];
      }
    }
    out->a[c] = minv;
  }
  return out;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(rmins) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, m0->n / cols, NULL, NULL);
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val;
    }
    out->a[r] = sum;
    tk_vec_base minv = m0->a[r * cols + 0];
    for (uint64_t c = 1; c < cols; c ++) {
      if (m0->a[r * cols + c] < minv) {
        minv = m0->a[r * cols + c];
      }
    }
    out->a[r] = minv;
  }
  return out;
}

static inline void tk_vec_pfx(scale) (
  tk_vec_pfx(t) *m0,
  tk_vec_base scale,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start >= end || start >= m0->n)
    return;
  if (end > m0->n)
    end = m0->n;
  for (size_t i = start; i < end; i ++)
    m0->a[i] *= scale;
}

static inline void tk_vec_pfx(add) (
  tk_vec_pfx(t) *m0,
  tk_vec_base add,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start >= end || start >= m0->n)
    return;
  if (end > m0->n)
    end = m0->n;
  for (size_t i = start; i < end; i ++) {
    m0->a[i] += add;
  }
}

static inline void tk_vec_pfx(add_scaled) (
  tk_vec_pfx(t) *m0,
  tk_vec_base add,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start >= end || start >= m0->n)
    return;
  if (end > m0->n)
    end = m0->n;
  for (size_t i = start; i < end; i ++) {
    m0->a[i] += add * (tk_vec_base) i;
  }
}

static inline void tk_vec_pfx(scalev) (
  tk_vec_pfx(t) *m0,
  tk_vec_pfx(t) *m1,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || !m1->n)
    return;
  if (start >= end)
    return;
  if (start >= m0->n || start >= m1->n)
    return;
  uint64_t lim0 = m0->n;
  uint64_t lim1 = m1->n;
  uint64_t lim = lim0 < lim1 ? lim0 : lim1;
  if (end > lim)
    end = lim;
  for (size_t i = start; i < end; i ++) {
    m0->a[i] *= m1->a[i];
  }
}

static inline void tk_vec_pfx(addv) (
  tk_vec_pfx(t) *m0,
  tk_vec_pfx(t) *m1,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || !m1->n)
    return;
  if (start >= end)
    return;
  if (start >= m0->n || start >= m1->n)
    return;
  uint64_t lim0 = m0->n;
  uint64_t lim1 = m1->n;
  uint64_t lim = lim0 < lim1 ? lim0 : lim1;
  if (end > lim)
    end = lim;
  for (size_t i = start; i < end; i ++) {
    m0->a[i] += m1->a[i];
  }
}

#ifdef tk_vec_abs
static inline void tk_vec_pfx(abs) (
  tk_vec_pfx(t) *m0,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start >= end || start >= m0->n)
    return;
  if (end > m0->n)
    end = m0->n;
  for (size_t i = start; i < end; i ++)
    m0->a[i] = tk_vec_abs(m0->a[i]);
}
#endif

static inline double tk_vec_pfx(dot) (tk_vec_pfx(t) *a, tk_vec_pfx(t) *b) {
  if (a->n != b->n)
    return 0;
  size_t n = a->n;
  double sum = 0.0;
  for (size_t i = 0; i < n; i ++)
    sum += a->a[i] * b->a[i];
  return sum;
}

static inline void tk_vec_pfx(fill) (tk_vec_pfx(t) *v, tk_vec_base x, uint64_t start, uint64_t end) {
  if (end <= start)
    return;
  if (end > v->n) {
    tk_vec_pfx(ensure)(v, end);
    v->n = end;
  }
  for (uint64_t i = start; i < end; i ++)
    v->a[i] = x;
}

static inline void tk_vec_pfx(multiply) (tk_vec_pfx(t) *a, tk_vec_pfx(t) *b, tk_vec_pfx(t) *c, uint64_t k, bool transpose_a, bool transpose_b) {
  size_t m = transpose_a ? k : a->n / k;
  size_t n = transpose_b ? k : b->n / k;
  tk_vec_pfx(ensure)(a, transpose_a ? k * m : m * k);
  tk_vec_pfx(ensure)(b, transpose_b ? k * n : k * n);
  tk_vec_pfx(ensure)(c, m * n);
  for (size_t i = 0; i < m; i ++) {
    for (size_t j = 0; j < n; j ++) {
      double sum = 0.0;
      for (size_t l = 0; l < k; l ++) {
        size_t a_idx = transpose_a ? l * m + i : i * k + l;
        size_t b_idx = transpose_b ? j * k + l : l * n + j;
        sum += a->a[a_idx] * b->a[b_idx];
      }
      c->a[i * n + j] = sum;
    }
  }
}

static inline void tk_vec_pfx(pow) (
  tk_vec_pfx(t) *m0,
  double e,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start > end || start >= m0->n)
    return;
  if (end > m0->n - 1)
    end = m0->n - 1;
  for (size_t i = start; i <= end; i ++)
    m0->a[i] = pow(m0->a[i], e);
}

static inline void tk_vec_pfx(log) (
  tk_vec_pfx(t) *m0,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start > end || start >= m0->n)
    return;
  if (end > m0->n - 1)
    end = m0->n - 1;
  for (size_t i = start; i <= end; i ++)
    m0->a[i] = log(m0->a[i]);
}

static inline void tk_vec_pfx(exp) (
  tk_vec_pfx(t) *m0,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start > end || start >= m0->n)
    return;
  if (end > m0->n - 1)
    end = m0->n - 1;
  for (size_t i = start; i <= end; i ++)
    m0->a[i] = exp(m0->a[i]);
}

static inline const char *tk_vec_pfx(raw) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  const char *fmt
) {
  if (fmt == NULL) {
    lua_pushlstring(L, (char *) m0->a, sizeof(tk_vec_base) * m0->n);
    return lua_tostring(L, -1);
  }
  void *out = NULL;
  size_t out_size = 0;
  if (!strcmp(fmt, "u32")) {
    out_size = sizeof(uint32_t);
    out = malloc(out_size * m0->n);
    if (!out) goto err_mem;
    for (size_t i = 0; i < m0->n; i ++) {
      tk_vec_base r = m0->a[i];
      if (r < 0)
        goto err_negative;
      if (r > (tk_vec_base) UINT32_MAX)
        goto err_toobig;
      ((uint32_t *) out)[i] = (uint32_t) r;
    }
  } else if (!strcmp(fmt, "u64")) {
    out_size = sizeof(uint64_t);
    out = malloc(out_size * m0->n);
    if (!out) goto err_mem;
    for (size_t i = 0; i < m0->n; i ++) {
      tk_vec_base r = m0->a[i];
      if (r < 0)
        goto err_negative;
      if (r > (tk_vec_base) UINT64_MAX)
        goto err_toobig;
      ((uint64_t *) out)[i] = (uint64_t) r;
    }
  } else if (!strcmp(fmt, "f32")) {
    out_size = sizeof(float);
    out = malloc(out_size * m0->n);
    if (!out)
      goto err_mem;
    for (size_t i = 0; i < m0->n; i ++)
      ((float *) out)[i] = (float) m0->a[i];
  } else if (!strcmp(fmt, "f64")) {
    out_size = sizeof(double);
    out = malloc(out_size * m0->n);
    if (!out)
      goto err_mem;
    for (size_t i = 0; i < m0->n; i ++)
      ((double *) out)[i] = (double) m0->a[i];
  } else {
    goto err_format;
  }
  lua_pushlstring(L, (char *) out, out_size * m0->n);
  free(out);
  return lua_tostring(L, -1);
err_format:
  tk_vec_err(L, raw, 2, "bad format", fmt);
  return NULL;
err_mem:
  tk_vec_err(L, raw, 2, "malloc error");
  return NULL;
err_negative:
  tk_vec_err(L, raw, 2, "unexpected negative");
  goto cleanup;
err_toobig:
  tk_vec_err(L, raw, 2, "number too large");
cleanup:
  free(out);
  return NULL;
}

static inline void tk_vec_pfx(fill_indices) (tk_vec_pfx(t) *v)
{
  for (uint64_t i = 0; i < v->n; i ++)
    v->a[i] = (tk_vec_base) i;
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
  size_t n;
  tk_lua_fread(L, &n, sizeof(size_t), 1, fh);
  tk_vec_pfx(t) *v = tk_vec_pfx(create)(L, n, 0, 0);
  tk_lua_fread(L, v->a, sizeof(tk_vec_base), n, fh);
  return v;
}

static inline int tk_vec_pfx(load_lua) (lua_State *L)
{
  lua_settop(L, 2);
  size_t len;
  const char *data = luaL_checklstring(L, 1, &len);
  bool isstr = lua_type(L, 2) == LUA_TBOOLEAN && lua_toboolean(L, 2);
  FILE *fh = isstr ? tk_lua_fmemopen(L, (char *) data, len, "r") : tk_lua_fopen(L, data, "r");
  tk_vec_pfx(load)(L, fh);
  tk_lua_fclose(L, fh);
  return 1;
}

static inline int tk_vec_pfx(create_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peekopt)(L, 1);
  uint64_t m = tk_lua_optunsigned(L, 2, "size", 0);
  if (m0 != NULL) {
    tk_vec_pfx(t) *m1 = tk_vec_pfx(create)(L, m > 0 ? m : m0->n, NULL, NULL);
    tk_vec_pfx(copy)(m1, m0, 0, (int64_t) m0->n, 0);
    return 1;
#ifndef tk_vec_limited
  } else if (lua_type(L, 1) == LUA_TTABLE) {
    tk_vec_pfx(t) *m0 = tk_vec_pfx(create)(L, m, NULL, NULL);
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
    tk_vec_pfx(create)(L, m, NULL, NULL);
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
  return 0;
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
  return 0;
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
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(clear)(m0);
  return 0;
}

static inline int tk_vec_pfx(zero_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(zero)(m0);
  return 0;
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
  lua_settop(L, 1);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(shuffle)(m0);
  return 1;
}

static inline int tk_vec_pfx(asc_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start = tk_lua_optunsigned(L, 2, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 3, "end", m0->n);
  tk_vec_pfx(asc)(m0, start, end);
  return 1;
}

static inline int tk_vec_pfx(desc_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start = tk_lua_optunsigned(L, 2, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 3, "end", m0->n);
  tk_vec_pfx(desc)(m0, start, end);
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

static inline int tk_vec_pfx(xasc_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start = tk_lua_optunsigned(L, 2, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 3, "end", m0->n);
  uint64_t end0 = tk_vec_pfx(xasc)(m0, start, end);
  if (start == 0 && end == m0->n)
    m0->n = end0;
  lua_pushinteger(L, (int64_t) end0);
  return 1;
}

static inline int tk_vec_pfx(xdesc_lua) (lua_State *L)
{
  lua_settop(L, 3);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t start = tk_lua_optunsigned(L, 2, "start", 0);
  uint64_t end = tk_lua_optunsigned(L, 3, "end", m0->n);
  uint64_t end0 = tk_vec_pfx(xdesc)(m0, start, end);
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
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
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
    scale = tk_lua_checkdouble(L, 2, "scale");
    start = 0;
    end = m0->n;
  } else if (t == 3) {
    scale = tk_lua_checkdouble(L, 2, "scale");
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, scale, 1, "expected either 2 or 4 arguments (vec, scale, or vec, scale, start, end)");
    return 0;
  }
  tk_vec_pfx(scale)(m0, scale, start, end);
  return 1;
}

static inline int tk_vec_pfx(add_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base add;
  uint64_t start, end;
  if (t == 2) {
    add = tk_lua_checkdouble(L, 2, "add");
    start = 0;
    end = m0->n;
  } else if (t == 4) {
    add = tk_lua_checkdouble(L, 2, "add");
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, add, 1, "expected either 2 or 4 arguments (vec, add or vec, add, start, end)");
    return 0;
  }
  tk_vec_pfx(add)(m0, add, start, end);
  return 1;
}

static inline int tk_vec_pfx(add_scaled_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base add;
  uint64_t start, end;
  if (t == 2) {
    add = tk_lua_checkdouble(L, 2, "add");
    start = 0;
    end = m0->n;
  } else if (t == 4) {
    add = tk_lua_checkdouble(L, 2, "add");
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, add_scaled, 1, "expected either 2 or 4 arguments (vec, add or vec, add, start, end)");
    return 0;
  }
  tk_vec_pfx(add_scaled)(m0, add, start, end);
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
  } else if (t == 3) {
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, scalev, 1, "expected either 2 or 4 arguments (vec, scale, or vec, scale, start, end)");
    return 0;
  }
  tk_vec_pfx(scalev)(m0, m1, start, end);
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
  return 1;
}

static inline int tk_vec_pfx(pow_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_base e;
  uint64_t start, end;
  if (t == 2) {
    e = tk_lua_checkdouble(L, 2, "e");
    start = 0;
    end = m0->n;
  } else if (t == 3) {
    e = tk_lua_checkdouble(L, 2, "e");
    start = tk_lua_checkunsigned(L, 3, "start");
    end = tk_lua_checkunsigned(L, 4, "end");
  } else {
    tk_vec_err(L, pow, 1, "expected either 2 or 4 arguments (vec, e or vec, e, start, end)");
    return 0;
  }
  tk_vec_pfx(pow)(m0, e, start, end);
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
  return 1;
}

static inline int tk_vec_pfx(fill_indices_lua) (lua_State *L)
{
  lua_settop(L, 1);
  tk_vec_pfx(t) *v = tk_vec_pfx(peek)(L, 1, "vector");
  tk_vec_pfx(fill_indices)(v);
  return 1;
}

static inline int tk_vec_pfx(raw_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  const char *fmt = tk_lua_optstring(L, 2, "fmt", NULL);
  tk_vec_pfx(raw)(L, m0, fmt);
  return 1;
}

#endif

static luaL_Reg tk_vec_pfx(lua_mt_fns)[] =
{
  { "copy", tk_vec_pfx(copy_lua) },
  { "destroy", tk_vec_pfx(destroy_lua) },
  { "size", tk_vec_pfx(size_lua) },
  { "capacity", tk_vec_pfx(capacity_lua) },
  { "resize", tk_vec_pfx(resize_lua) },
  { "setn", tk_vec_pfx(setn_lua) },
  { "ensure", tk_vec_pfx(ensure_lua) },
  { "shrink", tk_vec_pfx(shrink_lua) },
  { "clear", tk_vec_pfx(clear_lua) },
  { "zero", tk_vec_pfx(zero_lua) },
  { "transpose", tk_vec_pfx(transpose_lua) },

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
  { "xasc", tk_vec_pfx(xasc_lua) },
  { "xdesc", tk_vec_pfx(xdesc_lua) },
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
  { "multiply", tk_vec_pfx(multiply_lua) },
  { "dot", tk_vec_pfx(dot_lua) },
  { "magnitude", tk_vec_pfx(magnitude_lua) },
  { "sum", tk_vec_pfx(sum_lua) },
  { "csums", tk_vec_pfx(csums_lua) },
  { "rsums", tk_vec_pfx(rsums_lua) },
  { "max", tk_vec_pfx(max_lua) },
  { "cmaxs", tk_vec_pfx(rmaxs_lua) },
  { "rmaxs", tk_vec_pfx(cmaxs_lua) },
  { "min", tk_vec_pfx(min_lua) },
  { "cmins", tk_vec_pfx(rmins_lua) },
  { "rmins", tk_vec_pfx(cmins_lua) },
  { "fill", tk_vec_pfx(fill_lua) },
  { "fill_indices", tk_vec_pfx(fill_indices_lua) },
  { "raw", tk_vec_pfx(raw_lua) },

#endif

  { NULL, NULL }
};

static inline void tk_vec_pfx(suppress_unused_lua_mt_fns) (void)
  { (void) tk_vec_pfx(lua_mt_fns); }

static inline tk_vec_pfx(t) *tk_vec_pfx(create) (lua_State *L, size_t n, tk_vec_base *data, tk_vec_pfx(t) *v)
{
  tk_vec_pfx(t) v0;
  if (v == NULL) {
    v = L == NULL
      ? malloc(sizeof(tk_vec_pfx(t)))
      : tk_lua_newuserdata(L, tk_vec_pfx(t), tk_vec_mt, tk_vec_pfx(lua_mt_fns), tk_vec_pfx(gc_lua));
    if (v == NULL)
      return v;
    v0 = *v;
    bool lua_managed = L != NULL;
    kv_init(v0, lua_managed);
    kv_resize(tk_vec_base, v0, n);
    v0.n = n;
    *v = v0;
    return v;
  } else if (v->a && v->a != data) {
    free(v->a);
  }
  v->n = n;
  v->m = n;
  v->a = data;
  return v;
}

static inline tk_vec_pfx(t) *tk_vec_pfx(register) (lua_State *L, tk_vec_pfx(t) *v)
{
  tk_vec_pfx(t) *x = tk_lua_newuserdata(L, tk_vec_pfx(t), tk_vec_mt, tk_vec_pfx(lua_mt_fns), tk_vec_pfx(gc_lua));
  *x = *v;
  free(v);
  return x;
}

static luaL_Reg tk_vec_pfx(lua_fns)[] =
{
  { "create", tk_vec_pfx(create_lua) },
  { "load", tk_vec_pfx(load_lua) },
  { NULL, NULL }
};

static inline void tk_vec_pfx(suppress_unused_lua_fns) (void)
  { (void) tk_vec_pfx(lua_fns); }

#include <santoku/vec/undef.h>
