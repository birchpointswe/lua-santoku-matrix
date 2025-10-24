#include <omp.h>

#define tk_vec_pfx(name) tk_pp_strcat(tk_vec_name, name)

static inline void tk_vec_pfx(copy_indexed) (
  tk_vec_pfx(t) *m0,
  tk_vec_pfx(t) *m1,
  tk_ivec_t *indices
) {
  tk_vec_pfx(ensure)(m0, indices->n);
  #pragma omp parallel for
  for (uint64_t i = 0; i < indices->n; i++) {
    int64_t idx = indices->a[i];
    if (idx >= 0 && idx < (int64_t)m1->n) {
      m0->a[i] = m1->a[idx];
    }
  }
  if (m0->n < indices->n)
    m0->n = indices->n;
}

static inline int tk_vec_pfx(copy_indexed_lua) (lua_State *L)
{
  int t = lua_gettop(L);
  if (t >= 3 && lua_type(L, 3) == LUA_TUSERDATA) {
    tk_ivec_t *indices = tk_ivec_peekopt(L, 3);
    if (indices != NULL) {
      tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "dest");
      tk_vec_pfx(t) *m1 = tk_vec_pfx(peek)(L, 2, "source");
      tk_vec_pfx(copy_indexed)(m0, m1, indices);
      return 0;
    }
  }
  return tk_vec_pfx(copy_lua)(L);
}

#ifndef tk_vec_limited

static inline tk_ivec_t *tk_vec_pfx(rasc) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, m0->n, NULL, NULL);
  uint64_t rows = m0->n / cols;
  #pragma omp parallel
  {
    tk_rvec_t *ranks = tk_rvec_create(NULL, cols, NULL, NULL);
    #pragma omp for
    for (uint64_t r = 0; r < rows; r ++) {
      for (size_t c = 0; c < cols; c ++) {
        ranks->a[c].i = (int64_t) c;
        ranks->a[c].d = m0->a[r * cols + c];
      }
      tk_rvec_asc(ranks, 0, ranks->n);
      for (size_t c = 0; c < cols; c ++)
        out->a[r * cols + c] = ranks->a[c].i;
    }
    tk_rvec_destroy(ranks);
  }
  return out;
}

static inline tk_ivec_t *tk_vec_pfx(rdesc) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, m0->n, NULL, NULL);
  uint64_t rows = m0->n / cols;
  #pragma omp parallel
  {
    tk_rvec_t *ranks = tk_rvec_create(NULL, cols, NULL, NULL);
    #pragma omp for
    for (uint64_t r = 0; r < rows; r ++) {
      for (size_t c = 0; c < cols; c ++) {
        ranks->a[c].i = (int64_t) c;
        ranks->a[c].d = m0->a[r * cols + c];
      }
      tk_rvec_desc(ranks, 0, ranks->n);
      for (size_t c = 0; c < cols; c ++)
        out->a[r * cols + c] = ranks->a[c].i;
    }
    tk_rvec_destroy(ranks);
  }
  return out;
}

static inline tk_ivec_t *tk_vec_pfx(casc) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  uint64_t rows = m0->n / cols;
  tk_ivec_t *out = tk_ivec_create(L, m0->n, NULL, NULL);
  #pragma omp parallel
  {
    tk_rvec_t *ranks = tk_rvec_create(0, rows, NULL, NULL);
    #pragma omp for
    for (size_t c = 0; c < cols; c ++) {
      for (uint64_t r = 0; r < rows ; r ++) {
        ranks->a[r].i = (int64_t) r;
        ranks->a[r].d = m0->a[c * rows + r];
      }
      tk_rvec_asc(ranks, 0, ranks->n);
      for (size_t r = 0; r < rows; r ++)
        out->a[c * rows + r] = ranks->a[r].i;
    }
    tk_rvec_destroy(ranks);
  }
  return out;
}

static inline tk_ivec_t *tk_vec_pfx(cdesc) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  uint64_t rows = m0->n / cols;
  tk_ivec_t *out = tk_ivec_create(L, m0->n, NULL, NULL);
  #pragma omp parallel
  {
    tk_rvec_t *ranks = tk_rvec_create(0, rows, NULL, NULL);
    #pragma omp for
    for (size_t c = 0; c < cols; c ++) {
      for (uint64_t r = 0; r < rows ; r ++) {
        ranks->a[r].i = (int64_t) r;
        ranks->a[r].d = m0->a[c * rows + r];
      }
      tk_rvec_desc(ranks, 0, ranks->n);
      for (size_t r = 0; r < rows; r ++)
        out->a[c * rows + r] = ranks->a[r].i;
    }
    tk_rvec_destroy(ranks);
  }
  return out;
}

static inline tk_dvec_t *tk_vec_pfx(cmagnitudes) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_dvec_t *out = tk_dvec_create(L, cols, NULL, NULL);
  uint64_t rows = m0->n / cols;
  #pragma omp parallel for
  for (uint64_t c = 0; c < cols; c ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t r = 0; r < rows; r ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val * val;
    }
    out->a[c] = sqrt(sum);
  }
  return out;
}

static inline tk_dvec_t *tk_vec_pfx(rmagnitudes) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  uint64_t rows = m0->n / cols;
  tk_dvec_t *out = tk_dvec_create(L, rows, NULL, NULL);
  #pragma omp parallel for
  for (uint64_t r = 0; r < rows; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val * val;
    }
    out->a[r] = sqrt(sum);
  }
  return out;
}

static inline tk_ivec_t *tk_vec_pfx(cmaxargs) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, cols, NULL, NULL);
  uint64_t rows = m0->n / cols;
  #pragma omp parallel for
  for (uint64_t c = 0; c < cols; c ++) {
    uint64_t maxr = 0;
    tk_vec_base maxv = m0->a[0 * cols + c];
    for (size_t r = 1; r < rows; r ++) {
      if (m0->a[r * cols + c] > maxv) {
        maxr = r;
        maxv = m0->a[r * cols + c];
      }
    }
    out->a[c] = (int64_t) maxr;
  }
  return out;
}

static inline tk_ivec_t *tk_vec_pfx(rmaxargs) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  uint64_t rows = m0->n / cols;
  tk_ivec_t *out = tk_ivec_create(L, rows, NULL, NULL);
  #pragma omp parallel for
  for (uint64_t r = 0; r < rows; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val;
    }
    out->a[r] = sum;
    uint64_t maxc = 0;
    tk_vec_base maxv = m0->a[r * cols + 0];
    for (uint64_t c = 1; c < cols; c ++) {
      if (m0->a[r * cols + c] > maxv) {
        maxc = c;
        maxv = m0->a[r * cols + c];
      }
    }
    out->a[r] = (int64_t) maxc;
  }
  return out;
}

static inline tk_ivec_t *tk_vec_pfx(cminargs) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, cols, NULL, NULL);
  uint64_t rows = m0->n / cols;
  #pragma omp parallel for
  for (uint64_t c = 0; c < cols; c ++) {
    uint64_t minr = 0;
    tk_vec_base minv = m0->a[0 * cols + c];
    for (size_t r = 1; r < rows; r ++) {
      if (m0->a[r * cols + c] < minv) {
        minr = r;
        minv = m0->a[r * cols + c];
      }
    }
    out->a[c] = (int64_t) minr;
  }
  return out;
}

static inline tk_ivec_t *tk_vec_pfx(rminargs) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  uint64_t rows = m0->n / cols;
  tk_ivec_t *out = tk_ivec_create(L, rows, NULL, NULL);
  #pragma omp parallel for
  for (uint64_t r = 0; r < rows; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      sum += val;
    }
    out->a[r] = sum;
    uint64_t minc = 0;
    tk_vec_base minv = m0->a[r * cols + 0];
    for (uint64_t c = 1; c < cols; c ++) {
      if (m0->a[r * cols + c] < minv) {
        minc = c;
        minv = m0->a[r * cols + c];
      }
    }
    out->a[r] = (int64_t) minc;
  }
  return out;
}

static inline int tk_vec_pfx(rmagnitudes_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rmagnitudes)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(cmagnitudes_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cmagnitudes)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(rminargs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rminargs)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(cminargs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cminargs)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(rmaxargs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rmaxargs)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(cmaxargs_lua) (lua_State *L) {
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cmaxargs)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(rasc_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rasc)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(rdesc_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(rdesc)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(casc_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(casc)(L, m0, cols);
  return 1;
}

static inline int tk_vec_pfx(cdesc_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  uint64_t cols = tk_lua_checkunsigned(L, 2, "cols");
  tk_vec_pfx(cdesc)(L, m0, cols);
  return 1;
}

#endif

static inline void tk_vec_pfx(persist) (lua_State *L, tk_vec_pfx(t) *v, FILE *fh)
{
  tk_lua_fwrite(L, (char *) &v->n, sizeof(size_t), 1, fh);
  tk_lua_fwrite(L, (char *) v->a, sizeof(tk_vec_base) * v->n, 1, fh);
}

static inline int tk_vec_pfx(persist_lua) (lua_State *L)
{
  lua_settop(L, 2);
  tk_vec_pfx(t) *m0 = tk_vec_pfx(peek)(L, 1, "vector");
  bool tostr = lua_type(L, 2) == LUA_TNIL;
  FILE *fh;
  if (tostr)
    fh = tk_lua_tmpfile(L);
  else
    fh = tk_lua_fopen(L, luaL_checkstring(L, 2), "w");
  tk_vec_pfx(persist)(L, m0, fh);
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

static luaL_Reg tk_vec_pfx(lua_mt_ext_fns)[] =
{
  { "copy", tk_vec_pfx(copy_indexed_lua) },
  { "persist", tk_vec_pfx(persist_lua) },
#ifndef tk_vec_limited
  { "cmagnitudes", tk_vec_pfx(cmagnitudes_lua) },
  { "rmagnitudes", tk_vec_pfx(rmagnitudes_lua) },
  { "cmaxargs", tk_vec_pfx(rmaxargs_lua) },
  { "rmaxargs", tk_vec_pfx(cmaxargs_lua) },
  { "cminargs", tk_vec_pfx(rminargs_lua) },
  { "rminargs", tk_vec_pfx(cminargs_lua) },
  { "rasc", tk_vec_pfx(rasc_lua) },
  { "rdesc", tk_vec_pfx(rdesc_lua) },
  { "casc", tk_vec_pfx(casc_lua) },
  { "cdesc", tk_vec_pfx(cdesc_lua) },
#endif
  { NULL, NULL }
};

static inline void tk_vec_pfx(suppress_unused_lua_mt_ext_fns) (void)
  { (void) tk_vec_pfx(lua_mt_ext_fns); }

#include <santoku/vec/undef.h>
