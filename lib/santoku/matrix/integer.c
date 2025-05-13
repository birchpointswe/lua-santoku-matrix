#define tk_base_t int64_t
#define tk_sort(...) ks_introsort(int64_t, __VA_ARGS__)
#define tk_matrix_pushbase(...) lua_pushinteger(__VA_ARGS__)
#define tk_matrix_peekbase(...) luaL_checkinteger(__VA_ARGS__)
#define TK_MT "santoku_matrix_integer"
#define TK_OPEN luaopen_santoku_matrix_integer
#include "gen.h"

typedef struct { double chi2; uint64_t v; } tk_matrix_chi2_pair;
#define tk_matrix_chi2_pair_gt(a, b) ((a).chi2 > (b).chi2)
KSORT_INIT(chi2_desc, tk_matrix_chi2_pair, tk_matrix_chi2_pair_gt)

#include <errno.h>

static inline int tk_error (
  lua_State *L,
  const char *label,
  int err
) {
  lua_pushstring(L, label);
  lua_pushstring(L, strerror(err));
  tk_lua_callmod(L, 2, 0, "santoku.error", "error");
  return 1;
}

static inline void *tk_malloc (
  lua_State *L,
  size_t s
) {
  void *p = malloc(s);
  if (!p) {
    tk_error(L, "malloc failed", ENOMEM);
    return NULL;
  } else {
    return p;
  }
}

static inline tk_matrix_t *_tk_matrix_create (lua_State *, size_t, size_t, tk_base_t *, tk_matrix_t *);


static inline int tk_matrix_radd (lua_State *L)
{
  lua_settop(L, 5);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  luaL_checktype(L, 3, LUA_TNUMBER);
  luaL_checktype(L, 4, LUA_TNUMBER);
  lua_Integer rowstart = lua_tointeger(L, 2);
  lua_Integer rowend = lua_tointeger(L, 3);
  lua_Number add = lua_tonumber(L, 4);
  if (rowstart > rowend)
    luaL_error(L, "Error in radd: start row is greater than end row");
  size_t idxstart = tk_matrix_index(L, m0, rowstart, 1);
  size_t idxend = tk_matrix_index(L, m0, rowend, m0->columns);
  for (unsigned int i = idxstart; i <= idxend; i++)
    m0->data[i] += add;
  return 0;
}

static inline int tk_matrix_rmult (lua_State *L)
{
  lua_settop(L, 4);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  luaL_checktype(L, 3, LUA_TNUMBER);
  luaL_checktype(L, 4, LUA_TNUMBER);
  lua_Integer rowstart = lua_tointeger(L, 2);
  lua_Integer rowend = lua_tointeger(L, 3);
  lua_Number scal = lua_tonumber(L, 4);
  if (rowstart > rowend)
    luaL_error(L, "Error in rmult: start row is greater than end row");
  size_t idxstart = tk_matrix_index(L, m0, rowstart, 1);
  size_t idxend = tk_matrix_index(L, m0, rowend, m0->columns);
  tk_base_t *x = m0->data;
  for (unsigned int i = idxstart; i <= idxend; i++)
    x[i] *= scal;
  return 0;
}

static inline int tk_matrix_ramax (lua_State *L)
{
  lua_settop(L, 2);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  lua_Integer row = lua_tointeger(L, 2);
  size_t idx = tk_matrix_index(L, m0, row, 1);
  size_t idxval = idx;
  tk_base_t maxval = abs(m0->data[idx]);
  for (size_t i = idx + 1; i < idx + m0->columns; i++) {
    tk_base_t absval = abs(m0->data[i]);
    if (absval > maxval) {
      maxval = absval;
      idxval = i;
    }
  }
  lua_pushnumber(L, m0->data[idx + idxval]);
  lua_pushinteger(L, idxval);
  return 2;
}

static inline int tk_matrix_sums (lua_State *L)
{
  lua_settop(L, 3);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  tk_matrix_t *m1 = tk_matrix_peek(L, 2);
  luaL_checktype(L, 3, LUA_TNUMBER);
  if (m0->columns != m1->columns)
    luaL_error(L, "Error in sums: destination matrix columns don't match source matrix columns");
  lua_Integer rowdest = lua_tointeger(L, 3);
  size_t idxdest = tk_matrix_index(L, m1, rowdest, 1);
  size_t idxsrc = tk_matrix_index(L, m0, 1, 1);
  memcpy(&m1->data[idxdest], &m0->data[idxsrc], sizeof(tk_base_t) * m1->columns);
  for (size_t i = 2; i <= m0->rows; i ++) {
    idxsrc = tk_matrix_index(L, m0, i, 1);
    for (size_t i = 0; i < m0->columns; i++)
      m1->data[idxdest + i] += m0->data[idxsrc + i];
  }
  return 0;
}

static inline int tk_matrix_mmult (lua_State *L)
{
  lua_settop(L, 5);
  tk_matrix_t *a = tk_matrix_peek(L, 1);
  tk_matrix_t *b = tk_matrix_peek(L, 2);
  tk_matrix_t *c = tk_matrix_peek(L, 3);
  bool transpose_a = lua_toboolean(L, 4);
  bool transpose_b = lua_toboolean(L, 5);
  if (!transpose_a && !transpose_b) {
    if (a->columns != b->rows)
      luaL_error(L, "Error in mmult: columns of A don't match rows of B");
    if (a->rows != c->rows)
      luaL_error(L, "Error in mmult: rows of C don't match rows of A");
    if (b->columns != c->columns)
      luaL_error(L, "Error in mmult: columns of C don't match columns of B");
  } else if (transpose_a && !transpose_b) {
    if (a->rows != b->rows)
      luaL_error(L, "Error in mmult: rows of A don't match rows of B");
    if (a->columns != c->rows)
      luaL_error(L, "Error in mmult: rows of C don't match columns of A");
    if (b->columns != c->columns)
      luaL_error(L, "Error in mmult: columns of C don't match columns of B");
  } else if (!transpose_a && transpose_b) {
    if (a->columns != b->columns)
      luaL_error(L, "Error in mmult: columns of A don't match columns of B");
    if (a->rows != c->rows)
      luaL_error(L, "Error in mmult: rows of C don't match rows of A");
    if (b->rows != c->columns)
      luaL_error(L, "Error in mmult: columns of C don't match rows of B");
  } else if (transpose_a && transpose_b) {
    if (a->rows != b->columns)
      luaL_error(L, "Error in mmult: rows of A don't match columns of B");
    if (a->columns != c->columns)
      luaL_error(L, "Error in mmult: columns of C don't match columns of A");
    if (b->rows != c->rows)
      luaL_error(L, "Error in mmult: rows of C don't match rows of B");
  }
  for (size_t i = 0; i < c->rows; i++) {
    for (size_t j = 0; j < c->columns; j++) {
      tk_base_t sum = 0.0;
      for (size_t k = 0; k < a->columns; k++) {
        tk_base_t a_val = transpose_a ? a->data[k * a->columns + i] : a->data[i * a->columns + k];
        tk_base_t b_val = transpose_b ? b->data[j * b->columns + k] : b->data[k * b->columns + j];
        sum += a_val * b_val;
      }
      c->data[i * c->columns + j] = sum;
    }
  }
  return 0;
}

static inline int tk_matrix_magnitude (lua_State *L)
{
  lua_settop(L, 2);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  lua_Integer row = lua_tointeger(L, 2);
  size_t idx = tk_matrix_index(L, m0, row, 1);
  tk_base_t sum = 0.0;
  for (size_t i = 0; i < m0->columns; i++) {
    tk_base_t val = m0->data[idx + i];
    sum += val * val;
  }
  lua_pushnumber(L, sqrt(sum));
  return 1;
}

static inline int tk_matrix_flip_interleave (lua_State *L)
{
  lua_settop(L, 3);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  uint64_t n_samples = tk_lua_checkunsigned(L, 2);
  uint64_t n_features = tk_lua_checkunsigned(L, 3);
  lua_settop(L, 1);
  bool sorted = lua_toboolean(L, 4);
  if (!sorted) {
    lua_pushboolean(L, true);
    tk_matrix_sort(L);
    lua_settop(L, 1);
  }
  size_t values = m0->values;
  tk_base_t *data = m0->data;
  size_t total = n_samples * n_features;
  data = realloc(data, sizeof(tk_base_t) * total);
  if (data == NULL)
    luaL_error(L, "Error in realloc during flip_interleave");
  size_t write = values;
  size_t last = 0;
  tk_base_t s, k;
  for (size_t i_present = 0; i_present < values; i_present ++) {
    tk_base_t x = data[i_present];
    s = x / n_features;
    k = x % n_features;
    data[i_present] = (s * 2 * n_features) + k;
    for (size_t y = last; y < x; y ++) {
      tk_base_t s = y / n_features;
      tk_base_t k = y % n_features;
      data[write ++] = (s * 2 * n_features) + n_features + k;
    }
    last = x + 1;
  }
  for (size_t y = last; y < total; y ++) {
    tk_base_t s = y / n_features;
    tk_base_t k = y % n_features;
    data[write ++] = (s * 2 * n_features) + n_features + k;
  }
  m0->rows = 1;
  m0->columns = write;
  m0->values = write;
  m0->data = data;
  tk_matrix_sort(L);
  return 0;
}

static int tk_matrix_top_chi2 (lua_State *L)
{
  lua_settop(L, 6);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  int64_t *set_bits = m0->data;
  size_t n_set_bits = m0->values;
  const char *codes = luaL_checkstring(L, 2);
  uint64_t n_samples = tk_lua_checkunsigned(L, 3);
  uint64_t n_visible = tk_lua_checkunsigned(L, 4);
  uint64_t n_labels = tk_lua_checkunsigned(L, 5);
  uint64_t top_k = tk_lua_checkunsigned(L, 6);
  tk_matrix_chi2_pair *chis = malloc(n_visible * sizeof(tk_matrix_chi2_pair));
  for (unsigned int v = 0; v < n_visible; v ++)
    chis[v] = (tk_matrix_chi2_pair){ .chi2 = 0, .v = v };
  unsigned int *global_counts = tk_malloc(L, n_labels * sizeof(unsigned int));
  for (unsigned int y = 0; y < n_labels; y ++)
    global_counts[y] = 0;
  for (unsigned int s = 0; s < n_samples; s ++)
    for (unsigned int l = 0; l < n_labels; l ++)
      if (codes[s * (n_labels / CHAR_BIT) + (l / CHAR_BIT)] & (1 << (l % CHAR_BIT)))
        global_counts[l] ++;
  uint64_t *active_counts = tk_malloc(L, n_visible * n_labels * sizeof(uint64_t));
  for (uint64_t i = 0; i  < n_set_bits; i ++) {
    int64_t val = set_bits[i];
    unsigned int sample = val / n_visible;
    unsigned int feature = val % n_visible;
    for (uint64_t j = 0; j < n_labels; j ++)
      if (codes[sample * (n_labels / CHAR_BIT) + (j / CHAR_BIT)] & (1 << (j % CHAR_BIT)))
        active_counts[feature * n_labels + j] ++;
  }
  for (unsigned int v = 0; v < n_visible; v ++) {
    unsigned int active_total = 0;
    for (unsigned int y = 0; y < n_labels; y ++)
      active_total += active_counts[v * n_labels + y];
    unsigned int inactive_total = n_samples - active_total;
    if (active_total == 0)
      continue;
    double sum_chi2 = 0.0;
    for (unsigned int y = 0; y < n_labels; y ++) {
      unsigned int observed_active = active_counts[v * n_labels + y];
      unsigned int observed_inactive = global_counts[y] - observed_active;
      double expected_active = ((double)global_counts[y] * active_total) / n_samples;
      double expected_inactive = ((double)global_counts[y] * inactive_total) / n_samples;
      double chi2_y = 0.0;
      if (expected_active > 0.0)
        chi2_y += ((observed_active - expected_active) * (observed_active - expected_active)) / expected_active;
      if (expected_inactive > 0.0)
        chi2_y += ((observed_inactive - expected_inactive) * (observed_inactive - expected_inactive)) / expected_inactive;
      sum_chi2 += chi2_y;
    }
    chis[v].chi2 = sum_chi2 / n_labels;
  }
  free(global_counts);
  free(active_counts);
  ks_introsort(chi2_desc, n_visible, chis);
  unsigned int m = (unsigned int) top_k;
  if (m > n_visible)
    m = n_visible;
  tk_matrix_t *m1 = _tk_matrix_create(L, 1, m, NULL, NULL);
  for (unsigned int i = 0; i < m; i ++)
    m1->data[i] = (int64_t) chis[i].v;
  free(chis);
  return 1;
}

static int tk_matrix_filter (lua_State *L)
{
  lua_settop(L, 3);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  int64_t *set_bits = m0->data;
  size_t n_set_bits = m0->values;
  tk_matrix_t *m1 = tk_matrix_peek(L, 2);
  int64_t *top_v = m1->data;
  size_t n_top_v = m1->values;
  uint64_t n_visible = tk_lua_checkunsigned(L, 3);
  int64_t *vmap = tk_malloc(L, n_visible * sizeof(int64_t));
  for (unsigned int i = 0; i < n_visible; i ++)
    vmap[i] = -1;
  for (unsigned int i = 0; i < n_top_v; i ++)
    vmap[top_v[i]] = i;
  size_t write = 0;
  for (size_t i = 0; i < n_set_bits; i ++) {
    int64_t val = set_bits[i];
    uint64_t sample = val / n_visible;
    uint64_t feature = val % n_visible;
    int64_t new_feature = vmap[feature];
    if (new_feature == -1)
      continue;
    set_bits[write ++] = sample * n_top_v + new_feature;
  }
  m0->rows = 1;
  m0->columns = write;
  m0->values = write;
  free(vmap);
  return 0;
}

static int tk_matrix_raw_bitmap (lua_State *L)
{
  lua_settop(L, 3);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  int64_t *set_bits = m0->data;
  size_t n_set_bits = m0->values;
  uint64_t n_samples = tk_lua_checkunsigned(L, 2);
  uint64_t n_features = tk_lua_checkunsigned(L, 3);
  size_t len = (n_samples * n_features + CHAR_BIT - 1) / CHAR_BIT;
  char *out = tk_malloc(L, len);
  memset(out, 0, len);
  for (uint64_t i = 0; i < n_set_bits; i ++) {
    int64_t v = set_bits[i];
    if ((uint64_t) v >= n_samples * n_features)
      continue;
    out[v / CHAR_BIT] |= (1 << (v % CHAR_BIT));
  }
  lua_pushlstring(L, out, len);
  free(out);
  return 1;
}

static luaL_Reg tk_matrix_extra_fns[] =
{
  { "flip_interleave", tk_matrix_flip_interleave },
  { "top_chi2", tk_matrix_top_chi2 },
  { "filter", tk_matrix_filter },
  { "raw_bitmap", tk_matrix_raw_bitmap },
  { NULL, NULL }
};
