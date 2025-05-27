#include <santoku/threads.h>
#include <santoku/matrix/integer.conf.h>
#include <santoku/matrix/gen.h>
#include <santoku/matrix/integer.h>

static inline tk_matrix_t *_tk_matrix_create (lua_State *, size_t, size_t, tk_base_t *, tk_matrix_t *);

typedef enum {
  TK_MTX_CHI2,
  TK_MTX_MI,
} tk_mtx_stage_t;

typedef struct {
  uint64_t n_visible, n_hidden, n_samples;
  double *scores;
  uint64_t *counts, *feat_counts;
  int64_t *active_counts, *global_counts;
  int64_t *set_bits;
  uint64_t n_set_bits;
  int64_t *labels;
  char *codes;
} tk_mtx_ctx_t;

typedef struct {
  tk_mtx_ctx_t *state;
  uint64_t hfirst, hlast;
} tk_mtx_thread_t;

static void tk_eval_worker (void *dp, int sig)
{
  tk_mtx_thread_t *data = (tk_mtx_thread_t *) dp;
  tk_mtx_ctx_t *state = data->state;
  uint64_t n_visible = state->n_visible;
  uint64_t n_hidden = state->n_hidden;
  uint64_t n_samples = state->n_samples;
  double *scores = state->scores;
  uint64_t *counts = state->counts;
  int64_t *active_counts = state->active_counts;
  int64_t *global_counts = state->global_counts;
  uint64_t *feat_counts = state->feat_counts;
  switch ((tk_mtx_stage_t) sig) {

    case TK_MTX_CHI2:
      for (uint64_t b = data->hfirst; b <= data->hlast; b  ++) {
        double *scores_b = scores + b * n_visible;
        for (uint64_t f = 0; f < n_visible; f ++) {
          uint64_t A = active_counts[f * n_hidden + b];
          uint64_t G = global_counts[b];
          uint64_t C = feat_counts[f];
          if (C == 0 || G == 0 || C == n_samples || G == n_samples) {
            scores_b[f] = 0.0;
            continue;
          }
          uint64_t B = G - A;
          uint64_t C_ = C - A;
          uint64_t D = n_samples - C - B;
          double n = (double) n_samples;
          double E_A = ((double) C * (double) G) / n;
          double E_B = ((double)(n - C) * (double) G) / n;
          double E_C = ((double) C * (double)(n - G)) / n;
          double E_D = ((double)(n - C) * (double)(n - G)) / n;
          double chi2 = 0.0;
          if (E_A > 0)
            chi2 += ((A - E_A)*(A - E_A)) / E_A;
          if (E_B > 0)
            chi2 += ((B - E_B)*(B - E_B)) / E_B;
          if (E_C > 0)
            chi2 += ((C_ - E_C)*(C_ - E_C)) / E_C;
          if (E_D > 0)
            chi2 += ((D - E_D)*(D - E_D)) / E_D;
          scores_b[f] = chi2;
        }
      }
      break;

    case TK_MTX_MI:
      for (int64_t j = data->hfirst; j <= data->hlast; j ++) {
        double *scores_h = scores + j * n_visible;
        for (int64_t i = 0; i < n_visible; i ++) {
          uint64_t *c = counts + i * n_hidden * 4 + j * 4;
          for (int k = 0; k < 4; k ++)
            c[k] += 1;
          double total = c[0] + c[1] + c[2] + c[3];
          double mi = 0.0;
          if (total == 0.0)
            continue;
          for (unsigned int o = 0; o < 4; o ++) {
            if (c[o] == 0)
              continue;
            double p_fb = c[o] / total;
            unsigned int f = o >> 1;
            unsigned int b = o & 1;
            double pf = (c[2] + c[3]) / total;
            if (f == 0) pf = 1.0 - pf;
            double pb = (c[1] + c[3]) / total;
            if (b == 0) pb = 1.0 - pb;
            double d = pf * pb;
            if (d > 0) mi += p_fb * log2(p_fb / d);
          }
          scores_h[i] = mi;
        }
      }
      break;
  }
}

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
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "features");
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

static inline double *_tk_matrix_chi2_scores (
  lua_State *L,
  int64_t *set_bits,
  size_t n_set_bits,
  char *codes,
  int64_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  int64_t **activep,
  int64_t **globalp,
  unsigned int n_threads
) {
  tk_mtx_ctx_t ctx;
  tk_mtx_thread_t threads[n_threads];
  tk_threadpool_t *pool = tk_threads_create(L, n_threads, tk_eval_worker);
  ctx.set_bits = set_bits;
  ctx.n_set_bits = n_set_bits;
  ctx.codes = codes;
  ctx.labels = labels;
  ctx.n_visible = n_visible;
  ctx.n_hidden = n_hidden;
  ctx.n_samples = n_samples;
  ctx.active_counts = tk_malloc(L, ctx.n_visible * ctx.n_hidden * sizeof(int64_t));
  ctx.global_counts = tk_malloc(L, ctx.n_hidden * sizeof(int64_t));
  ctx.feat_counts  = tk_malloc(L, ctx.n_visible * sizeof(uint64_t));
  memset(ctx.active_counts, 0, ctx.n_visible * ctx.n_hidden * sizeof(int64_t));
  memset(ctx.global_counts, 0, ctx.n_hidden * sizeof(int64_t));
  memset(ctx.feat_counts,  0, ctx.n_visible * sizeof(uint64_t));
  for (unsigned int i = 0; i < n_threads; i ++) {
    tk_mtx_thread_t *data = threads + i;
    pool->threads[i].data = data;
    data->state = &ctx;
    tk_thread_range(i, n_threads, ctx.n_hidden, &data->hfirst, &data->hlast);
  }



  for (uint64_t s = 0; s < ctx.n_samples; s++) {
    const unsigned char *sample_bitmap
      = (const unsigned char *) (ctx.codes + s * (ctx.n_hidden / CHAR_BIT));
    for (uint64_t chunk = 0; chunk < (ctx.n_hidden / CHAR_BIT); chunk ++) {
      unsigned char byte = sample_bitmap[chunk];
      while (byte) {
        int bit = __builtin_ctz(byte);
        uint64_t b = chunk * CHAR_BIT + bit;
        if (b < ctx.n_hidden)
          ctx.global_counts[b] ++;
        byte &= byte - 1;
      }
    }
  }



  if (ctx.codes != NULL) {
    for (uint64_t i = 0; i < ctx.n_set_bits; i ++) {
      uint64_t s = ctx.set_bits[i] / ctx.n_visible;
      uint64_t f = ctx.set_bits[i] % ctx.n_visible;
      ctx.feat_counts[f] ++;
      const unsigned char *bitmap = (unsigned char *) (ctx.codes + s * (ctx.n_hidden / CHAR_BIT));
      for (uint64_t b = 0; b < ctx.n_hidden; b ++) {
        uint64_t chunk = b / CHAR_BIT, bit = b % CHAR_BIT;
        if (bitmap[chunk] & (1u << bit)) {
          ctx.active_counts[f * ctx.n_hidden + b] ++;
        }
      }
    }
  } else if (labels != NULL) {
    for (uint64_t i = 0; i < ctx.n_set_bits; i ++) {
      uint64_t s = ctx.set_bits[i] / ctx.n_visible;
      uint64_t f = ctx.set_bits[i] % ctx.n_visible;
      ctx.feat_counts[f] ++;
      int64_t label = ctx.labels[s];
      if ((size_t) label < ctx.n_hidden)
        ctx.active_counts[f * ctx.n_hidden + label] ++;
    }
    memset(ctx.global_counts, 0, ctx.n_hidden * sizeof(int64_t));
    for (uint64_t s = 0; s < ctx.n_samples; s ++) {
      int64_t label = ctx.labels[s];
      if ((size_t) label < ctx.n_hidden)
        ctx.global_counts[label] ++;
    }
  }



  ctx.scores = tk_malloc(L, ctx.n_hidden * ctx.n_visible * sizeof(double));
  memset(ctx.scores, 0, ctx.n_hidden * ctx.n_visible * sizeof(double));
  tk_threads_signal(pool, TK_MTX_CHI2);
  tk_threads_destroy(pool);

  if (activep != NULL)
    *activep = ctx.active_counts;
  else
    free(ctx.active_counts);

  if (globalp != NULL)
    *globalp = ctx.global_counts;
  else
    free(ctx.global_counts);

  free(ctx.feat_counts);

  return ctx.scores;
}

static inline double *_tk_matrix_mi_scores (
  lua_State *L,
  int64_t *set_bits,
  size_t n_set_bits,
  char *codes,
  int64_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  unsigned int n_threads
) {
  tk_mtx_ctx_t ctx;
  tk_mtx_thread_t threads[n_threads];
  tk_threadpool_t *pool = tk_threads_create(L, n_threads, tk_eval_worker);
  ctx.set_bits = set_bits;
  ctx.n_set_bits = n_set_bits;
  ctx.codes = codes;
  ctx.labels = labels;
  ctx.n_visible = n_visible;
  ctx.n_hidden = n_hidden;
  ctx.n_samples = n_samples;
  ctx.counts = tk_malloc(L, ctx.n_visible * ctx.n_hidden * 4 * sizeof(uint64_t));
  memset(ctx.counts, 0, ctx.n_visible * ctx.n_hidden * 4 * sizeof(uint64_t));
  for (unsigned int i = 0; i < n_threads; i ++) {
    tk_mtx_thread_t *data = threads + i;
    pool->threads[i].data = data;
    data->state = &ctx;
    tk_thread_range(i, n_threads, ctx.n_hidden, &data->hfirst, &data->hlast);
  }

  if (ctx.codes != NULL) {
    int64_t last_set_bit = -1;
    for (uint64_t i = 0; i < ctx.n_set_bits; i ++) {
      int64_t set_bit = ctx.set_bits[i];
      if (set_bit < 0) continue;
      for (uint64_t unset_bit = last_set_bit + 1; unset_bit < set_bit; unset_bit ++) {
        uint64_t sample = unset_bit / ctx.n_visible;
        uint64_t visible = unset_bit % ctx.n_visible;
        for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
          uint64_t chunk = j / CHAR_BIT;
          uint64_t bit = j % CHAR_BIT;
          bool hidden = (ctx.codes[sample * (ctx.n_hidden / CHAR_BIT) + chunk] & (1 << bit)) > 0;
          ctx.counts[visible * ctx.n_hidden * 4 + j * 4 + (hidden ? 1 : 0)] ++;
        }
      }
      uint64_t sample = set_bit / ctx.n_visible;
      uint64_t visible = set_bit % ctx.n_visible;
      for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
        uint64_t chunk = j / CHAR_BIT;
        uint64_t bit = j % CHAR_BIT;
        bool hidden = (ctx.codes[sample * (ctx.n_hidden / CHAR_BIT) + chunk] & (1 << bit)) > 0;
        ctx.counts[visible * ctx.n_hidden * 4 + j * 4 + 2 + (hidden ? 1 : 0)] ++;
      }
      last_set_bit = set_bit;
    }
    for (uint64_t unset_bit = last_set_bit + 1; unset_bit < n_samples * n_visible; unset_bit ++) {
      uint64_t sample = unset_bit / ctx.n_visible;
      uint64_t visible = unset_bit % ctx.n_visible;
      for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
        uint64_t chunk = j / CHAR_BIT;
        uint64_t bit = j % CHAR_BIT;
        bool hidden = (ctx.codes[sample * (ctx.n_hidden / CHAR_BIT) + chunk] & (1 << bit)) > 0;
        ctx.counts[visible * ctx.n_hidden * 4 + j * 4 + (hidden ? 1 : 0)] ++;
      }
    }
  } else if (ctx.labels != NULL) {
    int64_t last_set_bit = -1;
    for (uint64_t i = 0; i < n_set_bits; i ++) {
      int64_t set_bit = ctx.set_bits[i];
      if (set_bit < 0) continue;
      for (uint64_t unset_bit = last_set_bit + 1; unset_bit < set_bit; unset_bit ++) {
        uint64_t sample = unset_bit / ctx.n_visible;
        uint64_t visible = unset_bit % ctx.n_visible;
        int64_t label = ctx.labels[sample];
        for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
          bool hidden = ((size_t)j == (size_t)label);
          ctx.counts[visible * ctx.n_hidden * 4 + j * 4 + (hidden ? 1 : 0)] ++;
        }
      }
      uint64_t sample = set_bit / ctx.n_visible;
      uint64_t visible = set_bit % ctx.n_visible;
      int64_t label = ctx.labels[sample];
      for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
        bool hidden = ((size_t)j == (size_t)label);
        ctx.counts[visible * ctx.n_hidden * 4 + j * 4 + 2 + (hidden ? 1 : 0)] ++;
      }
      last_set_bit = set_bit;
    }
    for (uint64_t unset_bit = last_set_bit + 1; unset_bit < ctx.n_samples * ctx.n_visible; unset_bit ++) {
      uint64_t sample = unset_bit / ctx.n_visible;
      uint64_t visible = unset_bit % ctx.n_visible;
      int64_t label = ctx.labels[sample];
      for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
        bool hidden = ((size_t)j == (size_t)label);
        ctx.counts[visible * ctx.n_hidden * 4 + j * 4 + (hidden ? 1 : 0)] ++;
      }
    }
  }


  ctx.scores = tk_malloc(L, ctx.n_hidden * ctx.n_visible * sizeof(double));
  memset(ctx.scores, 0, ctx.n_hidden * ctx.n_visible * sizeof(double));
  tk_threads_signal(pool, TK_MTX_MI);
  tk_threads_destroy(pool);


  free(ctx.counts);

  return ctx.scores;
}

static inline int tk_matrix_score_chi2 (
  lua_State *L
) {
  lua_settop(L, 8);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  int64_t *set_bits = m0->data;
  size_t n_set_bits = m0->values;
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  bool ret_counts = lua_toboolean(L, 6);
  bool do_sums = lua_toboolean(L, 7);
  unsigned int n_threads = tk_threads_getn(L, 8, "threads", NULL);

  char *codes = NULL;
  int64_t *labels = NULL;

  if (lua_type(L, 2) == LUA_TSTRING) {
    codes = (char *) luaL_checkstring(L, 2);
  } else if (lua_type(L, 2) == LUA_TLIGHTUSERDATA) {
    codes = (char *) lua_touserdata(L, 2);
  } else {
    tk_matrix_t *m1 = tk_matrix_peek(L, 2);
    n_samples = m1->values < n_samples ? m1->values : n_samples;
    labels = m1->data;
  }

  int64_t *actives;
  int64_t *globals;
  double *scores = _tk_matrix_chi2_scores(L, set_bits, n_set_bits, codes, labels, n_samples, n_visible, n_hidden, &actives, &globals, n_threads);
  if (do_sums) {
    double *sums = tk_malloc(L, n_hidden * sizeof(double));
    memset(sums, 0, n_hidden * sizeof(double));
    for (uint64_t i = 0; i < n_hidden; i ++) {
      double sum = 0.0;
      for (uint64_t j = 0; j < n_visible; j ++)
        sum += scores[i * n_visible + j];
      sums[i] = sum;
    }
    free(scores);
    lua_pushlightuserdata(L, sums);
    lua_pushinteger(L, 1);
    lua_pushinteger(L, n_hidden);
    tk_lua_callmod(L, 3, 1, "santoku.matrix.number", "from_view");
  } else {
    lua_pushlightuserdata(L, scores);
    lua_pushinteger(L, n_hidden);
    lua_pushinteger(L, n_visible);
    tk_lua_callmod(L, 3, 1, "santoku.matrix.number", "from_view");
  }
  if (!ret_counts)
    return 1;
  lua_pushlightuserdata(L, actives);
  lua_pushinteger(L, n_hidden);
  lua_pushinteger(L, n_visible);
  tk_lua_callmod(L, 3, 1, "santoku.matrix.integer", "from_view");
  lua_pushlightuserdata(L, globals);
  lua_pushinteger(L, n_hidden);
  lua_pushinteger(L, n_visible);
  tk_lua_callmod(L, 3, 1, "santoku.matrix.integer", "from_view");
  return 3;
}

static inline int tk_matrix_score_mi (
  lua_State *L
) {
  lua_settop(L, 6);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  int64_t *set_bits = m0->data;
  size_t n_set_bits = m0->values;
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  unsigned int n_threads = tk_threads_getn(L, 6, "threads", NULL);

  char *codes = NULL;
  int64_t *labels = NULL;

  if (lua_type(L, 2) == LUA_TSTRING) {
    codes = (char *) luaL_checkstring(L, 2);
  } else if (lua_type(L, 2) == LUA_TLIGHTUSERDATA) {
    codes = (char *) lua_touserdata(L, 2);
  } else {
    tk_matrix_t *m1 = tk_matrix_peek(L, 2);
    n_samples = m1->values < n_samples ? m1->values : n_samples;
    labels = m1->data;
  }

  double *scores = _tk_matrix_mi_scores(L, set_bits, n_set_bits, codes, labels, n_samples, n_visible, n_hidden, n_threads);

  lua_pushlightuserdata(L, scores);
  lua_pushinteger(L, n_hidden);
  lua_pushinteger(L, n_visible);
  tk_lua_callmod(L, 3, 1, "santoku.matrix.number", "from_view");
  return 1;
}

static int tk_matrix_top_mi (lua_State *L)
{
  lua_settop(L, 7);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  int64_t *set_bits = m0->data;
  size_t n_set_bits = m0->values;
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  uint64_t top_k = tk_lua_checkunsigned(L, 6, "top_k");
  unsigned int n_threads = tk_threads_getn(L, 7, "threads", NULL);

  char *codes = NULL;
  int64_t *labels = NULL;

  if (lua_type(L, 2) == LUA_TSTRING) {
    codes = (char *) luaL_checkstring(L, 2);
  } else if (lua_type(L, 2) == LUA_TLIGHTUSERDATA) {
    codes = (char *) lua_touserdata(L, 2);
  } else {
    tk_matrix_t *m1 = tk_matrix_peek(L, 2);
    n_samples = m1->values < n_samples ? m1->values : n_samples;
    labels = m1->data;
  }


  ks_introsort(i64, n_set_bits, set_bits);


  double *scores = _tk_matrix_mi_scores(L, set_bits, n_set_bits, codes, labels, n_samples, n_visible, n_hidden, n_threads);
  tk_matrix_top_generic(L, scores, n_visible, n_hidden, top_k, 0);
  free(scores);
  return 1;
}


static int tk_matrix_top_chi2 (lua_State *L)
{
  lua_settop(L, 7);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  int64_t *set_bits = m0->data;
  size_t n_set_bits = m0->values;
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  uint64_t top_k = tk_lua_checkunsigned(L, 6, "top_k");
  unsigned int n_threads = tk_threads_getn(L, 7, "threads", NULL);

  char *codes = NULL;
  int64_t *labels = NULL;

  if (lua_type(L, 2) == LUA_TSTRING) {
    codes = (char *) luaL_checkstring(L, 2);
  } else if (lua_type(L, 2) == LUA_TLIGHTUSERDATA) {
    codes = (char *) lua_touserdata(L, 2);
  } else {
    tk_matrix_t *m1 = tk_matrix_peek(L, 2);
    n_samples = m1->values < n_samples ? m1->values : n_samples;
    labels = m1->data;
  }


  ks_introsort(i64, n_set_bits, set_bits);


  double *scores = _tk_matrix_chi2_scores(L, set_bits, n_set_bits, codes, labels, n_samples, n_visible, n_hidden, NULL, NULL, n_threads);
  tk_matrix_top_generic(L, scores, n_visible, n_hidden, top_k, 0);
  free(scores);
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
  uint64_t n_visible = tk_lua_checkunsigned(L, 3, "visible");
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
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "features");
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


static int tk_matrix_from_bitmap (lua_State *L)
{
  lua_settop(L, 3);
  const char *bm = tk_lua_checkustring(L, 1, "bitmap");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "features");
  kvec_t(int64_t) out;
  kv_init(out);
  for (uint64_t i = 0; i < n_samples; i ++)
    for (uint64_t j = 0; j < n_features; j ++) {
      uint64_t bit = i * n_features + j;
      uint64_t chunk = bit / CHAR_BIT;
      uint64_t pos = bit % CHAR_BIT;
      if (bm[chunk] & (1 << pos))
        kv_push(int64_t, out, (int64_t) bit);
    }
  kv_resize(int64_t, out, out.n);
  _tk_matrix_create(L, 1, out.n, out.a, NULL);
  return 1;
}


static int tk_matrix_extend_bits (lua_State *L)
{
  lua_settop(L, 4);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  int64_t *base_bits = m0->data;
  size_t n_base = m0->values;
  tk_matrix_t *m1 = tk_matrix_peek(L, 2);
  int64_t *ext_bits = m1->data;
  size_t n_ext = m1->values;
  uint64_t n_feat = tk_lua_checkunsigned(L, 3, "features");
  uint64_t n_extfeat = tk_lua_checkunsigned(L, 4, "extended");
  ks_introsort(int64_t, n_base, base_bits);
  ks_introsort(int64_t, n_ext, ext_bits);
  size_t total = n_base + n_ext;
  base_bits = tk_realloc(L, base_bits, total * sizeof(int64_t));
  for (size_t i = 0; i < n_base; i ++) {
    uint64_t bit = (uint64_t) base_bits[i];
    uint64_t sample = bit / n_feat;
    uint64_t old_off = sample * n_feat;
    uint64_t new_off = sample * (n_feat + n_extfeat);
    base_bits[i] = (int64_t) (bit - old_off + new_off);
  }
  for (size_t i = 0; i < n_ext; ++i) {
    uint64_t bit = (uint64_t) ext_bits[i];
    uint64_t sample = bit / n_extfeat;
    uint64_t old_off = sample * n_extfeat;
    uint64_t new_off = sample * (n_feat + n_extfeat);
    base_bits[n_base + i] = (int64_t) (bit - old_off + new_off);
  }
  ks_introsort(int64_t, total, base_bits);
  m0->data = base_bits;
  m0->rows = 1;
  m0->columns = total;
  m0->values = total;
  return 0;
}

static luaL_Reg tk_matrix_extra_fns[] =
{
  { "flip_interleave", tk_matrix_flip_interleave },
  { "top_chi2", tk_matrix_top_chi2 },
  { "top_mi", tk_matrix_top_mi },
  { "score_chi2", tk_matrix_score_chi2 },
  { "score_mi", tk_matrix_score_mi },
  { "filter", tk_matrix_filter },
  { "raw_bitmap", tk_matrix_raw_bitmap },
  { "from_bitmap", tk_matrix_from_bitmap },
  { "extend_bits", tk_matrix_extend_bits },
  { NULL, NULL }
};
