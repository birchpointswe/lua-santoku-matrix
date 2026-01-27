


#ifndef tk_parallel_sfx
#error "Must include santoku/parallel/tpl.h before this template"
#endif

static inline void tk_parallel_sfx(tk_vec_pfx(copy_indexed)) (
  tk_vec_pfx(t) *m0,
  tk_vec_pfx(t) *m1,
  tk_ivec_t *indices
) {
  tk_vec_pfx(ensure)(m0, indices->n);
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < indices->n; i++) {
    int64_t idx = indices->a[i];
    if (idx >= 0 && idx < (int64_t)m1->n) {
      m0->a[i] = m1->a[idx];
    }
  }
  if (m0->n < indices->n)
    m0->n = indices->n;
}

static inline void tk_parallel_sfx(tk_vec_pfx(scatter_indexed)) (
  tk_vec_pfx(t) *m0,
  tk_vec_pfx(t) *m1,
  tk_ivec_t *indices
) {
  int64_t max_idx = 0;
  for (uint64_t i = 0; i < indices->n; i++) {
    if (indices->a[i] > max_idx)
      max_idx = indices->a[i];
  }
  tk_vec_pfx(ensure)(m0, (uint64_t)(max_idx + 1));
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < indices->n; i++) {
    int64_t idx = indices->a[i];
    if (idx >= 0 && i < m1->n) {
      m0->a[idx] = m1->a[i];
    }
  }
  if (m0->n < (uint64_t)(max_idx + 1))
    m0->n = (uint64_t)(max_idx + 1);
}

#ifndef tk_vec_limited

static inline tk_ivec_t *tk_parallel_sfx(tk_vec_pfx(rasc)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, m0->n, NULL, NULL);
  uint64_t rows = m0->n / cols;
  TK_PARALLEL
  {
    tk_rvec_t *ranks = tk_rvec_create(NULL, cols, NULL, NULL);
    TK_FOR(schedule(static))
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

static inline tk_ivec_t *tk_parallel_sfx(tk_vec_pfx(rdesc)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, m0->n, NULL, NULL);
  uint64_t rows = m0->n / cols;
  TK_PARALLEL
  {
    tk_rvec_t *ranks = tk_rvec_create(NULL, cols, NULL, NULL);
    TK_FOR(schedule(static))
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

static inline tk_ivec_t *tk_parallel_sfx(tk_vec_pfx(casc)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  uint64_t rows = m0->n / cols;
  tk_ivec_t *out = tk_ivec_create(L, m0->n, NULL, NULL);
  TK_PARALLEL
  {
    tk_rvec_t *ranks = tk_rvec_create(0, rows, NULL, NULL);
    TK_FOR(schedule(static))
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

static inline tk_ivec_t *tk_parallel_sfx(tk_vec_pfx(cdesc)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  uint64_t rows = m0->n / cols;
  tk_ivec_t *out = tk_ivec_create(L, m0->n, NULL, NULL);
  TK_PARALLEL
  {
    tk_rvec_t *ranks = tk_rvec_create(0, rows, NULL, NULL);
    TK_FOR(schedule(static))
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

static inline tk_dvec_t *tk_parallel_sfx(tk_vec_pfx(cmagnitudes)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_dvec_t *out = tk_dvec_create(L, cols, NULL, NULL);
  uint64_t rows = m0->n / cols;
  TK_PARALLEL_FOR(schedule(static))
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

static inline tk_dvec_t *tk_parallel_sfx(tk_vec_pfx(rmagnitudes)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  uint64_t rows = m0->n / cols;
  tk_dvec_t *out = tk_dvec_create(L, rows, NULL, NULL);
  TK_PARALLEL_FOR(schedule(static))
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

static inline tk_ivec_t *tk_parallel_sfx(tk_vec_pfx(cmaxargs)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, cols, NULL, NULL);
  uint64_t rows = m0->n / cols;
  TK_PARALLEL_FOR(schedule(static))
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

static inline tk_ivec_t *tk_parallel_sfx(tk_vec_pfx(rmaxargs)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  uint64_t rows = m0->n / cols;
  tk_ivec_t *out = tk_ivec_create(L, rows, NULL, NULL);
  TK_PARALLEL_FOR(schedule(static))
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

static inline tk_ivec_t *tk_parallel_sfx(tk_vec_pfx(cminargs)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_ivec_t *out = tk_ivec_create(L, cols, NULL, NULL);
  uint64_t rows = m0->n / cols;
  TK_PARALLEL_FOR(schedule(static))
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

static inline tk_ivec_t *tk_parallel_sfx(tk_vec_pfx(rminargs)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  uint64_t rows = m0->n / cols;
  tk_ivec_t *out = tk_ivec_create(L, rows, NULL, NULL);
  TK_PARALLEL_FOR(schedule(static))
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

#endif
