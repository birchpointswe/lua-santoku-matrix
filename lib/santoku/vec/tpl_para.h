


#ifndef tk_parallel_sfx
#error "Must include santoku/parallel/tpl.h before this template"
#endif

static inline void tk_parallel_sfx(tk_vec_pfx(transpose)) (
  tk_vec_pfx(t) *m0,
  tk_vec_pfx(t) *m1,
  uint64_t cols
) {
  tk_vec_pfx(ensure)(m0, m1->n);
  uint64_t rows = m1->n / cols;
  TK_PARALLEL_FOR(collapse(2))
  for (uint64_t r = 0; r < rows; r ++)
    for (uint64_t c = 0; c < cols; c ++)
      m0->a[c * rows + r] = m1->a[r * cols + c];
}

#ifndef tk_vec_limited

static inline tk_vec_pfx(t) *tk_parallel_sfx(tk_vec_pfx(csums)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, cols, NULL, NULL);
  TK_PARALLEL_FOR(schedule(static))
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

static inline tk_vec_pfx(t) *tk_parallel_sfx(tk_vec_pfx(rsums)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, m0->n / cols, NULL, NULL);
TK_PARALLEL_FOR(schedule(static))
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

static inline tk_vec_pfx(t) *tk_parallel_sfx(tk_vec_pfx(cmaxs)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, cols, NULL, NULL);
TK_PARALLEL_FOR(schedule(static))
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

static inline tk_vec_pfx(t) *tk_parallel_sfx(tk_vec_pfx(rmaxs)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, m0->n / cols, NULL, NULL);
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    tk_vec_base sum = 0.0;
    for (uint64_t c = 0; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      if (val > sum)
        sum = val;
    }
    out->a[r] = sum;
  }
  return out;
}

static inline tk_vec_pfx(t) *tk_parallel_sfx(tk_vec_pfx(cmins)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, cols, NULL, NULL);
  TK_PARALLEL_FOR(schedule(static))
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

static inline tk_vec_pfx(t) *tk_parallel_sfx(tk_vec_pfx(rmins)) (
  lua_State *L,
  tk_vec_pfx(t) *m0,
  uint64_t cols
) {
  tk_vec_pfx(t) *out = tk_vec_pfx(create)(L, m0->n / cols, NULL, NULL);
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t r = 0; r < m0->n / cols; r ++) {
    tk_vec_base minv = m0->a[r * cols + 0];
    for (uint64_t c = 1; c < cols; c ++) {
      tk_vec_base val = m0->a[r * cols + c];
      if (val < minv)
        minv = val;
    }
    out->a[r] = minv;
  }
  return out;
}

static inline void tk_parallel_sfx(tk_vec_pfx(scale)) (
  tk_vec_pfx(t) *m0,
  tk_vec_base scale,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start >= end || start >= m0->n)
    return;
  if (end > m0->n)
    end = m0->n;
  TK_PARALLEL_FOR(schedule(static))
  for (size_t i = start; i < end; i ++)
    m0->a[i] *= scale;
}

static inline void tk_parallel_sfx(tk_vec_pfx(add)) (
  tk_vec_pfx(t) *m0,
  tk_vec_base add,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start >= end || start >= m0->n)
    return;
  if (end > m0->n)
    end = m0->n;
  TK_PARALLEL_FOR(schedule(static))
  for (size_t i = start; i < end; i ++) {
    m0->a[i] += add;
  }
}

static inline void tk_parallel_sfx(tk_vec_pfx(add_scaled)) (
  tk_vec_pfx(t) *m0,
  tk_vec_base add,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start >= end || start >= m0->n)
    return;
  if (end > m0->n)
    end = m0->n;
  TK_PARALLEL_FOR(schedule(static))
  for (size_t i = start; i < end; i ++) {
    m0->a[i] += add * (tk_vec_base) i;
  }
}

static inline void tk_parallel_sfx(tk_vec_pfx(scalev)) (
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
  TK_PARALLEL_FOR(schedule(static))
  for (size_t i = start; i < end; i ++) {
    m0->a[i] *= m1->a[i];
  }
}

static inline void tk_parallel_sfx(tk_vec_pfx(addv)) (
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
  TK_PARALLEL_FOR(schedule(static))
  for (size_t i = start; i < end; i ++) {
    m0->a[i] += m1->a[i];
  }
}

#ifdef tk_vec_abs
static inline void tk_parallel_sfx(tk_vec_pfx(abs)) (
  tk_vec_pfx(t) *m0,
  uint64_t start,
  uint64_t end
) {
  if (!m0->n || start >= end || start >= m0->n)
    return;
  if (end > m0->n)
    end = m0->n;
  TK_PARALLEL_FOR(schedule(static))
  for (size_t i = start; i < end; i ++)
    m0->a[i] = tk_vec_abs(m0->a[i]);
}
#endif

static inline double tk_parallel_sfx(tk_vec_pfx(dot)) (tk_vec_pfx(t) *a, tk_vec_pfx(t) *b) {
  if (a->n != b->n)
    return 0;
  size_t n = a->n;
  double sum = 0.0;
  TK_PARALLEL_FOR(reduction(+:sum))
  for (size_t i = 0; i < n; i ++)
    sum += a->a[i] * b->a[i];
  return sum;
}

static inline void tk_parallel_sfx(tk_vec_pfx(multiply)) (tk_vec_pfx(t) *a, tk_vec_pfx(t) *b, tk_vec_pfx(t) *c, uint64_t k, bool transpose_a, bool transpose_b) {
  size_t m = transpose_a ? k : a->n / k;
  size_t n = transpose_b ? k : b->n / k;
  tk_vec_pfx(ensure)(c, m * n);
  c->n = m * n;
  TK_PARALLEL_FOR(collapse(2))
  for (size_t i = 0; i < m; i ++) {
    for (size_t j = 0; j < n; j ++) {
      tk_vec_base sum = 0;
      for (size_t p = 0; p < k; p ++) {
        tk_vec_base a_val = transpose_a ? a->a[p * m + i] : a->a[i * k + p];
        tk_vec_base b_val = transpose_b ? b->a[j * k + p] : b->a[p * n + j];
        sum += a_val * b_val;
      }
      c->a[i * n + j] = sum;
    }
  }
}

static inline void tk_parallel_sfx(tk_vec_pfx(pow)) (
  tk_vec_pfx(t) *v,
  double exponent,
  uint64_t start,
  uint64_t end
) {
  if (end <= start)
    return;
  if (end > v->n) {
    tk_vec_pfx(ensure)(v, end);
    v->n = end;
  }
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = start; i < end; i ++)
    v->a[i] = pow(v->a[i], exponent);
}

static inline void tk_parallel_sfx(tk_vec_pfx(log)) (
  tk_vec_pfx(t) *v,
  uint64_t start,
  uint64_t end
) {
  if (end <= start)
    return;
  if (end > v->n) {
    tk_vec_pfx(ensure)(v, end);
    v->n = end;
  }
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = start; i < end; i ++)
    v->a[i] = log(v->a[i]);
}

static inline void tk_parallel_sfx(tk_vec_pfx(exp)) (
  tk_vec_pfx(t) *v,
  uint64_t start,
  uint64_t end
) {
  if (end <= start)
    return;
  if (end > v->n) {
    tk_vec_pfx(ensure)(v, end);
    v->n = end;
  }
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = start; i < end; i ++)
    v->a[i] = exp(v->a[i]);
}

static inline void tk_parallel_sfx(tk_vec_pfx(fill_indices)) (tk_vec_pfx(t) *v) {
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < v->n; i++)
    v->a[i] = (tk_vec_base) i;
}

static inline void tk_parallel_sfx(tk_vec_pfx(clamp)) (tk_vec_pfx(t) *v, tk_vec_base lo, tk_vec_base hi, uint64_t start, uint64_t end) {
  if (end <= start)
    return;
  if (end > v->n)
    end = v->n;
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = start; i < end; i ++) {
    if (tk_vec_lt(v->a[i], lo)) v->a[i] = lo;
    else if (tk_vec_gt(v->a[i], hi)) v->a[i] = hi;
  }
}

#endif

static inline void tk_parallel_sfx(tk_vec_pfx(fill)) (tk_vec_pfx(t) *v, tk_vec_base x, uint64_t start, uint64_t end) {
  if (end <= start)
    return;
  if (end > v->n) {
    tk_vec_pfx(ensure)(v, end);
    v->n = end;
  }
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = start; i < end; i ++)
    v->a[i] = x;
}
