#ifndef TK_PVEC_EXT_H
#define TK_PVEC_EXT_H

#include <math.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <santoku/dumap.h>

static inline tk_pvec_t *tk_pvec_from_ivec (
  lua_State *L,
  tk_ivec_t *ivec
) {
  tk_pvec_t *P = tk_pvec_create(L, ivec->n, 0, 0);
  for (int64_t i = 0; i < (int64_t) ivec->n; i ++)
    P->a[i] = tk_pair(i, ivec->a[i]);
  return P;
}

static inline tk_ivec_t *tk_pvec_keys (
  lua_State *L,
  tk_pvec_t *P,
  tk_ivec_t *out
) {
  tk_ivec_t *result = out ? out : tk_ivec_create(L, P->n, 0, 0);
  if (out)
    tk_ivec_ensure(result, P->n);
  for (uint64_t i = 0; i < P->n; i ++)
    result->a[i] = P->a[i].i;
  if (out)
    result->n = P->n;
  return result;
}

static inline tk_ivec_t *tk_pvec_values (
  lua_State *L,
  tk_pvec_t *P,
  tk_ivec_t *out
) {
  tk_ivec_t *result = out ? out : tk_ivec_create(L, P->n, 0, 0);
  if (out)
    tk_ivec_ensure(result, P->n);
  for (uint64_t i = 0; i < P->n; i ++)
    result->a[i] = P->a[i].p;
  if (out)
    result->n = P->n;
  return result;
}

static inline int tk_pvec_each_lua_iter (lua_State *L)
{
  lua_settop(L, 0);
  tk_pvec_t *m0 = tk_pvec_peek(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_pair_t v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, v.i);
  lua_pushinteger(L, v.p);
  return 2;
}

static inline int tk_pvec_ieach_lua_iter (lua_State *L)
{
  lua_settop(L, 0);
  tk_pvec_t *m0 = tk_pvec_peek(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_pair_t v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, (int64_t) n);
  lua_pushinteger(L, v.i);
  lua_pushinteger(L, v.p);
  return 3;
}

static inline int tk_pvec_each_lua (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_pvec_each_lua_iter, 2);
  return 1;
}

static inline int tk_pvec_ieach_lua (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_pvec_ieach_lua_iter, 2);
  return 1;
}

static inline void tk_pvec_ranks (
  tk_pvec_t *v,
  tk_dumap_t *r
) {
  tk_dumap_clear(r);
  int kha;
  uint32_t khi;
  for (size_t i = 0; i < v->n; i ++) {
    double rank = (double) i;
    size_t count = 1;
    while (i + 1 < v->n && v->a[i].p == v->a[i + 1].p) {
      count ++;
      i ++;
    }
    double average_rank = (rank + (rank + count - 1)) / 2.0;
    for (size_t j = 0; j < count; j ++) {
      khi = tk_dumap_put(r, v->a[i - j].i, &kha);
      tk_dumap_setval(r, khi, average_rank);
    }
  }
}




static inline size_t tk_pvec_scores_max_curvature (
  tk_pvec_t *v,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? v->a[0].p : 0;
    return n > 0 ? n - 1 : 0;
  }
  int64_t max_curv = 0;
  size_t max_idx = 1;
  for (size_t i = 1; i < n - 1; i++) {
    int64_t curv = llabs(v->a[i-1].p - 2 * v->a[i].p + v->a[i+1].p);
    if (curv > max_curv) {
      max_curv = curv;
      max_idx = i;
    }
  }

  if (max_curv == 0) {
    if (out_val) *out_val = v->a[n - 1].p;
    return n - 1;
  }
  if (out_val) *out_val = v->a[max_idx].p;
  return max_idx;
}

static inline size_t tk_pvec_scores_lmethod (
  tk_pvec_t *v,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? v->a[0].p : 0;
    return n > 0 ? n - 1 : 0;
  }

  int64_t min_val = v->a[0].p, max_val = v->a[0].p;
  for (size_t i = 1; i < n; i++) {
    if (v->a[i].p < min_val) min_val = v->a[i].p;
    if (v->a[i].p > max_val) max_val = v->a[i].p;
  }
  if (max_val == min_val) {
    if (out_val) *out_val = v->a[n - 1].p;
    return n - 1;
  }
  double best_rmse = DBL_MAX;
  size_t best_k = 1;
  for (size_t k = 1; k < n - 1; k++) {
    double sum_x1 = 0.0, sum_y1 = 0.0, sum_xy1 = 0.0, sum_xx1 = 0.0;
    for (size_t i = 0; i <= k; i++) {
      double x = (double)i;
      double y = (double)v->a[i].p;
      sum_x1 += x;
      sum_y1 += y;
      sum_xy1 += x * y;
      sum_xx1 += x * x;
    }
    size_t n1 = k + 1;
    double mean_x1 = sum_x1 / (double)n1;
    double mean_y1 = sum_y1 / (double)n1;
    double slope1 = (sum_xy1 - (double)n1 * mean_x1 * mean_y1) / (sum_xx1 - (double)n1 * mean_x1 * mean_x1 + 1e-10);
    double intercept1 = mean_y1 - slope1 * mean_x1;
    double sum_x2 = 0.0, sum_y2 = 0.0, sum_xy2 = 0.0, sum_xx2 = 0.0;
    for (size_t i = k + 1; i < n; i++) {
      double x = (double)i;
      double y = (double)v->a[i].p;
      sum_x2 += x;
      sum_y2 += y;
      sum_xy2 += x * y;
      sum_xx2 += x * x;
    }
    size_t n2 = n - k - 1;
    double mean_x2 = sum_x2 / (double)n2;
    double mean_y2 = sum_y2 / (double)n2;
    double slope2 = (sum_xy2 - (double)n2 * mean_x2 * mean_y2) / (sum_xx2 - (double)n2 * mean_x2 * mean_x2 + 1e-10);
    double intercept2 = mean_y2 - slope2 * mean_x2;
    double sse = 0.0;
    for (size_t i = 0; i <= k; i++) {
      double pred = slope1 * (double)i + intercept1;
      double err = (double)v->a[i].p - pred;
      sse += err * err;
    }
    for (size_t i = k + 1; i < n; i++) {
      double pred = slope2 * (double)i + intercept2;
      double err = (double)v->a[i].p - pred;
      sse += err * err;
    }
    double rmse = sqrt(sse / (double)n);
    if (rmse < best_rmse) {
      best_rmse = rmse;
      best_k = k;
    }
  }
  if (out_val) *out_val = v->a[best_k].p;
  return best_k;
}

static inline size_t tk_pvec_scores_max_gap (
  tk_pvec_t *v,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? v->a[0].p : 0;
    return n > 0 ? n - 1 : 0;
  }
  int64_t max_gap = 0;
  size_t max_idx = 0;
  for (size_t i = 0; i < n - 1; i++) {
    int64_t gap = llabs(v->a[i + 1].p - v->a[i].p);
    if (gap > max_gap) {
      max_gap = gap;
      max_idx = i;
    }
  }

  if (max_gap == 0) {
    if (out_val) *out_val = v->a[n - 1].p;
    return n - 1;
  }
  if (out_val) *out_val = v->a[max_idx].p;
  return max_idx;
}

static inline size_t tk_pvec_scores_plateau (
  tk_pvec_t *v,
  double tolerance,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n == 0) {
    if (out_val) *out_val = 0;
    return 0;
  }
  if (n == 1) {
    if (out_val) *out_val = v->a[0].p;
    return 0;
  }
  if (tolerance <= 0.0) tolerance = 0.01;

  int64_t min_score = v->a[0].p, max_score = v->a[0].p;
  for (size_t i = 1; i < n; i++) {
    if (v->a[i].p < min_score) min_score = v->a[i].p;
    if (v->a[i].p > max_score) max_score = v->a[i].p;
  }
  int64_t range = max_score - min_score;
  if (range <= 0) {
    if (out_val) *out_val = v->a[n - 1].p;
    return n - 1;
  }

  int64_t abs_tolerance = (int64_t)(tolerance * range);
  int64_t base = v->a[0].p;
  size_t end_idx = 0;
  for (size_t i = 1; i < n; i++) {
    if (llabs(v->a[i].p - base) <= abs_tolerance) {
      end_idx = i;
    } else {
      break;
    }
  }
  if (out_val) *out_val = v->a[end_idx].p;
  return end_idx;
}

static inline size_t tk_pvec_scores_kneedle (
  tk_pvec_t *v,
  double sensitivity,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? v->a[0].p : 0;
    return n > 0 ? n - 1 : 0;
  }
  if (sensitivity <= 0.0)
    sensitivity = 1.0;
  int64_t min_score = v->a[0].p;
  int64_t max_score = v->a[0].p;
  for (size_t i = 1; i < n; i++) {
    if (v->a[i].p < min_score) min_score = v->a[i].p;
    if (v->a[i].p > max_score) max_score = v->a[i].p;
  }
  double score_range = (double)(max_score - min_score);
  if (score_range < 1e-10) {
    if (out_val) *out_val = v->a[0].p;
    return n - 1;
  }
  double *normalized = (double *)malloc(n * sizeof(double));
  if (!normalized) {
    if (out_val) *out_val = v->a[0].p;
    return n - 1;
  }
  for (size_t i = 0; i < n; i++) {
    double x_norm = (double)i / (double)(n - 1);
    double y_norm = ((double)v->a[i].p - (double)min_score) / score_range;
    normalized[i] = y_norm - x_norm;
  }
  double *smoothed = (double *)malloc(n * sizeof(double));
  if (!smoothed) {
    free(normalized);
    if (out_val) *out_val = v->a[0].p;
    return n - 1;
  }
  smoothed[0] = normalized[0];
  smoothed[n - 1] = normalized[n - 1];
  for (size_t i = 1; i < n - 1; i++) {
    smoothed[i] = (normalized[i - 1] + normalized[i] + normalized[i + 1]) / 3.0;
  }
  double max_diff = -DBL_MAX;
  size_t knee_idx = 0;
  for (size_t i = 0; i < n; i++) {
    if (smoothed[i] > max_diff) {
      max_diff = smoothed[i];
      knee_idx = i;
    }
  }
  free(normalized);
  free(smoothed);
  double threshold = max_diff - sensitivity / (double)n;
  size_t final_knee = knee_idx;
  for (size_t i = 0; i < n; i++) {
    double x_norm = (double)i / (double)(n - 1);
    double y_norm = ((double)v->a[i].p - (double)min_score) / score_range;
    double diff = y_norm - x_norm;
    if (diff >= threshold) {
      final_knee = i;
    }
  }
  if (out_val) *out_val = v->a[final_knee].p;
  return final_knee;
}




static inline size_t tk_pvec_scores_first_gap (
  tk_pvec_t *v,
  double alpha,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? v->a[0].p : 0;
    return n > 0 ? n - 1 : 0;
  }
  if (alpha <= 0.0) alpha = 3.0;

  size_t n_gaps = n - 1;
  int64_t *gaps = (int64_t *)malloc(n_gaps * sizeof(int64_t));
  if (!gaps) {
    if (out_val) *out_val = v->a[n - 1].p;
    return n - 1;
  }

  for (size_t i = 0; i < n_gaps; i++) {
    gaps[i] = llabs(v->a[i + 1].p - v->a[i].p);
  }

  for (size_t i = 1; i < n_gaps; i++) {
    int64_t key = gaps[i];
    size_t j = i;
    while (j > 0 && gaps[j - 1] > key) {
      gaps[j] = gaps[j - 1];
      j--;
    }
    gaps[j] = key;
  }

  int64_t median_gap;
  if (n_gaps % 2 == 1) {
    median_gap = gaps[n_gaps / 2];
  } else {
    median_gap = (gaps[n_gaps / 2 - 1] + gaps[n_gaps / 2]) / 2;
  }
  free(gaps);

  if (median_gap == 0) {
    int64_t max_gap = 0;
    size_t max_idx = n - 1;
    for (size_t i = 0; i < n - 1; i++) {
      int64_t gap = llabs(v->a[i + 1].p - v->a[i].p);
      if (gap > max_gap) {
        max_gap = gap;
        max_idx = i;
      }
    }
    if (out_val) *out_val = v->a[max_idx].p;
    return max_idx;
  }

  int64_t threshold = (int64_t)(alpha * (double)median_gap);
  if (threshold < 1) threshold = 1;

  for (size_t i = 0; i < n - 1; i++) {
    int64_t gap = llabs(v->a[i + 1].p - v->a[i].p);
    if (gap > threshold) {
      if (out_val) *out_val = v->a[i].p;
      return i;
    }
  }

  if (out_val) *out_val = v->a[n - 1].p;
  return n - 1;
}




static inline size_t tk_pvec_scores_otsu (
  tk_pvec_t *v,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? v->a[0].p : 0;
    return n > 0 ? n - 1 : 0;
  }


  int64_t min_val = v->a[0].p, max_val = v->a[0].p;
  for (size_t i = 1; i < n; i++) {
    if (v->a[i].p < min_val) min_val = v->a[i].p;
    if (v->a[i].p > max_val) max_val = v->a[i].p;
  }
  if (max_val == min_val) {
    if (out_val) *out_val = v->a[n - 1].p;
    return n - 1;
  }


  double total_sum = 0.0;
  for (size_t i = 0; i < n; i++) {
    total_sum += (double)v->a[i].p;
  }




  double best_variance = -1.0;
  size_t best_k = 0;
  double sum0 = 0.0;

  for (size_t k = 0; k < n - 1; k++) {
    sum0 += (double)v->a[k].p;
    double sum1 = total_sum - sum0;

    size_t n0 = k + 1;
    size_t n1 = n - n0;

    double w0 = (double)n0 / (double)n;
    double w1 = (double)n1 / (double)n;

    double mean0 = sum0 / (double)n0;
    double mean1 = sum1 / (double)n1;


    double variance = w0 * w1 * (mean0 - mean1) * (mean0 - mean1);

    if (variance > best_variance) {
      best_variance = variance;
      best_k = k;
    }
  }

  if (out_val) *out_val = v->a[best_k].p;
  return best_k;
}

#endif
