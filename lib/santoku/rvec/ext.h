#ifndef TK_RVEC_EXT_H
#define TK_RVEC_EXT_H

#if defined(_OPENMP) && !defined(__EMSCRIPTEN__)
#include <omp.h>
#endif
#include <float.h>
#include <stdlib.h>
#include <santoku/dumap.h>
#include <santoku/iumap.h>
#include <santoku/pvec.h>
#include <santoku/evec.h>
#include <santoku/dvec/base.h>
#include <math.h>

static inline tk_rvec_t *tk_rvec_from_dvec (
  lua_State *L,
  tk_dvec_t *D
) {
  tk_rvec_t *R = tk_rvec_create(L, D->n);
  for (int64_t i = 0; i < (int64_t) D->n; i ++)
    R->a[i] = tk_rank(i, D->a[i]);
  return R;
}

static inline int tk_rvec_split (
  tk_rvec_t *P,
  tk_ivec_t *K,
  tk_dvec_t *V,
  bool append
) {
  if (!append) {
    tk_ivec_clear(K);
    tk_dvec_clear(V);
  }
  if (tk_ivec_ensure(K, K->n + P->n) != 0)
    return -1;
  if (tk_dvec_ensure(V, V->n + P->n) != 0)
    return -1;
  for (uint64_t i = 0; i < P->n; i ++) {
    tk_rank_t r = P->a[i];
    K->a[K->n ++] = r.i;
    V->a[V->n ++] = r.d;
  }
  return 0;
}

static inline tk_ivec_t *tk_rvec_keys (
  lua_State *L,
  tk_rvec_t *P,
  tk_ivec_t *out
) {
  tk_ivec_t *result = out ? out : tk_ivec_create(L, P->n);
  if (out)
    tk_ivec_ensure(result, P->n);
  for (uint64_t i = 0; i < P->n; i ++)
    result->a[i] = P->a[i].i;
  if (out)
    result->n = P->n;
  return result;
}

static inline tk_dvec_t *tk_rvec_values (
  lua_State *L,
  tk_rvec_t *P,
  tk_dvec_t *out
) {
  tk_dvec_t *result = out ? out : tk_dvec_create(L, P->n);
  if (out)
    tk_dvec_ensure(result, P->n);
  for (uint64_t i = 0; i < P->n; i ++)
    result->a[i] = P->a[i].d;
  if (out)
    result->n = P->n;
  return result;
}

static inline int tk_rvec_each_lua_iter (lua_State *L)
{
  lua_settop(L, 0);
  tk_rvec_t *m0 = tk_rvec_peek(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_rank_t v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, v.i);
  lua_pushnumber(L, v.d);
  return 2;
}

static inline int tk_rvec_ieach_lua_iter (lua_State *L)
{
  lua_settop(L, 0);
  tk_rvec_t *m0 = tk_rvec_peek(L, lua_upvalueindex(1), "vector");
  uint64_t n = tk_lua_checkunsigned(L, lua_upvalueindex(2), "idx");
  if (n >= m0->n)
    return 0;
  tk_rank_t v = m0->a[n];
  lua_pushinteger(L, (int64_t) n + 1);
  lua_replace(L, lua_upvalueindex(2));
  lua_pushinteger(L, (int64_t) n);
  lua_pushinteger(L, v.i);
  lua_pushnumber(L, v.d);
  return 3;
}

static inline int tk_rvec_each_lua (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_rvec_each_lua_iter, 2);
  return 1;
}

static inline int tk_rvec_ieach_lua (lua_State *L)
{
  lua_settop(L, 1);
  lua_pushinteger(L, 0);
  lua_pushcclosure(L, tk_rvec_ieach_lua_iter, 2);
  return 1;
}

static inline void tk_rvec_ranks (
  tk_rvec_t *v,
  tk_dumap_t *r
) {
  tk_dumap_clear(r);
  int kha;
  uint32_t khi;
  for (size_t i = 0; i < v->n; i ++) {
    double rank = (double) i;
    size_t count = 1;
    while (i + 1 < v->n && v->a[i].d == v->a[i + 1].d) {
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




static inline size_t tk_rvec_scores_max_curvature (
  tk_rvec_t *v,
  double *out_val
) {
  size_t n = v->n;
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? v->a[0].d : 0.0;
    return n > 0 ? n - 1 : 0;
  }
  double max_curv = 0.0;
  size_t max_idx = 1;
  for (size_t i = 1; i < n - 1; i++) {
    double curv = fabs(v->a[i-1].d - 2.0 * v->a[i].d + v->a[i+1].d);
    if (curv > max_curv) {
      max_curv = curv;
      max_idx = i;
    }
  }

  if (max_curv < 1e-10) {
    if (out_val) *out_val = v->a[n - 1].d;
    return n - 1;
  }
  if (out_val) *out_val = v->a[max_idx].d;
  return max_idx;
}

static inline size_t tk_rvec_scores_lmethod (
  tk_rvec_t *v,
  double *out_val
) {
  size_t n = v->n;
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? v->a[0].d : 0.0;
    return n > 0 ? n - 1 : 0;
  }

  double min_val = v->a[0].d, max_val = v->a[0].d;
  for (size_t i = 1; i < n; i++) {
    if (v->a[i].d < min_val) min_val = v->a[i].d;
    if (v->a[i].d > max_val) max_val = v->a[i].d;
  }
  if (max_val - min_val < 1e-10) {
    if (out_val) *out_val = v->a[n - 1].d;
    return n - 1;
  }
  double best_rmse = DBL_MAX;
  size_t best_k = 1;
  for (size_t k = 1; k < n - 1; k++) {
    double sum_x1 = 0.0, sum_y1 = 0.0, sum_xy1 = 0.0, sum_xx1 = 0.0;
    for (size_t i = 0; i <= k; i++) {
      double x = (double)i;
      double y = v->a[i].d;
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
      double y = v->a[i].d;
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
      double err = v->a[i].d - pred;
      sse += err * err;
    }
    for (size_t i = k + 1; i < n; i++) {
      double pred = slope2 * (double)i + intercept2;
      double err = v->a[i].d - pred;
      sse += err * err;
    }
    double rmse = sqrt(sse / (double)n);
    if (rmse < best_rmse) {
      best_rmse = rmse;
      best_k = k;
    }
  }
  if (out_val) *out_val = v->a[best_k].d;
  return best_k;
}

static inline size_t tk_rvec_scores_max_gap (
  tk_rvec_t *v,
  double *out_val
) {
  size_t n = v->n;
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? v->a[0].d : 0.0;
    return n > 0 ? n - 1 : 0;
  }
  double max_gap = 0.0;
  size_t max_idx = 0;
  for (size_t i = 0; i < n - 1; i++) {
    double gap = fabs(v->a[i + 1].d - v->a[i].d);
    if (gap > max_gap) {
      max_gap = gap;
      max_idx = i;
    }
  }

  if (max_gap < 1e-10) {
    if (out_val) *out_val = v->a[n - 1].d;
    return n - 1;
  }
  if (out_val) *out_val = v->a[max_idx].d;
  return max_idx;
}

static inline size_t tk_rvec_scores_plateau (
  tk_rvec_t *v,
  double tolerance,
  double *out_val
) {
  size_t n = v->n;
  if (n == 0) {
    if (out_val) *out_val = 0.0;
    return 0;
  }
  if (n == 1) {
    if (out_val) *out_val = v->a[0].d;
    return 0;
  }
  if (tolerance <= 0.0) tolerance = 0.01;

  double min_score = v->a[0].d, max_score = v->a[0].d;
  for (size_t i = 1; i < n; i++) {
    if (v->a[i].d < min_score) min_score = v->a[i].d;
    if (v->a[i].d > max_score) max_score = v->a[i].d;
  }
  double range = max_score - min_score;
  if (range <= 0.0) {
    if (out_val) *out_val = v->a[n - 1].d;
    return n - 1;
  }

  double abs_tolerance = tolerance * range;
  double base = v->a[0].d;
  size_t end_idx = 0;
  for (size_t i = 1; i < n; i++) {
    if (fabs(v->a[i].d - base) <= abs_tolerance) {
      end_idx = i;
    } else {
      break;
    }
  }
  if (out_val) *out_val = v->a[end_idx].d;
  return end_idx;
}

static inline size_t tk_rvec_scores_kneedle (
  tk_rvec_t *v,
  double sensitivity,
  double *out_val
) {
  size_t n = v->n;
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? v->a[0].d : 0.0;
    return n > 0 ? n - 1 : 0;
  }
  if (sensitivity <= 0.0)
    sensitivity = 1.0;
  double min_score = v->a[0].d;
  double max_score = v->a[0].d;
  for (size_t i = 1; i < n; i++) {
    if (v->a[i].d < min_score) min_score = v->a[i].d;
    if (v->a[i].d > max_score) max_score = v->a[i].d;
  }
  double score_range = max_score - min_score;
  if (score_range < 1e-10) {
    if (out_val) *out_val = v->a[0].d;
    return n - 1;
  }
  double *normalized = (double *)malloc(n * sizeof(double));
  if (!normalized) {
    if (out_val) *out_val = v->a[0].d;
    return n - 1;
  }
  for (size_t i = 0; i < n; i++) {
    double x_norm = (double)i / (double)(n - 1);
    double y_norm = (v->a[i].d - min_score) / score_range;
    normalized[i] = y_norm - x_norm;
  }
  double *smoothed = (double *)malloc(n * sizeof(double));
  if (!smoothed) {
    free(normalized);
    if (out_val) *out_val = v->a[0].d;
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
    double y_norm = (v->a[i].d - min_score) / score_range;
    double diff = y_norm - x_norm;
    if (diff >= threshold) {
      final_knee = i;
    }
  }
  if (out_val) *out_val = v->a[final_knee].d;
  return final_knee;
}




static inline size_t tk_rvec_scores_first_gap (
  tk_rvec_t *v,
  double alpha,
  double *out_val
) {
  size_t n = v->n;
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? v->a[0].d : 0.0;
    return n > 0 ? n - 1 : 0;
  }
  if (alpha <= 0.0) alpha = 3.0;

  size_t n_gaps = n - 1;
  double *gaps = (double *)malloc(n_gaps * sizeof(double));
  if (!gaps) {
    if (out_val) *out_val = v->a[n - 1].d;
    return n - 1;
  }

  for (size_t i = 0; i < n_gaps; i++) {
    gaps[i] = fabs(v->a[i + 1].d - v->a[i].d);
  }

  for (size_t i = 1; i < n_gaps; i++) {
    double key = gaps[i];
    size_t j = i;
    while (j > 0 && gaps[j - 1] > key) {
      gaps[j] = gaps[j - 1];
      j--;
    }
    gaps[j] = key;
  }

  double median_gap;
  if (n_gaps % 2 == 1) {
    median_gap = gaps[n_gaps / 2];
  } else {
    median_gap = (gaps[n_gaps / 2 - 1] + gaps[n_gaps / 2]) / 2.0;
  }
  free(gaps);

  if (median_gap <= 0.0) {
    double max_gap = 0.0;
    size_t max_idx = n - 1;
    for (size_t i = 0; i < n - 1; i++) {
      double gap = fabs(v->a[i + 1].d - v->a[i].d);
      if (gap > max_gap) {
        max_gap = gap;
        max_idx = i;
      }
    }
    if (out_val) *out_val = v->a[max_idx].d;
    return max_idx;
  }

  double threshold = alpha * median_gap;

  for (size_t i = 0; i < n - 1; i++) {
    double gap = fabs(v->a[i + 1].d - v->a[i].d);
    if (gap > threshold) {
      if (out_val) *out_val = v->a[i].d;
      return i;
    }
  }

  if (out_val) *out_val = v->a[n - 1].d;
  return n - 1;
}




static inline size_t tk_rvec_scores_otsu (
  tk_rvec_t *v,
  double *out_val
) {
  size_t n = v->n;
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? v->a[0].d : 0.0;
    return n > 0 ? n - 1 : 0;
  }


  double min_val = v->a[0].d, max_val = v->a[0].d;
  for (size_t i = 1; i < n; i++) {
    if (v->a[i].d < min_val) min_val = v->a[i].d;
    if (v->a[i].d > max_val) max_val = v->a[i].d;
  }
  if (max_val - min_val < 1e-10) {
    if (out_val) *out_val = v->a[n - 1].d;
    return n - 1;
  }


  double total_sum = 0.0;
  for (size_t i = 0; i < n; i++) {
    total_sum += v->a[i].d;
  }




  double best_variance = -1.0;
  size_t best_k = 0;
  double sum0 = 0.0;

  for (size_t k = 0; k < n - 1; k++) {
    sum0 += v->a[k].d;
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

  if (out_val) *out_val = v->a[best_k].d;
  return best_k;
}

#endif
