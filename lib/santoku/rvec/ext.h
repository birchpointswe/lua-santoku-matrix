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
  tk_rvec_t *R = tk_rvec_create(L, D->n, 0, 0);
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
  tk_ivec_t *result = out ? out : tk_ivec_create(L, P->n, 0, 0);
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
  tk_dvec_t *result = out ? out : tk_dvec_create(L, P->n, 0, 0);
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

static inline double tk_csr_pearson_distance(
  tk_ivec_t *expected_ids,
  tk_ivec_t *expected_neighbors,
  tk_dvec_t *expected_weights,
  int64_t exp_start,
  int64_t exp_end,
  tk_ivec_t *retrieved_ids,
  tk_pvec_t *bin_ranks,
  tk_dumap_t *rank_buffer_b
) {
  uint64_t m = (uint64_t)(exp_end - exp_start);
  if (m == 0 || !bin_ranks || bin_ranks->n == 0)
    return 0.0;

  tk_dumap_clear(rank_buffer_b);
  int kha;
  for (uint64_t i = 0; i < bin_ranks->n; i++) {
    int64_t neighbor_idx = bin_ranks->a[i].i;
    int64_t neighbor_id = retrieved_ids->a[neighbor_idx];
    double hamming = (double)bin_ranks->a[i].p;
    uint32_t khi = tk_dumap_put(rank_buffer_b, neighbor_id, &kha);
    tk_dumap_setval(rank_buffer_b, khi, hamming);
  }
  double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0, sum_y2 = 0.0;
  uint64_t n = 0;
  for (int64_t j = exp_start; j < exp_end; j++) {
    int64_t neighbor_idx = expected_neighbors->a[j];
    int64_t neighbor_id = expected_ids->a[neighbor_idx];
    uint32_t khi = tk_dumap_get(rank_buffer_b, neighbor_id);
    if (khi != tk_dumap_end(rank_buffer_b)) {
      double x = expected_weights->a[j];
      double y = tk_dumap_val(rank_buffer_b, khi);
      sum_x += x;
      sum_y += y;
      sum_xy += x * y;
      sum_x2 += x * x;
      sum_y2 += y * y;
      n++;
    }
  }
  if (n < 2)
    return 0.0;
  double n_d = (double)n;
  double numerator = n_d * sum_xy - sum_x * sum_y;
  double denom_x = n_d * sum_x2 - sum_x * sum_x;
  double denom_y = n_d * sum_y2 - sum_y * sum_y;
  if (denom_x < 1e-10 || denom_y < 1e-10)
    return 0.0;
  return numerator / sqrt(denom_x * denom_y);
}

static inline double tk_csr_spearman_distance(
  tk_ivec_t *expected_ids,
  tk_ivec_t *expected_neighbors,
  tk_dvec_t *expected_weights,
  int64_t exp_start,
  int64_t exp_end,
  tk_ivec_t *retrieved_ids,
  tk_pvec_t *sorted_bin_ranks,
  tk_rvec_t *weight_ranks_buffer,
  tk_dumap_t *weight_rank_map
) {
  uint64_t m = (uint64_t)(exp_end - exp_start);
  if (m == 0 || !sorted_bin_ranks || sorted_bin_ranks->n == 0)
    return 0.0;
  if (sorted_bin_ranks->n < 2)
    return 0.0;
  if (tk_rvec_ensure(weight_ranks_buffer, m) != 0)
    return 0.0;
  weight_ranks_buffer->n = m;
  uint64_t idx = 0;
  for (int64_t j = exp_start; j < exp_end; j++) {
    int64_t neighbor_idx = expected_neighbors->a[j];
    int64_t neighbor_id = expected_ids->a[neighbor_idx];
    double weight = expected_weights->a[j];
    weight_ranks_buffer->a[idx++] = (tk_rank_t) { neighbor_id, weight };
  }
  tk_rvec_desc(weight_ranks_buffer, 0, weight_ranks_buffer->n);
  tk_dumap_clear(weight_rank_map);
  for (uint64_t i = 0; i < weight_ranks_buffer->n; ) {
    double val = weight_ranks_buffer->a[i].d;
    uint64_t j = i;
    while (j < weight_ranks_buffer->n && weight_ranks_buffer->a[j].d == val)
      j++;
    double avg_rank = (double)(i + j - 1) / 2.0;
    for (uint64_t k = i; k < j; k++) {
      int64_t neighbor_id = weight_ranks_buffer->a[k].i;
      int kha;
      uint32_t khi = tk_dumap_put(weight_rank_map, neighbor_id, &kha);
      tk_dumap_setval(weight_rank_map, khi, avg_rank);
    }
    i = j;
  }
  double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0, sum_y2 = 0.0;
  uint64_t n = 0;
  for (uint64_t i = 0; i < sorted_bin_ranks->n; i++) {
    int64_t neighbor_idx = sorted_bin_ranks->a[i].i;
    int64_t neighbor_id = retrieved_ids->a[neighbor_idx];
    uint32_t khi = tk_dumap_get(weight_rank_map, neighbor_id);
    if (khi != tk_dumap_end(weight_rank_map)) {
      double x_rank = tk_dumap_val(weight_rank_map, khi);
      double y_rank = (double)i;
      sum_x += x_rank;
      sum_y += y_rank;
      sum_xy += x_rank * y_rank;
      sum_x2 += x_rank * x_rank;
      sum_y2 += y_rank * y_rank;
      n++;
    }
  }
  if (n < 2)
    return 0.0;
  double n_d = (double)n;
  double numerator = n_d * sum_xy - sum_x * sum_y;
  double denom_x = n_d * sum_x2 - sum_x * sum_x;
  double denom_y = n_d * sum_y2 - sum_y * sum_y;
  if (denom_x < 1e-10 || denom_y < 1e-10)
    return 0.0;
  return numerator / sqrt(denom_x * denom_y);
}

static inline double tk_csr_mean(
  tk_ivec_t *neighbors_a,
  tk_dvec_t *weights_a,
  int64_t start_a,
  int64_t end_a,
  tk_iuset_t *group_1
) {
  double sum = 0.0;
  uint64_t count = 0;
  for (int64_t j = start_a; j < end_a; j++) {
    int64_t neighbor = neighbors_a->a[j];
    if (tk_iuset_get(group_1, neighbor) != tk_iuset_end(group_1)) {
      sum += weights_a->a[j];
      count++;
    }
  }
  return count > 0 ? sum / count : 0.0;
}

static inline double tk_csr_min(
  tk_ivec_t *neighbors_a,
  tk_dvec_t *weights_a,
  int64_t start_a,
  int64_t end_a,
  tk_iuset_t *group_1
) {
  double min_val = 1.0;
  bool found = false;
  for (int64_t j = start_a; j < end_a; j++) {
    int64_t neighbor = neighbors_a->a[j];
    if (tk_iuset_get(group_1, neighbor) != tk_iuset_end(group_1)) {
      min_val = fmin(min_val, weights_a->a[j]);
      found = true;
    }
  }
  return found ? min_val : 0.0;
}

static inline double tk_csr_pearson(
  tk_ivec_t *neighbors_a,
  tk_dvec_t *weights_a,
  int64_t start_a,
  int64_t end_a,
  tk_iuset_t *group_1
) {
  uint64_t n_a = (uint64_t)(end_a - start_a);
  if (n_a < 2)
    return 0.0;
  double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0, sum_y2 = 0.0;
  for (int64_t j = start_a; j < end_a; j++) {
    int64_t neighbor = neighbors_a->a[j];
    double x = tk_iuset_get(group_1, neighbor) != tk_iuset_end(group_1) ? 1.0 : 0.0;
    double y = weights_a->a[j];
    sum_x += x;
    sum_y += y;
    sum_xy += x * y;
    sum_x2 += x * x;
    sum_y2 += y * y;
  }
  double numerator = n_a * sum_xy - sum_x * sum_y;
  double denom_x = n_a * sum_x2 - sum_x * sum_x;
  double denom_y = n_a * sum_y2 - sum_y * sum_y;

  if (denom_y < 1e-10)
    return 0.0;

  if (denom_x < 1e-10) {

    double prop_in_group = sum_x / n_a;
    return prop_in_group > 0.99 ? 1.0 : 0.0;
  }

  return numerator / sqrt(denom_x * denom_y);
}

static inline double tk_csr_variance_ratio(
  tk_ivec_t *neighbors_a,
  tk_dvec_t *weights_a,
  int64_t start_a,
  int64_t end_a,
  tk_iuset_t *group_1,
  tk_dumap_t *rank_buffer
) {
  uint64_t n_a = (uint64_t)(end_a - start_a);
  if (n_a == 0 || !group_1)
    return 0.0;
  tk_dumap_clear(rank_buffer);
  int kha;
  for (int64_t j = start_a; j < end_a; j++) {
    double rank = (double)(j - start_a);
    uint64_t count = 1;
    double weight = weights_a->a[j];
    while (j + 1 < end_a && weights_a->a[j + 1] == weight) {
      count++;
      j++;
    }
    double average_rank = (rank + (rank + count - 1)) / 2.0 + 1.0;
    for (uint64_t k = 0; k < count; k++) {
      int64_t neighbor = neighbors_a->a[(uint64_t)j - k];
      uint32_t khi = tk_dumap_put(rank_buffer, neighbor, &kha);
      tk_dumap_setval(rank_buffer, khi, average_rank);
    }
  }
  double rank_sum_0 = 0.0, rank_sum_1 = 0.0;
  uint64_t count_0 = 0, count_1 = 0;
  int64_t neighbor;
  double rank;
  tk_umap_foreach(rank_buffer, neighbor, rank, ({
    if (tk_iuset_get(group_1, neighbor) != tk_iuset_end(group_1)) {
      rank_sum_1 += rank;
      count_1++;
    } else {
      rank_sum_0 += rank;
      count_0++;
    }
  }))
  if (count_1 == 0)
    return 0.0;

  if (count_0 == 0)
    return 1.0;

  double mean_rank_0 = rank_sum_0 / count_0;
  double mean_rank_1 = rank_sum_1 / count_1;
  double overall_mean = (rank_sum_0 + rank_sum_1) / (count_0 + count_1);
  double ss_between = count_0 * (mean_rank_0 - overall_mean) * (mean_rank_0 - overall_mean)
                    + count_1 * (mean_rank_1 - overall_mean) * (mean_rank_1 - overall_mean);
  double ss_within = 0.0;
  tk_umap_foreach(rank_buffer, neighbor, rank, ({
    if (tk_iuset_get(group_1, neighbor) != tk_iuset_end(group_1)) {
      double diff = rank - mean_rank_1;
      ss_within += diff * diff;
    } else {
      double diff = rank - mean_rank_0;
      ss_within += diff * diff;
    }
  }))
  double ss_total = ss_between + ss_within;
  if (ss_total <= 1e-12) {
    return 0.0;
  }
  return ss_between / ss_total;
}

static inline double tk_csr_ndcg_distance(
  tk_ivec_t *expected_ids,
  tk_ivec_t *expected_neighbors,
  tk_dvec_t *expected_weights,
  int64_t exp_start,
  int64_t exp_end,
  tk_ivec_t *retrieved_ids,
  tk_pvec_t *sorted_bin_ranks,
  tk_dumap_t *weight_map
) {
  uint64_t m = (uint64_t)(exp_end - exp_start);
  if (m == 0 || !sorted_bin_ranks || sorted_bin_ranks->n == 0)
    return 0.0;

  tk_dumap_clear(weight_map);
  int kha;
  for (int64_t j = exp_start; j < exp_end; j++) {
    int64_t neighbor_idx = expected_neighbors->a[j];
    int64_t neighbor_id = expected_ids->a[neighbor_idx];
    double weight = expected_weights->a[j];
    uint32_t khi = tk_dumap_put(weight_map, neighbor_id, &kha);
    tk_dumap_setval(weight_map, khi, weight);
  }

  double dcg = 0.0;
  uint64_t n = sorted_bin_ranks->n;
  uint64_t i = 0;
  while (i < n) {
    int64_t current_hamming = sorted_bin_ranks->a[i].p;
    uint64_t tie_start = i;
    while (i < n && sorted_bin_ranks->a[i].p == current_hamming)
      i++;
    uint64_t tie_end = i;
    uint64_t tie_count = tie_end - tie_start;
    double discount_sum = 0.0;
    for (uint64_t pos = tie_start; pos < tie_end; pos++)
      discount_sum += log2((double)(pos + 2));
    double avg_discount = discount_sum / (double)tie_count;
    for (uint64_t j = tie_start; j < tie_end; j++) {
      int64_t neighbor_idx = sorted_bin_ranks->a[j].i;
      int64_t neighbor_id = retrieved_ids->a[neighbor_idx];
      uint32_t khi = tk_dumap_get(weight_map, neighbor_id);
      if (khi != tk_dumap_end(weight_map)) {
        double relevance = tk_dumap_val(weight_map, khi);
        dcg += relevance / avg_discount;
      }
    }
  }

  uint64_t k = n;
  tk_dvec_t *sorted_weights = tk_dvec_create(NULL, m, 0, 0);
  if (!sorted_weights)
    return 0.0;
  for (i = 0; i < m; i++)
    sorted_weights->a[i] = expected_weights->a[(uint64_t)exp_start + i];
  tk_dvec_desc(sorted_weights, 0, m);

  double idcg = 0.0;
  uint64_t idcg_count = (m < k) ? m : k;
  for (i = 0; i < idcg_count; i++) {
    double relevance = sorted_weights->a[i];
    double discount = log2((double)(i + 2));
    idcg += relevance / discount;
  }

  tk_dvec_destroy(sorted_weights);

  if (idcg < 1e-10)
    return 0.0;

  return dcg / idcg;
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
