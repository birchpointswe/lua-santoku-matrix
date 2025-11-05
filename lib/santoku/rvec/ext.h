#ifndef TK_RVEC_EXT_H
#define TK_RVEC_EXT_H

#include <omp.h>
#include <santoku/dumap.h>
#include <santoku/iumap.h>
#include <santoku/pvec.h>
#include <santoku/evec.h>
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

static inline double tk_csr_spearman(
  tk_ivec_t *neighbors_a,
  tk_dvec_t *weights_a,
  int64_t start_a,
  int64_t end_a,
  tk_pvec_t *bin_ranks,
  uint64_t max_hamming,
  tk_ivec_t *count_buffer,
  tk_dvec_t *avgrank_buffer,
  tk_dumap_t *rank_buffer_b
) {
  uint64_t m = (uint64_t)(end_a - start_a);
  if (m == 0 || !bin_ranks || bin_ranks->n == 0)
    return 0.0;
  for (uint64_t h = 0; h <= max_hamming && h < count_buffer->m; h++)
    count_buffer->a[h] = 0;
  for (uint64_t i = 0; i < bin_ranks->n; i++) {
    uint64_t hamming = (uint64_t)bin_ranks->a[i].p;
    if (hamming <= max_hamming && hamming < count_buffer->m) {
      count_buffer->a[hamming]++;
    }
  }
  uint64_t cumulative = 0;
  for (uint64_t h = 0; h <= max_hamming && h < avgrank_buffer->m; h++) {
    if (count_buffer->a[h] > 0) {
      uint64_t count = (uint64_t)count_buffer->a[h];
      double start_rank = (double)cumulative;
      double end_rank = (double)(cumulative + count - 1);
      avgrank_buffer->a[h] = (start_rank + end_rank) / 2.0;
      cumulative += count;
    }
  }
  tk_dumap_clear(rank_buffer_b);
  int kha;
  for (uint64_t i = 0; i < bin_ranks->n; i++) {
    int64_t neighbor_pos = bin_ranks->a[i].i;
    uint64_t hamming = (uint64_t)bin_ranks->a[i].p;
    if (hamming <= max_hamming && hamming < avgrank_buffer->m) {
      uint32_t khi = tk_dumap_put(rank_buffer_b, neighbor_pos, &kha);
      tk_dumap_setval(rank_buffer_b, khi, avgrank_buffer->a[hamming]);
    }
  }
  double sum_squared_diff = 0.0;
  uint64_t n = 0;
  int64_t j = start_a;
  while (j < end_a) {
    double weight = weights_a->a[j];
    int64_t tie_start = j;
    while (j + 1 < end_a && weights_a->a[j + 1] == weight)
      j++;
    int64_t tie_end = j;
    double avg_adj_rank = ((double)(tie_start - start_a) + (double)(tie_end - start_a)) / 2.0;
    for (int64_t t = tie_start; t <= tie_end; t++) {
      int64_t neighbor_pos = neighbors_a->a[t];
      uint32_t khi = tk_dumap_get(rank_buffer_b, neighbor_pos);
      if (khi != tk_dumap_end(rank_buffer_b)) {
        double ham_rank = tk_dumap_val(rank_buffer_b, khi);
        double diff = avg_adj_rank - ham_rank;
        sum_squared_diff += diff * diff;
        n++;
      }
    }
    j++;
  }
  if (n <= 1)
    return 1.0;
  double n_d = (double)n;
  return 1.0 - (6.0 * sum_squared_diff) / (n_d * (n_d * n_d - 1.0));
}

static inline double tk_csr_position(
  tk_ivec_t *neighbors_a,
  int64_t start_a,
  int64_t end_a,
  tk_pvec_t *bin_ranks,
  tk_dumap_t *pos_buffer
) {
  uint64_t m = (uint64_t)(end_a - start_a);
  if (m <= 1 || !bin_ranks || bin_ranks->n == 0)
    return 1.0;
  tk_dumap_clear(pos_buffer);
  int kha;
  for (int64_t j = start_a; j < end_a; j++) {
    int64_t neighbor_pos = neighbors_a->a[j];
    double position = (double)(j - start_a);
    uint32_t khi = tk_dumap_put(pos_buffer, neighbor_pos, &kha);
    tk_dumap_setval(pos_buffer, khi, position);
  }
  double sum_y = 0.0, sum_xy = 0.0, sum_y2 = 0.0;
  uint64_t n = 0;
  for (uint64_t i = 0; i < bin_ranks->n; i++) {
    int64_t neighbor_pos = bin_ranks->a[i].i;
    uint32_t khi = tk_dumap_get(pos_buffer, neighbor_pos);
    if (khi != tk_dumap_end(pos_buffer)) {
      double x = tk_dumap_val(pos_buffer, khi);
      double y = (double)bin_ranks->a[i].p;
      sum_y += y;
      sum_xy += x * y;
      sum_y2 += y * y;
      n++;
    }
  }
  if (n <= 1) return 1.0;
  double n_d = (double)n;
  double sum_x = n_d * (n_d - 1.0) / 2.0;
  double sum_x2 = n_d * (n_d - 1.0) * (2.0 * n_d - 1.0) / 6.0;
  double mean_x = sum_x / n_d;
  double mean_y = sum_y / n_d;
  double cov = (sum_xy / n_d) - (mean_x * mean_y);
  double var_x = (sum_x2 / n_d) - (mean_x * mean_x);
  double var_y = (sum_y2 / n_d) - (mean_y * mean_y);
  if (var_x <= 0.0 || var_y <= 0.0) return 0.0;
  return cov / sqrt(var_x * var_y);
}

static inline double tk_csr_biserial(
  tk_ivec_t *neighbors_a,
  tk_dvec_t *weights_a,
  int64_t start_a,
  int64_t end_a,
  tk_iuset_t *group_1
) {
  uint64_t n_a = (uint64_t)(end_a - start_a);
  if (n_a == 0 || !group_1)
    return 0.0;
  double rank_sum_1 = 0.0;
  uint64_t count_0 = 0, count_1 = 0;
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
      if (tk_iuset_get(group_1, neighbor) != tk_iuset_end(group_1)) {
        rank_sum_1 += average_rank;
        count_1++;
      } else {
        count_0++;
      }
    }
  }
  if (n_a == 0)
    return 0.0;
  if (count_0 == 0 && count_1 > 0)
    return 1.0;
  if (count_1 == 0 || count_0 == 0)
    return 0.0;
  double U1 = rank_sum_1 - (count_1 * (count_1 + 1)) / 2.0;
  return 1.0 - (2.0 * U1) / (count_0 * count_1);
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
  if (count_0 == 0 || count_1 == 0)
    return 0.0;
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

#endif
