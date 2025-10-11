#ifndef TK_RVEC_EXT_H
#define TK_RVEC_EXT_H

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
  tk_ivec_t *I = out ? out : tk_ivec_create(L, P->n, 0, 0);
  if (out)
    tk_ivec_ensure(I, P->n);
  for (uint64_t i = 0; i < P->n; i ++)
    I->a[i] = P->a[i].i;
  if (out)
    I->n = P->n;
  return I;
}

static inline tk_dvec_t *tk_rvec_values (
  lua_State *L,
  tk_rvec_t *P,
  tk_dvec_t *out
) {
  tk_dvec_t *I = out ? out : tk_dvec_create(L, P->n, 0, 0);
  if (out)
    tk_dvec_ensure(I, P->n);
  for (uint64_t i = 0; i < P->n; i ++)
    I->a[i] = P->a[i].d;
  if (out)
    I->n = P->n;
  return I;
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

static inline void tk_kendall_merge_sort_weighted(
  tk_evec_t *items,
  tk_edge_t *temp,
  uint64_t left,
  uint64_t right,
  double *discordant
) {
  if (right <= left + 1)
    return;

  uint64_t mid = left + (right - left) / 2;

  tk_kendall_merge_sort_weighted(items, temp, left, mid, discordant);
  tk_kendall_merge_sort_weighted(items, temp, mid, right, discordant);

  uint64_t i = left;
  uint64_t j = mid;
  uint64_t k = left;

  while (i < mid && j < right) {
    if (items->a[i].v <= items->a[j].v) {
      temp[k++] = items->a[i++];
    } else {
      for (uint64_t ii = i; ii < mid; ii++) {
        double w_ij = sqrt(fabs(items->a[ii].w * items->a[j].w));
        *discordant += w_ij;
      }
      temp[k++] = items->a[j++];
    }
  }

  while (i < mid) temp[k++] = items->a[i++];
  while (j < right) temp[k++] = items->a[j++];

  for (uint64_t idx = left; idx < right; idx++)
    items->a[idx] = temp[idx];
}

static inline void tk_kendall_merge_sort_unweighted(
  tk_evec_t *items,
  tk_edge_t *temp,
  uint64_t left,
  uint64_t right,
  uint64_t *discordant
) {
  if (right <= left + 1)
    return;

  uint64_t mid = left + (right - left) / 2;

  tk_kendall_merge_sort_unweighted(items, temp, left, mid, discordant);
  tk_kendall_merge_sort_unweighted(items, temp, mid, right, discordant);

  uint64_t i = left;
  uint64_t j = mid;
  uint64_t k = left;

  while (i < mid && j < right) {
    if (items->a[i].v <= items->a[j].v) {
      temp[k++] = items->a[i++];
    } else {
      *discordant += (mid - i);
      temp[k++] = items->a[j++];
    }
  }

  while (i < mid) temp[k++] = items->a[i++];
  while (j < right) temp[k++] = items->a[j++];

  for (uint64_t idx = left; idx < right; idx++)
    items->a[idx] = temp[idx];
}

static inline double tk_rvec_kendall(
  tk_rvec_t *ranks_a,
  tk_rvec_t *ranks_b
) {
  if (!ranks_a || !ranks_b || ranks_a->n == 0 || ranks_b->n == 0)
    return 0.0;

  tk_iumap_t *idx_b = tk_iumap_create(0, ranks_b->n);
  if (!idx_b)
    return 0.0;

  for (uint64_t i = 0; i < ranks_b->n; i++) {
    int kha;
    uint32_t khi = tk_iumap_put(idx_b, ranks_b->a[i].i, &kha);
    tk_iumap_setval(idx_b, khi, (int64_t)i);
  }

  tk_evec_t *items = tk_evec_create(0, ranks_a->n, 0, 0);
  if (!items) {
    tk_iumap_destroy(idx_b);
    return 0.0;
  }

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      tk_edge_t item = {.u = (int64_t)i, .v = tk_iumap_val(idx_b, khi), .w = 1.0};
      if (tk_evec_push(items, item) != 0) {
        tk_evec_destroy(items);
        tk_iumap_destroy(idx_b);
        return 0.0;
      }
    }
  }

  tk_iumap_destroy(idx_b);

  if (items->n <= 1) {
    tk_evec_destroy(items);
    return 1.0;
  }

  uint64_t n = items->n;
  uint64_t total_pairs = (n * (n - 1)) / 2;

  tk_edge_t *temp = malloc(n * sizeof(tk_edge_t));
  if (!temp) {
    tk_evec_destroy(items);
    return 0.0;
  }

  uint64_t discordant = 0;
  tk_kendall_merge_sort_unweighted(items, temp, 0, n, &discordant);

  tk_evec_destroy(items);
  free(temp);

  uint64_t concordant = total_pairs - discordant;
  return (double)(concordant - discordant) / (double)total_pairs;
}

static inline double tk_rvec_kendall_weighted(
  tk_rvec_t *ranks_a,
  tk_rvec_t *ranks_b
) {
  if (!ranks_a || !ranks_b || ranks_a->n == 0 || ranks_b->n == 0)
    return 0.0;

  tk_iumap_t *idx_b = tk_iumap_create(0, ranks_b->n);
  if (!idx_b)
    return 0.0;

  for (uint64_t i = 0; i < ranks_b->n; i++) {
    int kha;
    uint32_t khi = tk_iumap_put(idx_b, ranks_b->a[i].i, &kha);
    tk_iumap_setval(idx_b, khi, (int64_t)i);
  }

  tk_evec_t *items = tk_evec_create(0, ranks_a->n, 0, 0);
  if (!items) {
    tk_iumap_destroy(idx_b);
    return 0.0;
  }

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      tk_edge_t item = {.u = (int64_t)i, .v = tk_iumap_val(idx_b, khi), .w = ranks_a->a[i].d};
      if (tk_evec_push(items, item) != 0) {
        tk_evec_destroy(items);
        tk_iumap_destroy(idx_b);
        return 0.0;
      }
    }
  }

  tk_iumap_destroy(idx_b);

  if (items->n <= 1) {
    tk_evec_destroy(items);
    return 1.0;
  }

  double total_weight = 0.0;
  for (uint64_t i = 0; i < items->n; i++) {
    for (uint64_t j = i + 1; j < items->n; j++) {
      double w_ij = sqrt(fabs(items->a[i].w * items->a[j].w));
      total_weight += w_ij;
    }
  }

  if (total_weight <= 0.0) {
    tk_evec_destroy(items);
    return 0.0;
  }

  tk_edge_t *temp = malloc(items->n * sizeof(tk_edge_t));
  if (!temp) {
    tk_evec_destroy(items);
    return 0.0;
  }

  double discordant = 0.0;
  tk_kendall_merge_sort_weighted(items, temp, 0, items->n, &discordant);

  tk_evec_destroy(items);
  free(temp);

  double concordant = total_weight - discordant;
  return (concordant - discordant) / total_weight;
}

static inline double tk_rvec_spearman_weighted(
  tk_rvec_t *ranks_a,
  tk_rvec_t *ranks_b
) {
  if (!ranks_a || !ranks_b || ranks_a->n == 0 || ranks_b->n == 0)
    return 0.0;

  tk_iumap_t *idx_b = tk_iumap_create(0, ranks_b->n);
  if (!idx_b)
    return 0.0;

  for (uint64_t i = 0; i < ranks_b->n; i++) {
    int kha;
    uint32_t khi = tk_iumap_put(idx_b, ranks_b->a[i].i, &kha);
    tk_iumap_setval(idx_b, khi, (int64_t)i);
  }

  double sum_w = 0.0;
  double sum_w_ra = 0.0;
  double sum_w_rb = 0.0;
  uint64_t n = 0;

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      double rank_a = (double)i;
      double rank_b = (double)tk_iumap_val(idx_b, khi);
      double weight = ranks_a->a[i].d;
      sum_w += weight;
      sum_w_ra += weight * rank_a;
      sum_w_rb += weight * rank_b;
      n++;
    }
  }

  if (n <= 1 || sum_w <= 0.0) {
    tk_iumap_destroy(idx_b);
    return 1.0;
  }

  double mean_a = sum_w_ra / sum_w;
  double mean_b = sum_w_rb / sum_w;

  double cov = 0.0;
  double var_a = 0.0;
  double var_b = 0.0;

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      double rank_a = (double)i;
      double rank_b = (double)tk_iumap_val(idx_b, khi);
      double weight = ranks_a->a[i].d;
      double da = rank_a - mean_a;
      double db = rank_b - mean_b;
      cov += weight * da * db;
      var_a += weight * da * da;
      var_b += weight * db * db;
    }
  }

  tk_iumap_destroy(idx_b);

  cov /= sum_w;
  var_a /= sum_w;
  var_b /= sum_w;

  if (var_a <= 0.0 || var_b <= 0.0)
    return 1.0;

  return cov / sqrt(var_a * var_b);
}

static inline double tk_rvec_kendall_pvec(
  tk_rvec_t *ranks_a,
  tk_pvec_t *ranks_b
) {
  if (!ranks_a || !ranks_b || ranks_a->n == 0 || ranks_b->n == 0)
    return 0.0;

  tk_iumap_t *idx_b = tk_iumap_create(0, ranks_b->n);
  if (!idx_b)
    return 0.0;

  for (uint64_t i = 0; i < ranks_b->n; i++) {
    int kha;
    uint32_t khi = tk_iumap_put(idx_b, ranks_b->a[i].i, &kha);
    tk_iumap_setval(idx_b, khi, (int64_t)i);
  }

  tk_evec_t *items = tk_evec_create(0, ranks_a->n, 0, 0);
  if (!items) {
    tk_iumap_destroy(idx_b);
    return 0.0;
  }

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      tk_edge_t item = {
        .u = (int64_t)i,
        .v = tk_iumap_val(idx_b, khi),
        .w = 1.0
      };
      if (tk_evec_push(items, item) != 0) {
        tk_evec_destroy(items);
        tk_iumap_destroy(idx_b);
        return 0.0;
      }
    }
  }

  tk_iumap_destroy(idx_b);

  if (items->n <= 1) {
    tk_evec_destroy(items);
    return 1.0;
  }

  uint64_t n = items->n;
  uint64_t total_pairs = (n * (n - 1)) / 2;

  tk_edge_t *temp = malloc(n * sizeof(tk_edge_t));
  if (!temp) {
    tk_evec_destroy(items);
    return 0.0;
  }

  uint64_t discordant = 0;
  tk_kendall_merge_sort_unweighted(items, temp, 0, n, &discordant);

  tk_evec_destroy(items);
  free(temp);

  uint64_t concordant = total_pairs - discordant;
  return (double)(concordant - discordant) / (double)total_pairs;
}

static inline double tk_rvec_kendall_weighted_pvec(
  tk_rvec_t *ranks_a,
  tk_pvec_t *ranks_b
) {
  if (!ranks_a || !ranks_b || ranks_a->n == 0 || ranks_b->n == 0)
    return 0.0;

  tk_iumap_t *idx_b = tk_iumap_create(0, ranks_b->n);
  if (!idx_b)
    return 0.0;

  for (uint64_t i = 0; i < ranks_b->n; i++) {
    int kha;
    uint32_t khi = tk_iumap_put(idx_b, ranks_b->a[i].i, &kha);
    tk_iumap_setval(idx_b, khi, (int64_t)i);
  }

  tk_evec_t *items = tk_evec_create(0, ranks_a->n, 0, 0);
  if (!items) {
    tk_iumap_destroy(idx_b);
    return 0.0;
  }

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      tk_edge_t item = {
        .u = (int64_t)i,
        .v = tk_iumap_val(idx_b, khi),
        .w = ranks_a->a[i].d
      };
      if (tk_evec_push(items, item) != 0) {
        tk_evec_destroy(items);
        tk_iumap_destroy(idx_b);
        return 0.0;
      }
    }
  }

  tk_iumap_destroy(idx_b);

  if (items->n <= 1) {
    tk_evec_destroy(items);
    return 1.0;
  }

  double total_weight = 0.0;
  for (uint64_t i = 0; i < items->n; i++) {
    for (uint64_t j = i + 1; j < items->n; j++) {
      double w_ij = sqrt(fabs(items->a[i].w * items->a[j].w));
      total_weight += w_ij;
    }
  }

  if (total_weight <= 0.0) {
    tk_evec_destroy(items);
    return 0.0;
  }

  tk_edge_t *temp = malloc(items->n * sizeof(tk_edge_t));
  if (!temp) {
    tk_evec_destroy(items);
    return 0.0;
  }

  double discordant = 0.0;
  tk_kendall_merge_sort_weighted(items, temp, 0, items->n, &discordant);

  tk_evec_destroy(items);
  free(temp);

  double concordant = total_weight - discordant;
  return (concordant - discordant) / total_weight;
}

static inline double tk_rvec_spearman_weighted_pvec(
  tk_rvec_t *ranks_a,
  tk_pvec_t *ranks_b
) {
  if (!ranks_a || !ranks_b || ranks_a->n == 0 || ranks_b->n == 0)
    return 0.0;

  tk_iumap_t *idx_b = tk_iumap_create(0, ranks_b->n);
  if (!idx_b)
    return 0.0;

  for (uint64_t i = 0; i < ranks_b->n; i++) {
    int kha;
    uint32_t khi = tk_iumap_put(idx_b, ranks_b->a[i].i, &kha);
    tk_iumap_setval(idx_b, khi, (int64_t)i);
  }

  double sum_w = 0.0;
  double sum_w_ra = 0.0;
  double sum_w_rb = 0.0;
  uint64_t n = 0;

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      double rank_a = (double)i;
      double rank_b = (double)tk_iumap_val(idx_b, khi);
      double weight = ranks_a->a[i].d;
      sum_w += weight;
      sum_w_ra += weight * rank_a;
      sum_w_rb += weight * rank_b;
      n++;
    }
  }

  if (n <= 1 || sum_w <= 0.0) {
    tk_iumap_destroy(idx_b);
    return 1.0;
  }

  double mean_a = sum_w_ra / sum_w;
  double mean_b = sum_w_rb / sum_w;

  double cov = 0.0;
  double var_a = 0.0;
  double var_b = 0.0;

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      double rank_a = (double)i;
      double rank_b = (double)tk_iumap_val(idx_b, khi);
      double weight = ranks_a->a[i].d;
      double da = rank_a - mean_a;
      double db = rank_b - mean_b;
      cov += weight * da * db;
      var_a += weight * da * da;
      var_b += weight * db * db;
    }
  }

  tk_iumap_destroy(idx_b);

  cov /= sum_w;
  var_a /= sum_w;
  var_b /= sum_w;

  if (var_a <= 0.0 || var_b <= 0.0)
    return 1.0;

  return cov / sqrt(var_a * var_b);
}

static inline double tk_pvec_kendall(
  tk_pvec_t *ranks_a,
  tk_pvec_t *ranks_b
) {
  if (!ranks_a || !ranks_b || ranks_a->n == 0 || ranks_b->n == 0)
    return 0.0;

  tk_iumap_t *idx_b = tk_iumap_create(0, ranks_b->n);
  if (!idx_b)
    return 0.0;

  for (uint64_t i = 0; i < ranks_b->n; i++) {
    int kha;
    uint32_t khi = tk_iumap_put(idx_b, ranks_b->a[i].i, &kha);
    tk_iumap_setval(idx_b, khi, (int64_t)i);
  }

  tk_evec_t *items = tk_evec_create(0, ranks_a->n, 0, 0);
  if (!items) {
    tk_iumap_destroy(idx_b);
    return 0.0;
  }

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      tk_edge_t item = {.u = (int64_t)i, .v = tk_iumap_val(idx_b, khi), .w = 1.0};
      if (tk_evec_push(items, item) != 0) {
        tk_evec_destroy(items);
        tk_iumap_destroy(idx_b);
        return 0.0;
      }
    }
  }

  tk_iumap_destroy(idx_b);

  if (items->n <= 1) {
    tk_evec_destroy(items);
    return 1.0;
  }

  uint64_t n = items->n;
  uint64_t total_pairs = (n * (n - 1)) / 2;
  tk_edge_t *temp = malloc(n * sizeof(tk_edge_t));
  if (!temp) {
    tk_evec_destroy(items);
    return 0.0;
  }

  uint64_t discordant = 0;
  tk_kendall_merge_sort_unweighted(items, temp, 0, n, &discordant);

  tk_evec_destroy(items);
  free(temp);

  uint64_t concordant = total_pairs - discordant;
  return (double)(concordant - discordant) / (double)total_pairs;
}

static inline double tk_pvec_kendall_weighted(
  tk_pvec_t *ranks_a,
  tk_pvec_t *ranks_b
) {
  if (!ranks_a || !ranks_b || ranks_a->n == 0 || ranks_b->n == 0)
    return 0.0;

  tk_iumap_t *idx_b = tk_iumap_create(0, ranks_b->n);
  if (!idx_b)
    return 0.0;

  for (uint64_t i = 0; i < ranks_b->n; i++) {
    int kha;
    uint32_t khi = tk_iumap_put(idx_b, ranks_b->a[i].i, &kha);
    tk_iumap_setval(idx_b, khi, (int64_t)i);
  }

  tk_evec_t *items = tk_evec_create(0, ranks_a->n, 0, 0);
  if (!items) {
    tk_iumap_destroy(idx_b);
    return 0.0;
  }

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      tk_edge_t item = {.u = (int64_t)i, .v = tk_iumap_val(idx_b, khi), .w = (double)ranks_a->a[i].p};
      if (tk_evec_push(items, item) != 0) {
        tk_evec_destroy(items);
        tk_iumap_destroy(idx_b);
        return 0.0;
      }
    }
  }

  tk_iumap_destroy(idx_b);

  if (items->n <= 1) {
    tk_evec_destroy(items);
    return 1.0;
  }

  double total_weight = 0.0;
  for (uint64_t i = 0; i < items->n; i++) {
    for (uint64_t j = i + 1; j < items->n; j++) {
      double w_ij = sqrt(fabs(items->a[i].w * items->a[j].w));
      total_weight += w_ij;
    }
  }

  if (total_weight <= 0.0) {
    tk_evec_destroy(items);
    return 0.0;
  }

  tk_edge_t *temp = malloc(items->n * sizeof(tk_edge_t));
  if (!temp) {
    tk_evec_destroy(items);
    return 0.0;
  }

  double discordant = 0.0;
  tk_kendall_merge_sort_weighted(items, temp, 0, items->n, &discordant);

  tk_evec_destroy(items);
  free(temp);

  double concordant = total_weight - discordant;
  return (concordant - discordant) / total_weight;
}

static inline double tk_pvec_kendall_rvec(
  tk_pvec_t *ranks_a,
  tk_rvec_t *ranks_b
) {
  if (!ranks_a || !ranks_b || ranks_a->n == 0 || ranks_b->n == 0)
    return 0.0;

  tk_iumap_t *idx_b = tk_iumap_create(0, ranks_b->n);
  if (!idx_b)
    return 0.0;

  for (uint64_t i = 0; i < ranks_b->n; i++) {
    int kha;
    uint32_t khi = tk_iumap_put(idx_b, ranks_b->a[i].i, &kha);
    tk_iumap_setval(idx_b, khi, (int64_t)i);
  }

  tk_evec_t *items = tk_evec_create(0, ranks_a->n, 0, 0);
  if (!items) {
    tk_iumap_destroy(idx_b);
    return 0.0;
  }

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      tk_edge_t item = {.u = (int64_t)i, .v = tk_iumap_val(idx_b, khi), .w = 1.0};
      if (tk_evec_push(items, item) != 0) {
        tk_evec_destroy(items);
        tk_iumap_destroy(idx_b);
        return 0.0;
      }
    }
  }

  tk_iumap_destroy(idx_b);

  if (items->n <= 1) {
    tk_evec_destroy(items);
    return 1.0;
  }

  uint64_t n = items->n;
  uint64_t total_pairs = (n * (n - 1)) / 2;
  tk_edge_t *temp = malloc(n * sizeof(tk_edge_t));
  if (!temp) {
    tk_evec_destroy(items);
    return 0.0;
  }

  uint64_t discordant = 0;
  tk_kendall_merge_sort_unweighted(items, temp, 0, n, &discordant);

  tk_evec_destroy(items);
  free(temp);

  uint64_t concordant = total_pairs - discordant;
  return (double)(concordant - discordant) / (double)total_pairs;
}

static inline double tk_pvec_kendall_weighted_rvec(
  tk_pvec_t *ranks_a,
  tk_rvec_t *ranks_b
) {
  if (!ranks_a || !ranks_b || ranks_a->n == 0 || ranks_b->n == 0)
    return 0.0;

  tk_iumap_t *idx_b = tk_iumap_create(0, ranks_b->n);
  if (!idx_b)
    return 0.0;

  for (uint64_t i = 0; i < ranks_b->n; i++) {
    int kha;
    uint32_t khi = tk_iumap_put(idx_b, ranks_b->a[i].i, &kha);
    tk_iumap_setval(idx_b, khi, (int64_t)i);
  }

  tk_evec_t *items = tk_evec_create(0, ranks_a->n, 0, 0);
  if (!items) {
    tk_iumap_destroy(idx_b);
    return 0.0;
  }

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      tk_edge_t item = {.u = (int64_t)i, .v = tk_iumap_val(idx_b, khi), .w = (double)ranks_a->a[i].p};
      if (tk_evec_push(items, item) != 0) {
        tk_evec_destroy(items);
        tk_iumap_destroy(idx_b);
        return 0.0;
      }
    }
  }

  tk_iumap_destroy(idx_b);

  if (items->n <= 1) {
    tk_evec_destroy(items);
    return 1.0;
  }

  double total_weight = 0.0;
  for (uint64_t i = 0; i < items->n; i++) {
    for (uint64_t j = i + 1; j < items->n; j++) {
      double w_ij = sqrt(fabs(items->a[i].w * items->a[j].w));
      total_weight += w_ij;
    }
  }

  if (total_weight <= 0.0) {
    tk_evec_destroy(items);
    return 0.0;
  }

  tk_edge_t *temp = malloc(items->n * sizeof(tk_edge_t));
  if (!temp) {
    tk_evec_destroy(items);
    return 0.0;
  }

  double discordant = 0.0;
  tk_kendall_merge_sort_weighted(items, temp, 0, items->n, &discordant);

  tk_evec_destroy(items);
  free(temp);

  double concordant = total_weight - discordant;
  return (concordant - discordant) / total_weight;
}

static inline double tk_pvec_spearman_weighted(
  tk_pvec_t *ranks_a,
  tk_pvec_t *ranks_b
) {
  if (!ranks_a || !ranks_b || ranks_a->n == 0 || ranks_b->n == 0)
    return 0.0;

  tk_iumap_t *idx_b = tk_iumap_create(0, ranks_b->n);
  if (!idx_b)
    return 0.0;

  for (uint64_t i = 0; i < ranks_b->n; i++) {
    int kha;
    uint32_t khi = tk_iumap_put(idx_b, ranks_b->a[i].i, &kha);
    tk_iumap_setval(idx_b, khi, (int64_t)i);
  }

  double sum_w = 0.0;
  double sum_w_ra = 0.0;
  double sum_w_rb = 0.0;
  uint64_t n = 0;

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      double rank_a = (double)i;
      double rank_b = (double)tk_iumap_val(idx_b, khi);
      double weight = (double)ranks_a->a[i].p;
      sum_w += weight;
      sum_w_ra += weight * rank_a;
      sum_w_rb += weight * rank_b;
      n++;
    }
  }

  if (n <= 1 || sum_w <= 0.0) {
    tk_iumap_destroy(idx_b);
    return 1.0;
  }

  double mean_a = sum_w_ra / sum_w;
  double mean_b = sum_w_rb / sum_w;

  double cov = 0.0;
  double var_a = 0.0;
  double var_b = 0.0;

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      double rank_a = (double)i;
      double rank_b = (double)tk_iumap_val(idx_b, khi);
      double weight = (double)ranks_a->a[i].p;
      double da = rank_a - mean_a;
      double db = rank_b - mean_b;
      cov += weight * da * db;
      var_a += weight * da * da;
      var_b += weight * db * db;
    }
  }

  tk_iumap_destroy(idx_b);

  cov /= sum_w;
  var_a /= sum_w;
  var_b /= sum_w;

  if (var_a <= 0.0 || var_b <= 0.0)
    return 1.0;

  return cov / sqrt(var_a * var_b);
}

static inline double tk_pvec_spearman_weighted_rvec(
  tk_pvec_t *ranks_a,
  tk_rvec_t *ranks_b
) {
  if (!ranks_a || !ranks_b || ranks_a->n == 0 || ranks_b->n == 0)
    return 0.0;

  tk_iumap_t *idx_b = tk_iumap_create(0, ranks_b->n);
  if (!idx_b)
    return 0.0;

  for (uint64_t i = 0; i < ranks_b->n; i++) {
    int kha;
    uint32_t khi = tk_iumap_put(idx_b, ranks_b->a[i].i, &kha);
    tk_iumap_setval(idx_b, khi, (int64_t)i);
  }

  double sum_w = 0.0;
  double sum_w_ra = 0.0;
  double sum_w_rb = 0.0;
  uint64_t n = 0;

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      double rank_a = (double)i;
      double rank_b = (double)tk_iumap_val(idx_b, khi);
      double weight = (double)ranks_a->a[i].p;
      sum_w += weight;
      sum_w_ra += weight * rank_a;
      sum_w_rb += weight * rank_b;
      n++;
    }
  }

  if (n <= 1 || sum_w <= 0.0) {
    tk_iumap_destroy(idx_b);
    return 1.0;
  }

  double mean_a = sum_w_ra / sum_w;
  double mean_b = sum_w_rb / sum_w;

  double cov = 0.0;
  double var_a = 0.0;
  double var_b = 0.0;

  for (uint64_t i = 0; i < ranks_a->n; i++) {
    uint32_t khi = tk_iumap_get(idx_b, ranks_a->a[i].i);
    if (khi != tk_iumap_end(idx_b)) {
      double rank_a = (double)i;
      double rank_b = (double)tk_iumap_val(idx_b, khi);
      double weight = (double)ranks_a->a[i].p;
      double da = rank_a - mean_a;
      double db = rank_b - mean_b;
      cov += weight * da * db;
      var_a += weight * da * da;
      var_b += weight * db * db;
    }
  }

  tk_iumap_destroy(idx_b);

  cov /= sum_w;
  var_a /= sum_w;
  var_b /= sum_w;

  if (var_a <= 0.0 || var_b <= 0.0)
    return 1.0;

  return cov / sqrt(var_a * var_b);
}

static inline double tk_rvec_rank_biserial_binary(
  tk_rvec_t *ranks,
  tk_iuset_t *group_1,
  double *out_mean_rank_0,
  double *out_mean_rank_1
) {
  if (!ranks || !group_1 || ranks->n == 0)
    return 0.0;

  double rank_sum_0 = 0.0, rank_sum_1 = 0.0;
  uint64_t count_0 = 0, count_1 = 0;

  for (uint64_t i = 0; i < ranks->n; i++) {
    int64_t id = ranks->a[i].i;
    double rank = (double)(i + 1);

    if (tk_iuset_get(group_1, id) != tk_iuset_end(group_1)) {
      rank_sum_1 += rank;
      count_1++;
    } else {
      rank_sum_0 += rank;
      count_0++;
    }
  }

  if (count_0 == 0 || count_1 == 0 || ranks->n == 0)
    return 0.0;

  double mean_rank_0 = rank_sum_0 / count_0;
  double mean_rank_1 = rank_sum_1 / count_1;

  if (out_mean_rank_0) *out_mean_rank_0 = mean_rank_0;
  if (out_mean_rank_1) *out_mean_rank_1 = mean_rank_1;

  return (mean_rank_0 - mean_rank_1) / (ranks->n / 2.0);
}

static inline double tk_rvec_variance_ratio_binary(
  tk_rvec_t *values,
  tk_iuset_t *group_1,
  double *out_var_0,
  double *out_var_1
) {
  if (!values || !group_1 || values->n == 0)
    return 0.0;

  double sum_0 = 0.0, sum_1 = 0.0;
  uint64_t count_0 = 0, count_1 = 0;

  for (uint64_t i = 0; i < values->n; i++) {
    int64_t id = values->a[i].i;
    double value = values->a[i].d;

    if (tk_iuset_get(group_1, id) != tk_iuset_end(group_1)) {
      sum_1 += value;
      count_1++;
    } else {
      sum_0 += value;
      count_0++;
    }
  }

  if (count_0 < 2 || count_1 < 2)
    return 0.0;

  double mean_0 = sum_0 / count_0;
  double mean_1 = sum_1 / count_1;

  double var_0 = 0.0, var_1 = 0.0;

  for (uint64_t i = 0; i < values->n; i++) {
    int64_t id = values->a[i].i;
    double value = values->a[i].d;
    double diff;

    if (tk_iuset_get(group_1, id) != tk_iuset_end(group_1)) {
      diff = value - mean_1;
      var_1 += diff * diff;
    } else {
      diff = value - mean_0;
      var_0 += diff * diff;
    }
  }

  var_0 /= count_0;
  var_1 /= count_1;

  if (out_var_0) *out_var_0 = var_0;
  if (out_var_1) *out_var_1 = var_1;

  if (var_1 <= 0.0)
    return 0.0;

  return var_0 / var_1;
}

#endif
