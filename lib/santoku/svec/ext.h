#ifndef TK_SVEC_EXT_H
#define TK_SVEC_EXT_H

#if defined(_OPENMP) && !defined(__EMSCRIPTEN__)
#include <omp.h>
#endif
#include <santoku/rvec.h>
#include <santoku/dvec.h>
#include <santoku/iumap.h>
#include <stdatomic.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifndef TK_CVEC_BITS_BYTES
#define TK_CVEC_BITS_BYTES(n) (((n) + CHAR_BIT - 1) / CHAR_BIT)
#endif

static inline void tk_svec_copy_pkeys (tk_svec_t *m0, tk_pvec_t *m1, int64_t start, int64_t end, int64_t dest) {
  if (start < 0 || start >= end || start >= (int64_t) m1->n)
    return;
  if (end >= (int64_t) m1->n)
    end = (int64_t) m1->n;
  uint64_t m = (uint64_t) dest + (uint64_t) (end - start);
  tk_svec_ensure(m0, m);
  uint64_t write = m0->n;
  for (int64_t i = start; i < end; i ++)
    m0->a[write ++] = (int32_t) m1->a[i].i;
  if (m0->n < m)
    m0->n = m;
}

static inline void tk_svec_copy_pvalues (tk_svec_t *m0, tk_pvec_t *m1, int64_t start, int64_t end, int64_t dest) {
  if (start < 0 || start >= end || start >= (int64_t) m1->n)
    return;
  if (end >= (int64_t) m1->n)
    end = (int64_t) m1->n;
  uint64_t m = (uint64_t) dest + (uint64_t) (end - start);
  tk_svec_ensure(m0, m);
  uint64_t write = m0->n;
  for (int64_t i = start; i < end; i ++)
    m0->a[write ++] = (int32_t) m1->a[i].p;
  if (m0->n < m)
    m0->n = m;
}

static inline void tk_svec_copy_rkeys (tk_svec_t *m0, tk_rvec_t *m1, int64_t start, int64_t end, int64_t dest) {
  if (start < 0 || start >= end || start >= (int64_t) m1->n)
    return;
  if (end >= (int64_t) m1->n)
    end = (int64_t) m1->n;
  uint64_t m = (uint64_t) dest + (uint64_t) (end - start);
  tk_svec_ensure(m0, m);
  uint64_t write = m0->n;
  for (int64_t i = start; i < end; i ++)
    m0->a[write ++] = (int32_t) m1->a[i].i;
  if (m0->n < m)
    m0->n = m;
}

static inline void tk_svec_copy_rvalues (tk_svec_t *m0, tk_rvec_t *m1, int64_t start, int64_t end, int64_t dest) {
  if (start < 0 || start >= end || start >= (int64_t) m1->n)
    return;
  if (end >= (int64_t) m1->n)
    end = (int64_t) m1->n;
  uint64_t m = (uint64_t) dest + (uint64_t) (end - start);
  tk_svec_ensure(m0, m);
  uint64_t write = m0->n;
  for (int64_t i = start; i < end; i ++)
    m0->a[write ++] = (int32_t) m1->a[i].d;
  if (m0->n < m)
    m0->n = m;
}

static inline tk_ivec_t *tk_svec_to_ivec (lua_State *L, tk_svec_t *v) {
  tk_ivec_t *out = tk_ivec_create(L, v->n);
  for (uint64_t i = 0; i < v->n; i++)
    out->a[i] = (int64_t)v->a[i];
  return out;
}

static inline tk_dvec_t *tk_svec_to_dvec (lua_State *L, tk_svec_t *v) {
  tk_dvec_t *out = tk_dvec_create(L, v->n);
  for (uint64_t i = 0; i < v->n; i++)
    out->a[i] = (double)v->a[i];
  return out;
}

static inline void tk_svec_lookup (tk_svec_t *indices, tk_svec_t *source) {
  int32_t write_pos = 0;
  for (uint64_t i = 0; i < indices->n; i++) {
    int32_t idx = indices->a[i];
    if (idx >= 0 && idx < (int32_t) source->n) {
      indices->a[write_pos++] = source->a[idx];
    }
  }
  indices->n = (uint64_t) write_pos;
}

typedef enum {
  TK_SVEC_JACCARD,
  TK_SVEC_OVERLAP,
  TK_SVEC_DICE,
  TK_SVEC_TVERSKY,
  TK_SVEC_COSINE,
  TK_SVEC_MIN_KERNEL
} tk_svec_sim_type_t;

static inline double tk_svec_set_jaccard (double inter_w, double sum_a, double sum_b)
{
  double union_w = sum_a + sum_b - inter_w;
  return (union_w == 0.0) ? 0.0 : inter_w / union_w;
}

static inline double tk_svec_set_overlap (double inter_w, double sum_a, double sum_b)
{
  double min_w = (sum_a < sum_b) ? sum_a : sum_b;
  return (min_w == 0.0) ? 0.0 : inter_w / min_w;
}

static inline double tk_svec_set_dice (double inter_w, double sum_a, double sum_b)
{
  double denom = sum_a + sum_b;
  return (denom == 0.0) ? 0.0 : (2.0 * inter_w) / denom;
}

static inline double tk_svec_set_tversky (double inter_w, double sum_a, double sum_b, double alpha, double beta)
{
  double a_only = sum_a - inter_w;
  double b_only = sum_b - inter_w;
  if (a_only < 0.0)
    a_only = 0.0;
  if (b_only < 0.0)
    b_only = 0.0;
  double denom = inter_w + alpha * a_only + beta * b_only;
  return (denom == 0.0) ? 0.0 : inter_w / denom;
}

static inline void tk_svec_set_stats (
  int32_t *a, size_t alen,
  int32_t *b, size_t blen,
  tk_dvec_t *weights,
  double *inter_w,
  double *sum_a,
  double *sum_b
) {
  size_t i = 0, j = 0;
  double inter = 0.0, sa = 0.0, sb = 0.0;
  while (i < alen && j < blen) {
    int32_t ai = a[i], bj = b[j];
    if (ai == bj) {
      double w = (weights && ai >= 0 && ai < (int32_t)weights->n) ? weights->a[ai] : 1.0;
      inter += w;
      sa += w;
      sb += w;
      i++;
      j++;
    } else if (ai < bj) {
      double w = (weights && ai >= 0 && ai < (int32_t)weights->n) ? weights->a[ai] : 1.0;
      sa += w;
      i++;
    } else {
      double w = (weights && bj >= 0 && bj < (int32_t)weights->n) ? weights->a[bj] : 1.0;
      sb += w;
      j++;
    }
  }
  while (i < alen) {
    int32_t ai = a[i++];
    double w = (weights && ai >= 0 && ai < (int32_t)weights->n) ? weights->a[ai] : 1.0;
    sa += w;
  }
  while (j < blen) {
    int32_t bj = b[j++];
    double w = (weights && bj >= 0 && bj < (int32_t)weights->n) ? weights->a[bj] : 1.0;
    sb += w;
  }
  *inter_w = inter;
  *sum_a = sa;
  *sum_b = sb;
}

static inline double tk_svec_set_similarity (
  int32_t *a, size_t alen,
  int32_t *b, size_t blen,
  tk_dvec_t *weights,
  tk_svec_sim_type_t type,
  double tversky_alpha,
  double tversky_beta
) {
  double inter_w = 0.0, sum_a = 0.0, sum_b = 0.0;
  tk_svec_set_stats(a, alen, b, blen, weights, &inter_w, &sum_a, &sum_b);
  switch (type) {
    case TK_SVEC_JACCARD:
      return tk_svec_set_jaccard(inter_w, sum_a, sum_b);
    case TK_SVEC_OVERLAP:
      return tk_svec_set_overlap(inter_w, sum_a, sum_b);
    case TK_SVEC_DICE:
      return tk_svec_set_dice(inter_w, sum_a, sum_b);
    case TK_SVEC_TVERSKY:
      return tk_svec_set_tversky(inter_w, sum_a, sum_b, tversky_alpha, tversky_beta);
    default:
      return tk_svec_set_jaccard(inter_w, sum_a, sum_b);
  }
}

static inline double tk_svec_set_similarity_from_partial (
  double inter_w,
  double q_w,
  double e_w,
  tk_svec_sim_type_t type,
  double tversky_alpha,
  double tversky_beta
) {
  switch (type) {
    case TK_SVEC_JACCARD: {
      double uni_w = q_w + e_w - inter_w;
      return (uni_w == 0.0) ? 0.0 : inter_w / uni_w;
    }
    case TK_SVEC_OVERLAP: {
      double min_w = (q_w < e_w) ? q_w : e_w;
      return (min_w == 0.0) ? 0.0 : inter_w / min_w;
    }
    case TK_SVEC_TVERSKY: {
      double a_only = q_w - inter_w;
      double b_only = e_w - inter_w;
      if (a_only < 0.0) a_only = 0.0;
      if (b_only < 0.0) b_only = 0.0;
      double denom = inter_w + tversky_alpha * a_only + tversky_beta * b_only;
      return (denom == 0.0) ? 0.0 : inter_w / denom;
    }
    case TK_SVEC_DICE: {
      double denom = q_w + e_w;
      return (denom == 0.0) ? 0.0 : (2.0 * inter_w) / denom;
    }
    case TK_SVEC_COSINE: {
      double denom = sqrt(q_w * e_w);
      return (denom == 0.0) ? 0.0 : inter_w / denom;
    }
    case TK_SVEC_MIN_KERNEL: {
      return inter_w;
    }
    default: {
      double uni_w = q_w + e_w - inter_w;
      return (uni_w == 0.0) ? 0.0 : inter_w / uni_w;
    }
  }
}

static inline void tk_svec_set_weights_by_rank (
  int32_t *features,
  size_t n_features,
  tk_dvec_t *weights,
  tk_svec_t *ranks,
  uint64_t n_ranks,
  double *weights_by_rank
) {
  for (uint64_t r = 0; r < n_ranks; r ++)
    weights_by_rank[r] = 0.0;
  for (size_t i = 0; i < n_features; i ++) {
    int32_t fid = features[i];
    if (fid >= 0) {
      double w = (weights && fid < (int32_t)weights->n) ? weights->a[fid] : 1.0;
      int32_t rank = (ranks && fid < (int32_t)ranks->n) ? ranks->a[fid] : 0;
      if (rank >= 0 && rank < (int32_t)n_ranks) {
        weights_by_rank[rank] += w;
      }
    }
  }
}

static inline double tk_svec_set_similarity_by_rank (
  tk_dvec_t *wacc,
  int64_t vsid,
  double *q_weights_by_rank,
  double *e_weights_by_rank,
  uint64_t n_ranks,
  int64_t rank_decay_window,
  double rank_decay_sigma,
  double rank_decay_floor,
  tk_svec_sim_type_t type,
  double tversky_alpha,
  double tversky_beta
) {
  double q_total = 0.0, e_total = 0.0;
  for (uint64_t r = 0; r < n_ranks; r ++) {
    q_total += q_weights_by_rank[r];
    e_total += e_weights_by_rank[r];
  }
  if (q_total == 0.0 && e_total == 0.0)
    return 0.0;
  double total_weighted_sim = 0.0;
  double total_rank_weight = 0.0;
  for (uint64_t rank = 0; rank < n_ranks; rank ++) {
    double rank_weight = 1.0;
    if (rank_decay_window >= 0) {
      if (rank == 0) {
        rank_weight = 1.0;
      } else if (rank >= (uint64_t)rank_decay_window) {
        rank_weight = rank_decay_floor;
      } else {
        double t = (double)rank / (double)rank_decay_window;
        double sigmoid_arg = rank_decay_sigma * (t - 0.5);
        double sigmoid_val = 1.0 / (1.0 + exp(sigmoid_arg));
        rank_weight = rank_decay_floor + (1.0 - rank_decay_floor) * sigmoid_val;
      }
    }
    double inter_w = wacc->a[(int64_t)n_ranks * vsid + (int64_t)rank];
    double q_w = q_weights_by_rank[rank];
    double e_w = e_weights_by_rank[rank];
    double rank_sim;
    if (q_w > 0.0 || e_w > 0.0) {
      rank_sim = tk_svec_set_similarity_from_partial(inter_w, q_w, e_w, type, tversky_alpha, tversky_beta);
    } else {
      rank_sim = 0.0;
    }
    total_weighted_sim += rank_sim * rank_weight;
    total_rank_weight += rank_weight;
  }
  return (total_rank_weight > 0.0) ? total_weighted_sim / total_rank_weight : 0.0;
}

static inline int32_t tk_svec_set_find (int32_t *arr, int32_t start, int32_t end, int32_t value) {
  int32_t left = start;
  int32_t right = end - 1;
  while (left <= right) {
    int32_t mid = left + (right - left) / 2;
    if (arr[mid] == value)
      return mid;
    if (arr[mid] < value)
      left = mid + 1;
    else
      right = mid - 1;
  }
  return -(left + 1);
}

static inline void tk_svec_set_insert (tk_svec_t *vec, int32_t pos, int32_t value) {
  if (pos < 0) pos = -(pos + 1);
  if (pos < 0) pos = 0;
  if (pos > (int32_t)vec->n) pos = (int32_t)vec->n;
  tk_svec_ensure(vec, vec->n + 1);
  if (pos < (int32_t)vec->n) {
    memmove(vec->a + pos + 1, vec->a + pos, (vec->n - (size_t) pos) * sizeof(int32_t));
  }
  vec->a[pos] = value;
  vec->n++;
}

static inline tk_svec_t *tk_svec_set_intersect (lua_State *L, tk_svec_t *a, tk_svec_t *b, tk_svec_t *out) {
  if (out == NULL) {
    size_t min_size = a->n < b->n ? a->n : b->n;
    out = tk_svec_create(L, min_size);
  } else {
    tk_svec_clear(out);
    size_t min_size = a->n < b->n ? a->n : b->n;
    tk_svec_ensure(out, min_size);
  }
  size_t i = 0, j = 0;
  out->n = 0;
  while (i < a->n && j < b->n) {
    if (a->a[i] == b->a[j]) {
      out->a[out->n ++] = a->a[i];
      i ++;
      j ++;
    } else if (a->a[i] < b->a[j]) {
      i ++;
    } else {
      j ++;
    }
  }
  tk_svec_shrink(out);
  return out;
}

static inline tk_svec_t *tk_svec_set_union (lua_State *L, tk_svec_t *a, tk_svec_t *b, tk_svec_t *out) {
  if (out == NULL) {
    out = tk_svec_create(L, a->n + b->n);
  } else {
    tk_svec_clear(out);
    tk_svec_ensure(out, a->n + b->n);
  }
  size_t i = 0, j = 0;
  out->n = 0;
  while (i < a->n && j < b->n) {
    if (a->a[i] == b->a[j]) {
      out->a[out->n ++] = a->a[i];
      i ++;
      j ++;
    } else if (a->a[i] < b->a[j]) {
      out->a[out->n ++] = a->a[i];
      i ++;
    } else {
      out->a[out->n ++] = b->a[j];
      j ++;
    }
  }
  while (i < a->n) {
    out->a[out->n ++] = a->a[i];
    i ++;
  }
  while (j < b->n) {
    out->a[out->n ++] = b->a[j];
    j ++;
  }
  tk_svec_shrink(out);
  return out;
}

static inline size_t tk_svec_scores_max_curvature (
  tk_svec_t *v,
  int32_t *out_val
) {
  size_t n = v->n;
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? v->a[0] : 0;
    return n > 0 ? n - 1 : 0;
  }
  int32_t max_curv = 0;
  size_t max_idx = 1;
  for (size_t i = 1; i < n - 1; i++) {
    int32_t curv = abs(v->a[i-1] - 2 * v->a[i] + v->a[i+1]);
    if (curv > max_curv) {
      max_curv = curv;
      max_idx = i;
    }
  }
  if (max_curv == 0) {
    if (out_val) *out_val = v->a[n - 1];
    return n - 1;
  }
  if (out_val) *out_val = v->a[max_idx];
  return max_idx;
}

static inline size_t tk_svec_scores_lmethod (
  tk_svec_t *v,
  int32_t *out_val
) {
  size_t n = v->n;
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? v->a[0] : 0;
    return n > 0 ? n - 1 : 0;
  }
  int32_t min_val = v->a[0], max_val = v->a[0];
  for (size_t i = 1; i < n; i++) {
    if (v->a[i] < min_val) min_val = v->a[i];
    if (v->a[i] > max_val) max_val = v->a[i];
  }
  if (max_val == min_val) {
    if (out_val) *out_val = v->a[n - 1];
    return n - 1;
  }
  double best_rmse = DBL_MAX;
  size_t best_k = 1;
  for (size_t k = 1; k < n - 1; k++) {
    double sum_x1 = 0.0, sum_y1 = 0.0, sum_xy1 = 0.0, sum_xx1 = 0.0;
    for (size_t i = 0; i <= k; i++) {
      double x = (double)i;
      double y = (double)v->a[i];
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
      double y = (double)v->a[i];
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
      double err = (double)v->a[i] - pred;
      sse += err * err;
    }
    for (size_t i = k + 1; i < n; i++) {
      double pred = slope2 * (double)i + intercept2;
      double err = (double)v->a[i] - pred;
      sse += err * err;
    }
    double rmse = sqrt(sse / (double)n);
    if (rmse < best_rmse) {
      best_rmse = rmse;
      best_k = k;
    }
  }
  if (out_val) *out_val = v->a[best_k];
  return best_k;
}

static inline size_t tk_svec_scores_max_gap (
  tk_svec_t *v,
  int32_t *out_val
) {
  size_t n = v->n;
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? v->a[0] : 0;
    return n > 0 ? n - 1 : 0;
  }
  int32_t max_gap = 0;
  size_t max_idx = 0;
  for (size_t i = 0; i < n - 1; i++) {
    int32_t gap = abs(v->a[i + 1] - v->a[i]);
    if (gap > max_gap) {
      max_gap = gap;
      max_idx = i;
    }
  }
  if (max_gap == 0) {
    if (out_val) *out_val = v->a[n - 1];
    return n - 1;
  }
  if (out_val) *out_val = v->a[max_idx];
  return max_idx;
}

static inline size_t tk_svec_scores_plateau (
  tk_svec_t *v,
  double tolerance,
  int32_t *out_val
) {
  size_t n = v->n;
  if (n == 0) {
    if (out_val) *out_val = 0;
    return 0;
  }
  if (n == 1) {
    if (out_val) *out_val = v->a[0];
    return 0;
  }
  if (tolerance <= 0.0) tolerance = 0.01;

  int32_t min_score = v->a[0], max_score = v->a[0];
  for (size_t i = 1; i < n; i++) {
    if (v->a[i] < min_score) min_score = v->a[i];
    if (v->a[i] > max_score) max_score = v->a[i];
  }
  int32_t range = max_score - min_score;
  if (range <= 0) {
    if (out_val) *out_val = v->a[n - 1];
    return n - 1;
  }

  int32_t abs_tolerance = (int32_t)(tolerance * range);
  int32_t base = v->a[0];
  size_t end_idx = 0;
  for (size_t i = 1; i < n; i++) {
    if (abs(v->a[i] - base) <= abs_tolerance) {
      end_idx = i;
    } else {
      break;
    }
  }
  if (out_val) *out_val = v->a[end_idx];
  return end_idx;
}

static inline size_t tk_svec_scores_kneedle (
  tk_svec_t *v,
  double sensitivity,
  int32_t *out_val
) {
  size_t n = v->n;
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? v->a[0] : 0;
    return n > 0 ? n - 1 : 0;
  }
  if (sensitivity <= 0.0)
    sensitivity = 1.0;
  int32_t min_score = v->a[0];
  int32_t max_score = v->a[0];
  for (size_t i = 1; i < n; i++) {
    if (v->a[i] < min_score) min_score = v->a[i];
    if (v->a[i] > max_score) max_score = v->a[i];
  }
  double score_range = (double)(max_score - min_score);
  if (score_range < 1e-10) {
    if (out_val) *out_val = v->a[0];
    return n - 1;
  }
  double *normalized = (double *)malloc(n * sizeof(double));
  if (!normalized) {
    if (out_val) *out_val = v->a[0];
    return n - 1;
  }
  for (size_t i = 0; i < n; i++) {
    double x_norm = (double)i / (double)(n - 1);
    double y_norm = ((double)v->a[i] - (double)min_score) / score_range;
    normalized[i] = y_norm - x_norm;
  }
  double *smoothed = (double *)malloc(n * sizeof(double));
  if (!smoothed) {
    free(normalized);
    if (out_val) *out_val = v->a[0];
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
    double y_norm = ((double)v->a[i] - (double)min_score) / score_range;
    double diff = y_norm - x_norm;
    if (diff >= threshold) {
      final_knee = i;
    }
  }
  if (out_val) *out_val = v->a[final_knee];
  return final_knee;
}

static inline size_t tk_svec_scores_first_gap (
  tk_svec_t *v,
  double alpha,
  int32_t *out_val
) {
  size_t n = v->n;
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? v->a[0] : 0;
    return n > 0 ? n - 1 : 0;
  }
  if (alpha <= 0.0) alpha = 3.0;

  size_t n_gaps = n - 1;
  int32_t *gaps = (int32_t *)malloc(n_gaps * sizeof(int32_t));
  if (!gaps) {
    if (out_val) *out_val = v->a[n - 1];
    return n - 1;
  }

  for (size_t i = 0; i < n_gaps; i++) {
    gaps[i] = abs(v->a[i + 1] - v->a[i]);
  }

  for (size_t i = 1; i < n_gaps; i++) {
    int32_t key = gaps[i];
    size_t j = i;
    while (j > 0 && gaps[j - 1] > key) {
      gaps[j] = gaps[j - 1];
      j--;
    }
    gaps[j] = key;
  }

  int32_t median_gap;
  if (n_gaps % 2 == 1) {
    median_gap = gaps[n_gaps / 2];
  } else {
    median_gap = (gaps[n_gaps / 2 - 1] + gaps[n_gaps / 2]) / 2;
  }
  free(gaps);

  if (median_gap == 0) {
    int32_t max_gap = 0;
    size_t max_idx = n - 1;
    for (size_t i = 0; i < n - 1; i++) {
      int32_t gap = abs(v->a[i + 1] - v->a[i]);
      if (gap > max_gap) {
        max_gap = gap;
        max_idx = i;
      }
    }
    if (out_val) *out_val = v->a[max_idx];
    return max_idx;
  }

  int32_t threshold = (int32_t)(alpha * (double)median_gap);
  if (threshold < 1) threshold = 1;

  for (size_t i = 0; i < n - 1; i++) {
    int32_t gap = abs(v->a[i + 1] - v->a[i]);
    if (gap > threshold) {
      if (out_val) *out_val = v->a[i];
      return i;
    }
  }

  if (out_val) *out_val = v->a[n - 1];
  return n - 1;
}

static inline size_t tk_svec_scores_otsu (
  tk_svec_t *v,
  int32_t *out_val
) {
  size_t n = v->n;
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? v->a[0] : 0;
    return n > 0 ? n - 1 : 0;
  }

  int32_t min_val = v->a[0], max_val = v->a[0];
  for (size_t i = 1; i < n; i++) {
    if (v->a[i] < min_val) min_val = v->a[i];
    if (v->a[i] > max_val) max_val = v->a[i];
  }
  if (max_val == min_val) {
    if (out_val) *out_val = v->a[n - 1];
    return n - 1;
  }

  double total_sum = 0.0;
  for (size_t i = 0; i < n; i++) {
    total_sum += (double)v->a[i];
  }

  double best_variance = -1.0;
  size_t best_k = 0;
  double sum0 = 0.0;

  for (size_t k = 0; k < n - 1; k++) {
    sum0 += (double)v->a[k];
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

  if (out_val) *out_val = v->a[best_k];
  return best_k;
}

#endif
