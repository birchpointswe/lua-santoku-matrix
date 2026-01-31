#ifndef TK_IUSET_EXT_H
#define TK_IUSET_EXT_H

#if defined(_OPENMP) && !defined(__EMSCRIPTEN__)
#include <omp.h>
#endif
#include <santoku/klib.h>
#include <santoku/ivec.h>
#include <santoku/cvec.h>
#include <santoku/cvec/ext.h>
#include <santoku/rvec.h>
#include <santoku/rvec/ext.h>
#include <santoku/dvec.h>
#include <santoku/iumap.h>
#include <santoku/iuset/base.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include <santoku/iumap.h>

typedef enum {
  TK_POOL_SUM,
  TK_POOL_AVG,
  TK_POOL_MIN,
  TK_POOL_MAX
} tk_pool_t;

static inline tk_pool_t tk_pool_from_string (const char *s) {
  if (!s) return TK_POOL_MAX;
  if (strcmp(s, "sum") == 0) return TK_POOL_SUM;
  if (strcmp(s, "avg") == 0) return TK_POOL_AVG;
  if (strcmp(s, "min") == 0) return TK_POOL_MIN;
  return TK_POOL_MAX;
}

static inline double tk_probit (double p)
{
  if (p <= 0.0) return -1e10;
  if (p >= 1.0) return 1e10;
  static const double a[] = {
    -3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
     1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00
  };
  static const double b[] = {
    -5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
     6.680131188771972e+01, -1.328068155288572e+01
  };
  static const double c[] = {
    -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
    -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00
  };
  static const double d[] = {
    7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
    3.754408661907416e+00
  };
  double plow = 0.02425;
  double phigh = 1.0 - plow;
  double q, r;
  if (p < plow) {
    q = sqrt(-2.0 * log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
  } else if (p <= phigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
  } else {
    q = sqrt(-2.0 * log(1.0 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
  }
}

#define TK_SMOOTH_EPS_FACTOR 0.5
#define TK_SMOOTH_ADD 1.0

static inline double tk_bns_from_marginals (double N, double C, double P, double A)
{
  if (C <= 0 || C >= N || P <= 0 || P >= N)
    return 0.0;
  double tpr = (A + TK_SMOOTH_EPS_FACTOR) / (P + 2.0 * TK_SMOOTH_EPS_FACTOR);
  double fpr = (C - A + TK_SMOOTH_EPS_FACTOR) / (N - P + 2.0 * TK_SMOOTH_EPS_FACTOR);
  return fabs(tk_probit(tpr) - tk_probit(fpr));
}

static inline double tk_mi_from_marginals (double N, double C, double P, double A)
{
  if (C <= 0 || C >= N || P <= 0 || P >= N)
    return 0.0;
  double c00 = (N - C - P + A) + TK_SMOOTH_ADD;
  double c01 = (P - A) + TK_SMOOTH_ADD;
  double c10 = (C - A) + TK_SMOOTH_ADD;
  double c11 = A + TK_SMOOTH_ADD;
  double total = c00 + c01 + c10 + c11;
  if (total <= 0) return 0.0;
  double p_f1 = (c10 + c11) / total;
  double p_f0 = (c00 + c01) / total;
  double p_h1 = (c01 + c11) / total;
  double p_h0 = (c00 + c10) / total;
  double mi = 0.0;
  if (c00 > 0 && p_f0 > 0 && p_h0 > 0) {
    double p_joint = c00 / total;
    mi += p_joint * log2(p_joint / (p_f0 * p_h0));
  }
  if (c01 > 0 && p_f0 > 0 && p_h1 > 0) {
    double p_joint = c01 / total;
    mi += p_joint * log2(p_joint / (p_f0 * p_h1));
  }
  if (c10 > 0 && p_f1 > 0 && p_h0 > 0) {
    double p_joint = c10 / total;
    mi += p_joint * log2(p_joint / (p_f1 * p_h0));
  }
  if (c11 > 0 && p_f1 > 0 && p_h1 > 0) {
    double p_joint = c11 / total;
    mi += p_joint * log2(p_joint / (p_f1 * p_h1));
  }
  return mi;
}

static inline double tk_chi2_from_marginals (double N, double C, double P, double A)
{
  if (C <= 0 || C >= N || P <= 0 || P >= N)
    return 0.0;
  double N_s = N + 4.0;
  double C_s = C + 2.0;
  double P_s = P + 2.0;
  double o00 = N - C - P + A;
  double o01 = P - A;
  double o10 = C - A;
  double o11 = A;
  double e00 = (N_s - C_s) * (N_s - P_s) / N_s;
  double e01 = (N_s - C_s) * P_s / N_s;
  double e10 = C_s * (N_s - P_s) / N_s;
  double e11 = C_s * P_s / N_s;
  double chi2 = 0.0;
  if (e00 > 0) chi2 += (o00 - e00) * (o00 - e00) / e00;
  if (e01 > 0) chi2 += (o01 - e01) * (o01 - e01) / e01;
  if (e10 > 0) chi2 += (o10 - e10) * (o10 - e10) / e10;
  if (e11 > 0) chi2 += (o11 - e11) * (o11 - e11) / e11;
  return chi2;
}

static inline int64_t *tk_ivec_bits_sample_offsets (
  tk_ivec_t *bits,
  uint64_t n_samples,
  uint64_t n_features
) {
  int64_t *offsets = (int64_t *)calloc(n_samples + 1, sizeof(int64_t));
  if (!offsets) return NULL;
  for (uint64_t i = 0; i < bits->n; i++) {
    int64_t bit = bits->a[i];
    if (bit >= 0) {
      uint64_t s = (uint64_t)bit / n_features;
      if (s < n_samples)
        offsets[s + 1]++;
    }
  }
  for (uint64_t s = 1; s <= n_samples; s++)
    offsets[s] += offsets[s - 1];
  return offsets;
}

static inline void tk_build_pooled_union (
  lua_State *L,
  tk_ivec_t *class_offsets,
  tk_ivec_t *class_ids,
  tk_dvec_t *class_weights,
  uint64_t n_hidden,
  uint64_t union_top_k,
  tk_pool_t pool
) {
  tk_iumap_t *feature_counts = tk_iumap_create(0, 0);
  tk_dumap_t *feature_weights = tk_dumap_create(0, 0);
  for (uint64_t c = 0; c < n_hidden; c++) {
    int64_t start = class_offsets->a[c];
    int64_t end = class_offsets->a[c + 1];
    for (int64_t i = start; i < end; i++) {
      int64_t fid = class_ids->a[i];
      double w = class_weights->a[i];
      khint_t kit = tk_iumap_get(feature_counts, fid);
      if (kit != kh_end(feature_counts)) {
        kh_value(feature_counts, kit) += 1;
        khint_t wkit = tk_dumap_get(feature_weights, fid);
        double old_w = kh_value(feature_weights, wkit);
        switch (pool) {
          case TK_POOL_SUM: kh_value(feature_weights, wkit) = old_w + w; break;
          case TK_POOL_AVG: kh_value(feature_weights, wkit) = old_w + w; break;
          case TK_POOL_MIN: kh_value(feature_weights, wkit) = old_w < w ? old_w : w; break;
          case TK_POOL_MAX: kh_value(feature_weights, wkit) = old_w > w ? old_w : w; break;
        }
      } else {
        int absent;
        khint_t newk = tk_iumap_put(feature_counts, fid, &absent);
        kh_value(feature_counts, newk) = 1;
        khint_t wnewk = tk_dumap_put(feature_weights, fid, &absent);
        kh_value(feature_weights, wnewk) = w;
      }
    }
  }
  uint64_t n_union = (uint64_t)kh_size(feature_weights);
  if (pool == TK_POOL_AVG) {
    int64_t fid;
    tk_umap_foreach_keys(feature_weights, fid, ({
      khint_t wkit = tk_dumap_get(feature_weights, fid);
      khint_t ckit = tk_iumap_get(feature_counts, fid);
      kh_value(feature_weights, wkit) /= (double)kh_value(feature_counts, ckit);
    }));
  }
  tk_iumap_destroy(feature_counts);
  tk_ivec_t *union_ids = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *union_weights = tk_dvec_create(L, 0, 0, 0);
  if (union_top_k > 0 && union_top_k < n_union) {
    tk_rvec_t *heap = tk_rvec_create(NULL, 0, 0, 0);
    int64_t fid; double w;
    tk_umap_foreach(feature_weights, fid, w, ({
      tk_rank_t r = { fid, w };
      tk_rvec_hmin(heap, union_top_k, r);
    }));
    tk_rvec_desc(heap, 0, heap->n);
    for (uint64_t i = 0; i < heap->n; i++) {
      tk_ivec_push(union_ids, heap->a[i].i);
      tk_dvec_push(union_weights, heap->a[i].d);
    }
    tk_rvec_destroy(heap);
  } else {
    int64_t fid; double w;
    tk_rvec_t *sorted = tk_rvec_create(NULL, 0, 0, 0);
    tk_umap_foreach(feature_weights, fid, w, ({
      tk_rank_t r = { fid, w };
      tk_rvec_push(sorted, r);
    }));
    tk_rvec_desc(sorted, 0, sorted->n);
    for (uint64_t i = 0; i < sorted->n; i++) {
      tk_ivec_push(union_ids, sorted->a[i].i);
      tk_dvec_push(union_weights, sorted->a[i].d);
    }
    tk_rvec_destroy(sorted);
  }
  tk_dumap_destroy(feature_weights);
}

typedef double (*tk_score_from_marginals_fn)(double N, double C, double P, double A);

static inline tk_ivec_t *tk_ivec_bits_top_sparse_twophase (
  lua_State *L,
  tk_ivec_t *set_bits,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k,
  tk_pool_t pool,
  tk_score_from_marginals_fn score_fn
) {
  atomic_uint *feat_counts = (atomic_uint *)calloc(n_visible, sizeof(atomic_uint));
  atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
  if (!feat_counts || !label_counts) {
    free(feat_counts);
    free(label_counts);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  #pragma omp parallel for schedule(static)
  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit = set_bits->a[i];
    if (bit >= 0)
      atomic_fetch_add(&feat_counts[(uint64_t)bit % n_visible], 1);
  }
  #pragma omp parallel for schedule(static)
  for (uint64_t i = 0; i < labels->n; i++) {
    int64_t bit = labels->a[i];
    if (bit >= 0)
      atomic_fetch_add(&label_counts[(uint64_t)bit % n_hidden], 1);
  }

  int64_t *feat_sample_offsets = tk_ivec_bits_sample_offsets(set_bits, n_samples, n_visible);
  int64_t *label_sample_offsets = tk_ivec_bits_sample_offsets(labels, n_samples, n_hidden);
  if (!feat_sample_offsets || !label_sample_offsets) {
    free(feat_sample_offsets); free(label_sample_offsets);
    free(feat_counts); free(label_counts);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  int n_threads = 1;
  #pragma omp parallel
  { n_threads = omp_get_num_threads(); }

  tk_iumap_t **local_maps = (tk_iumap_t **)calloc((size_t)n_threads, sizeof(tk_iumap_t *));
  if (!local_maps) {
    free(feat_sample_offsets); free(label_sample_offsets);
    free(feat_counts); free(label_counts);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_maps[tid] = tk_iumap_create(0, 0);

    #pragma omp for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      int64_t feat_start = feat_sample_offsets[s];
      int64_t feat_end = feat_sample_offsets[s + 1];
      int64_t label_start = label_sample_offsets[s];
      int64_t label_end = label_sample_offsets[s + 1];
      for (int64_t fi = feat_start; fi < feat_end; fi++) {
        uint64_t f = (uint64_t)set_bits->a[fi] % n_visible;
        for (int64_t li = label_start; li < label_end; li++) {
          uint64_t h = (uint64_t)labels->a[li] % n_hidden;
          int64_t key = (int64_t)(f * n_hidden + h);
          tk_iumap_inc(local_maps[tid], key);
        }
      }
    }
  }
  free(feat_sample_offsets);
  free(label_sample_offsets);

  tk_iumap_t *active_counts = tk_iumap_create(0, 0);
  for (int t = 0; t < n_threads; t++) {
    if (local_maps[t]) {
      int64_t mk, mv;
      tk_umap_foreach(local_maps[t], mk, mv, ({
        khint_t kit = tk_iumap_get(active_counts, mk);
        if (kit != kh_end(active_counts)) {
          kh_value(active_counts, kit) += mv;
        } else {
          int absent;
          khint_t newk = tk_iumap_put(active_counts, mk, &absent);
          kh_value(active_counts, newk) = mv;
        }
      }));
      tk_iumap_destroy(local_maps[t]);
    }
  }
  free(local_maps);

  double *feat_max = (double *)malloc(n_visible * sizeof(double));
  double *feat_min = (double *)malloc(n_visible * sizeof(double));
  double *feat_sum = (double *)calloc(n_visible, sizeof(double));
  uint64_t *feat_count = (uint64_t *)calloc(n_visible, sizeof(uint64_t));
  if (!feat_max || !feat_min || !feat_sum || !feat_count) {
    free(feat_counts);
    free(label_counts);
    tk_iumap_destroy(active_counts);
    free(feat_max);
    free(feat_min);
    free(feat_sum);
    free(feat_count);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }
  #pragma omp parallel for schedule(static)
  for (uint64_t f = 0; f < n_visible; f++) {
    feat_max[f] = -DBL_MAX;
    feat_min[f] = DBL_MAX;
  }

  double N = (double)n_samples;
  int64_t k, v;
  tk_umap_foreach(active_counts, k, v, ({
    uint64_t f = (uint64_t)k / n_hidden;
    uint64_t h = (uint64_t)k % n_hidden;
    if (f >= n_visible || h >= n_hidden) continue;
    double C = (double)atomic_load(&feat_counts[f]);
    double P = (double)atomic_load(&label_counts[h]);
    double A = (double)v;
    double score = score_fn(N, C, P, A);
    feat_sum[f] += score;
    if (score > feat_max[f]) feat_max[f] = score;
    if (score < feat_min[f]) feat_min[f] = score;
    feat_count[f]++;
  }));

  int n_threads2 = 1;
  #pragma omp parallel
  { n_threads2 = omp_get_num_threads(); }
  tk_rvec_t **local_max_heaps = (tk_rvec_t **)calloc((size_t)n_threads2, sizeof(tk_rvec_t *));
  tk_rvec_t **local_pool_heaps = (tk_rvec_t **)calloc((size_t)n_threads2, sizeof(tk_rvec_t *));

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_max_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
    local_pool_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
    #pragma omp for schedule(static)
    for (uint64_t f = 0; f < n_visible; f++) {
      double C = (double)atomic_load(&feat_counts[f]);
      if (feat_count[f] > 0 && C > 0 && C < N) {
        tk_rank_t r = { (int64_t)f, feat_max[f] };
        tk_rvec_hmin(local_max_heaps[tid], top_k, r);
        if (pool != TK_POOL_MAX) {
          double score;
          switch (pool) {
            case TK_POOL_MIN: score = feat_min[f]; break;
            case TK_POOL_AVG: score = feat_sum[f] / (double)feat_count[f]; break;
            default:          score = feat_sum[f]; break;
          }
          tk_rank_t rp = { (int64_t)f, score };
          tk_rvec_hmin(local_pool_heaps[tid], top_k, rp);
        }
      }
    }
  }

  tk_rvec_t *max_heap = tk_rvec_create(NULL, 0, 0, 0);
  tk_rvec_t *pool_heap = tk_rvec_create(NULL, 0, 0, 0);
  for (int t = 0; t < n_threads2; t++) {
    for (uint64_t i = 0; i < local_max_heaps[t]->n; i++)
      tk_rvec_hmin(max_heap, top_k, local_max_heaps[t]->a[i]);
    for (uint64_t i = 0; i < local_pool_heaps[t]->n; i++)
      tk_rvec_hmin(pool_heap, top_k, local_pool_heaps[t]->a[i]);
    tk_rvec_destroy(local_max_heaps[t]);
    tk_rvec_destroy(local_pool_heaps[t]);
  }
  free(local_max_heaps);
  free(local_pool_heaps);

  tk_iuset_t *candidate_set = tk_iuset_create(0, 0);
  for (uint64_t i = 0; i < max_heap->n; i++) {
    int absent;
    tk_iuset_put(candidate_set, max_heap->a[i].i, &absent);
  }
  for (uint64_t i = 0; i < pool_heap->n; i++) {
    int absent;
    tk_iuset_put(candidate_set, pool_heap->a[i].i, &absent);
  }
  tk_rvec_destroy(max_heap);
  tk_rvec_destroy(pool_heap);

  uint64_t n_cand = tk_iuset_size(candidate_set);
  int64_t *cand_features = (int64_t *)malloc(n_cand * sizeof(int64_t));
  if (!cand_features) {
    tk_iuset_destroy(candidate_set);
    free(feat_counts);
    free(label_counts);
    tk_iumap_destroy(active_counts);
    free(feat_max);
    free(feat_min);
    free(feat_sum);
    free(feat_count);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }
  uint64_t ci = 0;
  int64_t ckey;
  tk_umap_foreach_keys(candidate_set, ckey, ({
    cand_features[ci++] = ckey;
  }));
  tk_iuset_destroy(candidate_set);

  free(feat_max);
  free(feat_min);
  free(feat_sum);
  free(feat_count);

  double *final_scores = (double *)calloc(n_cand, sizeof(double));
  if (!final_scores) {
    free(cand_features);
    free(feat_counts);
    free(label_counts);
    tk_iumap_destroy(active_counts);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  #pragma omp parallel for schedule(static)
  for (uint64_t cj = 0; cj < n_cand; cj++) {
    uint64_t f = (uint64_t)cand_features[cj];
    double C = (double)atomic_load(&feat_counts[f]);
    if (C <= 0 || C >= N) continue;
    double pool_sum = 0.0;
    double pool_max = -DBL_MAX;
    double pool_min = DBL_MAX;
    uint64_t pool_count = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      double P = (double)atomic_load(&label_counts[h]);
      if (P <= 0 || P >= N) continue;
      int64_t key = (int64_t)(f * n_hidden + h);
      khint_t kit = tk_iumap_get(active_counts, key);
      double A = (kit != kh_end(active_counts)) ? (double)kh_value(active_counts, kit) : 0.0;
      double score = score_fn(N, C, P, A);
      pool_sum += score;
      if (score > pool_max) pool_max = score;
      if (score < pool_min) pool_min = score;
      pool_count++;
    }
    if (pool_count > 0) {
      switch (pool) {
        case TK_POOL_MIN: final_scores[cj] = pool_min; break;
        case TK_POOL_MAX: final_scores[cj] = pool_max; break;
        case TK_POOL_AVG: final_scores[cj] = pool_sum / (double)pool_count; break;
        default: final_scores[cj] = pool_sum; break;
      }
    }
  }

  tk_rvec_t *top_heap = tk_rvec_create(NULL, 0, 0, 0);
  for (uint64_t cj = 0; cj < n_cand; cj++) {
    tk_rank_t r = { cand_features[cj], final_scores[cj] };
    tk_rvec_hmin(top_heap, top_k, r);
  }

  free(cand_features);
  free(final_scores);
  tk_iumap_destroy(active_counts);
  free(feat_counts);
  free(label_counts);

  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}

static inline tk_ivec_t *tk_ivec_bits_top_sparse_codes (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k,
  tk_pool_t pool,
  tk_score_from_marginals_fn score_fn
) {
  atomic_uint *feat_counts = (atomic_uint *)calloc(n_visible, sizeof(atomic_uint));
  atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
  if (!feat_counts || !label_counts) {
    free(feat_counts);
    free(label_counts);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  #pragma omp parallel for schedule(static)
  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit = set_bits->a[i];
    if (bit >= 0)
      atomic_fetch_add(&feat_counts[(uint64_t)bit % n_visible], 1);
  }
  #pragma omp parallel for schedule(static)
  for (uint64_t s = 0; s < n_samples; s++) {
    uint8_t *s_codes = (uint8_t *)(codes + s * TK_CVEC_BITS_BYTES(n_hidden));
    for (uint64_t b = 0; b < n_hidden; b++) {
      uint64_t byte_idx = b / CHAR_BIT;
      uint64_t bit_pos = b % CHAR_BIT;
      if (s_codes[byte_idx] & (1u << bit_pos))
        atomic_fetch_add(&label_counts[b], 1);
    }
  }

  int n_threads = 1;
  #pragma omp parallel
  { n_threads = omp_get_num_threads(); }

  tk_iumap_t **local_maps = (tk_iumap_t **)calloc((size_t)n_threads, sizeof(tk_iumap_t *));
  if (!local_maps) {
    free(feat_counts);
    free(label_counts);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_maps[tid] = tk_iumap_create(0, 0);

    #pragma omp for schedule(static)
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit_idx = set_bits->a[i];
      if (bit_idx < 0) continue;
      uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
      uint64_t feature_idx = (uint64_t)bit_idx % n_visible;
      if (sample_idx >= n_samples) continue;
      uint8_t *sample_codes = (uint8_t *)(codes + sample_idx * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint64_t bit_pos = b % CHAR_BIT;
        if (sample_codes[byte_idx] & (1u << bit_pos)) {
          int64_t key = (int64_t)(feature_idx * n_hidden + b);
          tk_iumap_inc(local_maps[tid], key);
        }
      }
    }
  }

  tk_iumap_t *active_counts = tk_iumap_create(0, 0);
  for (int t = 0; t < n_threads; t++) {
    if (local_maps[t]) {
      int64_t mk, mv;
      tk_umap_foreach(local_maps[t], mk, mv, ({
        khint_t kit = tk_iumap_get(active_counts, mk);
        if (kit != kh_end(active_counts)) {
          kh_value(active_counts, kit) += mv;
        } else {
          int absent;
          khint_t newk = tk_iumap_put(active_counts, mk, &absent);
          kh_value(active_counts, newk) = mv;
        }
      }));
      tk_iumap_destroy(local_maps[t]);
    }
  }
  free(local_maps);

  double *feat_max = (double *)malloc(n_visible * sizeof(double));
  double *feat_min = (double *)malloc(n_visible * sizeof(double));
  double *feat_sum = (double *)calloc(n_visible, sizeof(double));
  uint64_t *feat_count = (uint64_t *)calloc(n_visible, sizeof(uint64_t));
  if (!feat_max || !feat_min || !feat_sum || !feat_count) {
    free(feat_counts);
    free(label_counts);
    tk_iumap_destroy(active_counts);
    free(feat_max);
    free(feat_min);
    free(feat_sum);
    free(feat_count);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }
  #pragma omp parallel for schedule(static)
  for (uint64_t f = 0; f < n_visible; f++) {
    feat_max[f] = -DBL_MAX;
    feat_min[f] = DBL_MAX;
  }

  double N = (double)n_samples;
  int64_t k, v;
  tk_umap_foreach(active_counts, k, v, ({
    uint64_t f = (uint64_t)k / n_hidden;
    uint64_t h = (uint64_t)k % n_hidden;
    if (f >= n_visible || h >= n_hidden) continue;
    double C = (double)atomic_load(&feat_counts[f]);
    double P = (double)atomic_load(&label_counts[h]);
    double A = (double)v;
    double score = score_fn(N, C, P, A);
    feat_sum[f] += score;
    if (score > feat_max[f]) feat_max[f] = score;
    if (score < feat_min[f]) feat_min[f] = score;
    feat_count[f]++;
  }));

  int n_threads2 = 1;
  #pragma omp parallel
  { n_threads2 = omp_get_num_threads(); }
  tk_rvec_t **local_max_heaps = (tk_rvec_t **)calloc((size_t)n_threads2, sizeof(tk_rvec_t *));
  tk_rvec_t **local_pool_heaps = (tk_rvec_t **)calloc((size_t)n_threads2, sizeof(tk_rvec_t *));

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_max_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
    local_pool_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
    #pragma omp for schedule(static)
    for (uint64_t f = 0; f < n_visible; f++) {
      double C = (double)atomic_load(&feat_counts[f]);
      if (feat_count[f] > 0 && C > 0 && C < N) {
        tk_rank_t r = { (int64_t)f, feat_max[f] };
        tk_rvec_hmin(local_max_heaps[tid], top_k, r);
        if (pool != TK_POOL_MAX) {
          double score;
          switch (pool) {
            case TK_POOL_MIN: score = feat_min[f]; break;
            case TK_POOL_AVG: score = feat_sum[f] / (double)feat_count[f]; break;
            default:          score = feat_sum[f]; break;
          }
          tk_rank_t rp = { (int64_t)f, score };
          tk_rvec_hmin(local_pool_heaps[tid], top_k, rp);
        }
      }
    }
  }

  tk_rvec_t *max_heap = tk_rvec_create(NULL, 0, 0, 0);
  tk_rvec_t *pool_heap = tk_rvec_create(NULL, 0, 0, 0);
  for (int t = 0; t < n_threads2; t++) {
    for (uint64_t i = 0; i < local_max_heaps[t]->n; i++)
      tk_rvec_hmin(max_heap, top_k, local_max_heaps[t]->a[i]);
    for (uint64_t i = 0; i < local_pool_heaps[t]->n; i++)
      tk_rvec_hmin(pool_heap, top_k, local_pool_heaps[t]->a[i]);
    tk_rvec_destroy(local_max_heaps[t]);
    tk_rvec_destroy(local_pool_heaps[t]);
  }
  free(local_max_heaps);
  free(local_pool_heaps);

  tk_iuset_t *candidate_set = tk_iuset_create(0, 0);
  for (uint64_t i = 0; i < max_heap->n; i++) {
    int absent;
    tk_iuset_put(candidate_set, max_heap->a[i].i, &absent);
  }
  for (uint64_t i = 0; i < pool_heap->n; i++) {
    int absent;
    tk_iuset_put(candidate_set, pool_heap->a[i].i, &absent);
  }
  tk_rvec_destroy(max_heap);
  tk_rvec_destroy(pool_heap);

  uint64_t n_cand = tk_iuset_size(candidate_set);
  int64_t *cand_features = (int64_t *)malloc(n_cand * sizeof(int64_t));
  if (!cand_features) {
    tk_iuset_destroy(candidate_set);
    free(feat_counts);
    free(label_counts);
    tk_iumap_destroy(active_counts);
    free(feat_max);
    free(feat_min);
    free(feat_sum);
    free(feat_count);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }
  uint64_t ci = 0;
  int64_t ckey;
  tk_umap_foreach_keys(candidate_set, ckey, ({
    cand_features[ci++] = ckey;
  }));
  tk_iuset_destroy(candidate_set);

  free(feat_max);
  free(feat_min);
  free(feat_sum);
  free(feat_count);

  double *final_scores = (double *)calloc(n_cand, sizeof(double));
  if (!final_scores) {
    free(cand_features);
    free(feat_counts);
    free(label_counts);
    tk_iumap_destroy(active_counts);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  #pragma omp parallel for schedule(static)
  for (uint64_t cj = 0; cj < n_cand; cj++) {
    uint64_t f = (uint64_t)cand_features[cj];
    double C = (double)atomic_load(&feat_counts[f]);
    if (C <= 0 || C >= N) continue;
    double pool_sum = 0.0;
    double pool_max = -DBL_MAX;
    double pool_min = DBL_MAX;
    uint64_t pool_count = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      double P = (double)atomic_load(&label_counts[h]);
      if (P <= 0 || P >= N) continue;
      int64_t key = (int64_t)(f * n_hidden + h);
      khint_t kit = tk_iumap_get(active_counts, key);
      double A = (kit != kh_end(active_counts)) ? (double)kh_value(active_counts, kit) : 0.0;
      double score = score_fn(N, C, P, A);
      pool_sum += score;
      if (score > pool_max) pool_max = score;
      if (score < pool_min) pool_min = score;
      pool_count++;
    }
    if (pool_count > 0) {
      switch (pool) {
        case TK_POOL_MIN: final_scores[cj] = pool_min; break;
        case TK_POOL_MAX: final_scores[cj] = pool_max; break;
        case TK_POOL_AVG: final_scores[cj] = pool_sum / (double)pool_count; break;
        default: final_scores[cj] = pool_sum; break;
      }
    }
  }

  tk_rvec_t *top_heap = tk_rvec_create(NULL, 0, 0, 0);
  for (uint64_t cj = 0; cj < n_cand; cj++) {
    tk_rank_t r = { cand_features[cj], final_scores[cj] };
    tk_rvec_hmin(top_heap, top_k, r);
  }

  free(cand_features);
  free(final_scores);
  tk_iumap_destroy(active_counts);
  free(feat_counts);
  free(label_counts);

  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}

static inline void tk_iuset_union (tk_iuset_t *a, tk_iuset_t *b)
{
  int kha;
  int64_t x;
  tk_umap_foreach_keys(b, x, ({
    tk_iuset_put(a, x, &kha);
  }));
}

static inline double tk_iuset_jaccard (tk_iuset_t *a, tk_iuset_t *b)
{
  uint64_t intersection = tk_iuset_size(a);
  uint64_t union_count = 0;
  union_count = tk_iuset_size(a) + tk_iuset_size(b) - intersection;
  if (union_count == 0)
    return 0.0;
  return (double) intersection / (double) union_count;
}

static inline void tk_iuset_intersect (tk_iuset_t *a, tk_iuset_t *b)
{
  uint32_t i;
  int64_t k;
  tk_umap_foreach_iters(a, i, ({
    k = tk_iuset_key(a, i);
    if (!tk_iuset_contains(b, k))
      tk_iuset_del(a, i);
  }))
}

static inline void tk_iuset_subtract (tk_iuset_t *a, tk_iuset_t *b)
{
  uint32_t i;
  int64_t x;
  tk_umap_foreach_keys(b, x, ({
    i = tk_iuset_get(a, x);
    if (tk_iuset_exist(a, i))
      tk_iuset_del(a, i);
  }))
}

static inline void tk_iuset_union_iumap (tk_iuset_t *a, tk_iumap_t *b)
{
  int kha;
  int64_t x;
  tk_umap_foreach_keys(b, x, ({
    tk_iuset_put(a, x, &kha);
  }))
}

static inline void tk_iuset_intersect_iumap (tk_iuset_t *a, tk_iumap_t *b)
{
  uint32_t i;
  int64_t x;
  tk_umap_foreach_iters(a, i, ({
    x = tk_iuset_key(a, i);
    if (!tk_iumap_contains(b, x))
      tk_iuset_del(a, i);
  }))
}

static inline void tk_iuset_subtract_iumap (tk_iuset_t *a, tk_iumap_t *b)
{
  uint32_t i;
  int64_t x;
  tk_umap_foreach_keys(b, x, ({
    i = tk_iuset_get(a, x);
    if (tk_iuset_exist(a, i))
      tk_iuset_del(a, i);
  }))
}

static inline tk_iuset_t *tk_iuset_from_ivec (lua_State *L, tk_ivec_t *v)
{
  int kha;
  tk_iuset_t *s = tk_iuset_create(L, v->n);
  if (!s)
    return NULL;
  for (uint64_t i = 0; i < v->n; i ++) {
    tk_iuset_put(s, v->a[i], &kha);
    if (kha < 0) {
      tk_iuset_destroy(s);
      return NULL;
    }
  }
  return s;
}

static inline tk_ivec_t *tk_ivec_bits_select (
  tk_ivec_t *src_bits,
  tk_ivec_t *selected_features,
  tk_ivec_t *sample_ids,
  uint64_t n_features,
  tk_ivec_t *dest,
  uint64_t dest_sample,
  uint64_t dest_stride
) {
  if (dest != NULL && dest_sample == 0)
    tk_ivec_clear(dest);

  if (src_bits == NULL || src_bits->n == 0)
    return src_bits;

  if (dest == NULL && (selected_features == NULL || selected_features->n == 0) &&
      (sample_ids == NULL || sample_ids->n == 0))
    return src_bits;

  int64_t *feature_map = NULL;
  uint64_t n_new_features = n_features;
  if (selected_features != NULL && selected_features->n > 0) {
    feature_map = calloc(n_features, sizeof(int64_t));
    if (!feature_map) return NULL;
    for (uint64_t i = 0; i < n_features; i++)
      feature_map[i] = -1;
    for (uint64_t i = 0; i < selected_features->n; i++) {
      int64_t feat = selected_features->a[i];
      if (feat >= 0 && (uint64_t) feat < n_features)
        feature_map[feat] = (int64_t) i;
    }
    n_new_features = selected_features->n;
  }

  tk_iuset_t *sample_set = NULL;
  int64_t *sample_map = NULL;
  uint64_t max_sample = 0;
  if (sample_ids != NULL && sample_ids->n > 0) {
    for (uint64_t i = 0; i < src_bits->n; i++) {
      if (src_bits->a[i] >= 0) {
        uint64_t s = (uint64_t) src_bits->a[i] / n_features;
        if (s > max_sample) max_sample = s;
      }
    }
    sample_set = tk_iuset_from_ivec(NULL, sample_ids);
    if (!sample_set) {
      if (feature_map) free(feature_map);
      return NULL;
    }
    sample_map = calloc(max_sample + 1, sizeof(int64_t));
    if (!sample_map) {
      if (feature_map) free(feature_map);
      tk_iuset_destroy(sample_set);
      return NULL;
    }
    for (uint64_t i = 0; i <= max_sample; i++)
      sample_map[i] = -1;
    uint64_t new_idx = 0;
    for (uint64_t i = 0; i < sample_ids->n; i++) {
      int64_t sid = sample_ids->a[i];
      if (sid >= 0 && (uint64_t) sid <= max_sample)
        sample_map[sid] = (int64_t) new_idx++;
    }
  }

  uint64_t final_stride = (dest_stride > 0) ? dest_stride : n_new_features;

  if (dest != NULL) {
    for (size_t i = 0; i < src_bits->n; i++) {
      int64_t val = src_bits->a[i];
      if (val < 0)
        continue;
      uint64_t sample = (uint64_t) val / n_features;
      uint64_t feature = (uint64_t) val % n_features;
      int64_t new_sample = (int64_t) sample;
      if (sample_set != NULL) {
        if (!tk_iuset_contains(sample_set, (int64_t) sample))
          continue;
        if (sample <= max_sample) {
          new_sample = sample_map[sample];
          if (new_sample < 0)
            continue;
        } else {
          continue;
        }
      }
      int64_t new_feature = (int64_t) feature;
      if (feature_map != NULL) {
        new_feature = feature_map[feature];
        if (new_feature < 0)
          continue;
      }
      new_sample += (int64_t) dest_sample;
      int64_t new_val = new_sample * (int64_t) final_stride + new_feature;
      if (tk_ivec_push(dest, new_val) != 0) {
        if (feature_map) free(feature_map);
        if (sample_map) free(sample_map);
        if (sample_set) tk_iuset_destroy(sample_set);
        return NULL;
      }
    }
    tk_ivec_shrink(dest);
  } else {
    size_t write = 0;
    for (size_t i = 0; i < src_bits->n; i++) {
      int64_t val = src_bits->a[i];
      if (val < 0)
        continue;
      uint64_t sample = (uint64_t) val / n_features;
      uint64_t feature = (uint64_t) val % n_features;
      int64_t new_sample = (int64_t) sample;
      if (sample_set != NULL) {
        if (!tk_iuset_contains(sample_set, (int64_t) sample))
          continue;
        if (sample <= max_sample)
          new_sample = sample_map[sample];
        if (new_sample < 0)
          continue;
      }
      int64_t new_feature = (int64_t) feature;
      if (feature_map != NULL) {
        new_feature = feature_map[feature];
        if (new_feature < 0)
          continue;
      }
      src_bits->a[write++] = new_sample * (int64_t) final_stride + new_feature;
    }
    src_bits->n = write;
  }

  if (feature_map != NULL)
    free(feature_map);
  if (sample_map != NULL)
    free(sample_map);
  if (sample_set != NULL)
    tk_iuset_destroy(sample_set);
  return dest != NULL ? dest : src_bits;
}

static inline tk_cvec_t *tk_cvec_bits_select (
  tk_cvec_t *src_bitmap,
  tk_ivec_t *selected_features,
  tk_ivec_t *sample_ids,
  uint64_t n_features,
  tk_cvec_t *dest,
  uint64_t dest_sample,
  uint64_t dest_stride
) {
  if (src_bitmap == NULL)
    return NULL;

  uint64_t n_samples = src_bitmap->n / TK_CVEC_BITS_BYTES(n_features);

  if (dest == NULL && (selected_features == NULL || selected_features->n == 0) &&
      (sample_ids == NULL || sample_ids->n == 0))
    return src_bitmap;

  tk_iuset_t *sample_set = NULL;
  uint64_t n_output_samples = 0;

  if (sample_ids != NULL && sample_ids->n > 0) {
    sample_set = tk_iuset_from_ivec(NULL, sample_ids);
    if (!sample_set)
      return NULL;
    n_output_samples = sample_ids->n;
  } else {
    n_output_samples = n_samples;
  }

  uint64_t n_selected_features = (selected_features != NULL && selected_features->n > 0)
    ? selected_features->n : n_features;
  uint64_t in_bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t out_bytes_per_sample = TK_CVEC_BITS_BYTES(n_selected_features);
  uint8_t *src_data = (uint8_t *)src_bitmap->a;

  if (dest != NULL) {
    uint64_t final_stride = (dest_stride > 0) ? dest_stride : n_selected_features;
    uint64_t final_bytes_per_sample = TK_CVEC_BITS_BYTES(final_stride);

    if (dest_sample == 0) {
      tk_cvec_ensure(dest, n_output_samples * final_bytes_per_sample);
      dest->n = n_output_samples * final_bytes_per_sample;
    } else {
      tk_cvec_ensure(dest, (dest_sample + n_output_samples) * final_bytes_per_sample);
      dest->n = (dest_sample + n_output_samples) * final_bytes_per_sample;
    }

    uint8_t *dest_data = (uint8_t *)dest->a;
    uint64_t dest_idx = 0;

    for (uint64_t s = 0; s < n_samples; s++) {
      if (sample_set != NULL && !tk_iuset_contains(sample_set, (int64_t) s))
        continue;

      uint8_t *temp = calloc(final_bytes_per_sample, 1);
      if (!temp) {
        if (sample_set != NULL) tk_iuset_destroy(sample_set);
        return NULL;
      }

      uint64_t in_offset = s * in_bytes_per_sample;
      if (selected_features != NULL && selected_features->n > 0) {
        for (uint64_t i = 0; i < selected_features->n; i++) {
          int64_t src_bit = selected_features->a[i];
          if (src_bit >= 0 && (uint64_t) src_bit < n_features) {
            uint64_t src_byte = (uint64_t) src_bit / CHAR_BIT;
            uint8_t src_bit_pos = (uint64_t) src_bit % CHAR_BIT;
            if (src_data[in_offset + src_byte] & (1u << src_bit_pos)) {
              uint64_t dst_byte = i / CHAR_BIT;
              uint8_t dst_bit_pos = i % CHAR_BIT;
              temp[dst_byte] |= (1u << dst_bit_pos);
            }
          }
        }
      } else {
        memcpy(temp, src_data + in_offset, out_bytes_per_sample);
      }

      uint64_t out_offset = (dest_sample + dest_idx) * final_bytes_per_sample;
      memcpy(dest_data + out_offset, temp, out_bytes_per_sample);
      free(temp);
      dest_idx++;
    }
  } else {
    uint8_t *data = src_data;
    uint64_t write_sample = 0;

    for (uint64_t s = 0; s < n_samples; s++) {
      if (sample_set != NULL && !tk_iuset_contains(sample_set, (int64_t) s))
        continue;

      uint8_t *temp = malloc(out_bytes_per_sample);
      if (!temp) {
        if (sample_set != NULL) tk_iuset_destroy(sample_set);
        return NULL;
      }
      memset(temp, 0, out_bytes_per_sample);

      uint64_t in_offset = s * in_bytes_per_sample;
      if (selected_features != NULL && selected_features->n > 0) {
        for (uint64_t i = 0; i < selected_features->n; i++) {
          int64_t src_bit = selected_features->a[i];
          if (src_bit >= 0 && (uint64_t) src_bit < n_features) {
            uint64_t src_byte = (uint64_t) src_bit / CHAR_BIT;
            uint8_t src_bit_pos = (uint64_t) src_bit % CHAR_BIT;
            if (data[in_offset + src_byte] & (1u << src_bit_pos)) {
              uint64_t dst_byte = i / CHAR_BIT;
              uint8_t dst_bit_pos = i % CHAR_BIT;
              temp[dst_byte] |= (1u << dst_bit_pos);
            }
          }
        }
      } else {
        memcpy(temp, data + in_offset, out_bytes_per_sample);
      }

      uint64_t out_offset = write_sample * out_bytes_per_sample;
      memcpy(data + out_offset, temp, out_bytes_per_sample);
      free(temp);
      write_sample++;
    }

    src_bitmap->n = n_output_samples * out_bytes_per_sample;
  }

  if (sample_set != NULL)
    tk_iuset_destroy(sample_set);

  return dest != NULL ? dest : src_bitmap;
}

static inline tk_dvec_t *tk_dvec_mtx_select (
  tk_dvec_t *src_matrix,
  tk_ivec_t *selected_features,
  tk_ivec_t *sample_ids,
  uint64_t n_features,
  tk_dvec_t *dest,
  uint64_t dest_sample,
  uint64_t dest_stride
) {
  if (src_matrix == NULL)
    return NULL;

  uint64_t n_samples = src_matrix->n / n_features;

  if (dest == NULL && (selected_features == NULL || selected_features->n == 0) &&
      (sample_ids == NULL || sample_ids->n == 0))
    return src_matrix;

  tk_iuset_t *sample_set = NULL;
  uint64_t n_output_samples = 0;

  if (sample_ids != NULL && sample_ids->n > 0) {
    sample_set = tk_iuset_from_ivec(NULL, sample_ids);
    if (!sample_set)
      return NULL;
    n_output_samples = sample_ids->n;
  } else {
    n_output_samples = n_samples;
  }

  uint64_t n_selected_features = (selected_features != NULL && selected_features->n > 0)
    ? selected_features->n : n_features;
  double *src_data = src_matrix->a;

  if (dest != NULL) {
    uint64_t final_stride = (dest_stride > 0) ? dest_stride : n_selected_features;

    if (dest_sample == 0) {
      tk_dvec_ensure(dest, n_output_samples * final_stride);
      dest->n = n_output_samples * final_stride;
    } else {
      tk_dvec_ensure(dest, (dest_sample + n_output_samples) * final_stride);
      dest->n = (dest_sample + n_output_samples) * final_stride;
    }

    double *dest_data = dest->a;
    uint64_t dest_idx = 0;

    for (uint64_t s = 0; s < n_samples; s++) {
      if (sample_set != NULL && !tk_iuset_contains(sample_set, (int64_t) s))
        continue;

      uint64_t src_offset = s * n_features;
      uint64_t dest_offset = (dest_sample + dest_idx) * final_stride;

      if (selected_features != NULL && selected_features->n > 0) {
        for (uint64_t f = 0; f < selected_features->n; f++) {
          int64_t feat_idx = selected_features->a[f];
          if (feat_idx >= 0 && (uint64_t)feat_idx < n_features) {
            dest_data[dest_offset + f] = src_data[src_offset + (uint64_t)feat_idx];
          }
        }
      } else {
        memcpy(dest_data + dest_offset, src_data + src_offset, n_features * sizeof(double));
      }

      dest_idx++;
    }

    if (sample_set)
      tk_iuset_destroy(sample_set);

    return dest;
  } else {
    double *src_data = src_matrix->a;
    uint64_t write_idx = 0;

    for (uint64_t s = 0; s < n_samples; s++) {
      if (sample_set != NULL && !tk_iuset_contains(sample_set, (int64_t) s))
        continue;

      uint64_t src_offset = s * n_features;
      uint64_t dest_offset = write_idx * n_selected_features;

      if (selected_features != NULL && selected_features->n > 0) {
        for (uint64_t f = 0; f < selected_features->n; f++) {
          int64_t feat_idx = selected_features->a[f];
          if (feat_idx >= 0 && (uint64_t)feat_idx < n_features) {
            src_data[dest_offset + f] = src_data[src_offset + (uint64_t)feat_idx];
          }
        }
      } else {
        if (dest_offset != src_offset) {
          memmove(src_data + dest_offset, src_data + src_offset, n_features * sizeof(double));
        }
      }

      write_idx++;
    }

    src_matrix->n = n_output_samples * n_selected_features;

    if (sample_set)
      tk_iuset_destroy(sample_set);

    return src_matrix;
  }
}

static inline tk_ivec_t *tk_iuset_keys (lua_State *L, tk_iuset_t *S)
{
  tk_ivec_t *out = tk_ivec_create(L, tk_iuset_size(S), 0, 0);
  int64_t k;
  out->n = 0;
  tk_umap_foreach_keys(S, k, ({
    out->a[out->n ++] = k;
  }));
  return out;
}

static inline tk_ivec_t *tk_ivec_from_iuset (lua_State *L, tk_iuset_t *s)
{
  tk_ivec_t *v = tk_ivec_create(L, tk_iuset_size(s), 0, 0);
  int64_t x;
  v->n = 0;
  tk_umap_foreach_keys(s, x, ({
    v->a[v->n ++] = x;
  }))
  return v;
}

static inline void tk_ivec_bits_top_chi2 (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t per_class_k,
  uint64_t union_top_k,
  tk_pool_t pool
) {
  tk_ivec_t *offsets = tk_ivec_create(L, 0, 0, 0);
  tk_ivec_t *features = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);

  if (codes) {
    atomic_uint *feat_counts = (atomic_uint *)calloc(n_visible, sizeof(atomic_uint));
    atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
    if (!feat_counts || !label_counts) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit = set_bits->a[i];
      if (bit >= 0)
        atomic_fetch_add(&feat_counts[(uint64_t)bit % n_visible], 1);
    }
    #pragma omp parallel for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint64_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos))
          atomic_fetch_add(&label_counts[b], 1);
      }
    }

    int n_threads = 1;
    #pragma omp parallel
    { n_threads = omp_get_num_threads(); }
    tk_iumap_t **local_maps = (tk_iumap_t **)calloc((size_t)n_threads, sizeof(tk_iumap_t *));
    if (!local_maps) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      local_maps[tid] = tk_iumap_create(0, 0);
      #pragma omp for schedule(static)
      for (uint64_t i = 0; i < set_bits->n; i++) {
        int64_t bit_idx = set_bits->a[i];
        if (bit_idx < 0) continue;
        uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
        uint64_t feature_idx = (uint64_t)bit_idx % n_visible;
        if (sample_idx >= n_samples) continue;
        uint8_t *sample_codes = (uint8_t *)(codes + sample_idx * TK_CVEC_BITS_BYTES(n_hidden));
        for (uint64_t b = 0; b < n_hidden; b++) {
          uint64_t byte_idx = b / CHAR_BIT;
          uint64_t bit_pos = b % CHAR_BIT;
          if (sample_codes[byte_idx] & (1u << bit_pos))
            tk_iumap_inc(local_maps[tid], (int64_t)(feature_idx * n_hidden + b));
        }
      }
    }

    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    for (int t = 0; t < n_threads; t++) {
      if (local_maps[t]) {
        int64_t mk, mv;
        tk_umap_foreach(local_maps[t], mk, mv, ({
          khint_t kit = tk_iumap_get(active_counts, mk);
          if (kit != kh_end(active_counts))
            kh_value(active_counts, kit) += mv;
          else {
            int absent;
            khint_t newk = tk_iumap_put(active_counts, mk, &absent);
            kh_value(active_counts, newk) = mv;
          }
        }));
        tk_iumap_destroy(local_maps[t]);
      }
    }
    free(local_maps);

    double N = (double)n_samples;
    tk_rvec_t **class_heaps = (tk_rvec_t **)calloc(n_hidden, sizeof(tk_rvec_t *));
    #pragma omp parallel for schedule(static)
    for (uint64_t b = 0; b < n_hidden; b++) {
      class_heaps[b] = tk_rvec_create(NULL, 0, 0, 0);
      double P = (double)atomic_load(&label_counts[b]);
      if (P > 0 && P < N) {
        for (uint64_t f = 0; f < n_visible; f++) {
          double C = (double)atomic_load(&feat_counts[f]);
          if (C <= 0 || C >= N) continue;
          double A = 0.0;
          khint_t kit = tk_iumap_get(active_counts, (int64_t)(f * n_hidden + b));
          if (kit != kh_end(active_counts))
            A = (double)kh_value(active_counts, kit);
          if (A <= C * P / N) continue;
          double chi2 = tk_chi2_from_marginals(N, C, P, A);
          tk_rank_t r = { (int64_t)f, chi2 };
          tk_rvec_hmin(class_heaps[b], per_class_k, r);
        }
      }
      if (class_heaps[b]->n > 0)
        tk_rvec_desc(class_heaps[b], 0, class_heaps[b]->n);
    }
    tk_ivec_push(offsets, 0);
    for (uint64_t b = 0; b < n_hidden; b++) {
      for (uint64_t i = 0; i < class_heaps[b]->n; i++) {
        tk_ivec_push(features, class_heaps[b]->a[i].i);
        tk_dvec_push(weights, class_heaps[b]->a[i].d);
      }
      tk_ivec_push(offsets, (int64_t)features->n);
      tk_rvec_destroy(class_heaps[b]);
    }
    free(class_heaps);

    tk_iumap_destroy(active_counts);
    free(feat_counts); free(label_counts);

  } else if (labels) {
    atomic_uint *feat_counts = (atomic_uint *)calloc(n_visible, sizeof(atomic_uint));
    atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
    if (!feat_counts || !label_counts) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit = set_bits->a[i];
      if (bit >= 0)
        atomic_fetch_add(&feat_counts[(uint64_t)bit % n_visible], 1);
    }
    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0)
        atomic_fetch_add(&label_counts[(uint64_t)bit % n_hidden], 1);
    }

    int64_t *feat_sample_offsets = tk_ivec_bits_sample_offsets(set_bits, n_samples, n_visible);
    int64_t *label_sample_offsets = tk_ivec_bits_sample_offsets(labels, n_samples, n_hidden);
    if (!feat_sample_offsets || !label_sample_offsets) {
      free(feat_sample_offsets); free(label_sample_offsets);
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    int n_threads = 1;
    #pragma omp parallel
    { n_threads = omp_get_num_threads(); }
    tk_iumap_t **local_maps = (tk_iumap_t **)calloc((size_t)n_threads, sizeof(tk_iumap_t *));
    if (!local_maps) {
      free(feat_sample_offsets); free(label_sample_offsets);
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      local_maps[tid] = tk_iumap_create(0, 0);
      #pragma omp for schedule(static)
      for (uint64_t s = 0; s < n_samples; s++) {
        int64_t feat_start = feat_sample_offsets[s];
        int64_t feat_end = feat_sample_offsets[s + 1];
        int64_t label_start = label_sample_offsets[s];
        int64_t label_end = label_sample_offsets[s + 1];
        for (int64_t fi = feat_start; fi < feat_end; fi++) {
          uint64_t f = (uint64_t)set_bits->a[fi] % n_visible;
          for (int64_t li = label_start; li < label_end; li++) {
            uint64_t h = (uint64_t)labels->a[li] % n_hidden;
            tk_iumap_inc(local_maps[tid], (int64_t)(f * n_hidden + h));
          }
        }
      }
    }
    free(feat_sample_offsets);
    free(label_sample_offsets);

    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    for (int t = 0; t < n_threads; t++) {
      if (local_maps[t]) {
        int64_t mk, mv;
        tk_umap_foreach(local_maps[t], mk, mv, ({
          khint_t kit = tk_iumap_get(active_counts, mk);
          if (kit != kh_end(active_counts))
            kh_value(active_counts, kit) += mv;
          else {
            int absent;
            khint_t newk = tk_iumap_put(active_counts, mk, &absent);
            kh_value(active_counts, newk) = mv;
          }
        }));
        tk_iumap_destroy(local_maps[t]);
      }
    }
    free(local_maps);

    double N = (double)n_samples;
    tk_rvec_t **class_heaps = (tk_rvec_t **)calloc(n_hidden, sizeof(tk_rvec_t *));
    #pragma omp parallel for schedule(static)
    for (uint64_t b = 0; b < n_hidden; b++) {
      class_heaps[b] = tk_rvec_create(NULL, 0, 0, 0);
      double P = (double)atomic_load(&label_counts[b]);
      if (P > 0 && P < N) {
        for (uint64_t f = 0; f < n_visible; f++) {
          double C = (double)atomic_load(&feat_counts[f]);
          if (C <= 0 || C >= N) continue;
          double A = 0.0;
          khint_t kit = tk_iumap_get(active_counts, (int64_t)(f * n_hidden + b));
          if (kit != kh_end(active_counts))
            A = (double)kh_value(active_counts, kit);
          if (A <= C * P / N) continue;
          double chi2 = tk_chi2_from_marginals(N, C, P, A);
          tk_rank_t r = { (int64_t)f, chi2 };
          tk_rvec_hmin(class_heaps[b], per_class_k, r);
        }
      }
      if (class_heaps[b]->n > 0)
        tk_rvec_desc(class_heaps[b], 0, class_heaps[b]->n);
    }
    tk_ivec_push(offsets, 0);
    for (uint64_t b = 0; b < n_hidden; b++) {
      for (uint64_t i = 0; i < class_heaps[b]->n; i++) {
        tk_ivec_push(features, class_heaps[b]->a[i].i);
        tk_dvec_push(weights, class_heaps[b]->a[i].d);
      }
      tk_ivec_push(offsets, (int64_t)features->n);
      tk_rvec_destroy(class_heaps[b]);
    }
    free(class_heaps);

    tk_iumap_destroy(active_counts);
    free(feat_counts); free(label_counts);

  } else {
    for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
  }
  tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
  lua_insert(L, -5);
  lua_insert(L, -5);
}

static inline void tk_ivec_bits_top_mi (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t per_class_k,
  uint64_t union_top_k,
  tk_pool_t pool
) {
  tk_ivec_t *offsets = tk_ivec_create(L, 0, 0, 0);
  tk_ivec_t *features = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);

  if (codes) {
    atomic_uint *feat_counts = (atomic_uint *)calloc(n_visible, sizeof(atomic_uint));
    atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
    if (!feat_counts || !label_counts) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit = set_bits->a[i];
      if (bit >= 0)
        atomic_fetch_add(&feat_counts[(uint64_t)bit % n_visible], 1);
    }
    #pragma omp parallel for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint64_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos))
          atomic_fetch_add(&label_counts[b], 1);
      }
    }

    int n_threads = 1;
    #pragma omp parallel
    { n_threads = omp_get_num_threads(); }
    tk_iumap_t **local_maps = (tk_iumap_t **)calloc((size_t)n_threads, sizeof(tk_iumap_t *));
    if (!local_maps) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      local_maps[tid] = tk_iumap_create(0, 0);
      #pragma omp for schedule(static)
      for (uint64_t i = 0; i < set_bits->n; i++) {
        int64_t bit_idx = set_bits->a[i];
        if (bit_idx < 0) continue;
        uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
        uint64_t feature_idx = (uint64_t)bit_idx % n_visible;
        if (sample_idx >= n_samples) continue;
        uint8_t *sample_codes = (uint8_t *)(codes + sample_idx * TK_CVEC_BITS_BYTES(n_hidden));
        for (uint64_t b = 0; b < n_hidden; b++) {
          uint64_t byte_idx = b / CHAR_BIT;
          uint64_t bit_pos = b % CHAR_BIT;
          if (sample_codes[byte_idx] & (1u << bit_pos))
            tk_iumap_inc(local_maps[tid], (int64_t)(feature_idx * n_hidden + b));
        }
      }
    }

    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    for (int t = 0; t < n_threads; t++) {
      if (local_maps[t]) {
        int64_t mk, mv;
        tk_umap_foreach(local_maps[t], mk, mv, ({
          khint_t kit = tk_iumap_get(active_counts, mk);
          if (kit != kh_end(active_counts))
            kh_value(active_counts, kit) += mv;
          else {
            int absent;
            khint_t newk = tk_iumap_put(active_counts, mk, &absent);
            kh_value(active_counts, newk) = mv;
          }
        }));
        tk_iumap_destroy(local_maps[t]);
      }
    }
    free(local_maps);

    double N = (double)n_samples;
    tk_rvec_t **class_heaps = (tk_rvec_t **)calloc(n_hidden, sizeof(tk_rvec_t *));
    #pragma omp parallel for schedule(static)
    for (uint64_t b = 0; b < n_hidden; b++) {
      class_heaps[b] = tk_rvec_create(NULL, 0, 0, 0);
      double P = (double)atomic_load(&label_counts[b]);
      if (P > 0 && P < N) {
        for (uint64_t f = 0; f < n_visible; f++) {
          double C = (double)atomic_load(&feat_counts[f]);
          if (C <= 0 || C >= N) continue;
          double A = 0.0;
          khint_t kit = tk_iumap_get(active_counts, (int64_t)(f * n_hidden + b));
          if (kit != kh_end(active_counts))
            A = (double)kh_value(active_counts, kit);
          if (A <= C * P / N) continue;
          double mi = tk_mi_from_marginals(N, C, P, A);
          tk_rank_t r = { (int64_t)f, mi };
          tk_rvec_hmin(class_heaps[b], per_class_k, r);
        }
      }
      if (class_heaps[b]->n > 0)
        tk_rvec_desc(class_heaps[b], 0, class_heaps[b]->n);
    }
    tk_ivec_push(offsets, 0);
    for (uint64_t b = 0; b < n_hidden; b++) {
      for (uint64_t i = 0; i < class_heaps[b]->n; i++) {
        tk_ivec_push(features, class_heaps[b]->a[i].i);
        tk_dvec_push(weights, class_heaps[b]->a[i].d);
      }
      tk_ivec_push(offsets, (int64_t)features->n);
      tk_rvec_destroy(class_heaps[b]);
    }
    free(class_heaps);

    tk_iumap_destroy(active_counts);
    free(feat_counts); free(label_counts);

  } else if (labels) {
    atomic_uint *feat_counts = (atomic_uint *)calloc(n_visible, sizeof(atomic_uint));
    atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
    if (!feat_counts || !label_counts) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit = set_bits->a[i];
      if (bit >= 0)
        atomic_fetch_add(&feat_counts[(uint64_t)bit % n_visible], 1);
    }
    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0)
        atomic_fetch_add(&label_counts[(uint64_t)bit % n_hidden], 1);
    }

    int64_t *feat_sample_offsets = tk_ivec_bits_sample_offsets(set_bits, n_samples, n_visible);
    int64_t *label_sample_offsets = tk_ivec_bits_sample_offsets(labels, n_samples, n_hidden);
    if (!feat_sample_offsets || !label_sample_offsets) {
      free(feat_sample_offsets); free(label_sample_offsets);
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    int n_threads = 1;
    #pragma omp parallel
    { n_threads = omp_get_num_threads(); }
    tk_iumap_t **local_maps = (tk_iumap_t **)calloc((size_t)n_threads, sizeof(tk_iumap_t *));
    if (!local_maps) {
      free(feat_sample_offsets); free(label_sample_offsets);
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      local_maps[tid] = tk_iumap_create(0, 0);
      #pragma omp for schedule(static)
      for (uint64_t s = 0; s < n_samples; s++) {
        int64_t feat_start = feat_sample_offsets[s];
        int64_t feat_end = feat_sample_offsets[s + 1];
        int64_t label_start = label_sample_offsets[s];
        int64_t label_end = label_sample_offsets[s + 1];
        for (int64_t fi = feat_start; fi < feat_end; fi++) {
          uint64_t f = (uint64_t)set_bits->a[fi] % n_visible;
          for (int64_t li = label_start; li < label_end; li++) {
            uint64_t h = (uint64_t)labels->a[li] % n_hidden;
            tk_iumap_inc(local_maps[tid], (int64_t)(f * n_hidden + h));
          }
        }
      }
    }
    free(feat_sample_offsets);
    free(label_sample_offsets);

    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    for (int t = 0; t < n_threads; t++) {
      if (local_maps[t]) {
        int64_t mk, mv;
        tk_umap_foreach(local_maps[t], mk, mv, ({
          khint_t kit = tk_iumap_get(active_counts, mk);
          if (kit != kh_end(active_counts))
            kh_value(active_counts, kit) += mv;
          else {
            int absent;
            khint_t newk = tk_iumap_put(active_counts, mk, &absent);
            kh_value(active_counts, newk) = mv;
          }
        }));
        tk_iumap_destroy(local_maps[t]);
      }
    }
    free(local_maps);

    double N = (double)n_samples;
    tk_rvec_t **class_heaps = (tk_rvec_t **)calloc(n_hidden, sizeof(tk_rvec_t *));
    #pragma omp parallel for schedule(static)
    for (uint64_t b = 0; b < n_hidden; b++) {
      class_heaps[b] = tk_rvec_create(NULL, 0, 0, 0);
      double P = (double)atomic_load(&label_counts[b]);
      if (P > 0 && P < N) {
        for (uint64_t f = 0; f < n_visible; f++) {
          double C = (double)atomic_load(&feat_counts[f]);
          if (C <= 0 || C >= N) continue;
          double A = 0.0;
          khint_t kit = tk_iumap_get(active_counts, (int64_t)(f * n_hidden + b));
          if (kit != kh_end(active_counts))
            A = (double)kh_value(active_counts, kit);
          if (A <= C * P / N) continue;
          double mi = tk_mi_from_marginals(N, C, P, A);
          tk_rank_t r = { (int64_t)f, mi };
          tk_rvec_hmin(class_heaps[b], per_class_k, r);
        }
      }
      if (class_heaps[b]->n > 0)
        tk_rvec_desc(class_heaps[b], 0, class_heaps[b]->n);
    }
    tk_ivec_push(offsets, 0);
    for (uint64_t b = 0; b < n_hidden; b++) {
      for (uint64_t i = 0; i < class_heaps[b]->n; i++) {
        tk_ivec_push(features, class_heaps[b]->a[i].i);
        tk_dvec_push(weights, class_heaps[b]->a[i].d);
      }
      tk_ivec_push(offsets, (int64_t)features->n);
      tk_rvec_destroy(class_heaps[b]);
    }
    free(class_heaps);

    tk_iumap_destroy(active_counts);
    free(feat_counts); free(label_counts);

  } else {
    for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
  }
  tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
  lua_insert(L, -5);
  lua_insert(L, -5);
}

static inline void tk_ivec_bits_top_bns (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t per_class_k,
  uint64_t union_top_k,
  tk_pool_t pool
) {
  tk_ivec_t *offsets = tk_ivec_create(L, 0, 0, 0);
  tk_ivec_t *features = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);

  if (codes) {
    atomic_uint *feat_counts = (atomic_uint *)calloc(n_visible, sizeof(atomic_uint));
    atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
    if (!feat_counts || !label_counts) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit = set_bits->a[i];
      if (bit >= 0)
        atomic_fetch_add(&feat_counts[(uint64_t)bit % n_visible], 1);
    }
    #pragma omp parallel for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint64_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos))
          atomic_fetch_add(&label_counts[b], 1);
      }
    }

    int n_threads = 1;
    #pragma omp parallel
    { n_threads = omp_get_num_threads(); }
    tk_iumap_t **local_maps = (tk_iumap_t **)calloc((size_t)n_threads, sizeof(tk_iumap_t *));
    if (!local_maps) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      local_maps[tid] = tk_iumap_create(0, 0);
      #pragma omp for schedule(static)
      for (uint64_t i = 0; i < set_bits->n; i++) {
        int64_t bit_idx = set_bits->a[i];
        if (bit_idx < 0) continue;
        uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
        uint64_t feature_idx = (uint64_t)bit_idx % n_visible;
        if (sample_idx >= n_samples) continue;
        uint8_t *sample_codes = (uint8_t *)(codes + sample_idx * TK_CVEC_BITS_BYTES(n_hidden));
        for (uint64_t b = 0; b < n_hidden; b++) {
          uint64_t byte_idx = b / CHAR_BIT;
          uint64_t bit_pos = b % CHAR_BIT;
          if (sample_codes[byte_idx] & (1u << bit_pos))
            tk_iumap_inc(local_maps[tid], (int64_t)(feature_idx * n_hidden + b));
        }
      }
    }

    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    for (int t = 0; t < n_threads; t++) {
      if (local_maps[t]) {
        int64_t mk, mv;
        tk_umap_foreach(local_maps[t], mk, mv, ({
          khint_t kit = tk_iumap_get(active_counts, mk);
          if (kit != kh_end(active_counts))
            kh_value(active_counts, kit) += mv;
          else {
            int absent;
            khint_t newk = tk_iumap_put(active_counts, mk, &absent);
            kh_value(active_counts, newk) = mv;
          }
        }));
        tk_iumap_destroy(local_maps[t]);
      }
    }
    free(local_maps);

    double N = (double)n_samples;
    tk_rvec_t **class_heaps = (tk_rvec_t **)calloc(n_hidden, sizeof(tk_rvec_t *));
    #pragma omp parallel for schedule(static)
    for (uint64_t b = 0; b < n_hidden; b++) {
      class_heaps[b] = tk_rvec_create(NULL, 0, 0, 0);
      double P = (double)atomic_load(&label_counts[b]);
      if (P > 0 && P < N) {
        for (uint64_t f = 0; f < n_visible; f++) {
          double C = (double)atomic_load(&feat_counts[f]);
          if (C <= 0 || C >= N) continue;
          double A = 0.0;
          khint_t kit = tk_iumap_get(active_counts, (int64_t)(f * n_hidden + b));
          if (kit != kh_end(active_counts))
            A = (double)kh_value(active_counts, kit);
          if (A <= C * P / N) continue;
          double bns = tk_bns_from_marginals(N, C, P, A);
          tk_rank_t r = { (int64_t)f, bns };
          tk_rvec_hmin(class_heaps[b], per_class_k, r);
        }
      }
      if (class_heaps[b]->n > 0)
        tk_rvec_desc(class_heaps[b], 0, class_heaps[b]->n);
    }
    tk_ivec_push(offsets, 0);
    for (uint64_t b = 0; b < n_hidden; b++) {
      for (uint64_t i = 0; i < class_heaps[b]->n; i++) {
        tk_ivec_push(features, class_heaps[b]->a[i].i);
        tk_dvec_push(weights, class_heaps[b]->a[i].d);
      }
      tk_ivec_push(offsets, (int64_t)features->n);
      tk_rvec_destroy(class_heaps[b]);
    }
    free(class_heaps);

    tk_iumap_destroy(active_counts);
    free(feat_counts); free(label_counts);

  } else if (labels) {
    atomic_uint *feat_counts = (atomic_uint *)calloc(n_visible, sizeof(atomic_uint));
    atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
    if (!feat_counts || !label_counts) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit = set_bits->a[i];
      if (bit >= 0)
        atomic_fetch_add(&feat_counts[(uint64_t)bit % n_visible], 1);
    }
    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0)
        atomic_fetch_add(&label_counts[(uint64_t)bit % n_hidden], 1);
    }

    int64_t *feat_sample_offsets = tk_ivec_bits_sample_offsets(set_bits, n_samples, n_visible);
    int64_t *label_sample_offsets = tk_ivec_bits_sample_offsets(labels, n_samples, n_hidden);
    if (!feat_sample_offsets || !label_sample_offsets) {
      free(feat_sample_offsets); free(label_sample_offsets);
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    int n_threads = 1;
    #pragma omp parallel
    { n_threads = omp_get_num_threads(); }
    tk_iumap_t **local_maps = (tk_iumap_t **)calloc((size_t)n_threads, sizeof(tk_iumap_t *));
    if (!local_maps) {
      free(feat_sample_offsets); free(label_sample_offsets);
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
      lua_insert(L, -5);
      lua_insert(L, -5);
      return;
    }

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      local_maps[tid] = tk_iumap_create(0, 0);
      #pragma omp for schedule(static)
      for (uint64_t s = 0; s < n_samples; s++) {
        int64_t feat_start = feat_sample_offsets[s];
        int64_t feat_end = feat_sample_offsets[s + 1];
        int64_t label_start = label_sample_offsets[s];
        int64_t label_end = label_sample_offsets[s + 1];
        for (int64_t fi = feat_start; fi < feat_end; fi++) {
          uint64_t f = (uint64_t)set_bits->a[fi] % n_visible;
          for (int64_t li = label_start; li < label_end; li++) {
            uint64_t h = (uint64_t)labels->a[li] % n_hidden;
            tk_iumap_inc(local_maps[tid], (int64_t)(f * n_hidden + h));
          }
        }
      }
    }
    free(feat_sample_offsets);
    free(label_sample_offsets);

    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    for (int t = 0; t < n_threads; t++) {
      if (local_maps[t]) {
        int64_t mk, mv;
        tk_umap_foreach(local_maps[t], mk, mv, ({
          khint_t kit = tk_iumap_get(active_counts, mk);
          if (kit != kh_end(active_counts))
            kh_value(active_counts, kit) += mv;
          else {
            int absent;
            khint_t newk = tk_iumap_put(active_counts, mk, &absent);
            kh_value(active_counts, newk) = mv;
          }
        }));
        tk_iumap_destroy(local_maps[t]);
      }
    }
    free(local_maps);

    double N = (double)n_samples;
    tk_rvec_t **class_heaps = (tk_rvec_t **)calloc(n_hidden, sizeof(tk_rvec_t *));
    #pragma omp parallel for schedule(static)
    for (uint64_t b = 0; b < n_hidden; b++) {
      class_heaps[b] = tk_rvec_create(NULL, 0, 0, 0);
      double P = (double)atomic_load(&label_counts[b]);
      if (P > 0 && P < N) {
        for (uint64_t f = 0; f < n_visible; f++) {
          double C = (double)atomic_load(&feat_counts[f]);
          if (C <= 0 || C >= N) continue;
          double A = 0.0;
          khint_t kit = tk_iumap_get(active_counts, (int64_t)(f * n_hidden + b));
          if (kit != kh_end(active_counts))
            A = (double)kh_value(active_counts, kit);
          if (A <= C * P / N) continue;
          double bns = tk_bns_from_marginals(N, C, P, A);
          tk_rank_t r = { (int64_t)f, bns };
          tk_rvec_hmin(class_heaps[b], per_class_k, r);
        }
      }
      if (class_heaps[b]->n > 0)
        tk_rvec_desc(class_heaps[b], 0, class_heaps[b]->n);
    }
    tk_ivec_push(offsets, 0);
    for (uint64_t b = 0; b < n_hidden; b++) {
      for (uint64_t i = 0; i < class_heaps[b]->n; i++) {
        tk_ivec_push(features, class_heaps[b]->a[i].i);
        tk_dvec_push(weights, class_heaps[b]->a[i].d);
      }
      tk_ivec_push(offsets, (int64_t)features->n);
      tk_rvec_destroy(class_heaps[b]);
    }
    free(class_heaps);

    tk_iumap_destroy(active_counts);
    free(feat_counts); free(label_counts);

  } else {
    for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
  }
  tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
  lua_insert(L, -5);
  lua_insert(L, -5);
}

static inline tk_ivec_t *tk_ivec_bits_top_entropy (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_hidden,
  uint64_t top_k
) {
  uint64_t *bit_counts = (uint64_t *)calloc(n_hidden, sizeof(uint64_t));
  if (!bit_counts)
    return NULL;
  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit_idx = set_bits->a[i];
    if (bit_idx < 0)
      continue;
    uint64_t sample_idx = (uint64_t)bit_idx / n_hidden;
    uint64_t hidden_idx = (uint64_t)bit_idx % n_hidden;
    if (sample_idx >= n_samples || hidden_idx >= n_hidden)
      continue;
    bit_counts[hidden_idx]++;
  }
  int n_threads = 1;
  #pragma omp parallel
  { n_threads = omp_get_num_threads(); }
  tk_rvec_t **local_heaps = (tk_rvec_t **)calloc((size_t)n_threads, sizeof(tk_rvec_t *));
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
    #pragma omp for schedule(static)
    for (uint64_t h = 0; h < n_hidden; h++) {
      double p = (double)bit_counts[h] / (double)n_samples;
      double entropy = 0.0;
      if (p > 0.0 && p < 1.0)
        entropy = -(p * log2(p) + (1.0 - p) * log2(1.0 - p));
      tk_rank_t r = { (int64_t)h, entropy };
      tk_rvec_hmin(local_heaps[tid], top_k, r);
    }
  }
  tk_rvec_t *top_heap = tk_rvec_create(NULL, 0, 0, 0);
  for (int t = 0; t < n_threads; t++) {
    for (uint64_t i = 0; i < local_heaps[t]->n; i++)
      tk_rvec_hmin(top_heap, top_k, local_heaps[t]->a[i]);
    tk_rvec_destroy(local_heaps[t]);
  }
  free(local_heaps);
  free(bit_counts);
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}



static inline void tk_cvec_bits_top_chi2 (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  uint64_t per_class_k,
  uint64_t union_top_k,
  tk_pool_t pool
) {
  tk_ivec_t *offsets = tk_ivec_create(L, 0, 0, 0);
  tk_ivec_t *features = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);

  if (codes) {
    atomic_uint *active_counts = (atomic_uint *)calloc(n_features * n_hidden, sizeof(atomic_uint));
    atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
    atomic_uint *feat_counts = (atomic_uint *)calloc(n_features, sizeof(atomic_uint));
    if (!active_counts || !label_counts || !feat_counts) {
      free(active_counts); free(label_counts); free(feat_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      return;
    }

    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    #pragma omp parallel for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      uint8_t *sample_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
          atomic_fetch_add(&feat_counts[f], 1);
          for (uint64_t b = 0; b < n_hidden; b++) {
            uint64_t b_byte = b / CHAR_BIT;
            uint8_t b_bit = b % CHAR_BIT;
            if (sample_codes[b_byte] & (1u << b_bit))
              atomic_fetch_add(&active_counts[f * n_hidden + b], 1);
          }
        }
      }
    }

    #pragma omp parallel for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint8_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos))
          atomic_fetch_add(&label_counts[b], 1);
      }
    }

    double N = (double)n_samples;
    tk_ivec_push(offsets, 0);
    for (uint64_t b = 0; b < n_hidden; b++) {
      double P = (double)atomic_load(&label_counts[b]);
      tk_rvec_t *heap = tk_rvec_create(NULL, 0, 0, 0);
      if (P > 0 && P < N) {
        for (uint64_t f = 0; f < n_features; f++) {
          double C = (double)atomic_load(&feat_counts[f]);
          if (C <= 0 || C >= N) continue;
          double A = (double)atomic_load(&active_counts[f * n_hidden + b]);
          if (A <= C * P / N) continue;
          double chi2 = tk_chi2_from_marginals(N, C, P, A);
          tk_rank_t r = { (int64_t)f, chi2 };
          tk_rvec_hmin(heap, per_class_k, r);
        }
      }
      if (heap->n > 0)
        tk_rvec_desc(heap, 0, heap->n);
      for (uint64_t i = 0; i < heap->n; i++) {
        tk_ivec_push(features, heap->a[i].i);
        tk_dvec_push(weights, heap->a[i].d);
      }
      tk_ivec_push(offsets, (int64_t)features->n);
      tk_rvec_destroy(heap);
    }

    free(active_counts); free(feat_counts); free(label_counts);

  } else if (labels) {
    atomic_uint *feat_counts = (atomic_uint *)calloc(n_features, sizeof(atomic_uint));
    atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
    if (!feat_counts || !label_counts) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      return;
    }

    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    #pragma omp parallel for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx))
          atomic_fetch_add(&feat_counts[f], 1);
      }
    }
    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0)
        atomic_fetch_add(&label_counts[(uint64_t)bit % n_hidden], 1);
    }

    int64_t *label_sample_offsets = tk_ivec_bits_sample_offsets(labels, n_samples, n_hidden);
    if (!label_sample_offsets) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      return;
    }

    int n_threads = 1;
    #pragma omp parallel
    { n_threads = omp_get_num_threads(); }
    tk_iumap_t **local_maps = (tk_iumap_t **)calloc((size_t)n_threads, sizeof(tk_iumap_t *));
    if (!local_maps) {
      free(label_sample_offsets);
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      return;
    }

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      local_maps[tid] = tk_iumap_create(0, 0);
      #pragma omp for schedule(static)
      for (uint64_t s = 0; s < n_samples; s++) {
        uint64_t sample_offset = s * bytes_per_sample;
        int64_t label_start = label_sample_offsets[s];
        int64_t label_end = label_sample_offsets[s + 1];
        for (uint64_t f = 0; f < n_features; f++) {
          uint64_t byte_idx = f / CHAR_BIT;
          uint8_t bit_idx = f % CHAR_BIT;
          if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
            for (int64_t li = label_start; li < label_end; li++) {
              uint64_t h = (uint64_t)labels->a[li] % n_hidden;
              tk_iumap_inc(local_maps[tid], (int64_t)(f * n_hidden + h));
            }
          }
        }
      }
    }
    free(label_sample_offsets);

    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    for (int t = 0; t < n_threads; t++) {
      if (local_maps[t]) {
        int64_t mk, mv;
        tk_umap_foreach(local_maps[t], mk, mv, ({
          khint_t kit = tk_iumap_get(active_counts, mk);
          if (kit != kh_end(active_counts))
            kh_value(active_counts, kit) += mv;
          else {
            int absent;
            khint_t newk = tk_iumap_put(active_counts, mk, &absent);
            kh_value(active_counts, newk) = mv;
          }
        }));
        tk_iumap_destroy(local_maps[t]);
      }
    }
    free(local_maps);

    double N = (double)n_samples;
    tk_rvec_t **class_heaps = (tk_rvec_t **)calloc(n_hidden, sizeof(tk_rvec_t *));
    #pragma omp parallel for schedule(static)
    for (uint64_t b = 0; b < n_hidden; b++) {
      class_heaps[b] = tk_rvec_create(NULL, 0, 0, 0);
      double P = (double)atomic_load(&label_counts[b]);
      if (P > 0 && P < N) {
        for (uint64_t f = 0; f < n_features; f++) {
          double C = (double)atomic_load(&feat_counts[f]);
          if (C <= 0 || C >= N) continue;
          double A = 0.0;
          khint_t kit = tk_iumap_get(active_counts, (int64_t)(f * n_hidden + b));
          if (kit != kh_end(active_counts))
            A = (double)kh_value(active_counts, kit);
          if (A <= C * P / N) continue;
          double chi2 = tk_chi2_from_marginals(N, C, P, A);
          tk_rank_t r = { (int64_t)f, chi2 };
          tk_rvec_hmin(class_heaps[b], per_class_k, r);
        }
      }
      if (class_heaps[b]->n > 0)
        tk_rvec_desc(class_heaps[b], 0, class_heaps[b]->n);
    }
    tk_ivec_push(offsets, 0);
    for (uint64_t b = 0; b < n_hidden; b++) {
      for (uint64_t i = 0; i < class_heaps[b]->n; i++) {
        tk_ivec_push(features, class_heaps[b]->a[i].i);
        tk_dvec_push(weights, class_heaps[b]->a[i].d);
      }
      tk_ivec_push(offsets, (int64_t)features->n);
      tk_rvec_destroy(class_heaps[b]);
    }
    free(class_heaps);

    tk_iumap_destroy(active_counts);
    free(feat_counts); free(label_counts);

  } else {
    for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
  }
  tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
  lua_insert(L, -5);
  lua_insert(L, -5);
}

static inline void tk_cvec_bits_top_mi (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  uint64_t per_class_k,
  uint64_t union_top_k,
  tk_pool_t pool
) {
  tk_ivec_t *offsets = tk_ivec_create(L, 0, 0, 0);
  tk_ivec_t *features = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);

  if (codes) {
    atomic_uint *active_counts = (atomic_uint *)calloc(n_features * n_hidden, sizeof(atomic_uint));
    atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
    atomic_uint *feat_counts = (atomic_uint *)calloc(n_features, sizeof(atomic_uint));
    if (!active_counts || !label_counts || !feat_counts) {
      free(active_counts); free(label_counts); free(feat_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      return;
    }

    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    #pragma omp parallel for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      uint8_t *sample_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
          atomic_fetch_add(&feat_counts[f], 1);
          for (uint64_t b = 0; b < n_hidden; b++) {
            uint64_t b_byte = b / CHAR_BIT;
            uint8_t b_bit = b % CHAR_BIT;
            if (sample_codes[b_byte] & (1u << b_bit))
              atomic_fetch_add(&active_counts[f * n_hidden + b], 1);
          }
        }
      }
    }

    #pragma omp parallel for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint8_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos))
          atomic_fetch_add(&label_counts[b], 1);
      }
    }

    double N = (double)n_samples;
    tk_ivec_push(offsets, 0);
    for (uint64_t b = 0; b < n_hidden; b++) {
      double P = (double)atomic_load(&label_counts[b]);
      tk_rvec_t *heap = tk_rvec_create(NULL, 0, 0, 0);
      if (P > 0 && P < N) {
        for (uint64_t f = 0; f < n_features; f++) {
          double C = (double)atomic_load(&feat_counts[f]);
          if (C <= 0 || C >= N) continue;
          double A = (double)atomic_load(&active_counts[f * n_hidden + b]);
          if (A <= C * P / N) continue;
          double mi = tk_mi_from_marginals(N, C, P, A);
          tk_rank_t r = { (int64_t)f, mi };
          tk_rvec_hmin(heap, per_class_k, r);
        }
      }
      if (heap->n > 0)
        tk_rvec_desc(heap, 0, heap->n);
      for (uint64_t i = 0; i < heap->n; i++) {
        tk_ivec_push(features, heap->a[i].i);
        tk_dvec_push(weights, heap->a[i].d);
      }
      tk_ivec_push(offsets, (int64_t)features->n);
      tk_rvec_destroy(heap);
    }

    free(active_counts); free(feat_counts); free(label_counts);

  } else if (labels) {
    atomic_uint *feat_counts = (atomic_uint *)calloc(n_features, sizeof(atomic_uint));
    atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
    if (!feat_counts || !label_counts) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      return;
    }

    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    #pragma omp parallel for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx))
          atomic_fetch_add(&feat_counts[f], 1);
      }
    }
    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0)
        atomic_fetch_add(&label_counts[(uint64_t)bit % n_hidden], 1);
    }

    int64_t *label_sample_offsets = tk_ivec_bits_sample_offsets(labels, n_samples, n_hidden);
    if (!label_sample_offsets) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      return;
    }

    int n_threads = 1;
    #pragma omp parallel
    { n_threads = omp_get_num_threads(); }
    tk_iumap_t **local_maps = (tk_iumap_t **)calloc((size_t)n_threads, sizeof(tk_iumap_t *));
    if (!local_maps) {
      free(label_sample_offsets);
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      return;
    }

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      local_maps[tid] = tk_iumap_create(0, 0);
      #pragma omp for schedule(static)
      for (uint64_t s = 0; s < n_samples; s++) {
        uint64_t sample_offset = s * bytes_per_sample;
        int64_t label_start = label_sample_offsets[s];
        int64_t label_end = label_sample_offsets[s + 1];
        for (uint64_t f = 0; f < n_features; f++) {
          uint64_t byte_idx = f / CHAR_BIT;
          uint8_t bit_idx = f % CHAR_BIT;
          if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
            for (int64_t li = label_start; li < label_end; li++) {
              uint64_t h = (uint64_t)labels->a[li] % n_hidden;
              tk_iumap_inc(local_maps[tid], (int64_t)(f * n_hidden + h));
            }
          }
        }
      }
    }
    free(label_sample_offsets);

    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    for (int t = 0; t < n_threads; t++) {
      if (local_maps[t]) {
        int64_t mk, mv;
        tk_umap_foreach(local_maps[t], mk, mv, ({
          khint_t kit = tk_iumap_get(active_counts, mk);
          if (kit != kh_end(active_counts))
            kh_value(active_counts, kit) += mv;
          else {
            int absent;
            khint_t newk = tk_iumap_put(active_counts, mk, &absent);
            kh_value(active_counts, newk) = mv;
          }
        }));
        tk_iumap_destroy(local_maps[t]);
      }
    }
    free(local_maps);

    double N = (double)n_samples;
    tk_rvec_t **class_heaps = (tk_rvec_t **)calloc(n_hidden, sizeof(tk_rvec_t *));
    #pragma omp parallel for schedule(static)
    for (uint64_t b = 0; b < n_hidden; b++) {
      class_heaps[b] = tk_rvec_create(NULL, 0, 0, 0);
      double P = (double)atomic_load(&label_counts[b]);
      if (P > 0 && P < N) {
        for (uint64_t f = 0; f < n_features; f++) {
          double C = (double)atomic_load(&feat_counts[f]);
          if (C <= 0 || C >= N) continue;
          double A = 0.0;
          khint_t kit = tk_iumap_get(active_counts, (int64_t)(f * n_hidden + b));
          if (kit != kh_end(active_counts))
            A = (double)kh_value(active_counts, kit);
          if (A <= C * P / N) continue;
          double mi = tk_mi_from_marginals(N, C, P, A);
          tk_rank_t r = { (int64_t)f, mi };
          tk_rvec_hmin(class_heaps[b], per_class_k, r);
        }
      }
      if (class_heaps[b]->n > 0)
        tk_rvec_desc(class_heaps[b], 0, class_heaps[b]->n);
    }
    tk_ivec_push(offsets, 0);
    for (uint64_t b = 0; b < n_hidden; b++) {
      for (uint64_t i = 0; i < class_heaps[b]->n; i++) {
        tk_ivec_push(features, class_heaps[b]->a[i].i);
        tk_dvec_push(weights, class_heaps[b]->a[i].d);
      }
      tk_ivec_push(offsets, (int64_t)features->n);
      tk_rvec_destroy(class_heaps[b]);
    }
    free(class_heaps);

    tk_iumap_destroy(active_counts);
    free(feat_counts); free(label_counts);

  } else {
    for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
  }
  tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
  lua_insert(L, -5);
  lua_insert(L, -5);
}

static inline void tk_cvec_bits_top_bns (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  uint64_t per_class_k,
  uint64_t union_top_k,
  tk_pool_t pool
) {
  tk_ivec_t *offsets = tk_ivec_create(L, 0, 0, 0);
  tk_ivec_t *features = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);

  if (codes) {
    atomic_uint *active_counts = (atomic_uint *)calloc(n_features * n_hidden, sizeof(atomic_uint));
    atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
    atomic_uint *feat_counts = (atomic_uint *)calloc(n_features, sizeof(atomic_uint));
    if (!active_counts || !label_counts || !feat_counts) {
      free(active_counts); free(label_counts); free(feat_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      return;
    }

    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    #pragma omp parallel for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      uint8_t *sample_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
          atomic_fetch_add(&feat_counts[f], 1);
          for (uint64_t b = 0; b < n_hidden; b++) {
            uint64_t b_byte = b / CHAR_BIT;
            uint8_t b_bit = b % CHAR_BIT;
            if (sample_codes[b_byte] & (1u << b_bit))
              atomic_fetch_add(&active_counts[f * n_hidden + b], 1);
          }
        }
      }
    }

    #pragma omp parallel for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint8_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos))
          atomic_fetch_add(&label_counts[b], 1);
      }
    }

    double N = (double)n_samples;
    tk_ivec_push(offsets, 0);
    for (uint64_t b = 0; b < n_hidden; b++) {
      double P = (double)atomic_load(&label_counts[b]);
      tk_rvec_t *heap = tk_rvec_create(NULL, 0, 0, 0);
      if (P > 0 && P < N) {
        for (uint64_t f = 0; f < n_features; f++) {
          double C = (double)atomic_load(&feat_counts[f]);
          if (C <= 0 || C >= N) continue;
          double A = (double)atomic_load(&active_counts[f * n_hidden + b]);
          if (A <= C * P / N) continue;
          double bns = tk_bns_from_marginals(N, C, P, A);
          tk_rank_t r = { (int64_t)f, bns };
          tk_rvec_hmin(heap, per_class_k, r);
        }
      }
      if (heap->n > 0)
        tk_rvec_desc(heap, 0, heap->n);
      for (uint64_t i = 0; i < heap->n; i++) {
        tk_ivec_push(features, heap->a[i].i);
        tk_dvec_push(weights, heap->a[i].d);
      }
      tk_ivec_push(offsets, (int64_t)features->n);
      tk_rvec_destroy(heap);
    }

    free(active_counts); free(feat_counts); free(label_counts);

  } else if (labels) {
    atomic_uint *feat_counts = (atomic_uint *)calloc(n_features, sizeof(atomic_uint));
    atomic_uint *label_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
    if (!feat_counts || !label_counts) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      return;
    }

    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    #pragma omp parallel for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx))
          atomic_fetch_add(&feat_counts[f], 1);
      }
    }
    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0)
        atomic_fetch_add(&label_counts[(uint64_t)bit % n_hidden], 1);
    }

    int64_t *label_sample_offsets = tk_ivec_bits_sample_offsets(labels, n_samples, n_hidden);
    if (!label_sample_offsets) {
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      return;
    }

    int n_threads = 1;
    #pragma omp parallel
    { n_threads = omp_get_num_threads(); }
    tk_iumap_t **local_maps = (tk_iumap_t **)calloc((size_t)n_threads, sizeof(tk_iumap_t *));
    if (!local_maps) {
      free(label_sample_offsets);
      free(feat_counts); free(label_counts);
      for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
      return;
    }

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      local_maps[tid] = tk_iumap_create(0, 0);
      #pragma omp for schedule(static)
      for (uint64_t s = 0; s < n_samples; s++) {
        uint64_t sample_offset = s * bytes_per_sample;
        int64_t label_start = label_sample_offsets[s];
        int64_t label_end = label_sample_offsets[s + 1];
        for (uint64_t f = 0; f < n_features; f++) {
          uint64_t byte_idx = f / CHAR_BIT;
          uint8_t bit_idx = f % CHAR_BIT;
          if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
            for (int64_t li = label_start; li < label_end; li++) {
              uint64_t h = (uint64_t)labels->a[li] % n_hidden;
              tk_iumap_inc(local_maps[tid], (int64_t)(f * n_hidden + h));
            }
          }
        }
      }
    }
    free(label_sample_offsets);

    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    for (int t = 0; t < n_threads; t++) {
      if (local_maps[t]) {
        int64_t mk, mv;
        tk_umap_foreach(local_maps[t], mk, mv, ({
          khint_t kit = tk_iumap_get(active_counts, mk);
          if (kit != kh_end(active_counts))
            kh_value(active_counts, kit) += mv;
          else {
            int absent;
            khint_t newk = tk_iumap_put(active_counts, mk, &absent);
            kh_value(active_counts, newk) = mv;
          }
        }));
        tk_iumap_destroy(local_maps[t]);
      }
    }
    free(local_maps);

    double N = (double)n_samples;
    tk_rvec_t **class_heaps = (tk_rvec_t **)calloc(n_hidden, sizeof(tk_rvec_t *));
    #pragma omp parallel for schedule(static)
    for (uint64_t b = 0; b < n_hidden; b++) {
      class_heaps[b] = tk_rvec_create(NULL, 0, 0, 0);
      double P = (double)atomic_load(&label_counts[b]);
      if (P > 0 && P < N) {
        for (uint64_t f = 0; f < n_features; f++) {
          double C = (double)atomic_load(&feat_counts[f]);
          if (C <= 0 || C >= N) continue;
          double A = 0.0;
          khint_t kit = tk_iumap_get(active_counts, (int64_t)(f * n_hidden + b));
          if (kit != kh_end(active_counts))
            A = (double)kh_value(active_counts, kit);
          if (A <= C * P / N) continue;
          double bns = tk_bns_from_marginals(N, C, P, A);
          tk_rank_t r = { (int64_t)f, bns };
          tk_rvec_hmin(class_heaps[b], per_class_k, r);
        }
      }
      if (class_heaps[b]->n > 0)
        tk_rvec_desc(class_heaps[b], 0, class_heaps[b]->n);
    }
    tk_ivec_push(offsets, 0);
    for (uint64_t b = 0; b < n_hidden; b++) {
      for (uint64_t i = 0; i < class_heaps[b]->n; i++) {
        tk_ivec_push(features, class_heaps[b]->a[i].i);
        tk_dvec_push(weights, class_heaps[b]->a[i].d);
      }
      tk_ivec_push(offsets, (int64_t)features->n);
      tk_rvec_destroy(class_heaps[b]);
    }
    free(class_heaps);

    tk_iumap_destroy(active_counts);
    free(feat_counts); free(label_counts);

  } else {
    for (uint64_t b = 0; b <= n_hidden; b++) tk_ivec_push(offsets, 0);
  }
  tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
  lua_insert(L, -5);
  lua_insert(L, -5);
}

static inline tk_ivec_t *tk_cvec_bits_top_entropy (
  lua_State *L,
  tk_cvec_t *codes,
  uint64_t n_samples,
  uint64_t n_hidden,
  uint64_t top_k
) {
  atomic_uint *bit_counts = (atomic_uint *)calloc(n_hidden, sizeof(atomic_uint));
  if (!bit_counts)
    return NULL;
  #pragma omp parallel for schedule(static)
  for (uint64_t s = 0; s < n_samples; s++) {
    uint8_t *sample_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
    for (uint64_t h = 0; h < n_hidden; h++) {
      uint64_t byte_idx = h / CHAR_BIT;
      uint8_t bit_idx = h % CHAR_BIT;
      if (sample_codes[byte_idx] & (1u << bit_idx))
        atomic_fetch_add(&bit_counts[h], 1);
    }
  }
  int n_threads = 1;
  #pragma omp parallel
  { n_threads = omp_get_num_threads(); }
  tk_rvec_t **local_heaps = (tk_rvec_t **)calloc((size_t)n_threads, sizeof(tk_rvec_t *));
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
    #pragma omp for schedule(static)
    for (uint64_t h = 0; h < n_hidden; h++) {
      double p = (double)atomic_load(&bit_counts[h]) / (double)n_samples;
      double entropy = 0.0;
      if (p > 0.0 && p < 1.0)
        entropy = -(p * log2(p) + (1.0 - p) * log2(1.0 - p));
      tk_rank_t r = { (int64_t)h, entropy };
      tk_rvec_hmin(local_heaps[tid], top_k, r);
    }
  }
  tk_rvec_t *top_heap = tk_rvec_create(NULL, 0, 0, 0);
  for (int t = 0; t < n_threads; t++) {
    for (uint64_t i = 0; i < local_heaps[t]->n; i++)
      tk_rvec_hmin(top_heap, top_k, local_heaps[t]->a[i]);
    tk_rvec_destroy(local_heaps[t]);
  }
  free(local_heaps);
  free(bit_counts);
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}

static inline tk_ivec_t *tk_ivec_bits_top_df (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_visible,
  double min_df,
  double max_df,
  uint64_t top_k
) {
  atomic_uint *df_counts = (atomic_uint *)calloc(n_visible, sizeof(atomic_uint));
  if (!df_counts)
    return NULL;
  uint64_t *last_sample = (uint64_t *)malloc(n_visible * sizeof(uint64_t));
  if (!last_sample) {
    free(df_counts);
    return NULL;
  }
  for (uint64_t i = 0; i < n_visible; i++)
    last_sample[i] = UINT64_MAX;
  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit_idx = set_bits->a[i];
    if (bit_idx < 0)
      continue;
    uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
    uint64_t feature_idx = (uint64_t)bit_idx % n_visible;
    if (sample_idx < n_samples && feature_idx < n_visible) {
      if (last_sample[feature_idx] != sample_idx) {
        last_sample[feature_idx] = sample_idx;
        atomic_fetch_add(&df_counts[feature_idx], 1);
      }
    }
  }
  free(last_sample);
  double min_df_abs = min_df < 0 ? -min_df : min_df * n_samples;
  double max_df_abs = max_df < 0 ? -max_df : max_df * n_samples;
  int n_threads = 1;
  #pragma omp parallel
  { n_threads = omp_get_num_threads(); }
  tk_rvec_t **local_heaps = (tk_rvec_t **)calloc((size_t)n_threads, sizeof(tk_rvec_t *));
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
    #pragma omp for schedule(static)
    for (uint64_t i = 0; i < n_visible; i++) {
      double df_count = (double)atomic_load(&df_counts[i]);
      if (df_count <= 0 || df_count >= n_samples) continue;
      double idf = log((double)(n_samples + 1) / (df_count + 1));
      if (df_count >= min_df_abs && df_count <= max_df_abs) {
        tk_rank_t r = { (int64_t)i, idf };
        tk_rvec_hmin(local_heaps[tid], top_k, r);
      }
    }
  }
  tk_rvec_t *top_heap = tk_rvec_create(NULL, 0, 0, 0);
  for (int t = 0; t < n_threads; t++) {
    for (uint64_t i = 0; i < local_heaps[t]->n; i++)
      tk_rvec_hmin(top_heap, top_k, local_heaps[t]->a[i]);
    tk_rvec_destroy(local_heaps[t]);
  }
  free(local_heaps);
  free(df_counts);
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}


static inline tk_ivec_t *tk_cvec_bits_top_df (
  lua_State *L,
  tk_cvec_t *bitmap,
  uint64_t n_samples,
  uint64_t n_features,
  double min_df,
  double max_df,
  uint64_t top_k
) {
  atomic_uint *df_counts = (atomic_uint *)calloc(n_features, sizeof(atomic_uint));
  if (!df_counts)
    return NULL;
  uint8_t *data = (uint8_t *)bitmap->a;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  #pragma omp parallel for schedule(static)
  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t sample_offset = s * bytes_per_sample;
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t byte_idx = f / CHAR_BIT;
      uint8_t bit_idx = f % CHAR_BIT;
      if (data[sample_offset + byte_idx] & (1u << bit_idx))
        atomic_fetch_add(&df_counts[f], 1);
    }
  }
  double min_df_abs = min_df < 0 ? -min_df : min_df * n_samples;
  double max_df_abs = max_df < 0 ? -max_df : max_df * n_samples;
  int n_threads = 1;
  #pragma omp parallel
  { n_threads = omp_get_num_threads(); }
  tk_rvec_t **local_heaps = (tk_rvec_t **)calloc((size_t)n_threads, sizeof(tk_rvec_t *));
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
    #pragma omp for schedule(static)
    for (uint64_t i = 0; i < n_features; i++) {
      double df_count = (double)atomic_load(&df_counts[i]);
      if (df_count <= 0 || df_count >= n_samples) continue;
      double idf = log((double)(n_samples + 1) / (df_count + 1));
      if (df_count >= min_df_abs && df_count <= max_df_abs) {
        tk_rank_t r = { (int64_t)i, idf };
        tk_rvec_hmin(local_heaps[tid], top_k, r);
      }
    }
  }
  tk_rvec_t *top_heap = tk_rvec_create(NULL, 0, 0, 0);
  for (int t = 0; t < n_threads; t++) {
    for (uint64_t i = 0; i < local_heaps[t]->n; i++)
      tk_rvec_hmin(top_heap, top_k, local_heaps[t]->a[i]);
    tk_rvec_destroy(local_heaps[t]);
  }
  free(local_heaps);
  free(df_counts);
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}

static inline void tk_ivec_bits_top_reg_f (
  lua_State *L,
  tk_ivec_t *set_bits,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  uint64_t per_class_k,
  uint64_t union_top_k,
  tk_pool_t pool
) {
  tk_ivec_t *offsets = tk_ivec_create(L, 0, 0, 0);
  tk_ivec_t *features = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);

  atomic_uint *feat_counts = (atomic_uint *)calloc(n_features, sizeof(atomic_uint));
  tk_dvec_t **feat_sums = (tk_dvec_t **)calloc(n_hidden, sizeof(tk_dvec_t *));
  tk_dvec_t **feat_sum_sq = (tk_dvec_t **)calloc(n_hidden, sizeof(tk_dvec_t *));
  tk_dvec_t *overall_sums = tk_dvec_create(0, n_hidden, 0, 0);
  tk_dvec_t *total_sum_sq = tk_dvec_create(0, n_hidden, 0, 0);

  if (!feat_counts || !feat_sums || !feat_sum_sq) {
    free(feat_counts); free(feat_sums); free(feat_sum_sq);
    for (uint64_t h = 0; h <= n_hidden; h++) tk_ivec_push(offsets, 0);
    tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
    lua_insert(L, -5);
    lua_insert(L, -5);
    return;
  }

  for (uint64_t h = 0; h < n_hidden; h++) {
    feat_sums[h] = tk_dvec_create(0, n_features, 0, 0);
    feat_sum_sq[h] = tk_dvec_create(0, n_features, 0, 0);
    tk_dvec_zero(feat_sums[h]);
    tk_dvec_zero(feat_sum_sq[h]);
  }
  tk_dvec_zero(overall_sums);
  tk_dvec_zero(total_sum_sq);

  #pragma omp parallel
  {
    double *local_overall = (double *)calloc(n_hidden, sizeof(double));
    double *local_total_sq = (double *)calloc(n_hidden, sizeof(double));
    #pragma omp for schedule(static)
    for (uint64_t s = 0; s < n_samples; s++) {
      for (uint64_t h = 0; h < n_hidden; h++) {
        double y = targets->a[s * n_hidden + h];
        local_overall[h] += y;
        local_total_sq[h] += y * y;
      }
    }
    #pragma omp critical
    {
      for (uint64_t h = 0; h < n_hidden; h++) {
        overall_sums->a[h] += local_overall[h];
        total_sum_sq->a[h] += local_total_sq[h];
      }
    }
    free(local_overall);
    free(local_total_sq);
  }

  #pragma omp parallel for schedule(static)
  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit = set_bits->a[i];
    if (bit < 0) continue;
    uint64_t sample = (uint64_t)bit / n_features;
    uint64_t feature = (uint64_t)bit % n_features;
    if (sample >= n_samples || feature >= n_features) continue;
    atomic_fetch_add(&feat_counts[feature], 1);
    for (uint64_t h = 0; h < n_hidden; h++) {
      double y = targets->a[sample * n_hidden + h];
      #pragma omp atomic
      feat_sums[h]->a[feature] += y;
      #pragma omp atomic
      feat_sum_sq[h]->a[feature] += y * y;
    }
  }

  tk_ivec_push(offsets, 0);
  for (uint64_t h = 0; h < n_hidden; h++) {
    double overall_mean = overall_sums->a[h] / (double)n_samples;
    int n_threads = 1;
    #pragma omp parallel
    { n_threads = omp_get_num_threads(); }
    tk_rvec_t **local_heaps = (tk_rvec_t **)calloc((size_t)n_threads, sizeof(tk_rvec_t *));

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      local_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
      #pragma omp for schedule(static)
      for (uint64_t f = 0; f < n_features; f++) {
        int64_t n1 = (int64_t)atomic_load(&feat_counts[f]);
        int64_t n0 = (int64_t)n_samples - n1;
        if (n1 <= 1 || n0 <= 1) continue;

        double sum1 = feat_sums[h]->a[f];
        double sum0 = overall_sums->a[h] - sum1;
        double mean1 = sum1 / (double)n1;
        double mean0 = sum0 / (double)n0;

        double ssb = (double)n1 * (mean1 - overall_mean) * (mean1 - overall_mean) +
                     (double)n0 * (mean0 - overall_mean) * (mean0 - overall_mean);

        double sum_sq1 = feat_sum_sq[h]->a[f];
        double sum_sq0 = total_sum_sq->a[h] - sum_sq1;
        double ssw1 = sum_sq1 - (sum1 * sum1) / (double)n1;
        double ssw0 = sum_sq0 - (sum0 * sum0) / (double)n0;
        double ssw = ssw1 + ssw0;

        if (ssw < 1e-12) continue;
        double F = ssb / (ssw / (double)(n_samples - 2));

        tk_rank_t r = { (int64_t)f, F };
        tk_rvec_hmin(local_heaps[tid], per_class_k, r);
      }
    }

    tk_rvec_t *merged = tk_rvec_create(NULL, 0, 0, 0);
    for (int t = 0; t < n_threads; t++) {
      for (uint64_t j = 0; j < local_heaps[t]->n; j++)
        tk_rvec_hmin(merged, per_class_k, local_heaps[t]->a[j]);
      tk_rvec_destroy(local_heaps[t]);
    }
    free(local_heaps);

    tk_rvec_desc(merged, 0, merged->n);
    for (uint64_t j = 0; j < merged->n; j++) {
      tk_ivec_push(features, merged->a[j].i);
      tk_dvec_push(weights, merged->a[j].d);
    }
    tk_ivec_push(offsets, (int64_t)features->n);
    tk_rvec_destroy(merged);
  }

  free(feat_counts);
  for (uint64_t h = 0; h < n_hidden; h++) {
    tk_dvec_destroy(feat_sums[h]);
    tk_dvec_destroy(feat_sum_sq[h]);
  }
  free(feat_sums);
  free(feat_sum_sq);
  tk_dvec_destroy(overall_sums);
  tk_dvec_destroy(total_sum_sq);
  tk_build_pooled_union(L, offsets, features, weights, n_hidden, union_top_k, pool);
  lua_insert(L, -5);
  lua_insert(L, -5);
}

static inline tk_ivec_t *tk_ivec_bits_top_reg_pearson (
  lua_State *L,
  tk_ivec_t *set_bits,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t top_k
) {
  double overall_sum = 0.0, overall_sum_sq = 0.0;
  for (uint64_t i = 0; i < n_samples; i++) {
    double y = targets->a[i];
    overall_sum += y;
    overall_sum_sq += y * y;
  }
  double overall_mean = overall_sum / (double)n_samples;
  double overall_var = (overall_sum_sq / (double)n_samples) - (overall_mean * overall_mean);
  double overall_std = sqrt(overall_var);
  if (overall_std < 1e-12) {
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_dvec_t *feat_sums = tk_dvec_create(0, n_features, 0, 0);
  tk_ivec_zero(feat_counts);
  tk_dvec_zero(feat_sums);

  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit = set_bits->a[i];
    if (bit < 0) continue;
    uint64_t sample = (uint64_t)bit / n_features;
    uint64_t feature = (uint64_t)bit % n_features;
    if (sample >= n_samples || feature >= n_features) continue;
    feat_counts->a[feature]++;
    feat_sums->a[feature] += targets->a[sample];
  }

  int n_threads = 1;
  #pragma omp parallel
  { n_threads = omp_get_num_threads(); }
  tk_rvec_t **local_heaps = (tk_rvec_t **)calloc((size_t)n_threads, sizeof(tk_rvec_t *));
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
    #pragma omp for schedule(static)
    for (uint64_t f = 0; f < n_features; f++) {
      int64_t n1 = feat_counts->a[f];
      int64_t n0 = (int64_t)n_samples - n1;
      if (n1 == 0 || n0 == 0) continue;

      double sum1 = feat_sums->a[f];
      double sum0 = overall_sum - sum1;
      double mean1 = sum1 / (double)n1;
      double mean0 = sum0 / (double)n0;

      double r = (mean1 - mean0) / overall_std * sqrt((double)n1 * (double)n0 / ((double)n_samples * (double)n_samples));
      double score = r < 0 ? -r : r;

      tk_rank_t rk = { (int64_t)f, score };
      tk_rvec_hmin(local_heaps[tid], top_k, rk);
    }
  }

  tk_rvec_t *top_heap = tk_rvec_create(NULL, 0, 0, 0);
  for (int t = 0; t < n_threads; t++) {
    for (uint64_t i = 0; i < local_heaps[t]->n; i++)
      tk_rvec_hmin(top_heap, top_k, local_heaps[t]->a[i]);
    tk_rvec_destroy(local_heaps[t]);
  }
  free(local_heaps);

  tk_ivec_destroy(feat_counts);
  tk_dvec_destroy(feat_sums);

  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);

  return out;
}

static inline tk_ivec_t *tk_ivec_bits_top_reg_mi (
  lua_State *L,
  tk_ivec_t *set_bits,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t top_k,
  uint64_t n_bins
) {
  if (n_bins == 0) n_bins = 10;

  double y_min = targets->a[0], y_max = targets->a[0];
  for (uint64_t i = 1; i < n_samples; i++) {
    if (targets->a[i] < y_min) y_min = targets->a[i];
    if (targets->a[i] > y_max) y_max = targets->a[i];
  }
  double y_range = y_max - y_min;
  if (y_range < 1e-12) {
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  tk_ivec_t *bin_assignments = tk_ivec_create(0, n_samples, 0, 0);
  tk_ivec_t *bin_counts = tk_ivec_create(0, n_bins, 0, 0);
  tk_ivec_zero(bin_counts);
  for (uint64_t i = 0; i < n_samples; i++) {
    double normalized = (targets->a[i] - y_min) / y_range;
    uint64_t bin = (uint64_t)(normalized * (double)(n_bins - 1) + 0.5);
    if (bin >= n_bins) bin = n_bins - 1;
    bin_assignments->a[i] = (int64_t)bin;
    bin_counts->a[bin]++;
  }

  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_ivec_t *joint_counts = tk_ivec_create(0, n_features * n_bins, 0, 0);
  tk_ivec_zero(feat_counts);
  tk_ivec_zero(joint_counts);

  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit = set_bits->a[i];
    if (bit < 0) continue;
    uint64_t sample = (uint64_t)bit / n_features;
    uint64_t feature = (uint64_t)bit % n_features;
    if (sample >= n_samples || feature >= n_features) continue;
    feat_counts->a[feature]++;
    uint64_t bin = (uint64_t)bin_assignments->a[sample];
    joint_counts->a[feature * n_bins + bin]++;
  }

  int n_threads = 1;
  #pragma omp parallel
  { n_threads = omp_get_num_threads(); }
  tk_rvec_t **local_heaps = (tk_rvec_t **)calloc((size_t)n_threads, sizeof(tk_rvec_t *));
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
    #pragma omp for schedule(static)
    for (uint64_t f = 0; f < n_features; f++) {
      int64_t n1 = feat_counts->a[f];
      int64_t n0 = (int64_t)n_samples - n1;
      if (n1 == 0 || n0 == 0) continue;

      double mi = 0.0;
      double p1 = (double)n1 / (double)n_samples;
      double p0 = (double)n0 / (double)n_samples;

      for (uint64_t b = 0; b < n_bins; b++) {
        int64_t n_b = bin_counts->a[b];
        if (n_b == 0) continue;
        double p_b = (double)n_b / (double)n_samples;

        int64_t n_1b = joint_counts->a[f * n_bins + b];
        int64_t n_0b = n_b - n_1b;

        if (n_1b > 0) {
          double p_1b = (double)n_1b / (double)n_samples;
          mi += p_1b * (log(p_1b) - log(p1) - log(p_b));
        }
        if (n_0b > 0) {
          double p_0b = (double)n_0b / (double)n_samples;
          mi += p_0b * (log(p_0b) - log(p0) - log(p_b));
        }
      }

      tk_rank_t r = { (int64_t)f, mi };
      tk_rvec_hmin(local_heaps[tid], top_k, r);
    }
  }

  tk_rvec_t *top_heap = tk_rvec_create(NULL, 0, 0, 0);
  for (int t = 0; t < n_threads; t++) {
    for (uint64_t i = 0; i < local_heaps[t]->n; i++)
      tk_rvec_hmin(top_heap, top_k, local_heaps[t]->a[i]);
    tk_rvec_destroy(local_heaps[t]);
  }
  free(local_heaps);

  tk_ivec_destroy(bin_assignments);
  tk_ivec_destroy(bin_counts);
  tk_ivec_destroy(feat_counts);
  tk_ivec_destroy(joint_counts);

  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);

  return out;
}


static inline tk_ivec_t *tk_cvec_bits_top_reg_pearson (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t top_k
) {
  double overall_sum = 0.0, overall_sum_sq = 0.0;
  for (uint64_t i = 0; i < n_samples; i++) {
    double y = targets->a[i];
    overall_sum += y;
    overall_sum_sq += y * y;
  }
  double overall_mean = overall_sum / (double)n_samples;
  double overall_var = (overall_sum_sq / (double)n_samples) - (overall_mean * overall_mean);
  double overall_std = sqrt(overall_var);
  if (overall_std < 1e-12) {
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_dvec_t *feat_sums = tk_dvec_create(0, n_features, 0, 0);
  tk_ivec_zero(feat_counts);
  tk_dvec_zero(feat_sums);

  uint8_t *bitmap_data = (uint8_t *)bitmap->a;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);

  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t sample_offset = s * bytes_per_sample;
    double y = targets->a[s];
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t byte_idx = f / CHAR_BIT;
      uint8_t bit_idx = f % CHAR_BIT;
      if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
        feat_counts->a[f]++;
        feat_sums->a[f] += y;
      }
    }
  }

  int n_threads = 1;
  #pragma omp parallel
  { n_threads = omp_get_num_threads(); }
  tk_rvec_t **local_heaps = (tk_rvec_t **)calloc((size_t)n_threads, sizeof(tk_rvec_t *));
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
    #pragma omp for schedule(static)
    for (uint64_t f = 0; f < n_features; f++) {
      int64_t n1 = feat_counts->a[f];
      int64_t n0 = (int64_t)n_samples - n1;
      if (n1 == 0 || n0 == 0) continue;

      double sum1 = feat_sums->a[f];
      double sum0 = overall_sum - sum1;
      double mean1 = sum1 / (double)n1;
      double mean0 = sum0 / (double)n0;

      double r = (mean1 - mean0) / overall_std * sqrt((double)n1 * (double)n0 / ((double)n_samples * (double)n_samples));
      double score = r < 0 ? -r : r;

      tk_rank_t rk = { (int64_t)f, score };
      tk_rvec_hmin(local_heaps[tid], top_k, rk);
    }
  }

  tk_rvec_t *top_heap = tk_rvec_create(NULL, 0, 0, 0);
  for (int t = 0; t < n_threads; t++) {
    for (uint64_t i = 0; i < local_heaps[t]->n; i++)
      tk_rvec_hmin(top_heap, top_k, local_heaps[t]->a[i]);
    tk_rvec_destroy(local_heaps[t]);
  }
  free(local_heaps);

  tk_ivec_destroy(feat_counts);
  tk_dvec_destroy(feat_sums);

  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);

  return out;
}

static inline tk_ivec_t *tk_cvec_bits_top_reg_mi (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t top_k,
  uint64_t n_bins
) {
  if (n_bins == 0) n_bins = 10;

  double y_min = targets->a[0], y_max = targets->a[0];
  for (uint64_t i = 1; i < n_samples; i++) {
    if (targets->a[i] < y_min) y_min = targets->a[i];
    if (targets->a[i] > y_max) y_max = targets->a[i];
  }
  double y_range = y_max - y_min;
  if (y_range < 1e-12) {
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  tk_ivec_t *bin_assignments = tk_ivec_create(0, n_samples, 0, 0);
  tk_ivec_t *bin_counts = tk_ivec_create(0, n_bins, 0, 0);
  tk_ivec_zero(bin_counts);
  for (uint64_t i = 0; i < n_samples; i++) {
    double normalized = (targets->a[i] - y_min) / y_range;
    uint64_t bin = (uint64_t)(normalized * (double)(n_bins - 1) + 0.5);
    if (bin >= n_bins) bin = n_bins - 1;
    bin_assignments->a[i] = (int64_t)bin;
    bin_counts->a[bin]++;
  }

  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_ivec_t *joint_counts = tk_ivec_create(0, n_features * n_bins, 0, 0);
  tk_ivec_zero(feat_counts);
  tk_ivec_zero(joint_counts);

  uint8_t *bitmap_data = (uint8_t *)bitmap->a;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);

  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t sample_offset = s * bytes_per_sample;
    uint64_t bin = (uint64_t)bin_assignments->a[s];
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t byte_idx = f / CHAR_BIT;
      uint8_t bit_idx = f % CHAR_BIT;
      if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
        feat_counts->a[f]++;
        joint_counts->a[f * n_bins + bin]++;
      }
    }
  }

  int n_threads = 1;
  #pragma omp parallel
  { n_threads = omp_get_num_threads(); }
  tk_rvec_t **local_heaps = (tk_rvec_t **)calloc((size_t)n_threads, sizeof(tk_rvec_t *));
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    local_heaps[tid] = tk_rvec_create(NULL, 0, 0, 0);
    #pragma omp for schedule(static)
    for (uint64_t f = 0; f < n_features; f++) {
      int64_t n1 = feat_counts->a[f];
      int64_t n0 = (int64_t)n_samples - n1;
      if (n1 == 0 || n0 == 0) continue;

      double mi = 0.0;
      double p1 = (double)n1 / (double)n_samples;
      double p0 = (double)n0 / (double)n_samples;

      for (uint64_t b = 0; b < n_bins; b++) {
        int64_t n_b = bin_counts->a[b];
        if (n_b == 0) continue;
        double p_b = (double)n_b / (double)n_samples;

        int64_t n_1b = joint_counts->a[f * n_bins + b];
        int64_t n_0b = n_b - n_1b;

        if (n_1b > 0) {
          double p_1b = (double)n_1b / (double)n_samples;
          mi += p_1b * (log(p_1b) - log(p1) - log(p_b));
        }
        if (n_0b > 0) {
          double p_0b = (double)n_0b / (double)n_samples;
          mi += p_0b * (log(p_0b) - log(p0) - log(p_b));
        }
      }

      tk_rank_t r = { (int64_t)f, mi };
      tk_rvec_hmin(local_heaps[tid], top_k, r);
    }
  }

  tk_rvec_t *top_heap = tk_rvec_create(NULL, 0, 0, 0);
  for (int t = 0; t < n_threads; t++) {
    for (uint64_t i = 0; i < local_heaps[t]->n; i++)
      tk_rvec_hmin(top_heap, top_k, local_heaps[t]->a[i]);
    tk_rvec_destroy(local_heaps[t]);
  }
  free(local_heaps);

  tk_ivec_destroy(bin_assignments);
  tk_ivec_destroy(bin_counts);
  tk_ivec_destroy(feat_counts);
  tk_ivec_destroy(joint_counts);

  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);

  return out;
}


#endif
