#ifndef TK_IVEC_EXT_H
#define TK_IVEC_EXT_H

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
#include <float.h>

#ifndef TK_CVEC_BITS_BYTES
#define TK_CVEC_BITS_BYTES(n) (((n) + CHAR_BIT - 1) / CHAR_BIT)
#endif

static inline tk_ivec_t *tk_ivec_bits_from_cvec(lua_State *L, const char *bm, uint64_t n_samples, uint64_t n_features);
static inline tk_ivec_t *tk_ivec_bits_from_cvec_serial(lua_State *L, const char *bm, uint64_t n_samples, uint64_t n_features);
static inline tk_ivec_t *tk_ivec_bits_extend(tk_ivec_t *base, tk_ivec_t *ext, uint64_t n_base_features, uint64_t n_ext_features);
static inline tk_ivec_t *tk_ivec_bits_extend_serial(tk_ivec_t *base, tk_ivec_t *ext, uint64_t n_base_features, uint64_t n_ext_features);
static inline int tk_ivec_bits_extend_mapped(tk_ivec_t *base, tk_ivec_t *ext, tk_ivec_t *aids, tk_ivec_t *bids, uint64_t n_base_features, uint64_t n_ext_features, bool project);
static inline int tk_ivec_bits_extend_mapped_serial(tk_ivec_t *base, tk_ivec_t *ext, tk_ivec_t *aids, tk_ivec_t *bids, uint64_t n_base_features, uint64_t n_ext_features, bool project);
static inline tk_ivec_t *tk_ivec_bits_bipartite(lua_State *L, tk_ivec_t *src, uint64_t n_docs, uint64_t n_labels, uint64_t *out_n_edges);
static inline tk_ivec_t *tk_ivec_bits_bipartite_serial(lua_State *L, tk_ivec_t *src, uint64_t n_docs, uint64_t n_labels, uint64_t *out_n_edges);

static inline void tk_ivec_copy_pkeys (tk_ivec_t *m0, tk_pvec_t *m1, int64_t start, int64_t end, int64_t dest) {
  if (start < 0 || start >= end || start >= (int64_t) m1->n)
    return;
  if (end >= (int64_t) m1->n)
    end = (int64_t) m1->n;
  uint64_t m = (uint64_t) dest + (uint64_t) (end - start);
  tk_ivec_ensure(m0, m);
  uint64_t write = m0->n;
  for (int64_t i = start; i < end; i ++)
    m0->a[write ++] = m1->a[i].i;
  if (m0->n < m)
    m0->n = m;
}

static inline void tk_ivec_copy_pvalues (tk_ivec_t *m0, tk_pvec_t *m1, int64_t start, int64_t end, int64_t dest) {
  if (start < 0 || start >= end || start >= (int64_t) m1->n)
    return;
  if (end >= (int64_t) m1->n)
    end = (int64_t) m1->n;
  uint64_t m = (uint64_t) dest + (uint64_t) (end - start);
  tk_ivec_ensure(m0, m);
  uint64_t write = m0->n;
  for (int64_t i = start; i < end; i ++)
    m0->a[write ++] = m1->a[i].p;
  if (m0->n < m)
    m0->n = m;
}

static inline void tk_ivec_copy_rkeys (tk_ivec_t *m0, tk_rvec_t *m1, int64_t start, int64_t end, int64_t dest) {
  if (start < 0 || start >= end || start >= (int64_t) m1->n)
    return;
  if (end >= (int64_t) m1->n)
    end = (int64_t) m1->n;
  uint64_t m = (uint64_t) dest + (uint64_t) (end - start);
  tk_ivec_ensure(m0, m);
  uint64_t write = m0->n;
  for (int64_t i = start; i < end; i ++)
    m0->a[write ++] = m1->a[i].i;
  if (m0->n < m)
    m0->n = m;
}

static inline void tk_ivec_copy_rvalues (tk_ivec_t *m0, tk_rvec_t *m1, int64_t start, int64_t end, int64_t dest) {
  if (start < 0 || start >= end || start >= (int64_t) m1->n)
    return;
  if (end >= (int64_t) m1->n)
    end = (int64_t) m1->n;
  uint64_t m = (uint64_t) dest + (uint64_t) (end - start);
  tk_ivec_ensure(m0, m);
  uint64_t write = m0->n;
  for (int64_t i = start; i < end; i ++)
    m0->a[write ++] = m1->a[i].d;
  if (m0->n < m)
    m0->n = m;
}

static inline tk_ivec_t *tk_ivec_from_rvec (lua_State *L, tk_rvec_t *R) {
  tk_ivec_t *result = tk_ivec_create(L, R->n, 0, 0);
  for (int64_t i = 0; i < (int64_t) R->n; i ++)
    result->a[i] = R->a[i].i;
  return result;
}

static inline tk_dvec_t *tk_ivec_to_dvec (lua_State *L, tk_ivec_t *v) {
  tk_dvec_t *out = tk_dvec_create(L, v->n, NULL, NULL);
  for (uint64_t i = 0; i < v->n; i++)
    out->a[i] = (double)v->a[i];
  return out;
}

static inline void tk_ivec_lookup (tk_ivec_t *indices, tk_ivec_t *source) {
  int64_t write_pos = 0;
  for (uint64_t i = 0; i < indices->n; i++) {
    int64_t idx = indices->a[i];
    if (idx >= 0 && idx < (int64_t) source->n) {
      indices->a[write_pos++] = source->a[idx];
    }
  }
  indices->n = (uint64_t) write_pos;
}

static inline tk_rvec_t *tk_rvec_rankings (lua_State *L, tk_dvec_t *scores, uint64_t n_visible, uint64_t n_hidden) {
  tk_rvec_t *rankings = tk_rvec_create(L, n_hidden * n_visible, NULL, NULL);
  for (uint64_t h = 0; h < n_hidden; h ++)
    for (uint64_t v = 0; v < n_visible; v ++)
      rankings->a[h * n_visible + v] = tk_rank((int64_t) v, scores->a[h * n_visible + v]);
  for (uint64_t j = 0; j < n_hidden; j ++)
    tk_rvec_desc(rankings, j * n_visible, (j + 1) * n_visible);
  return rankings;
}

static inline tk_cvec_t *tk_ivec_bits_to_cvec (lua_State *L, tk_ivec_t *set_bits, uint64_t n_samples, uint64_t n_features, bool flip_interleave) {
  uint64_t output_features = flip_interleave ? (n_features * 2) : n_features;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(output_features);
  size_t len = n_samples * bytes_per_sample;
  tk_cvec_t *out_cvec = tk_cvec_create(L, len, NULL, NULL);
  uint8_t *out = (uint8_t *)out_cvec->a;
  if (flip_interleave) {
    uint64_t second_half_bit_offset = n_features;
    uint64_t second_half_byte_offset = second_half_bit_offset / CHAR_BIT;
    uint8_t second_half_bit_shift = second_half_bit_offset & (CHAR_BIT - 1);
    uint64_t remaining_bits = n_features % CHAR_BIT;
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      memset(out + sample_offset, 0, bytes_per_sample);
      if (second_half_bit_shift == 0) {
        uint64_t full_bytes = n_features / CHAR_BIT;
        memset(out + sample_offset + second_half_byte_offset, 0xFF, full_bytes);
        if (remaining_bits > 0) {
          uint8_t mask = (1u << remaining_bits) - 1;
          out[sample_offset + second_half_byte_offset + full_bytes] = mask;
        }
      } else {
        for (uint64_t k = 0; k < n_features; k++) {
          uint64_t bit_pos = n_features + k;
          uint64_t byte_off = sample_offset + (bit_pos / CHAR_BIT);
          uint8_t bit_off = bit_pos & (CHAR_BIT - 1);
          out[byte_off] |= (1u << bit_off);
        }
      }
    }
    for (uint64_t idx = 0; idx < set_bits->n; idx++) {
      int64_t v = set_bits->a[idx];
      if (v < 0) continue;
      uint64_t s = (uint64_t) v / n_features;
      uint64_t k = (uint64_t) v % n_features;
      uint64_t sample_offset = s * bytes_per_sample;
      uint64_t first_byte_off = sample_offset + (k / CHAR_BIT);
      uint8_t first_bit_off = k & (CHAR_BIT - 1);
      out[first_byte_off] |= (1u << first_bit_off);
      uint64_t second_bit_pos = n_features + k;
      uint64_t second_byte_off = sample_offset + (second_bit_pos / CHAR_BIT);
      uint8_t second_bit_off = second_bit_pos & (CHAR_BIT - 1);
      out[second_byte_off] &= ~(1u << second_bit_off);
    }
  } else {
    memset(out, 0, len);
    for (uint64_t idx = 0; idx < set_bits->n; idx ++) {
      int64_t v = set_bits->a[idx];
      if (v < 0) continue;
      uint64_t s = (uint64_t) v / n_features;
      uint64_t f = (uint64_t) v % n_features;
      uint64_t byte_off = s * bytes_per_sample + (f / CHAR_BIT);
      uint8_t bit_off = f & (CHAR_BIT - 1);
      out[byte_off] |= (uint8_t) (1u << bit_off);
    }
  }
  return out_cvec;
}

static inline uint64_t tk_ivec_bits_to_cvec_grouped (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_visible,
  tk_ivec_t *offsets,
  tk_ivec_t *features,
  bool flip_interleave
) {
  uint64_t n_classes = offsets->n > 0 ? offsets->n - 1 : 0;
  if (n_classes == 0) {
    tk_cvec_create(L, 0, NULL, NULL);
    return 0;
  }
  uint64_t max_k = 0;
  for (uint64_t c = 0; c < n_classes; c++) {
    uint64_t class_k = (uint64_t)(offsets->a[c + 1] - offsets->a[c]);
    if (class_k > max_k) max_k = class_k;
  }
  if (max_k == 0) {
    tk_cvec_create(L, 0, NULL, NULL);
    return 0;
  }
  uint64_t bytes_per_class = flip_interleave ? TK_CVEC_BITS_BYTES(max_k * 2) : TK_CVEC_BITS_BYTES(max_k);
  uint64_t bytes_per_sample = n_classes * bytes_per_class;
  uint64_t bits_per_class = bytes_per_class * CHAR_BIT;
  size_t len = n_samples * bytes_per_sample;
  tk_cvec_t *out_cvec = tk_cvec_create(L, len, NULL, NULL);
  uint8_t *out = (uint8_t *)out_cvec->a;

  tk_iumap_t *feat_map = tk_iumap_create(0, 0);
  if (!feat_map) return max_k;
  for (uint64_t c = 0; c < n_classes; c++) {
    uint64_t start = (uint64_t)offsets->a[c];
    uint64_t end = (uint64_t)offsets->a[c + 1];
    for (uint64_t i = start; i < end; i++) {
      int64_t fid = features->a[i];
      if (fid < 0 || (uint64_t)fid >= n_visible) continue;
      uint64_t local = i - start;
      int64_t key = (int64_t)((uint64_t)fid * n_classes + c);
      int absent;
      khint_t kit = tk_iumap_put(feat_map, key, &absent);
      kh_value(feat_map, kit) = (int64_t)local;
    }
  }

  if (flip_interleave) {
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      memset(out + sample_offset, 0, bytes_per_sample);
      for (uint64_t c = 0; c < n_classes; c++) {
        for (uint64_t local = 0; local < max_k; local++) {
          uint64_t neg_bit = c * bits_per_class + max_k + local;
          uint64_t byte_off = sample_offset + (neg_bit / CHAR_BIT);
          uint8_t bit_off = neg_bit & (CHAR_BIT - 1);
          out[byte_off] |= (1u << bit_off);
        }
      }
    }
    for (uint64_t idx = 0; idx < set_bits->n; idx++) {
      int64_t v = set_bits->a[idx];
      if (v < 0) continue;
      uint64_t s = (uint64_t)v / n_visible;
      uint64_t f = (uint64_t)v % n_visible;
      if (s >= n_samples) continue;
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t c = 0; c < n_classes; c++) {
        int64_t key = (int64_t)(f * n_classes + c);
        khint_t kit = tk_iumap_get(feat_map, key);
        if (kit == kh_end(feat_map)) continue;
        int64_t local = kh_value(feat_map, kit);
        uint64_t pos_bit = c * bits_per_class + (uint64_t)local;
        uint64_t pos_byte = sample_offset + (pos_bit / CHAR_BIT);
        uint8_t pos_bit_off = pos_bit & (CHAR_BIT - 1);
        out[pos_byte] |= (1u << pos_bit_off);
        uint64_t neg_bit = c * bits_per_class + max_k + (uint64_t)local;
        uint64_t neg_byte = sample_offset + (neg_bit / CHAR_BIT);
        uint8_t neg_bit_off = neg_bit & (CHAR_BIT - 1);
        out[neg_byte] &= ~(1u << neg_bit_off);
      }
    }
  } else {
    memset(out, 0, len);
    for (uint64_t idx = 0; idx < set_bits->n; idx++) {
      int64_t v = set_bits->a[idx];
      if (v < 0) continue;
      uint64_t s = (uint64_t)v / n_visible;
      uint64_t f = (uint64_t)v % n_visible;
      if (s >= n_samples) continue;
      for (uint64_t c = 0; c < n_classes; c++) {
        int64_t key = (int64_t)(f * n_classes + c);
        khint_t kit = tk_iumap_get(feat_map, key);
        if (kit == kh_end(feat_map)) continue;
        int64_t local = kh_value(feat_map, kit);
        uint64_t out_idx = c * max_k + (uint64_t)local;
        uint64_t byte_off = s * bytes_per_sample + (out_idx / CHAR_BIT);
        uint8_t bit_off = out_idx & (CHAR_BIT - 1);
        out[byte_off] |= (1u << bit_off);
      }
    }
  }

  tk_iumap_destroy(feat_map);
  return max_k;
}

typedef enum {
  TK_IVEC_JACCARD,
  TK_IVEC_OVERLAP,
  TK_IVEC_DICE,
  TK_IVEC_TVERSKY,
  TK_IVEC_COSINE,
  TK_IVEC_MIN_KERNEL
} tk_ivec_sim_type_t;

static inline double tk_ivec_set_jaccard (double inter_w, double sum_a, double sum_b)
{
  double union_w = sum_a + sum_b - inter_w;
  return (union_w == 0.0) ? 0.0 : inter_w / union_w;
}

static inline double tk_ivec_set_overlap (double inter_w, double sum_a, double sum_b)
{
  double min_w = (sum_a < sum_b) ? sum_a : sum_b;
  return (min_w == 0.0) ? 0.0 : inter_w / min_w;
}

static inline double tk_ivec_set_dice (double inter_w, double sum_a, double sum_b)
{
  double denom = sum_a + sum_b;
  return (denom == 0.0) ? 0.0 : (2.0 * inter_w) / denom;
}

static inline double tk_ivec_set_tversky (double inter_w, double sum_a, double sum_b, double alpha, double beta)
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

static inline void tk_ivec_set_stats (
  int64_t *a, size_t alen,
  int64_t *b, size_t blen,
  tk_dvec_t *weights,
  double *inter_w,
  double *sum_a,
  double *sum_b
) {
  size_t i = 0, j = 0;
  double inter = 0.0, sa = 0.0, sb = 0.0;
  while (i < alen && j < blen) {
    int64_t ai = a[i], bj = b[j];
    if (ai == bj) {
      double w = (weights && ai >= 0 && ai < (int64_t)weights->n) ? weights->a[ai] : 1.0;
      inter += w;
      sa += w;
      sb += w;
      i++;
      j++;
    } else if (ai < bj) {
      double w = (weights && ai >= 0 && ai < (int64_t)weights->n) ? weights->a[ai] : 1.0;
      sa += w;
      i++;
    } else {
      double w = (weights && bj >= 0 && bj < (int64_t)weights->n) ? weights->a[bj] : 1.0;
      sb += w;
      j++;
    }
  }
  while (i < alen) {
    int64_t ai = a[i++];
    double w = (weights && ai >= 0 && ai < (int64_t)weights->n) ? weights->a[ai] : 1.0;
    sa += w;
  }
  while (j < blen) {
    int64_t bj = b[j++];
    double w = (weights && bj >= 0 && bj < (int64_t)weights->n) ? weights->a[bj] : 1.0;
    sb += w;
  }
  *inter_w = inter;
  *sum_a = sa;
  *sum_b = sb;
}

static inline double tk_ivec_set_similarity (
  int64_t *a, size_t alen,
  int64_t *b, size_t blen,
  tk_dvec_t *weights,
  tk_ivec_sim_type_t type,
  double tversky_alpha,
  double tversky_beta
) {
  double inter_w = 0.0, sum_a = 0.0, sum_b = 0.0;
  tk_ivec_set_stats(a, alen, b, blen, weights, &inter_w, &sum_a, &sum_b);
  switch (type) {
    case TK_IVEC_JACCARD:
      return tk_ivec_set_jaccard(inter_w, sum_a, sum_b);
    case TK_IVEC_OVERLAP:
      return tk_ivec_set_overlap(inter_w, sum_a, sum_b);
    case TK_IVEC_DICE:
      return tk_ivec_set_dice(inter_w, sum_a, sum_b);
    case TK_IVEC_TVERSKY:
      return tk_ivec_set_tversky(inter_w, sum_a, sum_b, tversky_alpha, tversky_beta);
    default:
      return tk_ivec_set_jaccard(inter_w, sum_a, sum_b);
  }
}

static inline double tk_ivec_set_similarity_from_partial (
  double inter_w,
  double q_w,
  double e_w,
  tk_ivec_sim_type_t type,
  double tversky_alpha,
  double tversky_beta
) {
  switch (type) {
    case TK_IVEC_JACCARD: {
      double uni_w = q_w + e_w - inter_w;
      return (uni_w == 0.0) ? 0.0 : inter_w / uni_w;
    }
    case TK_IVEC_OVERLAP: {
      double min_w = (q_w < e_w) ? q_w : e_w;
      return (min_w == 0.0) ? 0.0 : inter_w / min_w;
    }
    case TK_IVEC_TVERSKY: {
      double a_only = q_w - inter_w;
      double b_only = e_w - inter_w;
      if (a_only < 0.0) a_only = 0.0;
      if (b_only < 0.0) b_only = 0.0;
      double denom = inter_w + tversky_alpha * a_only + tversky_beta * b_only;
      return (denom == 0.0) ? 0.0 : inter_w / denom;
    }
    case TK_IVEC_DICE: {
      double denom = q_w + e_w;
      return (denom == 0.0) ? 0.0 : (2.0 * inter_w) / denom;
    }
    case TK_IVEC_COSINE: {
      double denom = sqrt(q_w * e_w);
      return (denom == 0.0) ? 0.0 : inter_w / denom;
    }
    case TK_IVEC_MIN_KERNEL: {
      return inter_w;
    }
    default: {
      double uni_w = q_w + e_w - inter_w;
      return (uni_w == 0.0) ? 0.0 : inter_w / uni_w;
    }
  }
}

static inline void tk_ivec_set_weights_by_rank (
  int64_t *features,
  size_t n_features,
  tk_dvec_t *weights,
  tk_ivec_t *ranks,
  uint64_t n_ranks,
  double *weights_by_rank
) {
  for (uint64_t r = 0; r < n_ranks; r ++)
    weights_by_rank[r] = 0.0;
  for (size_t i = 0; i < n_features; i ++) {
    int64_t fid = features[i];
    if (fid >= 0) {
      double w = (weights && fid < (int64_t)weights->n) ? weights->a[fid] : 1.0;
      int64_t rank = (ranks && fid < (int64_t)ranks->n) ? ranks->a[fid] : 0;
      if (rank >= 0 && rank < (int64_t)n_ranks) {
        weights_by_rank[rank] += w;
      }
    }
  }
}

static inline double tk_ivec_set_similarity_by_rank (
  tk_dvec_t *wacc,
  int64_t vsid,
  double *q_weights_by_rank,
  double *e_weights_by_rank,
  uint64_t n_ranks,
  int64_t rank_decay_window,
  double rank_decay_sigma,
  double rank_decay_floor,
  tk_ivec_sim_type_t type,
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
      rank_sim = tk_ivec_set_similarity_from_partial(inter_w, q_w, e_w, type, tversky_alpha, tversky_beta);
    } else {
      rank_sim = 0.0;
    }
    total_weighted_sim += rank_sim * rank_weight;
    total_rank_weight += rank_weight;
  }
  return (total_rank_weight > 0.0) ? total_weighted_sim / total_rank_weight : 0.0;
}

static inline int64_t tk_ivec_set_find (int64_t *arr, int64_t start, int64_t end, int64_t value) {
  int64_t left = start;
  int64_t right = end - 1;
  while (left <= right) {
    int64_t mid = left + (right - left) / 2;
    if (arr[mid] == value)
      return mid;
    if (arr[mid] < value)
      left = mid + 1;
    else
      right = mid - 1;
  }
  return -(left + 1);
}
static inline void tk_ivec_set_insert (tk_ivec_t *vec, int64_t pos, int64_t value) {
  if (pos < 0) pos = -(pos + 1);
  if (pos < 0) pos = 0;
  if (pos > (int64_t)vec->n) pos = (int64_t)vec->n;
  tk_ivec_ensure(vec, vec->n + 1);
  if (pos < (int64_t)vec->n) {
    memmove(vec->a + pos + 1, vec->a + pos, (vec->n - (size_t) pos) * sizeof(int64_t));
  }
  vec->a[pos] = value;
  vec->n++;
}

static inline tk_ivec_t *tk_ivec_set_intersect (lua_State *L, tk_ivec_t *a, tk_ivec_t *b, tk_ivec_t *out) {
  if (out == NULL) {
    size_t min_size = a->n < b->n ? a->n : b->n;
    out = tk_ivec_create(L, min_size, 0, 0);
  } else {
    tk_ivec_clear(out);
    size_t min_size = a->n < b->n ? a->n : b->n;
    tk_ivec_ensure(out, min_size);
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
  tk_ivec_shrink(out);
  return out;
}

static inline tk_ivec_t *tk_ivec_set_union (lua_State *L, tk_ivec_t *a, tk_ivec_t *b, tk_ivec_t *out) {
  if (out == NULL) {
    out = tk_ivec_create(L, a->n + b->n, 0, 0);
  } else {
    tk_ivec_clear(out);
    tk_ivec_ensure(out, a->n + b->n);
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
  tk_ivec_shrink(out);
  return out;
}




static inline size_t tk_ivec_scores_max_curvature (
  tk_ivec_t *v,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? v->a[0] : 0;
    return n > 0 ? n - 1 : 0;
  }
  int64_t max_curv = 0;
  size_t max_idx = 1;
  for (size_t i = 1; i < n - 1; i++) {
    int64_t curv = llabs(v->a[i-1] - 2 * v->a[i] + v->a[i+1]);
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

static inline size_t tk_ivec_scores_lmethod (
  tk_ivec_t *v,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? v->a[0] : 0;
    return n > 0 ? n - 1 : 0;
  }

  int64_t min_val = v->a[0], max_val = v->a[0];
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

static inline size_t tk_ivec_scores_max_gap (
  tk_ivec_t *v,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? v->a[0] : 0;
    return n > 0 ? n - 1 : 0;
  }
  int64_t max_gap = 0;
  size_t max_idx = 0;
  for (size_t i = 0; i < n - 1; i++) {
    int64_t gap = llabs(v->a[i + 1] - v->a[i]);
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

static inline size_t tk_ivec_scores_plateau (
  tk_ivec_t *v,
  double tolerance,
  int64_t *out_val
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

  int64_t min_score = v->a[0], max_score = v->a[0];
  for (size_t i = 1; i < n; i++) {
    if (v->a[i] < min_score) min_score = v->a[i];
    if (v->a[i] > max_score) max_score = v->a[i];
  }
  int64_t range = max_score - min_score;
  if (range <= 0) {
    if (out_val) *out_val = v->a[n - 1];
    return n - 1;
  }

  int64_t abs_tolerance = (int64_t)(tolerance * range);
  int64_t base = v->a[0];
  size_t end_idx = 0;
  for (size_t i = 1; i < n; i++) {
    if (llabs(v->a[i] - base) <= abs_tolerance) {
      end_idx = i;
    } else {
      break;
    }
  }
  if (out_val) *out_val = v->a[end_idx];
  return end_idx;
}

static inline size_t tk_ivec_scores_kneedle (
  tk_ivec_t *v,
  double sensitivity,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? v->a[0] : 0;
    return n > 0 ? n - 1 : 0;
  }
  if (sensitivity <= 0.0)
    sensitivity = 1.0;
  int64_t min_score = v->a[0];
  int64_t max_score = v->a[0];
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



static inline size_t tk_ivec_scores_first_gap (
  tk_ivec_t *v,
  double alpha,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? v->a[0] : 0;
    return n > 0 ? n - 1 : 0;
  }
  if (alpha <= 0.0) alpha = 3.0;

  size_t n_gaps = n - 1;
  int64_t *gaps = (int64_t *)malloc(n_gaps * sizeof(int64_t));
  if (!gaps) {
    if (out_val) *out_val = v->a[n - 1];
    return n - 1;
  }

  for (size_t i = 0; i < n_gaps; i++) {
    gaps[i] = llabs(v->a[i + 1] - v->a[i]);
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
      int64_t gap = llabs(v->a[i + 1] - v->a[i]);
      if (gap > max_gap) {
        max_gap = gap;
        max_idx = i;
      }
    }
    if (out_val) *out_val = v->a[max_idx];
    return max_idx;
  }

  int64_t threshold = (int64_t)(alpha * (double)median_gap);
  if (threshold < 1) threshold = 1;

  for (size_t i = 0; i < n - 1; i++) {
    int64_t gap = llabs(v->a[i + 1] - v->a[i]);
    if (gap > threshold) {
      if (out_val) *out_val = v->a[i];
      return i;
    }
  }

  if (out_val) *out_val = v->a[n - 1];
  return n - 1;
}



static inline size_t tk_ivec_scores_otsu (
  tk_ivec_t *v,
  int64_t *out_val
) {
  size_t n = v->n;
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? v->a[0] : 0;
    return n > 0 ? n - 1 : 0;
  }


  int64_t min_val = v->a[0], max_val = v->a[0];
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


#define TK_GENERATE_SINGLE
#include <santoku/parallel/tpl.h>
#include <santoku/ivec/ext_tpl.h>
#undef TK_GENERATE_SINGLE

#include <santoku/parallel/tpl.h>
#include <santoku/ivec/ext_tpl.h>

#endif
