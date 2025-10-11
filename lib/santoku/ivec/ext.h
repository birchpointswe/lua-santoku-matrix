#ifndef TK_IVEC_EXT_H
#define TK_IVEC_EXT_H

#include <santoku/rvec.h>
#include <santoku/threads.h>
#include <santoku/dvec.h>
#include <santoku/iumap.h>
#include <stdatomic.h>
#include <math.h>
#include <limits.h>

#ifndef TK_CVEC_BITS_BYTES
#define TK_CVEC_BITS_BYTES(n) (((n) + CHAR_BIT - 1) / CHAR_BIT)
#endif

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
  tk_ivec_t *I = tk_ivec_create(L, R->n, 0, 0);
  for (int64_t i = 0; i < (int64_t) R->n; i ++)
    I->a[i] = R->a[i].i;
  return I;
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
static inline tk_ivec_t *tk_ivec_bits_from_cvec (lua_State *L, const char *bm, uint64_t n_samples, uint64_t n_features) {
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  for (uint64_t i = 0; i < n_samples; i ++)
    for (uint64_t j = 0; j < n_features; j ++) {
      uint64_t bit = i * n_features + j;
      uint64_t chunk = bit / CHAR_BIT;
      uint64_t pos = bit % CHAR_BIT;
      if (bm[chunk] & (1 << pos)) {
        if (tk_ivec_push(out, (int64_t) bit) != 0) {
          tk_ivec_destroy(out);
          return NULL;
        }
      }
    }
  tk_ivec_shrink(out);
  return out;
}

static inline void tk_ivec_bits_extend (tk_ivec_t *base, tk_ivec_t *ext, uint64_t n_feat, uint64_t n_extfeat) {
  tk_ivec_asc(base, 0, base->n);
  tk_ivec_asc(ext, 0, ext->n);
  size_t total = base->n + ext->n;
  tk_ivec_ensure(base, total);
  for (size_t i = 0; i < base->n; i ++) {
    uint64_t bit = (uint64_t) base->a[i];
    uint64_t sample = bit / n_feat;
    uint64_t old_off = sample * n_feat;
    uint64_t new_off = sample * (n_feat + n_extfeat);
    base->a[i] = (int64_t) (bit - old_off + new_off);
  }
  for (size_t i = 0; i < ext->n; i ++) {
    uint64_t bit = (uint64_t) ext->a[i];
    uint64_t sample = bit / n_extfeat;
    uint64_t old_off = sample * n_extfeat;
    uint64_t new_off = sample * (n_feat + n_extfeat);
    base->a[base->n + i] = (int64_t) (bit - old_off + new_off + n_feat);
  }
  base->n = total;
  tk_ivec_asc(base, 0, base->n);
}

static inline int tk_ivec_bits_extend_mapped (tk_ivec_t *base, tk_ivec_t *ext, tk_ivec_t *aids, tk_ivec_t *bids, uint64_t n_feat, uint64_t n_extfeat, bool project) {
  tk_ivec_asc(base, 0, base->n);
  tk_ivec_asc(ext, 0, ext->n);
  tk_iumap_t *a_id_to_pos = tk_iumap_from_ivec(0, aids);
  if (!a_id_to_pos)
    return -1;
  uint64_t n_only_b = 0;
  int64_t *b_to_final = (int64_t *)calloc(bids->n, sizeof(int64_t));
  int64_t *a_to_final = (int64_t *)calloc(aids->n, sizeof(int64_t));

  if (!b_to_final || !a_to_final) {
    free(b_to_final);
    free(a_to_final);
    tk_iumap_destroy(a_id_to_pos);
    return -1;
  }
  for (size_t i = 0; i < aids->n; i++)
    a_to_final[i] = (int64_t)i;
  int64_t next_pos = (int64_t)aids->n;
  for (size_t bi = 0; bi < bids->n; bi++) {
    khint_t khi = tk_iumap_get(a_id_to_pos, bids->a[bi]);
    if (khi != tk_iumap_end(a_id_to_pos)) {
      b_to_final[bi] = tk_iumap_val(a_id_to_pos, khi);
    } else {
      if (!project) {
        b_to_final[bi] = next_pos++;
        n_only_b++;
      } else {
        b_to_final[bi] = -1;
      }
    }
  }
  uint64_t final_n_samples = project ? aids->n : (aids->n + n_only_b);
  uint64_t n_total_feat = n_feat + n_extfeat;
  size_t old_aids_n = aids->n;
  if (!project) {
    if (tk_ivec_ensure(aids, final_n_samples) != 0) {
      free(b_to_final);
      free(a_to_final);
      tk_iumap_destroy(a_id_to_pos);
      return -1;
    }
    for (size_t i = 0; i < bids->n; i++) {
      if (b_to_final[i] >= (int64_t)old_aids_n)
        aids->a[aids->n++] = bids->a[i];
    }
  }
  size_t max_bits = 0;
  for (size_t i = 0; i < base->n; i++) {
    uint64_t bit = (uint64_t)base->a[i];
    uint64_t sample = bit / n_feat;
    if (sample < old_aids_n)
      max_bits++;
  }
  for (size_t i = 0; i < ext->n; i++) {
    uint64_t bit = (uint64_t)ext->a[i];
    uint64_t sample = bit / n_extfeat;
    if (sample < bids->n && b_to_final[sample] >= 0)
      max_bits++;
  }
  if (tk_ivec_ensure(base, max_bits) != 0) {
    free(b_to_final);
    free(a_to_final);
    tk_iumap_destroy(a_id_to_pos);
    return -1;
  }
  for (int64_t i = (int64_t)base->n - 1; i >= 0; i--) {
    uint64_t bit = (uint64_t)base->a[i];
    uint64_t sample = bit / n_feat;
    uint64_t feature = bit % n_feat;
    if (sample < old_aids_n) {
      uint64_t new_sample_pos = (uint64_t)a_to_final[sample];
      base->a[i] = (int64_t)(new_sample_pos * n_total_feat + feature);
    }
  }
  size_t insert_pos = base->n;
  for (size_t i = 0; i < ext->n; i++) {
    uint64_t bit = (uint64_t)ext->a[i];
    uint64_t b_sample = bit / n_extfeat;
    uint64_t feature = bit % n_extfeat;
    if (b_sample < bids->n && b_to_final[b_sample] >= 0) {
      uint64_t final_sample = (uint64_t)b_to_final[b_sample];
      base->a[insert_pos++] = (int64_t)(final_sample * n_total_feat + n_feat + feature);
    }
  }
  base->n = insert_pos;
  tk_ivec_asc(base, 0, base->n);
  free(b_to_final);
  free(a_to_final);
  tk_iumap_destroy(a_id_to_pos);
  return 0;
}
typedef enum { TK_IVEC_JACCARD, TK_IVEC_OVERLAP, TK_IVEC_DICE, TK_IVEC_TVERSKY } tk_ivec_sim_type_t;

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

#endif
