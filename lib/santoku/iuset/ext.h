#ifndef TK_IUSET_EXT_H
#define TK_IUSET_EXT_H

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

#include <santoku/iumap.h>

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
  uint64_t intersection = 0;
  uint64_t union_count = 0;
  int64_t x;
  tk_umap_foreach_keys(a, x, ({
    intersection ++;
  }));
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
  for (uint64_t i = 0; i < v->n; i ++)
    tk_iuset_put(s, v->a[i], &kha);
  return s;
}

static inline int tk_ivec_bits_filter (
  tk_ivec_t *set_bits,
  tk_ivec_t *top_v,
  tk_ivec_t *sample_ids,
  uint64_t n_visible
) {
  if ((top_v == NULL || top_v->n == 0) && (sample_ids == NULL || sample_ids->n == 0))
    return 0;

  int64_t *feature_map = NULL;
  uint64_t n_new_features = n_visible;
  if (top_v != NULL && top_v->n > 0) {
    feature_map = malloc(n_visible * sizeof(int64_t));
    if (!feature_map) return -1;
    for (uint64_t i = 0; i < n_visible; i ++)
      feature_map[i] = -1;
    for (uint64_t i = 0; i < top_v->n; i ++) {
      int64_t feat = top_v->a[i];
      if (feat >= 0 && (uint64_t) feat < n_visible)
        feature_map[feat] = (int64_t) i;
    }
    n_new_features = top_v->n;
  }

  tk_iuset_t *sample_set = NULL;
  int64_t *sample_map = NULL;
  uint64_t max_sample = 0;
  if (sample_ids != NULL && sample_ids->n > 0) {
    for (uint64_t i = 0; i < set_bits->n; i ++) {
      if (set_bits->a[i] >= 0) {
        uint64_t s = (uint64_t) set_bits->a[i] / n_visible;
        if (s > max_sample) max_sample = s;
      }
    }
    sample_set = tk_iuset_from_ivec(NULL, sample_ids);
    sample_map = malloc((max_sample + 1) * sizeof(int64_t));
    if (!sample_map) {
      if (feature_map) free(feature_map);
      tk_iuset_destroy(sample_set);
      return -1;
    }
    for (uint64_t i = 0; i <= max_sample; i ++)
      sample_map[i] = -1;
    uint64_t new_idx = 0;
    for (uint64_t i = 0; i < sample_ids->n; i ++) {
      int64_t sid = sample_ids->a[i];
      if (sid >= 0 && (uint64_t) sid <= max_sample)
        sample_map[sid] = (int64_t) new_idx ++;
    }
  }

  size_t write = 0;
  for (size_t i = 0; i < set_bits->n; i ++) {
    int64_t val = set_bits->a[i];
    if (val < 0)
      continue;
    uint64_t sample = (uint64_t) val / n_visible;
    uint64_t feature = (uint64_t) val % n_visible;
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
    set_bits->a[write ++] = new_sample * (int64_t) n_new_features + new_feature;
  }
  set_bits->n = write;

  if (feature_map != NULL)
    free(feature_map);
  if (sample_map != NULL)
    free(sample_map);
  if (sample_set != NULL)
    tk_iuset_destroy(sample_set);
  return 0;
}

static inline void tk_cvec_bits_filter (
  tk_cvec_t *bitmap,
  tk_ivec_t *selected_features,
  tk_ivec_t *sample_ids,
  uint64_t n_features
) {
  uint64_t n_samples = bitmap->n / TK_CVEC_BITS_BYTES(n_features);
  if ((selected_features == NULL || selected_features->n == 0) &&
    (sample_ids == NULL || sample_ids->n == 0))
    return;

  tk_iuset_t *sample_set = NULL;
  if (sample_ids != NULL && sample_ids->n > 0) {
    sample_set = tk_iuset_from_ivec(NULL, sample_ids);
  }

  uint64_t n_selected_features = (selected_features != NULL && selected_features->n > 0)
    ? selected_features->n : n_features;
  uint64_t in_bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t out_bytes_per_sample = TK_CVEC_BITS_BYTES(n_selected_features);
  uint8_t *data = (uint8_t *)bitmap->a;

  uint64_t n_output_samples = 0;
  if (sample_set != NULL) {
    for (uint64_t s = 0; s < n_samples; s ++) {
      if (tk_iuset_contains(sample_set, (int64_t) s))
        n_output_samples ++;
    }
  } else {
    n_output_samples = n_samples;
  }

  uint64_t write_sample = 0;
  for (uint64_t s = 0; s < n_samples; s ++) {
    if (sample_set != NULL && !tk_iuset_contains(sample_set, (int64_t) s))
      continue;
    uint8_t *temp = malloc(out_bytes_per_sample);
    if (!temp) {
      if (sample_set != NULL) tk_iuset_destroy(sample_set);
      return;
    }
    memset(temp, 0, out_bytes_per_sample);
    uint64_t in_offset = s * in_bytes_per_sample;
    if (selected_features != NULL && selected_features->n > 0) {
      for (uint64_t i = 0; i < selected_features->n; i ++) {
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
      memcpy(temp, data + in_offset, in_bytes_per_sample);
    }
    uint64_t out_offset = write_sample * out_bytes_per_sample;
    memcpy(data + out_offset, temp, out_bytes_per_sample);
    free(temp);
    write_sample ++;
  }

  bitmap->n = n_output_samples * out_bytes_per_sample;

  if (sample_set != NULL)
    tk_iuset_destroy(sample_set);
}

static inline void tk_ivec_bits_copy (
  tk_ivec_t *dest,
  tk_ivec_t *src_bits,
  tk_ivec_t *selected_features,
  tk_ivec_t *sample_ids,
  uint64_t n_visible,
  uint64_t dest_sample
) {
  if (dest_sample == 0)
    tk_ivec_clear(dest);
  if (src_bits == NULL || src_bits->n == 0)
    return;

  int64_t *feature_map = NULL;
  uint64_t n_new_features = n_visible;
  if (selected_features != NULL && selected_features->n > 0) {
    feature_map = malloc(n_visible * sizeof(int64_t));
    if (!feature_map)
      return;
    for (uint64_t i = 0; i < n_visible; i ++)
      feature_map[i] = -1;
    for (uint64_t i = 0; i < selected_features->n; i ++) {
      int64_t feat = selected_features->a[i];
      if (feat >= 0 && (uint64_t) feat < n_visible)
        feature_map[feat] = (int64_t) i;
    }
    n_new_features = selected_features->n;
  }

  tk_iuset_t *sample_set = NULL;
  int64_t *sample_map = NULL;
  uint64_t max_sample = 0;
  if (sample_ids != NULL && sample_ids->n > 0) {
    sample_set = tk_iuset_from_ivec(NULL, sample_ids);
    for (uint64_t i = 0; i < src_bits->n; i ++) {
      if (src_bits->a[i] >= 0) {
        uint64_t s = (uint64_t) src_bits->a[i] / n_visible;
        if (s > max_sample) max_sample = s;
      }
    }
    sample_map = malloc((max_sample + 1) * sizeof(int64_t));
    if (!sample_map) {
      if (feature_map) free(feature_map);
      tk_iuset_destroy(sample_set);
      return;
    }
    for (uint64_t i = 0; i <= max_sample; i ++)
      sample_map[i] = -1;
    uint64_t new_idx = 0;
    for (uint64_t i = 0; i < sample_ids->n; i ++) {
      int64_t sid = sample_ids->a[i];
      if (sid >= 0 && (uint64_t) sid <= max_sample)
        sample_map[sid] = (int64_t) new_idx ++;
    }
  }

  for (size_t i = 0; i < src_bits->n; i ++) {
    int64_t val = src_bits->a[i];
    if (val < 0)
      continue;
    uint64_t sample = (uint64_t) val / n_visible;
    uint64_t feature = (uint64_t) val % n_visible;
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
    int64_t new_val = new_sample * (int64_t) n_new_features + new_feature;
    tk_ivec_push(dest, new_val);
  }

  if (feature_map != NULL)
    free(feature_map);
  if (sample_map != NULL)
    free(sample_map);
  if (sample_set != NULL)
    tk_iuset_destroy(sample_set);

  tk_ivec_shrink(dest);
}

static inline void tk_cvec_bits_copy (
  tk_cvec_t *dest,
  tk_cvec_t *src_bitmap,
  tk_ivec_t *selected_features,
  tk_ivec_t *sample_ids,
  uint64_t n_features,
  uint64_t dest_sample,
  uint64_t dest_stride
) {
  uint64_t n_samples = src_bitmap->n / TK_CVEC_BITS_BYTES(n_features);
  tk_iuset_t *sample_set = NULL;
  if (sample_ids != NULL && sample_ids->n > 0)
    sample_set = tk_iuset_from_ivec(NULL, sample_ids);

  uint64_t n_selected_features = (selected_features != NULL && selected_features->n > 0)
    ? selected_features->n : n_features;
  uint64_t in_bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t out_bytes_per_sample = TK_CVEC_BITS_BYTES(n_selected_features);
  uint8_t *src_data = (uint8_t *)src_bitmap->a;

  uint64_t n_output_samples = 0;
  if (sample_set != NULL) {
    for (uint64_t s = 0; s < n_samples; s ++) {
      if (tk_iuset_contains(sample_set, (int64_t) s))
        n_output_samples ++;
    }
  } else {
    n_output_samples = n_samples;
  }

  uint64_t row_stride_bits;
  bool use_packed;

  if (dest_stride > 0) {
    use_packed = true;
    row_stride_bits = TK_CVEC_BITS_BYTES(dest_stride) * CHAR_BIT;
  } else {
    use_packed = false;
    row_stride_bits = out_bytes_per_sample * CHAR_BIT;
  }

  uint64_t final_samples = dest_sample + n_output_samples;
  uint64_t total_bytes;
  if (use_packed) {
    uint64_t total_bits = dest_sample * row_stride_bits + n_output_samples * n_selected_features;
    if (dest_stride > 0)
      total_bits = (dest_sample + 1) * row_stride_bits;
    total_bytes = TK_CVEC_BITS_BYTES(total_bits);
  } else {
    total_bytes = final_samples * out_bytes_per_sample;
  }

  tk_cvec_ensure(dest, total_bytes);
  uint8_t *dest_data = (uint8_t *)dest->a;
  if (dest_sample == 0) {
    dest->n = total_bytes;
    memset(dest_data, 0, total_bytes);
  } else {
    uint64_t old_size = dest->n;
    dest->n = total_bytes;
    if (total_bytes > old_size)
      memset(dest_data + old_size, 0, total_bytes - old_size);
  }

  uint64_t write_sample = dest_sample;
  uint64_t write_bit_offset = 0;
  if (use_packed)
    write_bit_offset = dest_sample * row_stride_bits;

  for (uint64_t s = 0; s < n_samples; s ++) {
    if (sample_set != NULL && !tk_iuset_contains(sample_set, (int64_t) s))
      continue;
    uint8_t *temp = malloc(out_bytes_per_sample);
    if (!temp) {
      if (sample_set != NULL)
        tk_iuset_destroy(sample_set);
      return;
    }
    memset(temp, 0, out_bytes_per_sample);
    uint64_t in_offset = s * in_bytes_per_sample;
    if (selected_features != NULL && selected_features->n > 0) {
      for (uint64_t i = 0; i < selected_features->n; i ++) {
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
      memcpy(temp, src_data + in_offset, in_bytes_per_sample);
    }
    if (use_packed) {
      for (uint64_t i = 0; i < n_selected_features; i++) {
        uint64_t src_byte = i / CHAR_BIT;
        uint8_t src_bit_pos = i % CHAR_BIT;
        if (temp[src_byte] & (1u << src_bit_pos)) {
          uint64_t dst_byte = write_bit_offset / CHAR_BIT;
          uint8_t dst_bit_pos = write_bit_offset % CHAR_BIT;
          dest_data[dst_byte] |= (1u << dst_bit_pos);
        }
        write_bit_offset++;
      }
    } else {
      uint64_t out_offset = write_sample * out_bytes_per_sample;
      memcpy(dest_data + out_offset, temp, out_bytes_per_sample);
      write_sample++;
    }
    free(temp);
  }

  if (sample_set != NULL)
    tk_iuset_destroy(sample_set);
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

static inline tk_ivec_t *tk_ivec_bits_top_mi (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);

  if (codes) {

    tk_ivec_t *counts = tk_ivec_create(0, n_visible * n_hidden * 4, 0, 0);
    tk_ivec_zero(counts);
    uint64_t prev_sample = UINT64_MAX;
    uint8_t *sample_codes = NULL;
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit_idx = set_bits->a[i];
      if (bit_idx < 0)
        continue;
      uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
      uint64_t feature_idx = (uint64_t)bit_idx % n_visible;
      if (sample_idx >= n_samples || feature_idx >= n_visible)
        continue;
      if (sample_idx != prev_sample) {
        prev_sample = sample_idx;
        sample_codes = (uint8_t *)(codes + sample_idx * TK_CVEC_BITS_BYTES(n_hidden));
      }
      for (uint64_t j = 0; j < n_hidden; j++) {
        uint64_t byte_idx = j / CHAR_BIT;
        uint64_t bit_pos = j % CHAR_BIT;
        bool hidden = (sample_codes[byte_idx] & (1u << bit_pos)) != 0;
        counts->a[feature_idx * n_hidden * 4 + j * 4 + 2 + (hidden ? 1 : 0)]++;
      }
    }

    tk_iuset_t *sample_features = tk_iuset_create(0, 0);
    for (uint64_t s = 0; s < n_samples; s++) {
      tk_iuset_clear(sample_features);
      for (uint64_t i = 0; i < set_bits->n; i++) {
        int64_t bit_idx = set_bits->a[i];
        if (bit_idx < 0)
          continue;
        uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
        uint64_t feature_idx = (uint64_t)bit_idx % n_visible;
        if (sample_idx == s) {
          int absent;
          tk_iuset_put(sample_features, (int64_t)feature_idx, &absent);
        }
      }
      uint8_t *s_codes = (uint8_t *)(codes + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t f = 0; f < n_visible; f++) {
        if (!tk_iuset_contains(sample_features, (int64_t)f)) {
          for (uint64_t j = 0; j < n_hidden; j++) {
            uint64_t byte_idx = j / CHAR_BIT;
            uint64_t bit_pos = j % CHAR_BIT;
            bool hidden = (s_codes[byte_idx] & (1u << bit_pos)) != 0;
            counts->a[f * n_hidden * 4 + j * 4 + (hidden ? 1 : 0)]++;
          }
        }
      }
    }
    tk_iuset_destroy(sample_features);

    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    for (uint64_t f = 0; f < n_visible; f++) {
      double max_mi = 0.0;
      for (uint64_t j = 0; j < n_hidden; j++) {
        int64_t c[4];
        int64_t *counts_ptr = counts->a + f * n_hidden * 4 + j * 4;
        for (int k = 0; k < 4; k++)
          c[k] = counts_ptr[k] + 1;
        double total = c[0] + c[1] + c[2] + c[3];
        double mi = 0.0;
        if (total > 0.0) {
          for (unsigned int o = 0; o < 4; o++) {
            if (c[o] == 0)
              continue;
            double p_fb = c[o] / total;
            unsigned int feat = o >> 1;
            unsigned int hid = o & 1;
            double pf = (c[2] + c[3]) / total;
            if (feat == 0) pf = 1.0 - pf;
            double ph = (c[1] + c[3]) / total;
            if (hid == 0) ph = 1.0 - ph;
            double d = pf * ph;
            if (d > 0) mi += p_fb * log2(p_fb / d);
          }
        }
        if (mi > max_mi) max_mi = mi;
      }
      tk_rank_t r = { (int64_t)f, max_mi };
      tk_rvec_hmin(top_heap, top_k, r);
    }

    tk_ivec_destroy(counts);

    tk_rvec_desc(top_heap, 0, top_heap->n);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
    tk_rvec_keys(L, top_heap, out);
    tk_rvec_values(L, top_heap, weights);
    tk_rvec_destroy(top_heap);
    return out;

  } else if (labels) {

    tk_ivec_asc(labels, 0, labels->n);
    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    tk_ivec_t *feat_counts = tk_ivec_create(0, n_visible, 0, 0);
    tk_ivec_t *label_counts = tk_ivec_create(0, n_hidden, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_zero(label_counts);
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit = set_bits->a[i];
      if (bit >= 0) {
        uint64_t f = (uint64_t)bit % n_visible;
        feat_counts->a[f]++;
      }
    }

    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0) {
        uint64_t h = (uint64_t)bit % n_hidden;
        label_counts->a[h]++;
      }
    }

    size_t si = 0, li = 0;
    while (si < set_bits->n) {
      if (set_bits->a[si] < 0) {
        si++;
        continue;
      }
      uint64_t s_sample = (uint64_t)set_bits->a[si] / n_visible;
      uint64_t f = (uint64_t)set_bits->a[si] % n_visible;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden < s_sample) {
        li++;
      }
      if (li >= labels->n || labels->a[li] < 0 ||
        (uint64_t)labels->a[li] / n_hidden > s_sample) {
        si++;
        continue;
      }
      size_t li_start = li;
      while (li < labels->n && labels->a[li] >= 0 &&
        (uint64_t)labels->a[li] / n_hidden == s_sample) {
        uint64_t h = (uint64_t)labels->a[li] % n_hidden;
        int64_t key = (int64_t)(f * n_hidden + h);
        tk_iumap_inc(active_counts, key);
        li++;
      }
      si++;
      if (si < set_bits->n && set_bits->a[si] >= 0 && (uint64_t)set_bits->a[si] / n_visible == s_sample) {
        li = li_start;
      }
    }

    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    tk_dvec_t *feat_max_mi = tk_dvec_create(0, n_visible, 0, 0);
    tk_dvec_zero(feat_max_mi);

    int64_t k, v;
    tk_umap_foreach(active_counts, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t h = (uint64_t)k % n_hidden;
      if (f >= n_visible || h >= n_hidden)
        continue;
      int64_t count_11 = v;
      int64_t feat_total = feat_counts->a[f];
      int64_t label_total = label_counts->a[h];
      int64_t count_10 = feat_total - count_11;
      int64_t count_01 = label_total - count_11;
      int64_t count_00 = (int64_t)n_samples - count_11 - count_10 - count_01;
      double observed_density = (double)(set_bits->n + labels->n) / (double)(n_samples * (n_visible + n_hidden));
      double smooth = observed_density * n_samples;
      if (smooth < 0.001)
        smooth = 0.001;
      double expected = (double)(feat_total * label_total) / n_samples;
      if (expected > 0 && count_11 < 0.5 * sqrt(expected))
        continue;
      double c[4];
      c[0] = count_00 + smooth;
      c[1] = count_01 + smooth;
      c[2] = count_10 + smooth;
      c[3] = count_11 + smooth;
      double total = c[0] + c[1] + c[2] + c[3];
      double mi = 0.0;
      if (total > 0.0) {
        for (unsigned int o = 0; o < 4; o++) {
          if (c[o] == 0)
            continue;
          double p_fb = c[o] / total;
          unsigned int feat = o >> 1;
          unsigned int hid = o & 1;
          double pf = (c[2] + c[3]) / total;
          if (feat == 0) pf = 1.0 - pf;
          double ph = (c[1] + c[3]) / total;
          if (hid == 0) ph = 1.0 - ph;
          double d = pf * ph;
          if (d > 0) mi += p_fb * log2(p_fb / d);
        }
      }
      if (mi > feat_max_mi->a[f]) {
        feat_max_mi->a[f] = mi;
      }
    }));

    for (uint64_t f = 0; f < n_visible; f++) {
      if (feat_max_mi->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_max_mi->a[f] };
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_max_mi);

    tk_rvec_desc(top_heap, 0, top_heap->n);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
    tk_rvec_keys(L, top_heap, out);
    tk_rvec_values(L, top_heap, weights);
    tk_rvec_destroy(top_heap);
    return out;

  } else {

    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;

  }
}

static inline tk_ivec_t *tk_ivec_bits_top_chi2 (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);

  if (codes) {

    tk_ivec_t *active_counts = tk_ivec_create(0, n_visible * n_hidden, 0, 0);
    tk_ivec_t *label_counts = tk_ivec_create(0, n_hidden, 0, 0);
    tk_ivec_t *feat_counts = tk_ivec_create(0, n_visible, 0, 0);
    tk_ivec_zero(active_counts);
    tk_ivec_zero(label_counts);
    tk_ivec_zero(feat_counts);

    uint64_t prev_sample = UINT64_MAX;
    uint8_t *sample_codes = NULL;
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit_idx = set_bits->a[i];
      if (bit_idx < 0)
        continue;
      uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
      uint64_t feature_idx = (uint64_t)bit_idx % n_visible;
      if (sample_idx >= n_samples || feature_idx >= n_visible)
        continue;
      feat_counts->a[feature_idx]++;
      if (sample_idx != prev_sample) {
        prev_sample = sample_idx;
        sample_codes = (uint8_t *)(codes + sample_idx * TK_CVEC_BITS_BYTES(n_hidden));
      }
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint64_t bit_pos = b % CHAR_BIT;
        if (sample_codes[byte_idx] & (1u << bit_pos)) {
          active_counts->a[feature_idx * n_hidden + b]++;
        }
      }
    }

    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint64_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos)) {
          label_counts->a[b]++;
        }
      }
    }

    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    for (uint64_t f = 0; f < n_visible; f++) {
      double max_chi2 = 0.0;
      for (uint64_t b = 0; b < n_hidden; b++) {
        int64_t A = active_counts->a[f * n_hidden + b];
        int64_t G = label_counts->a[b];
        int64_t C = feat_counts->a[f];
        if (C == 0 || G == 0 || C == (int64_t)n_samples || G == (int64_t)n_samples)
          continue;
        int64_t B = G - A;
        int64_t C_ = C - A;
        int64_t D = (int64_t)n_samples - C - B;
        double n = (double)n_samples;
        double E_A = ((double)C * (double)G) / n;
        double E_B = ((double)(n - C) * (double)G) / n;
        double E_C = ((double)C * (double)(n - G)) / n;
        double E_D = ((double)(n - C) * (double)(n - G)) / n;
        double chi2 = 0.0;
        if (E_A > 0) chi2 += ((A - E_A) * (A - E_A)) / E_A;
        if (E_B > 0) chi2 += ((B - E_B) * (B - E_B)) / E_B;
        if (E_C > 0) chi2 += ((C_ - E_C) * (C_ - E_C)) / E_C;
        if (E_D > 0) chi2 += ((D - E_D) * (D - E_D)) / E_D;
        if (chi2 > max_chi2) max_chi2 = chi2;
      }
      tk_rank_t r = { (int64_t)f, max_chi2 };
      tk_rvec_hmin(top_heap, top_k, r);
    }

    tk_ivec_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    tk_rvec_desc(top_heap, 0, top_heap->n);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
    tk_rvec_keys(L, top_heap, out);
    tk_rvec_values(L, top_heap, weights);
    tk_rvec_destroy(top_heap);
    return out;

  } else if (labels) {

    tk_ivec_asc(labels, 0, labels->n);
    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    tk_ivec_t *feat_counts = tk_ivec_create(0, n_visible, 0, 0);
    tk_ivec_t *label_counts = tk_ivec_create(0, n_hidden, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_zero(label_counts);
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit = set_bits->a[i];
      if (bit >= 0) {
        uint64_t f = (uint64_t)bit % n_visible;
        feat_counts->a[f]++;
      }
    }

    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0) {
        uint64_t h = (uint64_t)bit % n_hidden;
        label_counts->a[h]++;
      }
    }

    size_t si = 0, li = 0;
    while (si < set_bits->n) {
      if (set_bits->a[si] < 0) {
        si++;
        continue;
      }
      uint64_t s_sample = (uint64_t)set_bits->a[si] / n_visible;
      uint64_t f = (uint64_t)set_bits->a[si] % n_visible;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden < s_sample) {
        li++;
      }
      if (li >= labels->n || labels->a[li] < 0 || (uint64_t)labels->a[li] / n_hidden > s_sample) {
        si++;
        continue;
      }
      size_t li_start = li;
      while (li < labels->n && labels->a[li] >= 0 &&
        (uint64_t)labels->a[li] / n_hidden == s_sample) {
        uint64_t h = (uint64_t)labels->a[li] % n_hidden;
        int64_t key = (int64_t)(f * n_hidden + h);
        tk_iumap_inc(active_counts, key);
        li++;
      }
      si++;
      if (si < set_bits->n && set_bits->a[si] >= 0 &&
        (uint64_t)set_bits->a[si] / n_visible == s_sample) {
        li = li_start;
      }
    }

    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    tk_dvec_t *feat_max_chi2 = tk_dvec_create(0, n_visible, 0, 0);
    tk_dvec_zero(feat_max_chi2);

    int64_t k, v;
    tk_umap_foreach(active_counts, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t b = (uint64_t)k % n_hidden;
      if (f >= n_visible || b >= n_hidden)
        continue;
      int64_t A = v;
      int64_t C = feat_counts->a[f];
      int64_t G = label_counts->a[b];
      if (C == 0 || G == 0 || C == (int64_t)n_samples || G == (int64_t)n_samples)
        continue;
      double expected = (double)(C * G) / n_samples;
      if (expected > 0 && A < 0.5 * sqrt(expected))
        continue;
      double observed_density = (double)(set_bits->n + labels->n) / (double)(n_samples * (n_visible + n_hidden));
      double smooth = observed_density * n_samples;
      if (smooth < 0.001)
        smooth = 0.001;
      int64_t B = G - A;
      int64_t C_ = C - A;
      int64_t D = (int64_t)n_samples - C - B;
      double n = (double)(n_samples) + 4 * smooth;
      double C_smooth = C + 2 * smooth;
      double G_smooth = G + 2 * smooth;
      double E_A = (C_smooth * G_smooth) / n;
      double E_B = ((n - C_smooth) * G_smooth) / n;
      double E_C = (C_smooth * (n - G_smooth)) / n;
      double E_D = ((n - C_smooth) * (n - G_smooth)) / n;
      double chi2 = 0.0;
      if (E_A > 0) chi2 += ((A - E_A) * (A - E_A)) / E_A;
      if (E_B > 0) chi2 += ((B - E_B) * (B - E_B)) / E_B;
      if (E_C > 0) chi2 += ((C_ - E_C) * (C_ - E_C)) / E_C;
      if (E_D > 0) chi2 += ((D - E_D) * (D - E_D)) / E_D;
      if (chi2 > feat_max_chi2->a[f])
        feat_max_chi2->a[f] = chi2;
    }));

    for (uint64_t f = 0; f < n_visible; f++) {
      if (feat_max_chi2->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_max_chi2->a[f] };
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_max_chi2);

    tk_rvec_desc(top_heap, 0, top_heap->n);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
    tk_rvec_keys(L, top_heap, out);
    tk_rvec_values(L, top_heap, weights);
    tk_rvec_destroy(top_heap);
    return out;

  } else {

    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;

  }
}

static inline tk_ivec_t *tk_ivec_bits_top_entropy (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_hidden,
  uint64_t top_k
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);
  tk_ivec_t *bit_counts = tk_ivec_create(0, n_hidden, 0, 0);
  tk_ivec_zero(bit_counts);
  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit_idx = set_bits->a[i];
    if (bit_idx < 0)
      continue;
    uint64_t sample_idx = (uint64_t)bit_idx / n_hidden;
    uint64_t hidden_idx = (uint64_t)bit_idx % n_hidden;
    if (sample_idx >= n_samples || hidden_idx >= n_hidden)
      continue;
    bit_counts->a[hidden_idx]++;
  }
  tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
  for (uint64_t h = 0; h < n_hidden; h++) {
    double p = (double)bit_counts->a[h] / (double)n_samples;
    double entropy = 0.0;
    if (p > 0.0 && p < 1.0)
      entropy = -(p * log2(p) + (1.0 - p) * log2(1.0 - p));
    tk_rank_t r = { (int64_t)h, entropy };
    tk_rvec_hmin(top_heap, top_k, r);
  }
  tk_ivec_destroy(bit_counts);
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}

static inline tk_ivec_t *tk_cvec_bits_top_mi (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  uint64_t top_k
) {

  if (codes) {

    tk_ivec_t *counts = tk_ivec_create(0, n_features * n_hidden * 4, 0, 0);
    tk_ivec_zero(counts);
    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      uint8_t *sample_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        bool visible = (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) != 0;
        for (uint64_t j = 0; j < n_hidden; j++) {
          uint64_t j_byte = j / CHAR_BIT;
          uint8_t j_bit = j % CHAR_BIT;
          bool hidden = (sample_codes[j_byte] & (1u << j_bit)) != 0;
          counts->a[f * n_hidden * 4 + j * 4 + (visible ? 2 : 0) + (hidden ? 1 : 0)]++;
        }
      }
    }
    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    for (uint64_t f = 0; f < n_features; f++) {
      double max_mi = 0.0;
      for (uint64_t j = 0; j < n_hidden; j++) {
        int64_t c[4];
        int64_t *counts_ptr = counts->a + f * n_hidden * 4 + j * 4;
        for (int k = 0; k < 4; k++) {
          c[k] = counts_ptr[k] + 1;
        }
        double total = c[0] + c[1] + c[2] + c[3];
        double mi = 0.0;
        if (total > 0.0) {
          for (unsigned int o = 0; o < 4; o++) {
            if (c[o] == 0)
              continue;
            double p_fb = c[o] / total;
            unsigned int feat = o >> 1;
            unsigned int hid = o & 1;
            double pf = (c[2] + c[3]) / total;
            if (feat == 0) pf = 1.0 - pf;
            double ph = (c[1] + c[3]) / total;
            if (hid == 0) ph = 1.0 - ph;
            double d = pf * ph;
            if (d > 0) mi += p_fb * log2(p_fb / d);
          }
        }
        if (mi > max_mi) max_mi = mi;
      }
      tk_rank_t r = { (int64_t)f, max_mi };
      tk_rvec_hmin(top_heap, top_k, r);
    }

    tk_ivec_destroy(counts);
    tk_rvec_desc(top_heap, 0, top_heap->n);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
    tk_rvec_keys(L, top_heap, out);
    tk_rvec_values(L, top_heap, weights);
    tk_rvec_destroy(top_heap);
    return out;

  } else if (labels) {

    tk_ivec_asc(labels, 0, labels->n);
    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
    tk_ivec_t *label_counts = tk_ivec_create(0, n_hidden, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_zero(label_counts);
    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    uint64_t total_feature_bits = 0;
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
          feat_counts->a[f]++;
          total_feature_bits++;
        }
      }
    }

    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0) {
        uint64_t h = (uint64_t)bit % n_hidden;
        label_counts->a[h]++;
      }
    }

    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      int64_t label_start = tk_ivec_set_find(labels->a, 0, (int64_t) labels->n, (int64_t)(s * n_hidden));
      if (label_start < 0)
        label_start = -(label_start + 1);
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (!(bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)))
          continue;
        for (int64_t li = label_start;
        li < (int64_t) labels->n && labels->a[li] >= 0 &&
        (uint64_t)labels->a[li] / n_hidden == s;
        li++) {
          uint64_t h = (uint64_t)labels->a[li] % n_hidden;
          int64_t key = (int64_t)(f * n_hidden + h);
          tk_iumap_inc(active_counts, key);
        }
      }
    }

    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    tk_dvec_t *feat_max_mi = tk_dvec_create(0, n_features, 0, 0);
    tk_dvec_zero(feat_max_mi);

    double observed_density = (double)(total_feature_bits + labels->n) / (double)(n_samples * (n_features + n_hidden));
    double smooth = observed_density * n_samples;
    if (smooth < 0.001)
      smooth = 0.001;
    int64_t k, v;
    tk_umap_foreach(active_counts, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t h = (uint64_t)k % n_hidden;
      if (f >= n_features || h >= n_hidden)
        continue;
      int64_t count_11 = v;
      int64_t feat_total = feat_counts->a[f];
      int64_t label_total = label_counts->a[h];
      int64_t count_10 = feat_total - count_11;
      int64_t count_01 = label_total - count_11;
      int64_t count_00 = (int64_t)n_samples - count_11 - count_10 - count_01;
      double expected = (double)(feat_total * label_total) / n_samples;
      if (expected > 0 && count_11 < 0.5 * sqrt(expected))
        continue;
      double c[4];
      c[0] = count_00 + smooth;
      c[1] = count_01 + smooth;
      c[2] = count_10 + smooth;
      c[3] = count_11 + smooth;
      double total = c[0] + c[1] + c[2] + c[3];
      double mi = 0.0;
      if (total > 0.0) {
        for (unsigned int o = 0; o < 4; o++) {
          if (c[o] == 0)
            continue;
          double p_fb = c[o] / total;
          unsigned int feat = o >> 1;
          unsigned int hid = o & 1;
          double pf = (c[2] + c[3]) / total;
          if (feat == 0) pf = 1.0 - pf;
          double ph = (c[1] + c[3]) / total;
          if (hid == 0) ph = 1.0 - ph;
          double d = pf * ph;
          if (d > 0) mi += p_fb * log2(p_fb / d);
        }
      }
      if (mi > feat_max_mi->a[f]) {
        feat_max_mi->a[f] = mi;
      }
    }));

    for (uint64_t f = 0; f < n_features; f++) {
      if (feat_max_mi->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_max_mi->a[f] };
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_max_mi);
    tk_iumap_destroy(active_counts);
    tk_rvec_desc(top_heap, 0, top_heap->n);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
    tk_rvec_keys(L, top_heap, out);
    tk_rvec_values(L, top_heap, weights);
    tk_rvec_destroy(top_heap);
    return out;


  } else {

    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;

  }
}

static inline tk_ivec_t *tk_cvec_bits_top_chi2 (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  uint64_t top_k
) {

  if (codes) {

    tk_ivec_t *active_counts = tk_ivec_create(0, n_features * n_hidden, 0, 0);
    tk_ivec_t *label_counts = tk_ivec_create(0, n_hidden, 0, 0);
    tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
    tk_ivec_zero(active_counts);
    tk_ivec_zero(label_counts);
    tk_ivec_zero(feat_counts);

    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      uint8_t *sample_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
          feat_counts->a[f]++;
          for (uint64_t b = 0; b < n_hidden; b++) {
            uint64_t b_byte = b / CHAR_BIT;
            uint8_t b_bit = b % CHAR_BIT;
            if (sample_codes[b_byte] & (1u << b_bit)) {
              active_counts->a[f * n_hidden + b]++;
            }
          }
        }
      }
    }

    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint8_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos)) {
          label_counts->a[b]++;
        }
      }
    }

    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    for (uint64_t f = 0; f < n_features; f++) {
      double max_chi2 = 0.0;
      for (uint64_t b = 0; b < n_hidden; b++) {
        int64_t A = active_counts->a[f * n_hidden + b];
        int64_t G = label_counts->a[b];
        int64_t C = feat_counts->a[f];
        if (C == 0 || G == 0 || C == (int64_t)n_samples || G == (int64_t)n_samples)
          continue;
        int64_t B = G - A;
        int64_t C_ = C - A;
        int64_t D = (int64_t)n_samples - C - B;
        double n = (double)n_samples;
        double E_A = ((double)C * (double)G) / n;
        double E_B = ((double)(n - C) * (double)G) / n;
        double E_C = ((double)C * (double)(n - G)) / n;
        double E_D = ((double)(n - C) * (double)(n - G)) / n;
        double chi2 = 0.0;
        if (E_A > 0) chi2 += ((A - E_A) * (A - E_A)) / E_A;
        if (E_B > 0) chi2 += ((B - E_B) * (B - E_B)) / E_B;
        if (E_C > 0) chi2 += ((C_ - E_C) * (C_ - E_C)) / E_C;
        if (E_D > 0) chi2 += ((D - E_D) * (D - E_D)) / E_D;
        if (chi2 > max_chi2) max_chi2 = chi2;
      }
      tk_rank_t r = { (int64_t)f, max_chi2 };
      tk_rvec_hmin(top_heap, top_k, r);
    }

    tk_ivec_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    tk_rvec_desc(top_heap, 0, top_heap->n);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
    tk_rvec_keys(L, top_heap, out);
    tk_rvec_values(L, top_heap, weights);
    tk_rvec_destroy(top_heap);
    return out;

  } else if (labels) {

    tk_ivec_asc(labels, 0, labels->n);
    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
    tk_ivec_t *label_counts = tk_ivec_create(0, n_hidden, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_zero(label_counts);

    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    uint64_t total_feature_bits = 0;
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
          feat_counts->a[f]++;
          total_feature_bits++;
        }
      }
    }

    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0) {
        uint64_t h = (uint64_t)bit % n_hidden;
        label_counts->a[h]++;
      }
    }

    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      int64_t label_start = tk_ivec_set_find(
        labels->a, 0, (int64_t) labels->n, (int64_t)(s * n_hidden)
      );
      if (label_start < 0) {
        label_start = -(label_start + 1);
      }
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (!(bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)))
          continue;
        for (int64_t li = label_start;
        li < (int64_t) labels->n && labels->a[li] >= 0 &&
        (uint64_t)labels->a[li] / n_hidden == s;
        li++) {
          uint64_t h = (uint64_t)labels->a[li] % n_hidden;
          int64_t key = (int64_t)(f * n_hidden + h);
          tk_iumap_inc(active_counts, key);
        }
      }
    }

    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    tk_dvec_t *feat_max_chi2 = tk_dvec_create(0, n_features, 0, 0);
    tk_dvec_zero(feat_max_chi2);

    double observed_density = (double)(total_feature_bits + labels->n) / (double)(n_samples * (n_features + n_hidden));
    double smooth = observed_density * n_samples;
    if (smooth < 0.001)
      smooth = 0.001;

    int64_t k, v;
    tk_umap_foreach(active_counts, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t b = (uint64_t)k % n_hidden;
      if (f >= n_features || b >= n_hidden)
        continue;
      int64_t A = v;
      int64_t C = feat_counts->a[f];
      int64_t G = label_counts->a[b];
      if (C == 0 || G == 0 || C == (int64_t)n_samples || G == (int64_t)n_samples)
        continue;
      double expected = (double)(C * G) / n_samples;
      if (expected > 0 && A < 0.5 * sqrt(expected))
        continue;
      int64_t B = G - A;
      int64_t C_ = C - A;
      int64_t D = (int64_t)n_samples - C - B;
      double n = (double)(n_samples) + 4 * smooth;
      double C_smooth = C + 2 * smooth;
      double G_smooth = G + 2 * smooth;
      double E_A = (C_smooth * G_smooth) / n;
      double E_B = ((n - C_smooth) * G_smooth) / n;
      double E_C = (C_smooth * (n - G_smooth)) / n;
      double E_D = ((n - C_smooth) * (n - G_smooth)) / n;
      double chi2 = 0.0;
      if (E_A > 0) chi2 += ((A - E_A) * (A - E_A)) / E_A;
      if (E_B > 0) chi2 += ((B - E_B) * (B - E_B)) / E_B;
      if (E_C > 0) chi2 += ((C_ - E_C) * (C_ - E_C)) / E_C;
      if (E_D > 0) chi2 += ((D - E_D) * (D - E_D)) / E_D;
      if (chi2 > feat_max_chi2->a[f])
        feat_max_chi2->a[f] = chi2;
    }));

    for (uint64_t f = 0; f < n_features; f++) {
      if (feat_max_chi2->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_max_chi2->a[f] };
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_max_chi2);

    tk_rvec_desc(top_heap, 0, top_heap->n);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
    tk_rvec_keys(L, top_heap, out);
    tk_rvec_values(L, top_heap, weights);
    tk_rvec_destroy(top_heap);
    return out;

  } else {

    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;

  }
}

static inline tk_ivec_t *tk_ivec_bits_top_lift (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);

  if (codes) {

    tk_ivec_t *feat_counts = tk_ivec_create(0, n_visible, 0, 0);
    tk_ivec_t *label_counts = tk_ivec_create(0, n_hidden, 0, 0);
    tk_ivec_t *cooccur_counts = tk_ivec_create(0, n_visible * n_hidden, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_zero(label_counts);
    tk_ivec_zero(cooccur_counts);

    uint64_t prev_sample = UINT64_MAX;
    uint8_t *sample_codes = NULL;
    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit_idx = set_bits->a[i];
      if (bit_idx < 0)
        continue;
      uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
      uint64_t feature_idx = (uint64_t)bit_idx % n_visible;
      if (sample_idx >= n_samples || feature_idx >= n_visible)
        continue;
      feat_counts->a[feature_idx]++;
      if (sample_idx != prev_sample) {
        prev_sample = sample_idx;
        sample_codes = (uint8_t *)(codes + sample_idx * TK_CVEC_BITS_BYTES(n_hidden));
      }
      for (uint64_t h = 0; h < n_hidden; h++) {
        uint64_t byte_idx = h / CHAR_BIT;
        uint64_t bit_pos = h % CHAR_BIT;
        if (sample_codes[byte_idx] & (1u << bit_pos))
          cooccur_counts->a[feature_idx * n_hidden + h]++;
      }
    }

    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t h = 0; h < n_hidden; h++) {
        uint64_t byte_idx = h / CHAR_BIT;
        uint64_t bit_pos = h % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos))
          label_counts->a[h]++;
      }
    }

    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    for (uint64_t f = 0; f < n_visible; f++) {
      double max_score = 0.0;
      int64_t feat_total = feat_counts->a[f];
      if (feat_total == 0)
        continue;
      for (uint64_t h = 0; h < n_hidden; h++) {
        int64_t label_total = label_counts->a[h];
        if (label_total == 0)
          continue;
        int64_t cooccur = cooccur_counts->a[f * n_hidden + h];
        if (cooccur == 0)
          continue;
        double p_label_given_feature = (double)cooccur / feat_total;
        double p_label = (double)label_total / n_samples;
        double lift = p_label_given_feature / (p_label + 1e-10);
        double expected = (double)(feat_total * label_total) / n_samples;
        if (expected > 1.0 && cooccur < sqrt(expected))
          continue;
        double weighted_lift = lift * log2(1 + cooccur);
        if (weighted_lift > max_score)
          max_score = weighted_lift;
      }
      if (max_score > 0) {
        tk_rank_t r = { (int64_t)f, max_score };
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_ivec_destroy(cooccur_counts);

    tk_rvec_desc(top_heap, 0, top_heap->n);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
    tk_rvec_keys(L, top_heap, out);
    tk_rvec_values(L, top_heap, weights);
    tk_rvec_destroy(top_heap);
    return out;

  } else if (labels) {

    tk_ivec_asc(labels, 0, labels->n);
    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    tk_ivec_t *feat_counts = tk_ivec_create(0, n_visible, 0, 0);
    tk_ivec_t *label_counts = tk_ivec_create(0, n_hidden, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_zero(label_counts);

    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit = set_bits->a[i];
      if (bit >= 0) {
        uint64_t f = (uint64_t)bit % n_visible;
        feat_counts->a[f]++;
      }
    }

    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0) {
        uint64_t h = (uint64_t)bit % n_hidden;
        label_counts->a[h]++;
      }
    }

    size_t si = 0, li = 0;
    while (si < set_bits->n) {
      if (set_bits->a[si] < 0) {
        si++;
        continue;
      }
      uint64_t s_sample = (uint64_t)set_bits->a[si] / n_visible;
      uint64_t f = (uint64_t)set_bits->a[si] % n_visible;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden < s_sample) {
        li++;
      }
      if (li >= labels->n || labels->a[li] < 0 ||
        (uint64_t)labels->a[li] / n_hidden > s_sample) {
        si++;
        continue;
      }
      size_t li_start = li;
      while (li < labels->n && labels->a[li] >= 0 &&
        (uint64_t)labels->a[li] / n_hidden == s_sample) {
        uint64_t h = (uint64_t)labels->a[li] % n_hidden;
        int64_t key = (int64_t)(f * n_hidden + h);
        tk_iumap_inc(active_counts, key);
        li++;
      }
      si++;
      if (si < set_bits->n && set_bits->a[si] >= 0 && (uint64_t)set_bits->a[si] / n_visible == s_sample)
        li = li_start;
    }

    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    tk_dvec_t *feat_max_score = tk_dvec_create(0, n_visible, 0, 0);
    tk_dvec_zero(feat_max_score);

    int64_t k, v;
    tk_umap_foreach(active_counts, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t h = (uint64_t)k % n_hidden;
      if (f >= n_visible || h >= n_hidden)
        continue;
      int64_t cooccur = v;
      int64_t feat_total = feat_counts->a[f];
      int64_t label_total = label_counts->a[h];
      double p_label_given_feature = (double)cooccur / feat_total;
      double p_label = (double)label_total / n_samples;
      double lift = p_label_given_feature / (p_label + 1e-10);
      double expected = (double)(feat_total * label_total) / n_samples;
      if (expected > 1.0 && cooccur < sqrt(expected))
        continue;
      double weighted_lift = lift * log2(1 + cooccur);
      if (weighted_lift > feat_max_score->a[f])
        feat_max_score->a[f] = weighted_lift;
    }));
    for (uint64_t f = 0; f < n_visible; f++) {
      if (feat_max_score->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_max_score->a[f] };
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    tk_rvec_desc(top_heap, 0, top_heap->n);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
    tk_rvec_keys(L, top_heap, out);
    tk_rvec_values(L, top_heap, weights);
    tk_rvec_destroy(top_heap);
    return out;

  } else {

    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;

  }
}

static inline tk_ivec_t *tk_cvec_bits_top_entropy (
  lua_State *L,
  tk_cvec_t *codes,
  uint64_t n_samples,
  uint64_t n_hidden,
  uint64_t top_k
) {
  tk_ivec_t *bit_counts = tk_ivec_create(0, n_hidden, 0, 0);
  tk_ivec_zero(bit_counts);
  for (uint64_t s = 0; s < n_samples; s++) {
    uint8_t *sample_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
    for (uint64_t h = 0; h < n_hidden; h++) {
      uint64_t byte_idx = h / CHAR_BIT;
      uint8_t bit_idx = h % CHAR_BIT;
      if (sample_codes[byte_idx] & (1u << bit_idx))
        bit_counts->a[h]++;
    }
  }
  tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
  for (uint64_t h = 0; h < n_hidden; h++) {
    double p = (double)bit_counts->a[h] / (double)n_samples;
    double entropy = 0.0;
    if (p > 0.0 && p < 1.0)
      entropy = -(p * log2(p) + (1.0 - p) * log2(1.0 - p));
    tk_rank_t r = { (int64_t)h, entropy };
    tk_rvec_hmin(top_heap, top_k, r);
  }
  tk_ivec_destroy(bit_counts);
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
  uint64_t n_features,
  double min_df,
  double max_df,
  uint64_t top_k
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);
  tk_iuset_t **feature_docs = (tk_iuset_t **)calloc(n_features, sizeof(tk_iuset_t *));
  for (uint64_t i = 0; i < n_features; i++)
    feature_docs[i] = tk_iuset_create(0, 0);
  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit_idx = set_bits->a[i];
    if (bit_idx < 0)
      continue;
    uint64_t sample_idx = (uint64_t)bit_idx / n_features;
    uint64_t feature_idx = (uint64_t)bit_idx % n_features;
    if (sample_idx < n_samples && feature_idx < n_features) {
      int absent;
      tk_iuset_put(feature_docs[feature_idx], (int64_t)sample_idx, &absent);
    }
  }
  tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
  double min_df_abs = min_df < 0 ? -min_df : min_df * n_samples;
  double max_df_abs = max_df < 0 ? -max_df : max_df * n_samples;
  for (uint64_t i = 0; i < n_features; i++) {
    double df_count = (double)tk_iuset_size(feature_docs[i]);
    double idf = log((double)(n_samples + 1) / (df_count + 1));
    if (df_count >= min_df_abs && df_count <= max_df_abs) {
      tk_rank_t r = { (int64_t)i, idf };
      tk_rvec_hmin(top_heap, top_k, r);
    }
  }
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  for (uint64_t i = 0; i < n_features; i++)
    tk_iuset_destroy(feature_docs[i]);
  free(feature_docs);
  tk_rvec_destroy(top_heap);
  return out;
}

static inline tk_ivec_t *tk_cvec_bits_top_lift (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  uint64_t top_k
) {

  if (codes) {

    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
    tk_ivec_t *label_counts = tk_ivec_create(0, n_hidden, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_zero(label_counts);

    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint8_t *codes_data = (uint8_t *)codes->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    uint64_t codes_bytes_per_sample = TK_CVEC_BITS_BYTES(n_hidden);

    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx))
          feat_counts->a[f]++;
      }
    }

    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t codes_offset = s * codes_bytes_per_sample;
      for (uint64_t h = 0; h < n_hidden; h++) {
        uint64_t byte_idx = h / CHAR_BIT;
        uint8_t bit_idx = h % CHAR_BIT;
        if (codes_data[codes_offset + byte_idx] & (1u << bit_idx))
          label_counts->a[h]++;
      }
    }

    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      uint64_t codes_offset = s * codes_bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t f_byte_idx = f / CHAR_BIT;
        uint8_t f_bit_idx = f % CHAR_BIT;
        if (!(bitmap_data[sample_offset + f_byte_idx] & (1u << f_bit_idx)))
          continue;
        for (uint64_t h = 0; h < n_hidden; h++) {
          uint64_t h_byte_idx = h / CHAR_BIT;
          uint8_t h_bit_idx = h % CHAR_BIT;
          if (codes_data[codes_offset + h_byte_idx] & (1u << h_bit_idx)) {
            int64_t key = (int64_t)(f * n_hidden + h);
            tk_iumap_inc(active_counts, key);
          }
        }
      }
    }

    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    tk_dvec_t *feat_max_score = tk_dvec_create(0, n_features, 0, 0);
    tk_dvec_zero(feat_max_score);

    int64_t k, v;
    tk_umap_foreach(active_counts, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t h = (uint64_t)k % n_hidden;
      if (f >= n_features || h >= n_hidden)
        continue;
      int64_t cooccur = v;
      int64_t feat_total = feat_counts->a[f];
      int64_t label_total = label_counts->a[h];
      double p_label_given_feature = (double)cooccur / feat_total;
      double p_label = (double)label_total / n_samples;
      double lift = p_label_given_feature / (p_label + 1e-10);
      double expected = (double)(feat_total * label_total) / n_samples;
      if (expected > 1.0 && cooccur < sqrt(expected))
        continue;
      double weighted_lift = lift * log2(1 + cooccur);
      if (weighted_lift > feat_max_score->a[f])
        feat_max_score->a[f] = weighted_lift;
    }));

    for (uint64_t f = 0; f < n_features; f++) {
      if (feat_max_score->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_max_score->a[f] };
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_max_score);

    tk_rvec_desc(top_heap, 0, top_heap->n);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
    tk_rvec_keys(L, top_heap, out);
    tk_rvec_values(L, top_heap, weights);
    tk_rvec_destroy(top_heap);
    return out;

  } else if (labels) {

    tk_ivec_asc(labels, 0, labels->n);
    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
    tk_ivec_t *label_counts = tk_ivec_create(0, n_hidden, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_zero(label_counts);

    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx))
          feat_counts->a[f]++;
      }
    }

    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0) {
        uint64_t h = (uint64_t)bit % n_hidden;
        label_counts->a[h]++;
      }
    }

    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      int64_t label_start = tk_ivec_set_find(labels->a, 0, (int64_t) labels->n, (int64_t)(s * n_hidden));
      if (label_start < 0)
        label_start = -(label_start + 1);
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (!(bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)))
          continue;
        for (int64_t li = label_start;
        li < (int64_t) labels->n && labels->a[li] >= 0 &&
        (uint64_t)labels->a[li] / n_hidden == s;
        li++) {
          uint64_t h = (uint64_t)labels->a[li] % n_hidden;
          int64_t key = (int64_t)(f * n_hidden + h);
          tk_iumap_inc(active_counts, key);
        }
      }
    }

    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    tk_dvec_t *feat_max_score = tk_dvec_create(0, n_features, 0, 0);
    tk_dvec_zero(feat_max_score);

    int64_t k, v;
    tk_umap_foreach(active_counts, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t h = (uint64_t)k % n_hidden;
      if (f >= n_features || h >= n_hidden)
        continue;
      int64_t cooccur = v;
      int64_t feat_total = feat_counts->a[f];
      int64_t label_total = label_counts->a[h];
      double p_label_given_feature = (double)cooccur / feat_total;
      double p_label = (double)label_total / n_samples;
      double lift = p_label_given_feature / (p_label + 1e-10);
      double expected = (double)(feat_total * label_total) / n_samples;
      if (expected > 1.0 && cooccur < sqrt(expected))
        continue;
      double weighted_lift = lift * log2(1 + cooccur);
      if (weighted_lift > feat_max_score->a[f])
        feat_max_score->a[f] = weighted_lift;
    }));

    for (uint64_t f = 0; f < n_features; f++) {
      if (feat_max_score->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_max_score->a[f] };
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_max_score);

    tk_rvec_desc(top_heap, 0, top_heap->n);
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
    tk_rvec_keys(L, top_heap, out);
    tk_rvec_values(L, top_heap, weights);
    tk_rvec_destroy(top_heap);
    return out;

    return out;

  } else {

    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;

  }
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
  tk_iuset_t **feature_docs = (tk_iuset_t **)calloc(n_features, sizeof(tk_iuset_t *));
  for (uint64_t i = 0; i < n_features; i++)
    feature_docs[i] = tk_iuset_create(0, 0);
  uint8_t *data = (uint8_t *)bitmap->a;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t sample_offset = s * bytes_per_sample;
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t byte_idx = f / CHAR_BIT;
      uint8_t bit_idx = f % CHAR_BIT;
      if (data[sample_offset + byte_idx] & (1u << bit_idx)) {
        int absent;
        tk_iuset_put(feature_docs[f], (int64_t)s, &absent);
      }
    }
  }
  tk_rvec_t *top_heap = tk_rvec_create(L, 0, 0, 0);
  double min_df_abs = min_df < 0 ? -min_df : min_df * n_samples;
  double max_df_abs = max_df < 0 ? -max_df : max_df * n_samples;
  for (uint64_t i = 0; i < n_features; i++) {
    double df_count = (double)tk_iuset_size(feature_docs[i]);
    double idf = log((double)(n_samples + 1) / (df_count + 1));
    if (df_count >= min_df_abs && df_count <= max_df_abs) {
      tk_rank_t r = { (int64_t)i, idf };
      tk_rvec_hmin(top_heap, top_k, r);
    }
  }
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  for (uint64_t i = 0; i < n_features; i++)
    tk_iuset_destroy(feature_docs[i]);
  free(feature_docs);
  tk_rvec_destroy(top_heap);
  return out;
}

#endif
