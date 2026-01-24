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
    #pragma omp parallel for
    for (uint64_t f = 0; f < n_visible; f++) {
      double sum_mi = 0.0;
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
        sum_mi += mi;
      }
      tk_rank_t r = { (int64_t)f, sum_mi };
      #pragma omp critical
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
    tk_dvec_t *feat_sum_mi = tk_dvec_create(0, n_visible, 0, 0);
    tk_dvec_zero(feat_sum_mi);

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
      double c[4];
      c[0] = count_00 + 1;
      c[1] = count_01 + 1;
      c[2] = count_10 + 1;
      c[3] = count_11 + 1;
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
      feat_sum_mi->a[f] += mi;
    }));

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_visible; f++) {
      if (feat_sum_mi->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_sum_mi->a[f] };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_sum_mi);

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
    #pragma omp parallel for
    for (uint64_t f = 0; f < n_visible; f++) {
      double sum_chi2 = 0.0;
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
        sum_chi2 += chi2;
      }
      tk_rank_t r = { (int64_t)f, sum_chi2 };
      #pragma omp critical
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
    tk_dvec_t *feat_sum_chi2 = tk_dvec_create(0, n_visible, 0, 0);
    tk_dvec_zero(feat_sum_chi2);

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
      int64_t B = G - A;
      int64_t C_ = C - A;
      int64_t D = (int64_t)n_samples - C - B;
      double n = (double)(n_samples) + 4;
      double C_smooth = C + 2;
      double G_smooth = G + 2;
      double E_A = (C_smooth * G_smooth) / n;
      double E_B = ((n - C_smooth) * G_smooth) / n;
      double E_C = (C_smooth * (n - G_smooth)) / n;
      double E_D = ((n - C_smooth) * (n - G_smooth)) / n;
      double chi2 = 0.0;
      if (E_A > 0) chi2 += ((A - E_A) * (A - E_A)) / E_A;
      if (E_B > 0) chi2 += ((B - E_B) * (B - E_B)) / E_B;
      if (E_C > 0) chi2 += ((C_ - E_C) * (C_ - E_C)) / E_C;
      if (E_D > 0) chi2 += ((D - E_D) * (D - E_D)) / E_D;
      feat_sum_chi2->a[f] += chi2;
    }));

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_visible; f++) {
      if (feat_sum_chi2->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_sum_chi2->a[f] };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_sum_chi2);

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






static inline tk_ivec_t *tk_ivec_bits_top_coherence (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k,
  double lambda
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);

  if (!codes) {
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  tk_ivec_t *active_counts = tk_ivec_create(0, n_visible * n_hidden, 0, 0);
  tk_ivec_t *feat_counts = tk_ivec_create(0, n_visible, 0, 0);
  tk_ivec_zero(active_counts);
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

  tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
  #pragma omp parallel for
  for (uint64_t f = 0; f < n_visible; f++) {
    int64_t k_f = feat_counts->a[f];
    if (k_f < 2)
      continue;
    double sum_disagreements = 0.0;
    for (uint64_t b = 0; b < n_hidden; b++) {
      int64_t n1 = active_counts->a[f * n_hidden + b];
      int64_t n0 = k_f - n1;
      sum_disagreements += (double)n1 * (double)n0;
    }
    double num_pairs = (double)k_f * (double)(k_f - 1) / 2.0;
    double mean_hamming = sum_disagreements / num_pairs;
    double normalized = mean_hamming / (double)n_hidden;
    double coherence = 1.0 - normalized;
    double penalty = lambda / sqrt(num_pairs);
    double score = coherence - penalty;
    if (score <= 0.0)
      continue;
    tk_rank_t r = { (int64_t)f, score };
    #pragma omp critical
    tk_rvec_hmin(top_heap, top_k, r);
  }

  tk_ivec_destroy(active_counts);
  tk_ivec_destroy(feat_counts);

  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
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
  #pragma omp parallel for
  for (uint64_t h = 0; h < n_hidden; h++) {
    double p = (double)bit_counts->a[h] / (double)n_samples;
    double entropy = 0.0;
    if (p > 0.0 && p < 1.0)
      entropy = -(p * log2(p) + (1.0 - p) * log2(1.0 - p));
    tk_rank_t r = { (int64_t)h, entropy };
    #pragma omp critical
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
          counts->a[f * n_hidden * 4 + j * 4 + (visible ? 2u : 0u) + (hidden ? 1u : 0u)]++;
        }
      }
    }
    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    #pragma omp parallel for
    for (uint64_t f = 0; f < n_features; f++) {
      double sum_mi = 0.0;
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
        sum_mi += mi;
      }
      tk_rank_t r = { (int64_t)f, sum_mi };
      #pragma omp critical
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
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
          feat_counts->a[f]++;
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
    tk_dvec_t *feat_sum_mi = tk_dvec_create(0, n_features, 0, 0);
    tk_dvec_zero(feat_sum_mi);

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
      double c[4];
      c[0] = count_00 + 1;
      c[1] = count_01 + 1;
      c[2] = count_10 + 1;
      c[3] = count_11 + 1;
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
      feat_sum_mi->a[f] += mi;
    }));

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_features; f++) {
      if (feat_sum_mi->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_sum_mi->a[f] };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_sum_mi);
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
    #pragma omp parallel for
    for (uint64_t f = 0; f < n_features; f++) {
      double sum_chi2 = 0.0;
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
        sum_chi2 += chi2;
      }
      tk_rank_t r = { (int64_t)f, sum_chi2 };
      #pragma omp critical
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
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
          feat_counts->a[f]++;
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
    tk_dvec_t *feat_sum_chi2 = tk_dvec_create(0, n_features, 0, 0);
    tk_dvec_zero(feat_sum_chi2);

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
      int64_t B = G - A;
      int64_t C_ = C - A;
      int64_t D = (int64_t)n_samples - C - B;
      double n = (double)(n_samples) + 4;
      double C_smooth = C + 2;
      double G_smooth = G + 2;
      double E_A = (C_smooth * G_smooth) / n;
      double E_B = ((n - C_smooth) * G_smooth) / n;
      double E_C = (C_smooth * (n - G_smooth)) / n;
      double E_D = ((n - C_smooth) * (n - G_smooth)) / n;
      double chi2 = 0.0;
      if (E_A > 0) chi2 += ((A - E_A) * (A - E_A)) / E_A;
      if (E_B > 0) chi2 += ((B - E_B) * (B - E_B)) / E_B;
      if (E_C > 0) chi2 += ((C_ - E_C) * (C_ - E_C)) / E_C;
      if (E_D > 0) chi2 += ((D - E_D) * (D - E_D)) / E_D;
      feat_sum_chi2->a[f] += chi2;
    }));

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_features; f++) {
      if (feat_sum_chi2->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_sum_chi2->a[f] };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_sum_chi2);

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
    #pragma omp parallel for
    for (uint64_t f = 0; f < n_visible; f++) {
      double sum_score = 0.0;
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
        double weighted_lift = lift * log2(1 + cooccur);
        sum_score += weighted_lift;
      }
      if (sum_score > 0) {
        tk_rank_t r = { (int64_t)f, sum_score };
        #pragma omp critical
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
    tk_dvec_t *feat_sum_score = tk_dvec_create(0, n_visible, 0, 0);
    tk_dvec_zero(feat_sum_score);

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
      double weighted_lift = lift * log2(1 + cooccur);
      feat_sum_score->a[f] += weighted_lift;
    }));
    #pragma omp parallel for
    for (uint64_t f = 0; f < n_visible; f++) {
      if (feat_sum_score->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_sum_score->a[f] };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_sum_score);

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



static inline tk_ivec_t *tk_cvec_bits_top_coherence (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  uint64_t top_k,
  double lambda
) {
  if (!codes) {
    tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return out;
  }

  tk_ivec_t *active_counts = tk_ivec_create(0, n_features * n_hidden, 0, 0);
  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_ivec_zero(active_counts);
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

  tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    int64_t k_f = feat_counts->a[f];
    if (k_f < 2)
      continue;
    double sum_disagreements = 0.0;
    for (uint64_t b = 0; b < n_hidden; b++) {
      int64_t n1 = active_counts->a[f * n_hidden + b];
      int64_t n0 = k_f - n1;
      sum_disagreements += (double)n1 * (double)n0;
    }
    double num_pairs = (double)k_f * (double)(k_f - 1) / 2.0;
    double mean_hamming = sum_disagreements / num_pairs;
    double normalized = mean_hamming / (double)n_hidden;
    double coherence = 1.0 - normalized;
    double penalty = lambda / sqrt(num_pairs);
    double score = coherence - penalty;
    if (score <= 0.0)
      continue;
    tk_rank_t r = { (int64_t)f, score };
    #pragma omp critical
    tk_rvec_hmin(top_heap, top_k, r);
  }

  tk_ivec_destroy(active_counts);
  tk_ivec_destroy(feat_counts);

  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
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
  #pragma omp parallel for
  for (uint64_t h = 0; h < n_hidden; h++) {
    double p = (double)bit_counts->a[h] / (double)n_samples;
    double entropy = 0.0;
    if (p > 0.0 && p < 1.0)
      entropy = -(p * log2(p) + (1.0 - p) * log2(1.0 - p));
    tk_rank_t r = { (int64_t)h, entropy };
    #pragma omp critical
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
  uint64_t n_visible,
  double min_df,
  double max_df,
  uint64_t top_k
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);
  tk_iuset_t **feature_docs = (tk_iuset_t **)calloc(n_visible, sizeof(tk_iuset_t *));
  if (!feature_docs)
    return NULL;
  for (uint64_t i = 0; i < n_visible; i++)
    feature_docs[i] = tk_iuset_create(0, 0);
  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit_idx = set_bits->a[i];
    if (bit_idx < 0)
      continue;
    uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
    uint64_t feature_idx = (uint64_t)bit_idx % n_visible;
    if (sample_idx < n_samples && feature_idx < n_visible) {
      int absent;
      tk_iuset_put(feature_docs[feature_idx], (int64_t)sample_idx, &absent);
    }
  }
  tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
  double min_df_abs = min_df < 0 ? -min_df : min_df * n_samples;
  double max_df_abs = max_df < 0 ? -max_df : max_df * n_samples;
  #pragma omp parallel for
  for (uint64_t i = 0; i < n_visible; i++) {
    double df_count = (double)tk_iuset_size(feature_docs[i]);
    double idf = log((double)(n_samples + 1) / (df_count + 1));
    if (df_count >= min_df_abs && df_count <= max_df_abs) {
      tk_rank_t r = { (int64_t)i, idf };
      #pragma omp critical
      tk_rvec_hmin(top_heap, top_k, r);
    }
  }
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  for (uint64_t i = 0; i < n_visible; i++)
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
    tk_dvec_t *feat_sum_score = tk_dvec_create(0, n_features, 0, 0);
    tk_dvec_zero(feat_sum_score);

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
      double weighted_lift = lift * log2(1 + cooccur);
      feat_sum_score->a[f] += weighted_lift;
    }));

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_features; f++) {
      if (feat_sum_score->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_sum_score->a[f] };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_sum_score);

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
    tk_dvec_t *feat_sum_score = tk_dvec_create(0, n_features, 0, 0);
    tk_dvec_zero(feat_sum_score);

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
      double weighted_lift = lift * log2(1 + cooccur);
      feat_sum_score->a[f] += weighted_lift;
    }));

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_features; f++) {
      if (feat_sum_score->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_sum_score->a[f] };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_sum_score);

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
  if (!feature_docs)
    return NULL;
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
  #pragma omp parallel for
  for (uint64_t i = 0; i < n_features; i++) {
    double df_count = (double)tk_iuset_size(feature_docs[i]);
    double idf = log((double)(n_samples + 1) / (df_count + 1));
    if (df_count >= min_df_abs && df_count <= max_df_abs) {
      tk_rank_t r = { (int64_t)i, idf };
      #pragma omp critical
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

static inline void tk_ivec_bits_top_chi2_ind (
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

    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_visible; f++) {
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
        if (chi2 > 0) {
          tk_rank_t r = { (int64_t)f, chi2 };
          #pragma omp critical
          tk_rvec_hmin(per_dim_heaps[b], top_k, r);
        }
      }
    }

    tk_ivec_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int absent;
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
      }
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({
      ids_union->a[ids_union->n++] = key;
    }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

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

    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    tk_dvec_t **feat_chi2_per_dim = (tk_dvec_t **)malloc(n_hidden * sizeof(tk_dvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++) {
      feat_chi2_per_dim[h] = tk_dvec_create(0, n_visible, 0, 0);
      tk_dvec_zero(feat_chi2_per_dim[h]);
    }

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
      int64_t B = G - A;
      int64_t C_ = C - A;
      int64_t D = (int64_t)n_samples - C - B;
      double n = (double)(n_samples) + 4;
      double C_smooth = C + 2;
      double G_smooth = G + 2;
      double E_A = (C_smooth * G_smooth) / n;
      double E_B = ((n - C_smooth) * G_smooth) / n;
      double E_C = (C_smooth * (n - G_smooth)) / n;
      double E_D = ((n - C_smooth) * (n - G_smooth)) / n;
      double chi2 = 0.0;
      if (E_A > 0) chi2 += ((A - E_A) * (A - E_A)) / E_A;
      if (E_B > 0) chi2 += ((B - E_B) * (B - E_B)) / E_B;
      if (E_C > 0) chi2 += ((C_ - E_C) * (C_ - E_C)) / E_C;
      if (E_D > 0) chi2 += ((D - E_D) * (D - E_D)) / E_D;
      feat_chi2_per_dim[b]->a[f] = chi2;
    }));

    #pragma omp parallel for
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t f = 0; f < n_visible; f++) {
        if (feat_chi2_per_dim[h]->a[f] > 0) {
          tk_rank_t r = { (int64_t)f, feat_chi2_per_dim[h]->a[f] };
          #pragma omp critical
          tk_rvec_hmin(per_dim_heaps[h], top_k, r);
        }
      }
    }

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_dvec_destroy(feat_chi2_per_dim[h]);
    free(feat_chi2_per_dim);

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int absent;
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
      }
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({
      ids_union->a[ids_union->n++] = key;
    }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

  } else {

    tk_ivec_create(L, 0, 0, 0);
    tk_ivec_create(L, 1, 0, 0)->a[0] = 0;
    tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);

  }
}

static inline void tk_ivec_bits_top_mi_ind (
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

    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_visible; f++) {
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
        if (mi > 0) {
          tk_rank_t r = { (int64_t)f, mi };
          #pragma omp critical
          tk_rvec_hmin(per_dim_heaps[j], top_k, r);
        }
      }
    }

    tk_ivec_destroy(counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int absent;
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
      }
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({
      ids_union->a[ids_union->n++] = key;
    }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

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

    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    tk_dvec_t **feat_mi_per_dim = (tk_dvec_t **)malloc(n_hidden * sizeof(tk_dvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++) {
      feat_mi_per_dim[h] = tk_dvec_create(0, n_visible, 0, 0);
      tk_dvec_zero(feat_mi_per_dim[h]);
    }

    int64_t k, v;
    tk_umap_foreach(active_counts, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t b = (uint64_t)k % n_hidden;
      if (f >= n_visible || b >= n_hidden)
        continue;
      int64_t c[4];
      c[3] = v;
      c[2] = feat_counts->a[f] - v;
      c[1] = label_counts->a[b] - v;
      c[0] = (int64_t)n_samples - c[1] - c[2] - c[3];
      for (int i = 0; i < 4; i++)
        c[i] += 1;
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
      feat_mi_per_dim[b]->a[f] = mi;
    }));

    #pragma omp parallel for
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t f = 0; f < n_visible; f++) {
        if (feat_mi_per_dim[h]->a[f] > 0) {
          tk_rank_t r = { (int64_t)f, feat_mi_per_dim[h]->a[f] };
          #pragma omp critical
          tk_rvec_hmin(per_dim_heaps[h], top_k, r);
        }
      }
    }

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_dvec_destroy(feat_mi_per_dim[h]);
    free(feat_mi_per_dim);

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int absent;
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
      }
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({
      ids_union->a[ids_union->n++] = key;
    }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

  } else {

    tk_ivec_create(L, 0, 0, 0);
    tk_ivec_create(L, 1, 0, 0)->a[0] = 0;
    tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);

  }
}

static inline void tk_cvec_bits_top_chi2_ind (
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
    uint64_t codes_bytes_per_sample = TK_CVEC_BITS_BYTES(n_hidden);
    uint8_t *codes_data = (uint8_t *)codes->a;

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
            active_counts->a[f * n_hidden + h]++;
          }
        }
      }
    }

    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_features; f++) {
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
        if (chi2 > 0) {
          tk_rank_t r = { (int64_t)f, chi2 };
          #pragma omp critical
          tk_rvec_hmin(per_dim_heaps[b], top_k, r);
        }
      }
    }

    tk_ivec_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int absent;
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
      }
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({
      ids_union->a[ids_union->n++] = key;
    }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

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
      size_t li = 0;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden < s) {
        li++;
      }
      if (li >= labels->n || labels->a[li] < 0 || (uint64_t)labels->a[li] / n_hidden != s)
        continue;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (!(bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)))
          continue;
        size_t li_inner = li;
        while (li_inner < labels->n && labels->a[li_inner] >= 0 &&
          (uint64_t)labels->a[li_inner] / n_hidden == s) {
          uint64_t h = (uint64_t)labels->a[li_inner] % n_hidden;
          int64_t key = (int64_t)(f * n_hidden + h);
          tk_iumap_inc(active_counts, key);
          li_inner++;
        }
      }
    }

    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    tk_dvec_t **feat_chi2_per_dim = (tk_dvec_t **)malloc(n_hidden * sizeof(tk_dvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++) {
      feat_chi2_per_dim[h] = tk_dvec_create(0, n_features, 0, 0);
      tk_dvec_zero(feat_chi2_per_dim[h]);
    }

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
      int64_t B = G - A;
      int64_t C_ = C - A;
      int64_t D = (int64_t)n_samples - C - B;
      double n = (double)(n_samples) + 4;
      double C_smooth = C + 2;
      double G_smooth = G + 2;
      double E_A = (C_smooth * G_smooth) / n;
      double E_B = ((n - C_smooth) * G_smooth) / n;
      double E_C = (C_smooth * (n - G_smooth)) / n;
      double E_D = ((n - C_smooth) * (n - G_smooth)) / n;
      double chi2 = 0.0;
      if (E_A > 0) chi2 += ((A - E_A) * (A - E_A)) / E_A;
      if (E_B > 0) chi2 += ((B - E_B) * (B - E_B)) / E_B;
      if (E_C > 0) chi2 += ((C_ - E_C) * (C_ - E_C)) / E_C;
      if (E_D > 0) chi2 += ((D - E_D) * (D - E_D)) / E_D;
      feat_chi2_per_dim[b]->a[f] = chi2;
    }));

    #pragma omp parallel for
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t f = 0; f < n_features; f++) {
        if (feat_chi2_per_dim[h]->a[f] > 0) {
          tk_rank_t r = { (int64_t)f, feat_chi2_per_dim[h]->a[f] };
          #pragma omp critical
          tk_rvec_hmin(per_dim_heaps[h], top_k, r);
        }
      }
    }

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_dvec_destroy(feat_chi2_per_dim[h]);
    free(feat_chi2_per_dim);

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int absent;
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
      }
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({
      ids_union->a[ids_union->n++] = key;
    }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

  } else {

    tk_ivec_create(L, 0, 0, 0);
    tk_ivec_create(L, 1, 0, 0)->a[0] = 0;
    tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);

  }
}

static inline void tk_ivec_bits_individualize (
  lua_State *L,
  tk_ivec_t *toks,
  tk_ivec_t *offsets,
  tk_ivec_t *ids,
  uint64_t union_size,
  uint64_t n_hidden
) {
  uint64_t n_samples = 0;
  for (uint64_t i = 0; i < toks->n; i++) {
    if (toks->a[i] >= 0) {
      uint64_t s = (uint64_t)toks->a[i] / union_size;
      if (s + 1 > n_samples)
        n_samples = s + 1;
    }
  }

  if (n_samples == 0 || n_hidden == 0 || union_size == 0) {
    tk_ivec_t *empty_toks = tk_ivec_create(L, 0, 0, 0);
    (void)empty_toks;
    tk_ivec_t *empty_offsets = tk_ivec_create(L, 1, 0, 0);
    empty_offsets->a[0] = 0;
    empty_offsets->n = 1;
    return;
  }

  tk_iumap_t **dim_maps = (tk_iumap_t **)malloc(n_hidden * sizeof(tk_iumap_t *));
  for (uint64_t h = 0; h < n_hidden; h++) {
    dim_maps[h] = tk_iumap_create(0, 0);
    int64_t start = offsets->a[h];
    int64_t end = offsets->a[h + 1];
    for (int64_t i = start; i < end; i++) {
      int64_t union_idx = ids->a[i];
      int absent;
      khint_t k = tk_iumap_put(dim_maps[h], union_idx, &absent);
      kh_value(dim_maps[h], k) = i - start;
    }
  }

  uint64_t total_slots = n_samples * n_hidden;
  int64_t *counts = (int64_t *)calloc(total_slots, sizeof(int64_t));

  for (uint64_t i = 0; i < toks->n; i++) {
    if (toks->a[i] < 0)
      continue;
    uint64_t sample = (uint64_t)toks->a[i] / union_size;
    int64_t union_idx = toks->a[i] % (int64_t)union_size;
    for (uint64_t h = 0; h < n_hidden; h++) {
      khint_t k = tk_iumap_get(dim_maps[h], union_idx);
      if (k != kh_end(dim_maps[h])) {
        counts[sample * n_hidden + h]++;
      }
    }
  }

  int64_t *temp_offsets = (int64_t *)malloc((total_slots + 1) * sizeof(int64_t));
  temp_offsets[0] = 0;
  for (uint64_t i = 0; i < total_slots; i++) {
    temp_offsets[i + 1] = temp_offsets[i] + counts[i];
  }

  uint64_t total_output = (uint64_t)temp_offsets[total_slots];
  tk_ivec_t *ind_toks = tk_ivec_create(L, total_output, 0, 0);
  ind_toks->n = total_output;

  tk_ivec_t *ind_offsets = tk_ivec_create(L, total_slots + 1, 0, 0);
  for (uint64_t i = 0; i <= total_slots; i++)
    ind_offsets->a[i] = temp_offsets[i];
  ind_offsets->n = total_slots + 1;

  int64_t *cursors = (int64_t *)malloc(total_slots * sizeof(int64_t));
  for (uint64_t i = 0; i < total_slots; i++)
    cursors[i] = temp_offsets[i];

  for (uint64_t i = 0; i < toks->n; i++) {
    if (toks->a[i] < 0)
      continue;
    uint64_t sample = (uint64_t)toks->a[i] / union_size;
    int64_t union_idx = toks->a[i] % (int64_t)union_size;
    for (uint64_t h = 0; h < n_hidden; h++) {
      khint_t k = tk_iumap_get(dim_maps[h], union_idx);
      if (k != kh_end(dim_maps[h])) {
        int64_t local_idx = kh_value(dim_maps[h], k);
        uint64_t slot = sample * n_hidden + h;
        ind_toks->a[cursors[slot]++] = local_idx;
      }
    }
  }

  free(temp_offsets);
  free(cursors);
  free(counts);
  for (uint64_t h = 0; h < n_hidden; h++)
    tk_iumap_destroy(dim_maps[h]);
  free(dim_maps);
}

static inline void tk_cvec_bits_from_ind (
  lua_State *L,
  tk_ivec_t *ind_toks,
  tk_ivec_t *ind_offsets,
  tk_ivec_t *feat_offsets,
  uint64_t n_samples,
  bool flip_interleave
) {
  uint64_t n_hidden = feat_offsets->n - 1;

  if (n_hidden == 0 || n_samples == 0) {
    tk_cvec_create(L, 0, 0, 0);
    tk_ivec_t *dim_offsets = tk_ivec_create(L, 1, 0, 0);
    dim_offsets->a[0] = 0;
    dim_offsets->n = 1;
    return;
  }

  uint64_t total_bytes = 0;
  for (uint64_t h = 0; h < n_hidden; h++) {
    uint64_t k_h = (uint64_t)(feat_offsets->a[h + 1] - feat_offsets->a[h]);
    uint64_t output_bits = flip_interleave ? (2 * k_h) : k_h;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(output_bits);
    total_bytes += n_samples * bytes_per_sample;
  }

  tk_cvec_t *bitmap = tk_cvec_create(L, total_bytes, 0, 0);
  memset(bitmap->a, 0, total_bytes);

  tk_ivec_t *dim_offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
  dim_offsets->n = n_hidden + 1;
  dim_offsets->a[0] = 0;

  uint64_t byte_offset = 0;
  for (uint64_t h = 0; h < n_hidden; h++) {
    uint64_t k_h = (uint64_t)(feat_offsets->a[h + 1] - feat_offsets->a[h]);
    uint64_t output_bits = flip_interleave ? (2 * k_h) : k_h;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(output_bits);

    if (flip_interleave) {
      for (uint64_t s = 0; s < n_samples; s++) {
        uint8_t *row = (uint8_t *)(bitmap->a + byte_offset + s * bytes_per_sample);
        for (uint64_t k = 0; k < k_h; k++) {
          uint64_t absent_bit = k_h + k;
          uint64_t byte_idx = absent_bit / CHAR_BIT;
          uint8_t bit_idx = absent_bit % CHAR_BIT;
          row[byte_idx] |= (1u << bit_idx);
        }
      }
    }

    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t ind_slot = s * n_hidden + h;
      int64_t start = ind_offsets->a[ind_slot];
      int64_t end = ind_offsets->a[ind_slot + 1];
      uint8_t *row = (uint8_t *)(bitmap->a + byte_offset + s * bytes_per_sample);

      for (int64_t i = start; i < end; i++) {
        int64_t local_idx = ind_toks->a[i];
        if (local_idx >= 0 && (uint64_t)local_idx < k_h) {
          uint64_t byte_idx = (uint64_t)local_idx / CHAR_BIT;
          uint8_t bit_idx = (uint8_t)(local_idx % CHAR_BIT);
          row[byte_idx] |= (1u << bit_idx);
          if (flip_interleave) {
            uint64_t absent_bit = k_h + (uint64_t)local_idx;
            uint64_t absent_byte = absent_bit / CHAR_BIT;
            uint8_t absent_bit_idx = absent_bit % CHAR_BIT;
            row[absent_byte] &= ~(1u << absent_bit_idx);
          }
        }
      }
    }

    byte_offset += n_samples * bytes_per_sample;
    dim_offsets->a[h + 1] = (int64_t)byte_offset;
  }
}

static inline int tk_ivec_bits_extend_ind (
  lua_State *L,
  tk_ivec_t *base_toks,
  tk_ivec_t *base_offsets,
  tk_ivec_t *base_ids,
  tk_ivec_t *ext_toks,
  tk_ivec_t *ext_offsets,
  tk_ivec_t *ext_ids,
  tk_ivec_t *base_feat_offsets,
  tk_ivec_t *ext_feat_offsets,
  uint64_t n_hidden,
  bool project
) {
  if (!base_toks || !base_offsets || !base_ids || !ext_toks || !ext_offsets || !ext_ids || !base_feat_offsets || !ext_feat_offsets)
    return -1;

  if (n_hidden == 0 || base_feat_offsets->n != n_hidden + 1 || ext_feat_offsets->n != n_hidden + 1) {
    tk_ivec_t *merged_offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    for (uint64_t i = 0; i <= n_hidden; i++)
      merged_offsets->a[i] = 0;
    merged_offsets->n = n_hidden + 1;
    return 0;
  }

  uint64_t n_base_samples = 0;
  if (base_offsets->n > 0)
    n_base_samples = (base_offsets->n - 1) / n_hidden;
  uint64_t n_ext_samples = 0;
  if (ext_offsets->n > 0)
    n_ext_samples = (ext_offsets->n - 1) / n_hidden;

  tk_iumap_t *base_id_to_pos = tk_iumap_from_ivec(0, base_ids);
  if (!base_id_to_pos)
    return -1;

  tk_iumap_t *ext_id_to_pos = tk_iumap_from_ivec(0, ext_ids);
  if (!ext_id_to_pos) {
    tk_iumap_destroy(base_id_to_pos);
    return -1;
  }

  uint64_t final_n_samples = 0;
  int64_t *base_to_final = (int64_t *)calloc(n_base_samples, sizeof(int64_t));
  int64_t *ext_to_final = (int64_t *)calloc(n_ext_samples, sizeof(int64_t));
  if (!base_to_final || !ext_to_final) {
    free(base_to_final);
    free(ext_to_final);
    tk_iumap_destroy(base_id_to_pos);
    tk_iumap_destroy(ext_id_to_pos);
    return -1;
  }

  for (uint64_t i = 0; i < n_base_samples; i++)
    base_to_final[i] = -1;
  for (uint64_t i = 0; i < n_ext_samples; i++)
    ext_to_final[i] = -1;

  if (project) {
    for (uint64_t ei = 0; ei < n_ext_samples; ei++) {
      if (ei >= ext_ids->n)
        break;
      int64_t ext_id = ext_ids->a[ei];
      khint_t k = tk_iumap_get(base_id_to_pos, ext_id);
      if (k != tk_iumap_end(base_id_to_pos)) {
        uint64_t base_pos = (uint64_t)tk_iumap_val(base_id_to_pos, k);
        if (base_to_final[base_pos] < 0) {
          base_to_final[base_pos] = (int64_t)final_n_samples;
          ext_to_final[ei] = (int64_t)final_n_samples;
          final_n_samples++;
        }
      }
    }
  } else {
    for (uint64_t bi = 0; bi < n_base_samples; bi++) {
      if (bi >= base_ids->n)
        break;
      base_to_final[bi] = (int64_t)final_n_samples++;
    }
    for (uint64_t ei = 0; ei < n_ext_samples; ei++) {
      if (ei >= ext_ids->n)
        break;
      int64_t ext_id = ext_ids->a[ei];
      khint_t k = tk_iumap_get(base_id_to_pos, ext_id);
      if (k != tk_iumap_end(base_id_to_pos)) {
        uint64_t base_pos = (uint64_t)tk_iumap_val(base_id_to_pos, k);
        ext_to_final[ei] = base_to_final[base_pos];
      } else {
        ext_to_final[ei] = (int64_t)final_n_samples++;
      }
    }
  }

  uint64_t total_slots = final_n_samples * n_hidden;
  int64_t *counts = (int64_t *)calloc(total_slots, sizeof(int64_t));
  if (!counts) {
    free(base_to_final);
    free(ext_to_final);
    tk_iumap_destroy(base_id_to_pos);
    tk_iumap_destroy(ext_id_to_pos);
    return -1;
  }

  for (uint64_t bi = 0; bi < n_base_samples; bi++) {
    if (base_to_final[bi] < 0)
      continue;
    uint64_t final_pos = (uint64_t)base_to_final[bi];
    for (uint64_t h = 0; h < n_hidden; h++) {
      uint64_t base_slot = bi * n_hidden + h;
      if (base_slot + 1 >= base_offsets->n)
        continue;
      int64_t start = base_offsets->a[base_slot];
      int64_t end = base_offsets->a[base_slot + 1];
      counts[final_pos * n_hidden + h] += (end - start);
    }
  }

  for (uint64_t ei = 0; ei < n_ext_samples; ei++) {
    if (ext_to_final[ei] < 0)
      continue;
    uint64_t final_pos = (uint64_t)ext_to_final[ei];
    for (uint64_t h = 0; h < n_hidden; h++) {
      uint64_t ext_slot = ei * n_hidden + h;
      if (ext_slot + 1 >= ext_offsets->n)
        continue;
      int64_t start = ext_offsets->a[ext_slot];
      int64_t end = ext_offsets->a[ext_slot + 1];
      counts[final_pos * n_hidden + h] += (end - start);
    }
  }

  int64_t *temp_offsets = (int64_t *)malloc((total_slots + 1) * sizeof(int64_t));
  if (!temp_offsets) {
    free(counts);
    free(base_to_final);
    free(ext_to_final);
    tk_iumap_destroy(base_id_to_pos);
    tk_iumap_destroy(ext_id_to_pos);
    return -1;
  }

  temp_offsets[0] = 0;
  for (uint64_t i = 0; i < total_slots; i++)
    temp_offsets[i + 1] = temp_offsets[i] + counts[i];

  uint64_t total_output = (uint64_t)temp_offsets[total_slots];
  tk_ivec_t *merged_toks = tk_ivec_create(L, total_output, 0, 0);
  merged_toks->n = total_output;

  tk_ivec_t *merged_offsets = tk_ivec_create(L, total_slots + 1, 0, 0);
  for (uint64_t i = 0; i <= total_slots; i++)
    merged_offsets->a[i] = temp_offsets[i];
  merged_offsets->n = total_slots + 1;

  int64_t *cursors = (int64_t *)malloc(total_slots * sizeof(int64_t));
  if (!cursors) {
    free(temp_offsets);
    free(counts);
    free(base_to_final);
    free(ext_to_final);
    tk_iumap_destroy(base_id_to_pos);
    tk_iumap_destroy(ext_id_to_pos);
    return -1;
  }

  for (uint64_t i = 0; i < total_slots; i++)
    cursors[i] = temp_offsets[i];

  for (uint64_t bi = 0; bi < n_base_samples; bi++) {
    if (base_to_final[bi] < 0)
      continue;
    uint64_t final_pos = (uint64_t)base_to_final[bi];
    for (uint64_t h = 0; h < n_hidden; h++) {
      uint64_t base_slot = bi * n_hidden + h;
      if (base_slot + 1 >= base_offsets->n)
        continue;
      int64_t start = base_offsets->a[base_slot];
      int64_t end = base_offsets->a[base_slot + 1];
      uint64_t final_slot = final_pos * n_hidden + h;
      for (int64_t i = start; i < end; i++) {
        merged_toks->a[cursors[final_slot]++] = base_toks->a[i];
      }
    }
  }

  for (uint64_t ei = 0; ei < n_ext_samples; ei++) {
    if (ext_to_final[ei] < 0)
      continue;
    uint64_t final_pos = (uint64_t)ext_to_final[ei];
    for (uint64_t h = 0; h < n_hidden; h++) {
      uint64_t ext_slot = ei * n_hidden + h;
      if (ext_slot + 1 >= ext_offsets->n)
        continue;
      int64_t start = ext_offsets->a[ext_slot];
      int64_t end = ext_offsets->a[ext_slot + 1];
      uint64_t final_slot = final_pos * n_hidden + h;
      int64_t base_feat_offset = base_feat_offsets->a[h];
      int64_t ext_feat_offset = ext_feat_offsets->a[h];
      int64_t base_n_feats = base_feat_offsets->a[h + 1] - base_feat_offset;
      for (int64_t i = start; i < end; i++) {
        int64_t local_idx = ext_toks->a[i] - ext_feat_offset;
        merged_toks->a[cursors[final_slot]++] = base_n_feats + local_idx;
      }
    }
  }

  tk_ivec_t *merged_feat_offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
  merged_feat_offsets->a[0] = 0;
  for (uint64_t h = 0; h < n_hidden; h++) {
    int64_t base_n = base_feat_offsets->a[h + 1] - base_feat_offsets->a[h];
    int64_t ext_n = ext_feat_offsets->a[h + 1] - ext_feat_offsets->a[h];
    merged_feat_offsets->a[h + 1] = merged_feat_offsets->a[h] + base_n + ext_n;
  }
  merged_feat_offsets->n = n_hidden + 1;

  base_toks->n = 0;
  base_offsets->n = 0;
  for (uint64_t i = 0; i < merged_toks->n; i++) {
    if (tk_ivec_ensure(base_toks, i + 1) != 0) {
      free(cursors);
      free(temp_offsets);
      free(counts);
      free(base_to_final);
      free(ext_to_final);
      tk_iumap_destroy(base_id_to_pos);
      tk_iumap_destroy(ext_id_to_pos);
      return -1;
    }
    base_toks->a[i] = merged_toks->a[i];
    base_toks->n = i + 1;
  }

  for (uint64_t i = 0; i < merged_offsets->n; i++) {
    if (tk_ivec_ensure(base_offsets, i + 1) != 0) {
      free(cursors);
      free(temp_offsets);
      free(counts);
      free(base_to_final);
      free(ext_to_final);
      tk_iumap_destroy(base_id_to_pos);
      tk_iumap_destroy(ext_id_to_pos);
      return -1;
    }
    base_offsets->a[i] = merged_offsets->a[i];
    base_offsets->n = i + 1;
  }

  if (project) {
    base_ids->n = 0;
    for (uint64_t bi = 0; bi < n_base_samples; bi++) {
      if (base_to_final[bi] >= 0) {
        if (tk_ivec_ensure(base_ids, base_ids->n + 1) != 0) {
          free(cursors);
          free(temp_offsets);
          free(counts);
          free(base_to_final);
          free(ext_to_final);
          tk_iumap_destroy(base_id_to_pos);
          tk_iumap_destroy(ext_id_to_pos);
          return -1;
        }
        base_ids->a[base_ids->n++] = base_ids->a[bi];
      }
    }
  } else {
    uint64_t orig_base_n = base_ids->n;
    for (uint64_t ei = 0; ei < n_ext_samples; ei++) {
      if (ext_to_final[ei] >= 0 && (uint64_t)ext_to_final[ei] >= orig_base_n) {
        if (tk_ivec_ensure(base_ids, base_ids->n + 1) != 0) {
          free(cursors);
          free(temp_offsets);
          free(counts);
          free(base_to_final);
          free(ext_to_final);
          tk_iumap_destroy(base_id_to_pos);
          tk_iumap_destroy(ext_id_to_pos);
          return -1;
        }
        if (ei < ext_ids->n)
          base_ids->a[base_ids->n++] = ext_ids->a[ei];
      }
    }
  }

  for (uint64_t h = 0; h < n_hidden + 1; h++) {
    if (tk_ivec_ensure(base_feat_offsets, h + 1) != 0) {
      free(cursors);
      free(temp_offsets);
      free(counts);
      free(base_to_final);
      free(ext_to_final);
      tk_iumap_destroy(base_id_to_pos);
      tk_iumap_destroy(ext_id_to_pos);
      return -1;
    }
    base_feat_offsets->a[h] = merged_feat_offsets->a[h];
  }
  base_feat_offsets->n = n_hidden + 1;

  free(cursors);
  free(temp_offsets);
  free(counts);
  free(base_to_final);
  free(ext_to_final);
  tk_iumap_destroy(base_id_to_pos);
  tk_iumap_destroy(ext_id_to_pos);

  return 0;
}








static inline tk_ivec_t *tk_ivec_bits_top_reg_f (
  lua_State *L,
  tk_ivec_t *set_bits,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t top_k
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);

  double overall_sum = 0.0;
  for (uint64_t i = 0; i < n_samples; i++)
    overall_sum += targets->a[i];
  double overall_mean = overall_sum / (double)n_samples;

  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_dvec_t *feat_sums = tk_dvec_create(0, n_features, 0, 0);
  tk_dvec_t *feat_sum_sq = tk_dvec_create(0, n_features, 0, 0);
  tk_ivec_zero(feat_counts);
  tk_dvec_zero(feat_sums);
  tk_dvec_zero(feat_sum_sq);

  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit = set_bits->a[i];
    if (bit < 0) continue;
    uint64_t sample = (uint64_t)bit / n_features;
    uint64_t feature = (uint64_t)bit % n_features;
    if (sample >= n_samples || feature >= n_features) continue;
    double y = targets->a[sample];
    feat_counts->a[feature]++;
    feat_sums->a[feature] += y;
    feat_sum_sq->a[feature] += y * y;
  }

  double total_sum_sq = 0.0;
  for (uint64_t i = 0; i < n_samples; i++) {
    double y = targets->a[i];
    total_sum_sq += y * y;
  }

  tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    int64_t n1 = feat_counts->a[f];
    int64_t n0 = (int64_t)n_samples - n1;
    if (n1 <= 1 || n0 <= 1) continue;

    double sum1 = feat_sums->a[f];
    double sum0 = overall_sum - sum1;
    double mean1 = sum1 / (double)n1;
    double mean0 = sum0 / (double)n0;

    double ssb = (double)n1 * (mean1 - overall_mean) * (mean1 - overall_mean) +
                 (double)n0 * (mean0 - overall_mean) * (mean0 - overall_mean);

    double sum_sq1 = feat_sum_sq->a[f];
    double sum_sq0 = total_sum_sq - sum_sq1;
    double ssw1 = sum_sq1 - (sum1 * sum1) / (double)n1;
    double ssw0 = sum_sq0 - (sum0 * sum0) / (double)n0;
    double ssw = ssw1 + ssw0;

    if (ssw < 1e-12) continue;
    double F = ssb / (ssw / (double)(n_samples - 2));

    tk_rank_t r = { (int64_t)f, F };
    #pragma omp critical
    tk_rvec_hmin(top_heap, top_k, r);
  }

  tk_ivec_destroy(feat_counts);
  tk_dvec_destroy(feat_sums);
  tk_dvec_destroy(feat_sum_sq);

  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0, 0, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}



static inline tk_ivec_t *tk_ivec_bits_top_reg_pearson (
  lua_State *L,
  tk_ivec_t *set_bits,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t top_k
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);

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

  tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
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
    #pragma omp critical
    tk_rvec_hmin(top_heap, top_k, rk);
  }

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
  tk_ivec_asc(set_bits, 0, set_bits->n);

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

  tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
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
    #pragma omp critical
    tk_rvec_hmin(top_heap, top_k, r);
  }

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





static inline tk_ivec_t *tk_cvec_bits_top_reg_f (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t top_k
) {
  double overall_sum = 0.0;
  for (uint64_t i = 0; i < n_samples; i++)
    overall_sum += targets->a[i];
  double overall_mean = overall_sum / (double)n_samples;

  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_dvec_t *feat_sums = tk_dvec_create(0, n_features, 0, 0);
  tk_dvec_t *feat_sum_sq = tk_dvec_create(0, n_features, 0, 0);
  tk_ivec_zero(feat_counts);
  tk_dvec_zero(feat_sums);
  tk_dvec_zero(feat_sum_sq);

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
        feat_sum_sq->a[f] += y * y;
      }
    }
  }

  double total_sum_sq = 0.0;
  for (uint64_t i = 0; i < n_samples; i++) {
    double y = targets->a[i];
    total_sum_sq += y * y;
  }

  tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    int64_t n1 = feat_counts->a[f];
    int64_t n0 = (int64_t)n_samples - n1;
    if (n1 <= 1 || n0 <= 1) continue;

    double sum1 = feat_sums->a[f];
    double sum0 = overall_sum - sum1;
    double mean1 = sum1 / (double)n1;
    double mean0 = sum0 / (double)n0;

    double ssb = (double)n1 * (mean1 - overall_mean) * (mean1 - overall_mean) +
                 (double)n0 * (mean0 - overall_mean) * (mean0 - overall_mean);

    double sum_sq1 = feat_sum_sq->a[f];
    double sum_sq0 = total_sum_sq - sum_sq1;
    double ssw1 = sum_sq1 - (sum1 * sum1) / (double)n1;
    double ssw0 = sum_sq0 - (sum0 * sum0) / (double)n0;
    double ssw = ssw1 + ssw0;

    if (ssw < 1e-12) continue;
    double F = ssb / (ssw / (double)(n_samples - 2));

    tk_rank_t r = { (int64_t)f, F };
    #pragma omp critical
    tk_rvec_hmin(top_heap, top_k, r);
  }

  tk_ivec_destroy(feat_counts);
  tk_dvec_destroy(feat_sums);
  tk_dvec_destroy(feat_sum_sq);

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

  tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
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
    #pragma omp critical
    tk_rvec_hmin(top_heap, top_k, rk);
  }

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

  tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
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
    #pragma omp critical
    tk_rvec_hmin(top_heap, top_k, r);
  }

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

static inline void tk_cvec_bits_top_mi_ind (
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
    uint64_t codes_bytes_per_sample = TK_CVEC_BITS_BYTES(n_hidden);
    uint8_t *codes_data = (uint8_t *)codes->a;

    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      uint64_t codes_offset = s * codes_bytes_per_sample;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t f_byte_idx = f / CHAR_BIT;
        uint8_t f_bit_idx = f % CHAR_BIT;
        bool visible = (bitmap_data[sample_offset + f_byte_idx] & (1u << f_bit_idx)) != 0;
        for (uint64_t h = 0; h < n_hidden; h++) {
          uint64_t h_byte_idx = h / CHAR_BIT;
          uint8_t h_bit_idx = h % CHAR_BIT;
          bool hidden = (codes_data[codes_offset + h_byte_idx] & (1u << h_bit_idx)) != 0;
          unsigned int cell = (visible ? 2u : 0u) + (hidden ? 1u : 0u);
          counts->a[f * n_hidden * 4 + h * 4 + cell]++;
        }
      }
    }

    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_features; f++) {
      for (uint64_t h = 0; h < n_hidden; h++) {
        int64_t c[4];
        int64_t *counts_ptr = counts->a + f * n_hidden * 4 + h * 4;
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
        if (mi > 0) {
          tk_rank_t r = { (int64_t)f, mi };
          #pragma omp critical
          tk_rvec_hmin(per_dim_heaps[h], top_k, r);
        }
      }
    }

    tk_ivec_destroy(counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int absent;
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
      }
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({
      ids_union->a[ids_union->n++] = key;
    }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

  } else if (labels) {

    tk_ivec_asc(labels, 0, labels->n);
    tk_iumap_t *counts_map = tk_iumap_create(0, 0);
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

    size_t li = 0;
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden < s)
        li++;
      if (li >= labels->n || labels->a[li] < 0 || (uint64_t)labels->a[li] / n_hidden != s)
        continue;
      size_t li_start = li;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (!(bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)))
          continue;
        li = li_start;
        while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden == s) {
          uint64_t h = (uint64_t)labels->a[li] % n_hidden;
          int64_t key = (int64_t)(f * n_hidden + h);
          tk_iumap_inc(counts_map, key);
          li++;
        }
      }
    }

    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    for (uint64_t f = 0; f < n_features; f++) {
      int64_t n_f = feat_counts->a[f];
      if (n_f == 0) continue;
      for (uint64_t h = 0; h < n_hidden; h++) {
        int64_t n_h = label_counts->a[h];
        if (n_h == 0) continue;
        int64_t key = (int64_t)(f * n_hidden + h);
        khint_t k = tk_iumap_get(counts_map, key);
        int64_t n_fh = (k != kh_end(counts_map)) ? kh_value(counts_map, k) : 0;
        int64_t c[4];
        c[0] = (int64_t)n_samples - n_f - n_h + n_fh + 1;
        c[1] = n_h - n_fh + 1;
        c[2] = n_f - n_fh + 1;
        c[3] = n_fh + 1;
        double total = c[0] + c[1] + c[2] + c[3];
        double mi = 0.0;
        for (unsigned int o = 0; o < 4; o++) {
          if (c[o] == 0) continue;
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
        if (mi > 0) {
          tk_rank_t r = { (int64_t)f, mi };
          tk_rvec_hmin(per_dim_heaps[h], top_k, r);
        }
      }
    }

    tk_iumap_destroy(counts_map);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int absent;
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
      }
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t ukey;
    tk_umap_foreach_keys(union_set, ukey, ({
      ids_union->a[ids_union->n++] = ukey;
    }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

  } else {

    tk_ivec_create(L, 0, 0, 0);
    tk_ivec_create(L, 1, 0, 0);
    tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);

  }
}

static inline void tk_ivec_bits_top_coherence_ind (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k,
  double lambda
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);

  if (!codes) {
    tk_ivec_create(L, 0, 0, 0);
    tk_ivec_create(L, 1, 0, 0);
    tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return;
  }

  tk_ivec_t *active_counts = tk_ivec_create(0, n_visible * n_hidden, 0, 0);
  tk_ivec_t *feat_counts = tk_ivec_create(0, n_visible, 0, 0);
  tk_ivec_zero(active_counts);
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
      if (sample_codes[byte_idx] & (1u << bit_pos))
        active_counts->a[feature_idx * n_hidden + b]++;
    }
  }

  tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
  for (uint64_t h = 0; h < n_hidden; h++)
    per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
  for (uint64_t f = 0; f < n_visible; f++) {
    int64_t k_f = feat_counts->a[f];
    if (k_f < 2)
      continue;
    double num_pairs = (double)k_f * (double)(k_f - 1) / 2.0;
    double penalty = lambda / sqrt(num_pairs);
    for (uint64_t b = 0; b < n_hidden; b++) {
      int64_t n1 = active_counts->a[f * n_hidden + b];
      int64_t n0 = k_f - n1;
      double disagreements = (double)n1 * (double)n0;
      double normalized = disagreements / num_pairs;
      double coherence = 1.0 - normalized;
      double score = coherence - penalty;
      if (score <= 0.0)
        continue;
      tk_rank_t r = { (int64_t)f, score };
      #pragma omp critical
      tk_rvec_hmin(per_dim_heaps[b], top_k, r);
    }
  }

  tk_ivec_destroy(active_counts);
  tk_ivec_destroy(feat_counts);

  for (uint64_t h = 0; h < n_hidden; h++)
    tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

  tk_iuset_t *union_set = tk_iuset_create(0, 0);
  for (uint64_t h = 0; h < n_hidden; h++) {
    for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
      int absent;
      tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
    }
  }

  tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
  ids_union->n = 0;
  int64_t key;
  tk_umap_foreach_keys(union_set, key, ({
    ids_union->a[ids_union->n++] = key;
  }));
  tk_ivec_asc(ids_union, 0, ids_union->n);

  tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
  for (uint64_t i = 0; i < ids_union->n; i++) {
    int absent;
    khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
    kh_value(id_to_union_idx, k) = (int64_t)i;
  }

  uint64_t total_ids = 0;
  for (uint64_t h = 0; h < n_hidden; h++)
    total_ids += per_dim_heaps[h]->n;

  tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
  tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

  offsets->a[0] = 0;
  ids->n = 0;
  weights->n = 0;
  for (uint64_t h = 0; h < n_hidden; h++) {
    for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
      int64_t feat_id = per_dim_heaps[h]->a[i].i;
      khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
      ids->a[ids->n++] = kh_value(id_to_union_idx, k);
      weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
    }
    offsets->a[h + 1] = (int64_t)ids->n;
  }
  offsets->n = n_hidden + 1;

  tk_iumap_destroy(id_to_union_idx);
  tk_iuset_destroy(union_set);
  for (uint64_t h = 0; h < n_hidden; h++)
    tk_rvec_destroy(per_dim_heaps[h]);
  free(per_dim_heaps);
}

static inline void tk_cvec_bits_top_coherence_ind (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  uint64_t top_k,
  double lambda
) {
  if (!codes) {
    tk_ivec_create(L, 0, 0, 0);
    tk_ivec_create(L, 1, 0, 0);
    tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return;
  }

  tk_ivec_t *active_counts = tk_ivec_create(0, n_features * n_hidden, 0, 0);
  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_ivec_zero(active_counts);
  tk_ivec_zero(feat_counts);

  uint8_t *bitmap_data = (uint8_t *)bitmap->a;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t codes_bytes_per_sample = TK_CVEC_BITS_BYTES(n_hidden);
  uint8_t *codes_data = (uint8_t *)codes->a;

  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t sample_offset = s * bytes_per_sample;
    uint64_t codes_offset = s * codes_bytes_per_sample;
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t f_byte_idx = f / CHAR_BIT;
      uint8_t f_bit_idx = f % CHAR_BIT;
      if (!(bitmap_data[sample_offset + f_byte_idx] & (1u << f_bit_idx)))
        continue;
      feat_counts->a[f]++;
      for (uint64_t h = 0; h < n_hidden; h++) {
        uint64_t h_byte_idx = h / CHAR_BIT;
        uint8_t h_bit_idx = h % CHAR_BIT;
        if (codes_data[codes_offset + h_byte_idx] & (1u << h_bit_idx))
          active_counts->a[f * n_hidden + h]++;
      }
    }
  }

  tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
  for (uint64_t h = 0; h < n_hidden; h++)
    per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    int64_t k_f = feat_counts->a[f];
    if (k_f < 2)
      continue;
    double num_pairs = (double)k_f * (double)(k_f - 1) / 2.0;
    double penalty = lambda / sqrt(num_pairs);
    for (uint64_t b = 0; b < n_hidden; b++) {
      int64_t n1 = active_counts->a[f * n_hidden + b];
      int64_t n0 = k_f - n1;
      double disagreements = (double)n1 * (double)n0;
      double normalized = disagreements / num_pairs;
      double coherence = 1.0 - normalized;
      double score = coherence - penalty;
      if (score <= 0.0)
        continue;
      tk_rank_t r = { (int64_t)f, score };
      #pragma omp critical
      tk_rvec_hmin(per_dim_heaps[b], top_k, r);
    }
  }

  tk_ivec_destroy(active_counts);
  tk_ivec_destroy(feat_counts);

  for (uint64_t h = 0; h < n_hidden; h++)
    tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

  tk_iuset_t *union_set = tk_iuset_create(0, 0);
  for (uint64_t h = 0; h < n_hidden; h++) {
    for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
      int absent;
      tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
    }
  }

  tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
  ids_union->n = 0;
  int64_t key;
  tk_umap_foreach_keys(union_set, key, ({
    ids_union->a[ids_union->n++] = key;
  }));
  tk_ivec_asc(ids_union, 0, ids_union->n);

  tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
  for (uint64_t i = 0; i < ids_union->n; i++) {
    int absent;
    khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
    kh_value(id_to_union_idx, k) = (int64_t)i;
  }

  uint64_t total_ids = 0;
  for (uint64_t h = 0; h < n_hidden; h++)
    total_ids += per_dim_heaps[h]->n;

  tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
  tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

  offsets->a[0] = 0;
  ids->n = 0;
  weights->n = 0;
  for (uint64_t h = 0; h < n_hidden; h++) {
    for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
      int64_t feat_id = per_dim_heaps[h]->a[i].i;
      khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
      ids->a[ids->n++] = kh_value(id_to_union_idx, k);
      weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
    }
    offsets->a[h + 1] = (int64_t)ids->n;
  }
  offsets->n = n_hidden + 1;

  tk_iumap_destroy(id_to_union_idx);
  tk_iuset_destroy(union_set);
  for (uint64_t h = 0; h < n_hidden; h++)
    tk_rvec_destroy(per_dim_heaps[h]);
  free(per_dim_heaps);
}

static inline void tk_ivec_bits_top_lift_ind (
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

    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_visible; f++) {
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
        double weighted_lift = lift * log2(1 + cooccur);
        tk_rank_t r = { (int64_t)f, weighted_lift };
        #pragma omp critical
        tk_rvec_hmin(per_dim_heaps[h], top_k, r);
      }
    }

    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_ivec_destroy(cooccur_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int absent;
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
      }
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({
      ids_union->a[ids_union->n++] = key;
    }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

  } else if (labels) {

    tk_ivec_asc(labels, 0, labels->n);
    tk_iumap_t *cooccur_map = tk_iumap_create(0, 0);
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
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden < s_sample)
        li++;
      if (li >= labels->n || labels->a[li] < 0 || (uint64_t)labels->a[li] / n_hidden > s_sample) {
        si++;
        continue;
      }
      size_t li_start = li;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden == s_sample) {
        uint64_t h = (uint64_t)labels->a[li] % n_hidden;
        int64_t key = (int64_t)(f * n_hidden + h);
        tk_iumap_inc(cooccur_map, key);
        li++;
      }
      si++;
      if (si < set_bits->n && set_bits->a[si] >= 0 && (uint64_t)set_bits->a[si] / n_visible == s_sample)
        li = li_start;
    }

    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    int64_t k, v;
    tk_umap_foreach(cooccur_map, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t h = (uint64_t)k % n_hidden;
      if (f >= n_visible || h >= n_hidden)
        continue;
      int64_t feat_total = feat_counts->a[f];
      int64_t label_total = label_counts->a[h];
      if (feat_total == 0 || label_total == 0)
        continue;
      double p_label_given_feature = (double)v / feat_total;
      double p_label = (double)label_total / n_samples;
      double lift = p_label_given_feature / (p_label + 1e-10);
      double weighted_lift = lift * log2(1 + v);
      tk_rank_t r = { (int64_t)f, weighted_lift };
      tk_rvec_hmin(per_dim_heaps[h], top_k, r);
    }));

    tk_iumap_destroy(cooccur_map);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int absent;
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
      }
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t ukey;
    tk_umap_foreach_keys(union_set, ukey, ({
      ids_union->a[ids_union->n++] = ukey;
    }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t ki = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, ki) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t ki = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, ki);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

  } else {

    tk_ivec_create(L, 0, 0, 0);
    tk_ivec_create(L, 1, 0, 0);
    tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);

  }
}

static inline void tk_cvec_bits_top_lift_ind (
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

    tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
    tk_ivec_t *label_counts = tk_ivec_create(0, n_hidden, 0, 0);
    tk_ivec_t *cooccur_counts = tk_ivec_create(0, n_features * n_hidden, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_zero(label_counts);
    tk_ivec_zero(cooccur_counts);

    uint8_t *bitmap_data = (uint8_t *)bitmap->a;
    uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
    uint64_t codes_bytes_per_sample = TK_CVEC_BITS_BYTES(n_hidden);
    uint8_t *codes_data = (uint8_t *)codes->a;

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
          if (codes_data[codes_offset + h_byte_idx] & (1u << h_bit_idx))
            cooccur_counts->a[f * n_hidden + h]++;
        }
      }
    }

    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_features; f++) {
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
        double weighted_lift = lift * log2(1 + cooccur);
        tk_rank_t r = { (int64_t)f, weighted_lift };
        #pragma omp critical
        tk_rvec_hmin(per_dim_heaps[h], top_k, r);
      }
    }

    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_ivec_destroy(cooccur_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int absent;
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
      }
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({
      ids_union->a[ids_union->n++] = key;
    }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

  } else if (labels) {

    tk_ivec_asc(labels, 0, labels->n);
    tk_iumap_t *cooccur_map = tk_iumap_create(0, 0);
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

    size_t li = 0;
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden < s)
        li++;
      if (li >= labels->n || labels->a[li] < 0 || (uint64_t)labels->a[li] / n_hidden != s)
        continue;
      size_t li_start = li;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (!(bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)))
          continue;
        li = li_start;
        while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden == s) {
          uint64_t h = (uint64_t)labels->a[li] % n_hidden;
          int64_t key = (int64_t)(f * n_hidden + h);
          tk_iumap_inc(cooccur_map, key);
          li++;
        }
      }
    }

    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    int64_t k, v;
    tk_umap_foreach(cooccur_map, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t h = (uint64_t)k % n_hidden;
      if (f >= n_features || h >= n_hidden)
        continue;
      int64_t feat_total = feat_counts->a[f];
      int64_t label_total = label_counts->a[h];
      if (feat_total == 0 || label_total == 0)
        continue;
      double p_label_given_feature = (double)v / feat_total;
      double p_label = (double)label_total / n_samples;
      double lift = p_label_given_feature / (p_label + 1e-10);
      double weighted_lift = lift * log2(1 + v);
      tk_rank_t r = { (int64_t)f, weighted_lift };
      tk_rvec_hmin(per_dim_heaps[h], top_k, r);
    }));

    tk_iumap_destroy(cooccur_map);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int absent;
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
      }
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t ukey;
    tk_umap_foreach_keys(union_set, ukey, ({
      ids_union->a[ids_union->n++] = ukey;
    }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t ki = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, ki) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t ki = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, ki);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

  } else {

    tk_ivec_create(L, 0, 0, 0);
    tk_ivec_create(L, 1, 0, 0);
    tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);

  }
}

static inline void tk_ivec_bits_top_reg_f_ind (
  lua_State *L,
  tk_ivec_t *set_bits,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_targets,
  uint64_t top_k
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);

  tk_dvec_t *overall_sums = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_t *overall_means = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_t *total_sum_sqs = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_zero(overall_sums);
  tk_dvec_zero(total_sum_sqs);

  for (uint64_t i = 0; i < n_samples; i++) {
    for (uint64_t t = 0; t < n_targets; t++) {
      double y = targets->a[i * n_targets + t];
      overall_sums->a[t] += y;
      total_sum_sqs->a[t] += y * y;
    }
  }
  for (uint64_t t = 0; t < n_targets; t++)
    overall_means->a[t] = overall_sums->a[t] / (double)n_samples;

  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_dvec_t *feat_sums = tk_dvec_create(0, n_features * n_targets, 0, 0);
  tk_dvec_t *feat_sum_sq = tk_dvec_create(0, n_features * n_targets, 0, 0);
  tk_ivec_zero(feat_counts);
  tk_dvec_zero(feat_sums);
  tk_dvec_zero(feat_sum_sq);

  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit = set_bits->a[i];
    if (bit < 0) continue;
    uint64_t sample = (uint64_t)bit / n_features;
    uint64_t feature = (uint64_t)bit % n_features;
    if (sample >= n_samples || feature >= n_features) continue;
    feat_counts->a[feature]++;
    for (uint64_t t = 0; t < n_targets; t++) {
      double y = targets->a[sample * n_targets + t];
      feat_sums->a[feature * n_targets + t] += y;
      feat_sum_sq->a[feature * n_targets + t] += y * y;
    }
  }

  tk_rvec_t **per_target_heaps = (tk_rvec_t **)malloc(n_targets * sizeof(tk_rvec_t *));
  for (uint64_t t = 0; t < n_targets; t++)
    per_target_heaps[t] = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    int64_t n1 = feat_counts->a[f];
    int64_t n0 = (int64_t)n_samples - n1;
    if (n1 <= 1 || n0 <= 1) continue;
    for (uint64_t t = 0; t < n_targets; t++) {
      double sum1 = feat_sums->a[f * n_targets + t];
      double sum0 = overall_sums->a[t] - sum1;
      double mean1 = sum1 / (double)n1;
      double mean0 = sum0 / (double)n0;
      double overall_mean = overall_means->a[t];
      double ssb = (double)n1 * (mean1 - overall_mean) * (mean1 - overall_mean) +
                   (double)n0 * (mean0 - overall_mean) * (mean0 - overall_mean);
      double sum_sq1 = feat_sum_sq->a[f * n_targets + t];
      double sum_sq0 = total_sum_sqs->a[t] - sum_sq1;
      double ssw1 = sum_sq1 - (sum1 * sum1) / (double)n1;
      double ssw0 = sum_sq0 - (sum0 * sum0) / (double)n0;
      double ssw = ssw1 + ssw0;
      if (ssw < 1e-12) continue;
      double F = ssb / (ssw / (double)(n_samples - 2));
      tk_rank_t r = { (int64_t)f, F };
      #pragma omp critical
      tk_rvec_hmin(per_target_heaps[t], top_k, r);
    }
  }

  tk_ivec_destroy(feat_counts);
  tk_dvec_destroy(feat_sums);
  tk_dvec_destroy(feat_sum_sq);
  tk_dvec_destroy(overall_sums);
  tk_dvec_destroy(overall_means);
  tk_dvec_destroy(total_sum_sqs);

  for (uint64_t t = 0; t < n_targets; t++)
    tk_rvec_desc(per_target_heaps[t], 0, per_target_heaps[t]->n);

  tk_iuset_t *union_set = tk_iuset_create(0, 0);
  for (uint64_t t = 0; t < n_targets; t++) {
    for (uint64_t i = 0; i < per_target_heaps[t]->n; i++) {
      int absent;
      tk_iuset_put(union_set, per_target_heaps[t]->a[i].i, &absent);
    }
  }

  tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
  ids_union->n = 0;
  int64_t key;
  tk_umap_foreach_keys(union_set, key, ({
    ids_union->a[ids_union->n++] = key;
  }));
  tk_ivec_asc(ids_union, 0, ids_union->n);

  tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
  for (uint64_t i = 0; i < ids_union->n; i++) {
    int absent;
    khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
    kh_value(id_to_union_idx, k) = (int64_t)i;
  }

  uint64_t total_ids = 0;
  for (uint64_t t = 0; t < n_targets; t++)
    total_ids += per_target_heaps[t]->n;

  tk_ivec_t *offsets = tk_ivec_create(L, n_targets + 1, 0, 0);
  tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

  offsets->a[0] = 0;
  ids->n = 0;
  weights->n = 0;
  for (uint64_t t = 0; t < n_targets; t++) {
    for (uint64_t i = 0; i < per_target_heaps[t]->n; i++) {
      int64_t feat_id = per_target_heaps[t]->a[i].i;
      khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
      ids->a[ids->n++] = kh_value(id_to_union_idx, k);
      weights->a[weights->n++] = per_target_heaps[t]->a[i].d;
    }
    offsets->a[t + 1] = (int64_t)ids->n;
  }
  offsets->n = n_targets + 1;

  tk_iumap_destroy(id_to_union_idx);
  tk_iuset_destroy(union_set);
  for (uint64_t t = 0; t < n_targets; t++)
    tk_rvec_destroy(per_target_heaps[t]);
  free(per_target_heaps);
}

static inline void tk_cvec_bits_top_reg_f_ind (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_targets,
  uint64_t top_k
) {
  tk_dvec_t *overall_sums = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_t *overall_means = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_t *total_sum_sqs = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_zero(overall_sums);
  tk_dvec_zero(total_sum_sqs);

  for (uint64_t i = 0; i < n_samples; i++) {
    for (uint64_t t = 0; t < n_targets; t++) {
      double y = targets->a[i * n_targets + t];
      overall_sums->a[t] += y;
      total_sum_sqs->a[t] += y * y;
    }
  }
  for (uint64_t t = 0; t < n_targets; t++)
    overall_means->a[t] = overall_sums->a[t] / (double)n_samples;

  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_dvec_t *feat_sums = tk_dvec_create(0, n_features * n_targets, 0, 0);
  tk_dvec_t *feat_sum_sq = tk_dvec_create(0, n_features * n_targets, 0, 0);
  tk_ivec_zero(feat_counts);
  tk_dvec_zero(feat_sums);
  tk_dvec_zero(feat_sum_sq);

  uint8_t *bitmap_data = (uint8_t *)bitmap->a;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);

  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t sample_offset = s * bytes_per_sample;
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t byte_idx = f / CHAR_BIT;
      uint8_t bit_idx = f % CHAR_BIT;
      if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
        feat_counts->a[f]++;
        for (uint64_t t = 0; t < n_targets; t++) {
          double y = targets->a[s * n_targets + t];
          feat_sums->a[f * n_targets + t] += y;
          feat_sum_sq->a[f * n_targets + t] += y * y;
        }
      }
    }
  }

  tk_rvec_t **per_target_heaps = (tk_rvec_t **)malloc(n_targets * sizeof(tk_rvec_t *));
  for (uint64_t t = 0; t < n_targets; t++)
    per_target_heaps[t] = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    int64_t n1 = feat_counts->a[f];
    int64_t n0 = (int64_t)n_samples - n1;
    if (n1 <= 1 || n0 <= 1) continue;
    for (uint64_t t = 0; t < n_targets; t++) {
      double sum1 = feat_sums->a[f * n_targets + t];
      double sum0 = overall_sums->a[t] - sum1;
      double mean1 = sum1 / (double)n1;
      double mean0 = sum0 / (double)n0;
      double overall_mean = overall_means->a[t];
      double ssb = (double)n1 * (mean1 - overall_mean) * (mean1 - overall_mean) +
                   (double)n0 * (mean0 - overall_mean) * (mean0 - overall_mean);
      double sum_sq1 = feat_sum_sq->a[f * n_targets + t];
      double sum_sq0 = total_sum_sqs->a[t] - sum_sq1;
      double ssw1 = sum_sq1 - (sum1 * sum1) / (double)n1;
      double ssw0 = sum_sq0 - (sum0 * sum0) / (double)n0;
      double ssw = ssw1 + ssw0;
      if (ssw < 1e-12) continue;
      double F = ssb / (ssw / (double)(n_samples - 2));
      tk_rank_t r = { (int64_t)f, F };
      #pragma omp critical
      tk_rvec_hmin(per_target_heaps[t], top_k, r);
    }
  }

  tk_ivec_destroy(feat_counts);
  tk_dvec_destroy(feat_sums);
  tk_dvec_destroy(feat_sum_sq);
  tk_dvec_destroy(overall_sums);
  tk_dvec_destroy(overall_means);
  tk_dvec_destroy(total_sum_sqs);

  for (uint64_t t = 0; t < n_targets; t++)
    tk_rvec_desc(per_target_heaps[t], 0, per_target_heaps[t]->n);

  tk_iuset_t *union_set = tk_iuset_create(0, 0);
  for (uint64_t t = 0; t < n_targets; t++) {
    for (uint64_t i = 0; i < per_target_heaps[t]->n; i++) {
      int absent;
      tk_iuset_put(union_set, per_target_heaps[t]->a[i].i, &absent);
    }
  }

  tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
  ids_union->n = 0;
  int64_t key;
  tk_umap_foreach_keys(union_set, key, ({
    ids_union->a[ids_union->n++] = key;
  }));
  tk_ivec_asc(ids_union, 0, ids_union->n);

  tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
  for (uint64_t i = 0; i < ids_union->n; i++) {
    int absent;
    khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
    kh_value(id_to_union_idx, k) = (int64_t)i;
  }

  uint64_t total_ids = 0;
  for (uint64_t t = 0; t < n_targets; t++)
    total_ids += per_target_heaps[t]->n;

  tk_ivec_t *offsets = tk_ivec_create(L, n_targets + 1, 0, 0);
  tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

  offsets->a[0] = 0;
  ids->n = 0;
  weights->n = 0;
  for (uint64_t t = 0; t < n_targets; t++) {
    for (uint64_t i = 0; i < per_target_heaps[t]->n; i++) {
      int64_t feat_id = per_target_heaps[t]->a[i].i;
      khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
      ids->a[ids->n++] = kh_value(id_to_union_idx, k);
      weights->a[weights->n++] = per_target_heaps[t]->a[i].d;
    }
    offsets->a[t + 1] = (int64_t)ids->n;
  }
  offsets->n = n_targets + 1;

  tk_iumap_destroy(id_to_union_idx);
  tk_iuset_destroy(union_set);
  for (uint64_t t = 0; t < n_targets; t++)
    tk_rvec_destroy(per_target_heaps[t]);
  free(per_target_heaps);
}

static inline void tk_ivec_bits_top_reg_pearson_ind (
  lua_State *L,
  tk_ivec_t *set_bits,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_targets,
  uint64_t top_k
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);

  tk_dvec_t *overall_sums = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_t *overall_sum_sqs = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_t *overall_means = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_t *overall_stds = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_zero(overall_sums);
  tk_dvec_zero(overall_sum_sqs);

  for (uint64_t i = 0; i < n_samples; i++) {
    for (uint64_t t = 0; t < n_targets; t++) {
      double y = targets->a[i * n_targets + t];
      overall_sums->a[t] += y;
      overall_sum_sqs->a[t] += y * y;
    }
  }
  bool any_valid = false;
  for (uint64_t t = 0; t < n_targets; t++) {
    overall_means->a[t] = overall_sums->a[t] / (double)n_samples;
    double var = (overall_sum_sqs->a[t] / (double)n_samples) - (overall_means->a[t] * overall_means->a[t]);
    overall_stds->a[t] = sqrt(var);
    if (overall_stds->a[t] >= 1e-12) any_valid = true;
  }
  if (!any_valid) {
    tk_dvec_destroy(overall_sums);
    tk_dvec_destroy(overall_sum_sqs);
    tk_dvec_destroy(overall_means);
    tk_dvec_destroy(overall_stds);
    tk_ivec_create(L, 0, 0, 0);
    tk_ivec_create(L, n_targets + 1, 0, 0);
    tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return;
  }

  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_dvec_t *feat_sums = tk_dvec_create(0, n_features * n_targets, 0, 0);
  tk_ivec_zero(feat_counts);
  tk_dvec_zero(feat_sums);

  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit = set_bits->a[i];
    if (bit < 0) continue;
    uint64_t sample = (uint64_t)bit / n_features;
    uint64_t feature = (uint64_t)bit % n_features;
    if (sample >= n_samples || feature >= n_features) continue;
    feat_counts->a[feature]++;
    for (uint64_t t = 0; t < n_targets; t++)
      feat_sums->a[feature * n_targets + t] += targets->a[sample * n_targets + t];
  }

  tk_rvec_t **per_target_heaps = (tk_rvec_t **)malloc(n_targets * sizeof(tk_rvec_t *));
  for (uint64_t t = 0; t < n_targets; t++)
    per_target_heaps[t] = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    int64_t n1 = feat_counts->a[f];
    int64_t n0 = (int64_t)n_samples - n1;
    if (n1 == 0 || n0 == 0) continue;
    double p = (double)n1 / (double)n_samples;
    double feat_std = sqrt(p * (1.0 - p));
    if (feat_std < 1e-12) continue;
    for (uint64_t t = 0; t < n_targets; t++) {
      if (overall_stds->a[t] < 1e-12) continue;
      double mean1 = feat_sums->a[f * n_targets + t] / (double)n1;
      double r_pb = (mean1 - overall_means->a[t]) * sqrt(p * (1.0 - p)) / overall_stds->a[t];
      double score = fabs(r_pb);
      tk_rank_t r = { (int64_t)f, score };
      #pragma omp critical
      tk_rvec_hmin(per_target_heaps[t], top_k, r);
    }
  }

  tk_ivec_destroy(feat_counts);
  tk_dvec_destroy(feat_sums);
  tk_dvec_destroy(overall_sums);
  tk_dvec_destroy(overall_sum_sqs);
  tk_dvec_destroy(overall_means);
  tk_dvec_destroy(overall_stds);

  for (uint64_t t = 0; t < n_targets; t++)
    tk_rvec_desc(per_target_heaps[t], 0, per_target_heaps[t]->n);

  tk_iuset_t *union_set = tk_iuset_create(0, 0);
  for (uint64_t t = 0; t < n_targets; t++) {
    for (uint64_t i = 0; i < per_target_heaps[t]->n; i++) {
      int absent;
      tk_iuset_put(union_set, per_target_heaps[t]->a[i].i, &absent);
    }
  }

  tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
  ids_union->n = 0;
  int64_t key;
  tk_umap_foreach_keys(union_set, key, ({
    ids_union->a[ids_union->n++] = key;
  }));
  tk_ivec_asc(ids_union, 0, ids_union->n);

  tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
  for (uint64_t i = 0; i < ids_union->n; i++) {
    int absent;
    khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
    kh_value(id_to_union_idx, k) = (int64_t)i;
  }

  uint64_t total_ids = 0;
  for (uint64_t t = 0; t < n_targets; t++)
    total_ids += per_target_heaps[t]->n;

  tk_ivec_t *offsets = tk_ivec_create(L, n_targets + 1, 0, 0);
  tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

  offsets->a[0] = 0;
  ids->n = 0;
  weights->n = 0;
  for (uint64_t t = 0; t < n_targets; t++) {
    for (uint64_t i = 0; i < per_target_heaps[t]->n; i++) {
      int64_t feat_id = per_target_heaps[t]->a[i].i;
      khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
      ids->a[ids->n++] = kh_value(id_to_union_idx, k);
      weights->a[weights->n++] = per_target_heaps[t]->a[i].d;
    }
    offsets->a[t + 1] = (int64_t)ids->n;
  }
  offsets->n = n_targets + 1;

  tk_iumap_destroy(id_to_union_idx);
  tk_iuset_destroy(union_set);
  for (uint64_t t = 0; t < n_targets; t++)
    tk_rvec_destroy(per_target_heaps[t]);
  free(per_target_heaps);
}

static inline void tk_cvec_bits_top_reg_pearson_ind (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_targets,
  uint64_t top_k
) {
  tk_dvec_t *overall_sums = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_t *overall_sum_sqs = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_t *overall_means = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_t *overall_stds = tk_dvec_create(0, n_targets, 0, 0);
  tk_dvec_zero(overall_sums);
  tk_dvec_zero(overall_sum_sqs);

  for (uint64_t i = 0; i < n_samples; i++) {
    for (uint64_t t = 0; t < n_targets; t++) {
      double y = targets->a[i * n_targets + t];
      overall_sums->a[t] += y;
      overall_sum_sqs->a[t] += y * y;
    }
  }
  bool any_valid = false;
  for (uint64_t t = 0; t < n_targets; t++) {
    overall_means->a[t] = overall_sums->a[t] / (double)n_samples;
    double var = (overall_sum_sqs->a[t] / (double)n_samples) - (overall_means->a[t] * overall_means->a[t]);
    overall_stds->a[t] = sqrt(var);
    if (overall_stds->a[t] >= 1e-12) any_valid = true;
  }
  if (!any_valid) {
    tk_dvec_destroy(overall_sums);
    tk_dvec_destroy(overall_sum_sqs);
    tk_dvec_destroy(overall_means);
    tk_dvec_destroy(overall_stds);
    tk_ivec_create(L, 0, 0, 0);
    tk_ivec_create(L, n_targets + 1, 0, 0);
    tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
    return;
  }

  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_dvec_t *feat_sums = tk_dvec_create(0, n_features * n_targets, 0, 0);
  tk_ivec_zero(feat_counts);
  tk_dvec_zero(feat_sums);

  uint8_t *bitmap_data = (uint8_t *)bitmap->a;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);

  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t sample_offset = s * bytes_per_sample;
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t byte_idx = f / CHAR_BIT;
      uint8_t bit_idx = f % CHAR_BIT;
      if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
        feat_counts->a[f]++;
        for (uint64_t t = 0; t < n_targets; t++)
          feat_sums->a[f * n_targets + t] += targets->a[s * n_targets + t];
      }
    }
  }

  tk_rvec_t **per_target_heaps = (tk_rvec_t **)malloc(n_targets * sizeof(tk_rvec_t *));
  for (uint64_t t = 0; t < n_targets; t++)
    per_target_heaps[t] = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    int64_t n1 = feat_counts->a[f];
    int64_t n0 = (int64_t)n_samples - n1;
    if (n1 == 0 || n0 == 0) continue;
    double p = (double)n1 / (double)n_samples;
    double feat_std = sqrt(p * (1.0 - p));
    if (feat_std < 1e-12) continue;
    for (uint64_t t = 0; t < n_targets; t++) {
      if (overall_stds->a[t] < 1e-12) continue;
      double mean1 = feat_sums->a[f * n_targets + t] / (double)n1;
      double r_pb = (mean1 - overall_means->a[t]) * sqrt(p * (1.0 - p)) / overall_stds->a[t];
      double score = fabs(r_pb);
      tk_rank_t r = { (int64_t)f, score };
      #pragma omp critical
      tk_rvec_hmin(per_target_heaps[t], top_k, r);
    }
  }

  tk_ivec_destroy(feat_counts);
  tk_dvec_destroy(feat_sums);
  tk_dvec_destroy(overall_sums);
  tk_dvec_destroy(overall_sum_sqs);
  tk_dvec_destroy(overall_means);
  tk_dvec_destroy(overall_stds);

  for (uint64_t t = 0; t < n_targets; t++)
    tk_rvec_desc(per_target_heaps[t], 0, per_target_heaps[t]->n);

  tk_iuset_t *union_set = tk_iuset_create(0, 0);
  for (uint64_t t = 0; t < n_targets; t++) {
    for (uint64_t i = 0; i < per_target_heaps[t]->n; i++) {
      int absent;
      tk_iuset_put(union_set, per_target_heaps[t]->a[i].i, &absent);
    }
  }

  tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
  ids_union->n = 0;
  int64_t key;
  tk_umap_foreach_keys(union_set, key, ({
    ids_union->a[ids_union->n++] = key;
  }));
  tk_ivec_asc(ids_union, 0, ids_union->n);

  tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
  for (uint64_t i = 0; i < ids_union->n; i++) {
    int absent;
    khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
    kh_value(id_to_union_idx, k) = (int64_t)i;
  }

  uint64_t total_ids = 0;
  for (uint64_t t = 0; t < n_targets; t++)
    total_ids += per_target_heaps[t]->n;

  tk_ivec_t *offsets = tk_ivec_create(L, n_targets + 1, 0, 0);
  tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

  offsets->a[0] = 0;
  ids->n = 0;
  weights->n = 0;
  for (uint64_t t = 0; t < n_targets; t++) {
    for (uint64_t i = 0; i < per_target_heaps[t]->n; i++) {
      int64_t feat_id = per_target_heaps[t]->a[i].i;
      khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
      ids->a[ids->n++] = kh_value(id_to_union_idx, k);
      weights->a[weights->n++] = per_target_heaps[t]->a[i].d;
    }
    offsets->a[t + 1] = (int64_t)ids->n;
  }
  offsets->n = n_targets + 1;

  tk_iumap_destroy(id_to_union_idx);
  tk_iuset_destroy(union_set);
  for (uint64_t t = 0; t < n_targets; t++)
    tk_rvec_destroy(per_target_heaps[t]);
  free(per_target_heaps);
}

static inline void tk_ivec_bits_top_reg_mi_ind (
  lua_State *L,
  tk_ivec_t *set_bits,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_targets,
  uint64_t top_k,
  uint64_t n_bins
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);
  if (n_bins == 0) n_bins = 10;

  tk_ivec_t *bin_assignments = tk_ivec_create(0, n_samples * n_targets, 0, 0);
  tk_ivec_t *bin_counts = tk_ivec_create(0, n_targets * n_bins, 0, 0);
  tk_ivec_zero(bin_counts);

  for (uint64_t t = 0; t < n_targets; t++) {
    double y_min = targets->a[t], y_max = targets->a[t];
    for (uint64_t i = 1; i < n_samples; i++) {
      double y = targets->a[i * n_targets + t];
      if (y < y_min) y_min = y;
      if (y > y_max) y_max = y;
    }
    double y_range = y_max - y_min;
    for (uint64_t i = 0; i < n_samples; i++) {
      uint64_t bin;
      if (y_range < 1e-12) {
        bin = 0;
      } else {
        double normalized = (targets->a[i * n_targets + t] - y_min) / y_range;
        bin = (uint64_t)(normalized * (double)(n_bins - 1) + 0.5);
        if (bin >= n_bins) bin = n_bins - 1;
      }
      bin_assignments->a[i * n_targets + t] = (int64_t)bin;
      bin_counts->a[t * n_bins + bin]++;
    }
  }

  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_ivec_t *joint_counts = tk_ivec_create(0, n_features * n_targets * n_bins, 0, 0);
  tk_ivec_zero(feat_counts);
  tk_ivec_zero(joint_counts);

  for (uint64_t i = 0; i < set_bits->n; i++) {
    int64_t bit = set_bits->a[i];
    if (bit < 0) continue;
    uint64_t sample = (uint64_t)bit / n_features;
    uint64_t feature = (uint64_t)bit % n_features;
    if (sample >= n_samples || feature >= n_features) continue;
    feat_counts->a[feature]++;
    for (uint64_t t = 0; t < n_targets; t++) {
      uint64_t bin = (uint64_t)bin_assignments->a[sample * n_targets + t];
      joint_counts->a[feature * n_targets * n_bins + t * n_bins + bin]++;
    }
  }

  tk_rvec_t **per_target_heaps = (tk_rvec_t **)malloc(n_targets * sizeof(tk_rvec_t *));
  for (uint64_t t = 0; t < n_targets; t++)
    per_target_heaps[t] = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    int64_t n1 = feat_counts->a[f];
    int64_t n0 = (int64_t)n_samples - n1;
    if (n1 == 0 || n0 == 0) continue;
    double p1 = (double)n1 / (double)n_samples;
    double p0 = (double)n0 / (double)n_samples;
    for (uint64_t t = 0; t < n_targets; t++) {
      double mi = 0.0;
      for (uint64_t b = 0; b < n_bins; b++) {
        int64_t n_b = bin_counts->a[t * n_bins + b];
        if (n_b == 0) continue;
        double p_b = (double)n_b / (double)n_samples;
        int64_t n_1b = joint_counts->a[f * n_targets * n_bins + t * n_bins + b];
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
      #pragma omp critical
      tk_rvec_hmin(per_target_heaps[t], top_k, r);
    }
  }

  tk_ivec_destroy(bin_assignments);
  tk_ivec_destroy(bin_counts);
  tk_ivec_destroy(feat_counts);
  tk_ivec_destroy(joint_counts);

  for (uint64_t t = 0; t < n_targets; t++)
    tk_rvec_desc(per_target_heaps[t], 0, per_target_heaps[t]->n);

  tk_iuset_t *union_set = tk_iuset_create(0, 0);
  for (uint64_t t = 0; t < n_targets; t++) {
    for (uint64_t i = 0; i < per_target_heaps[t]->n; i++) {
      int absent;
      tk_iuset_put(union_set, per_target_heaps[t]->a[i].i, &absent);
    }
  }

  tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
  ids_union->n = 0;
  int64_t key;
  tk_umap_foreach_keys(union_set, key, ({
    ids_union->a[ids_union->n++] = key;
  }));
  tk_ivec_asc(ids_union, 0, ids_union->n);

  tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
  for (uint64_t i = 0; i < ids_union->n; i++) {
    int absent;
    khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
    kh_value(id_to_union_idx, k) = (int64_t)i;
  }

  uint64_t total_ids = 0;
  for (uint64_t t = 0; t < n_targets; t++)
    total_ids += per_target_heaps[t]->n;

  tk_ivec_t *offsets = tk_ivec_create(L, n_targets + 1, 0, 0);
  tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

  offsets->a[0] = 0;
  ids->n = 0;
  weights->n = 0;
  for (uint64_t t = 0; t < n_targets; t++) {
    for (uint64_t i = 0; i < per_target_heaps[t]->n; i++) {
      int64_t feat_id = per_target_heaps[t]->a[i].i;
      khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
      ids->a[ids->n++] = kh_value(id_to_union_idx, k);
      weights->a[weights->n++] = per_target_heaps[t]->a[i].d;
    }
    offsets->a[t + 1] = (int64_t)ids->n;
  }
  offsets->n = n_targets + 1;

  tk_iumap_destroy(id_to_union_idx);
  tk_iuset_destroy(union_set);
  for (uint64_t t = 0; t < n_targets; t++)
    tk_rvec_destroy(per_target_heaps[t]);
  free(per_target_heaps);
}

static inline void tk_cvec_bits_top_reg_mi_ind (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_dvec_t *targets,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_targets,
  uint64_t top_k,
  uint64_t n_bins
) {
  if (n_bins == 0) n_bins = 10;

  tk_ivec_t *bin_assignments = tk_ivec_create(0, n_samples * n_targets, 0, 0);
  tk_ivec_t *bin_counts = tk_ivec_create(0, n_targets * n_bins, 0, 0);
  tk_ivec_zero(bin_counts);

  for (uint64_t t = 0; t < n_targets; t++) {
    double y_min = targets->a[t], y_max = targets->a[t];
    for (uint64_t i = 1; i < n_samples; i++) {
      double y = targets->a[i * n_targets + t];
      if (y < y_min) y_min = y;
      if (y > y_max) y_max = y;
    }
    double y_range = y_max - y_min;
    for (uint64_t i = 0; i < n_samples; i++) {
      uint64_t bin;
      if (y_range < 1e-12) {
        bin = 0;
      } else {
        double normalized = (targets->a[i * n_targets + t] - y_min) / y_range;
        bin = (uint64_t)(normalized * (double)(n_bins - 1) + 0.5);
        if (bin >= n_bins) bin = n_bins - 1;
      }
      bin_assignments->a[i * n_targets + t] = (int64_t)bin;
      bin_counts->a[t * n_bins + bin]++;
    }
  }

  tk_ivec_t *feat_counts = tk_ivec_create(0, n_features, 0, 0);
  tk_ivec_t *joint_counts = tk_ivec_create(0, n_features * n_targets * n_bins, 0, 0);
  tk_ivec_zero(feat_counts);
  tk_ivec_zero(joint_counts);

  uint8_t *bitmap_data = (uint8_t *)bitmap->a;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);

  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t sample_offset = s * bytes_per_sample;
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t byte_idx = f / CHAR_BIT;
      uint8_t bit_idx = f % CHAR_BIT;
      if (bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)) {
        feat_counts->a[f]++;
        for (uint64_t t = 0; t < n_targets; t++) {
          uint64_t bin = (uint64_t)bin_assignments->a[s * n_targets + t];
          joint_counts->a[f * n_targets * n_bins + t * n_bins + bin]++;
        }
      }
    }
  }

  tk_rvec_t **per_target_heaps = (tk_rvec_t **)malloc(n_targets * sizeof(tk_rvec_t *));
  for (uint64_t t = 0; t < n_targets; t++)
    per_target_heaps[t] = tk_rvec_create(0, 0, 0, 0);

  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    int64_t n1 = feat_counts->a[f];
    int64_t n0 = (int64_t)n_samples - n1;
    if (n1 == 0 || n0 == 0) continue;
    double p1 = (double)n1 / (double)n_samples;
    double p0 = (double)n0 / (double)n_samples;
    for (uint64_t t = 0; t < n_targets; t++) {
      double mi = 0.0;
      for (uint64_t b = 0; b < n_bins; b++) {
        int64_t n_b = bin_counts->a[t * n_bins + b];
        if (n_b == 0) continue;
        double p_b = (double)n_b / (double)n_samples;
        int64_t n_1b = joint_counts->a[f * n_targets * n_bins + t * n_bins + b];
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
      #pragma omp critical
      tk_rvec_hmin(per_target_heaps[t], top_k, r);
    }
  }

  tk_ivec_destroy(bin_assignments);
  tk_ivec_destroy(bin_counts);
  tk_ivec_destroy(feat_counts);
  tk_ivec_destroy(joint_counts);

  for (uint64_t t = 0; t < n_targets; t++)
    tk_rvec_desc(per_target_heaps[t], 0, per_target_heaps[t]->n);

  tk_iuset_t *union_set = tk_iuset_create(0, 0);
  for (uint64_t t = 0; t < n_targets; t++) {
    for (uint64_t i = 0; i < per_target_heaps[t]->n; i++) {
      int absent;
      tk_iuset_put(union_set, per_target_heaps[t]->a[i].i, &absent);
    }
  }

  tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
  ids_union->n = 0;
  int64_t key;
  tk_umap_foreach_keys(union_set, key, ({
    ids_union->a[ids_union->n++] = key;
  }));
  tk_ivec_asc(ids_union, 0, ids_union->n);

  tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
  for (uint64_t i = 0; i < ids_union->n; i++) {
    int absent;
    khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
    kh_value(id_to_union_idx, k) = (int64_t)i;
  }

  uint64_t total_ids = 0;
  for (uint64_t t = 0; t < n_targets; t++)
    total_ids += per_target_heaps[t]->n;

  tk_ivec_t *offsets = tk_ivec_create(L, n_targets + 1, 0, 0);
  tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
  tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

  offsets->a[0] = 0;
  ids->n = 0;
  weights->n = 0;
  for (uint64_t t = 0; t < n_targets; t++) {
    for (uint64_t i = 0; i < per_target_heaps[t]->n; i++) {
      int64_t feat_id = per_target_heaps[t]->a[i].i;
      khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
      ids->a[ids->n++] = kh_value(id_to_union_idx, k);
      weights->a[weights->n++] = per_target_heaps[t]->a[i].d;
    }
    offsets->a[t + 1] = (int64_t)ids->n;
  }
  offsets->n = n_targets + 1;

  tk_iumap_destroy(id_to_union_idx);
  tk_iuset_destroy(union_set);
  for (uint64_t t = 0; t < n_targets; t++)
    tk_rvec_destroy(per_target_heaps[t]);
  free(per_target_heaps);
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

static inline tk_ivec_t *tk_ivec_bits_top_bns (
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

    double N = (double)n_samples;
    double eps = 0.5 / N;
    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_visible; f++) {
      double C = (double)feat_counts->a[f];
      if (C == 0 || C == N) continue;
      double sum_bns = 0.0;
      for (uint64_t b = 0; b < n_hidden; b++) {
        double P = (double)label_counts->a[b];
        if (P == 0 || P == N) continue;
        double a = (double)active_counts->a[f * n_hidden + b];
        double tpr_raw = a / P;
        double fpr_raw = (C - a) / (N - P);
        double tpr = tpr_raw < eps ? eps : (tpr_raw > 1.0 - eps ? 1.0 - eps : tpr_raw);
        double fpr = fpr_raw < eps ? eps : (fpr_raw > 1.0 - eps ? 1.0 - eps : fpr_raw);
        double bns = fabs(tk_probit(tpr) - tk_probit(fpr));
        sum_bns += bns;
      }
      if (sum_bns > 0) {
        tk_rank_t r = { (int64_t)f, sum_bns };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
      }
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
      if (set_bits->a[si] < 0) { si++; continue; }
      uint64_t s_sample = (uint64_t)set_bits->a[si] / n_visible;
      uint64_t f = (uint64_t)set_bits->a[si] % n_visible;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden < s_sample)
        li++;
      if (li >= labels->n || labels->a[li] < 0 || (uint64_t)labels->a[li] / n_hidden > s_sample) {
        si++; continue;
      }
      size_t li_start = li;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden == s_sample) {
        uint64_t h = (uint64_t)labels->a[li] % n_hidden;
        int64_t key = (int64_t)(f * n_hidden + h);
        tk_iumap_inc(active_counts, key);
        li++;
      }
      si++;
      if (si < set_bits->n && set_bits->a[si] >= 0 && (uint64_t)set_bits->a[si] / n_visible == s_sample)
        li = li_start;
    }

    double N = (double)n_samples;
    double eps = 0.5 / N;
    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    tk_dvec_t *feat_sum_bns = tk_dvec_create(0, n_visible, 0, 0);
    tk_dvec_zero(feat_sum_bns);

    int64_t k, v;
    tk_umap_foreach(active_counts, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t b = (uint64_t)k % n_hidden;
      if (f >= n_visible || b >= n_hidden) continue;
      double C = (double)feat_counts->a[f];
      double P = (double)label_counts->a[b];
      if (C == 0 || P == 0 || C == N || P == N) continue;
      double a = (double)v;
      double tpr_raw = a / P;
      double fpr_raw = (C - a) / (N - P);
      double tpr = tpr_raw < eps ? eps : (tpr_raw > 1.0 - eps ? 1.0 - eps : tpr_raw);
      double fpr = fpr_raw < eps ? eps : (fpr_raw > 1.0 - eps ? 1.0 - eps : fpr_raw);
      double bns = fabs(tk_probit(tpr) - tk_probit(fpr));
      feat_sum_bns->a[f] += bns;
    }));

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_visible; f++) {
      if (feat_sum_bns->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_sum_bns->a[f] };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_sum_bns);

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

static inline void tk_ivec_bits_top_bns_ind (
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
      if (bit_idx < 0) continue;
      uint64_t sample_idx = (uint64_t)bit_idx / n_visible;
      uint64_t feature_idx = (uint64_t)bit_idx % n_visible;
      if (sample_idx >= n_samples || feature_idx >= n_visible) continue;
      feat_counts->a[feature_idx]++;
      if (sample_idx != prev_sample) {
        prev_sample = sample_idx;
        sample_codes = (uint8_t *)(codes + sample_idx * TK_CVEC_BITS_BYTES(n_hidden));
      }
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint64_t bit_pos = b % CHAR_BIT;
        if (sample_codes[byte_idx] & (1u << bit_pos))
          active_counts->a[feature_idx * n_hidden + b]++;
      }
    }

    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint64_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos))
          label_counts->a[b]++;
      }
    }

    double N = (double)n_samples;
    double eps = 0.5 / N;
    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_visible; f++) {
      double C = (double)feat_counts->a[f];
      if (C == 0 || C == N) continue;
      for (uint64_t b = 0; b < n_hidden; b++) {
        double P = (double)label_counts->a[b];
        if (P == 0 || P == N) continue;
        double a = (double)active_counts->a[f * n_hidden + b];
        double tpr_raw = a / P;
        double fpr_raw = (C - a) / (N - P);
        double tpr = tpr_raw < eps ? eps : (tpr_raw > 1.0 - eps ? 1.0 - eps : tpr_raw);
        double fpr = fpr_raw < eps ? eps : (fpr_raw > 1.0 - eps ? 1.0 - eps : fpr_raw);
        double bns = fabs(tk_probit(tpr) - tk_probit(fpr));
        if (bns > 0) {
          tk_rank_t r = { (int64_t)f, bns };
          #pragma omp critical
          tk_rvec_hmin(per_dim_heaps[b], top_k, r);
        }
      }
    }

    tk_ivec_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    int absent;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++)
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({ ids_union->a[ids_union->n++] = key; }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

  } else if (labels) {

    tk_ivec_asc(labels, 0, labels->n);
    tk_iumap_t *active_counts = tk_iumap_create(0, 0);
    tk_ivec_t *feat_counts = tk_ivec_create(0, n_visible, 0, 0);
    tk_ivec_t *label_counts = tk_ivec_create(0, n_hidden, 0, 0);
    tk_ivec_zero(feat_counts);
    tk_ivec_zero(label_counts);

    for (uint64_t i = 0; i < set_bits->n; i++) {
      int64_t bit = set_bits->a[i];
      if (bit >= 0) feat_counts->a[(uint64_t)bit % n_visible]++;
    }

    for (uint64_t i = 0; i < labels->n; i++) {
      int64_t bit = labels->a[i];
      if (bit >= 0) label_counts->a[(uint64_t)bit % n_hidden]++;
    }

    size_t si = 0, li = 0;
    while (si < set_bits->n) {
      if (set_bits->a[si] < 0) { si++; continue; }
      uint64_t s_sample = (uint64_t)set_bits->a[si] / n_visible;
      uint64_t f = (uint64_t)set_bits->a[si] % n_visible;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden < s_sample)
        li++;
      if (li >= labels->n || labels->a[li] < 0 || (uint64_t)labels->a[li] / n_hidden > s_sample) {
        si++; continue;
      }
      size_t li_start = li;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden == s_sample) {
        uint64_t h = (uint64_t)labels->a[li] % n_hidden;
        tk_iumap_inc(active_counts, (int64_t)(f * n_hidden + h));
        li++;
      }
      si++;
      if (si < set_bits->n && set_bits->a[si] >= 0 && (uint64_t)set_bits->a[si] / n_visible == s_sample)
        li = li_start;
    }

    double N = (double)n_samples;
    double eps = 0.5 / N;
    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    int64_t k, v;
    tk_umap_foreach(active_counts, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t b = (uint64_t)k % n_hidden;
      if (f >= n_visible || b >= n_hidden) continue;
      double C = (double)feat_counts->a[f];
      double P = (double)label_counts->a[b];
      if (C == 0 || P == 0 || C == N || P == N) continue;
      double a = (double)v;
      double tpr_raw = a / P;
      double fpr_raw = (C - a) / (N - P);
      double tpr = tpr_raw < eps ? eps : (tpr_raw > 1.0 - eps ? 1.0 - eps : tpr_raw);
      double fpr = fpr_raw < eps ? eps : (fpr_raw > 1.0 - eps ? 1.0 - eps : fpr_raw);
      double bns = fabs(tk_probit(tpr) - tk_probit(fpr));
      if (bns > 0) {
        tk_rank_t r = { (int64_t)f, bns };
        tk_rvec_hmin(per_dim_heaps[b], top_k, r);
      }
    }));

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    int absent;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++)
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({ ids_union->a[ids_union->n++] = key; }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

  } else {
    tk_ivec_create(L, 0, 0, 0);
    tk_ivec_create(L, 1, 0, 0);
    tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
  }
}

static inline tk_ivec_t *tk_cvec_bits_top_bns (
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
            if (sample_codes[b_byte] & (1u << b_bit))
              active_counts->a[f * n_hidden + b]++;
          }
        }
      }
    }

    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint8_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos))
          label_counts->a[b]++;
      }
    }

    double N = (double)n_samples;
    double eps = 0.5 / N;
    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_features; f++) {
      double C = (double)feat_counts->a[f];
      if (C == 0 || C == N) continue;
      double sum_bns = 0.0;
      for (uint64_t b = 0; b < n_hidden; b++) {
        double P = (double)label_counts->a[b];
        if (P == 0 || P == N) continue;
        double a = (double)active_counts->a[f * n_hidden + b];
        double tpr_raw = a / P;
        double fpr_raw = (C - a) / (N - P);
        double tpr = tpr_raw < eps ? eps : (tpr_raw > 1.0 - eps ? 1.0 - eps : tpr_raw);
        double fpr = fpr_raw < eps ? eps : (fpr_raw > 1.0 - eps ? 1.0 - eps : fpr_raw);
        double bns = fabs(tk_probit(tpr) - tk_probit(fpr));
        sum_bns += bns;
      }
      if (sum_bns > 0) {
        tk_rank_t r = { (int64_t)f, sum_bns };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
      }
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
      if (bit >= 0)
        label_counts->a[(uint64_t)bit % n_hidden]++;
    }

    size_t li = 0;
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden < s)
        li++;
      if (li >= labels->n || labels->a[li] < 0 || (uint64_t)labels->a[li] / n_hidden != s)
        continue;
      size_t li_start = li;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (!(bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)))
          continue;
        li = li_start;
        while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden == s) {
          uint64_t h = (uint64_t)labels->a[li] % n_hidden;
          tk_iumap_inc(active_counts, (int64_t)(f * n_hidden + h));
          li++;
        }
      }
    }

    double N = (double)n_samples;
    double eps = 0.5 / N;
    tk_rvec_t *top_heap = tk_rvec_create(0, 0, 0, 0);
    tk_dvec_t *feat_sum_bns = tk_dvec_create(0, n_features, 0, 0);
    tk_dvec_zero(feat_sum_bns);

    int64_t k, v;
    tk_umap_foreach(active_counts, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t b = (uint64_t)k % n_hidden;
      if (f >= n_features || b >= n_hidden) continue;
      double C = (double)feat_counts->a[f];
      double P = (double)label_counts->a[b];
      if (C == 0 || P == 0 || C == N || P == N) continue;
      double a = (double)v;
      double tpr_raw = a / P;
      double fpr_raw = (C - a) / (N - P);
      double tpr = tpr_raw < eps ? eps : (tpr_raw > 1.0 - eps ? 1.0 - eps : tpr_raw);
      double fpr = fpr_raw < eps ? eps : (fpr_raw > 1.0 - eps ? 1.0 - eps : fpr_raw);
      double bns = fabs(tk_probit(tpr) - tk_probit(fpr));
      feat_sum_bns->a[f] += bns;
    }));

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_features; f++) {
      if (feat_sum_bns->a[f] > 0) {
        tk_rank_t r = { (int64_t)f, feat_sum_bns->a[f] };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
      }
    }

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);
    tk_dvec_destroy(feat_sum_bns);

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

static inline void tk_cvec_bits_top_bns_ind (
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
            if (sample_codes[b_byte] & (1u << b_bit))
              active_counts->a[f * n_hidden + b]++;
          }
        }
      }
    }

    for (uint64_t s = 0; s < n_samples; s++) {
      uint8_t *s_codes = (uint8_t *)(codes->a + s * TK_CVEC_BITS_BYTES(n_hidden));
      for (uint64_t b = 0; b < n_hidden; b++) {
        uint64_t byte_idx = b / CHAR_BIT;
        uint8_t bit_pos = b % CHAR_BIT;
        if (s_codes[byte_idx] & (1u << bit_pos))
          label_counts->a[b]++;
      }
    }

    double N = (double)n_samples;
    double eps = 0.5 / N;
    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    #pragma omp parallel for
    for (uint64_t f = 0; f < n_features; f++) {
      double C = (double)feat_counts->a[f];
      if (C == 0 || C == N) continue;
      for (uint64_t b = 0; b < n_hidden; b++) {
        double P = (double)label_counts->a[b];
        if (P == 0 || P == N) continue;
        double a = (double)active_counts->a[f * n_hidden + b];
        double tpr_raw = a / P;
        double fpr_raw = (C - a) / (N - P);
        double tpr = tpr_raw < eps ? eps : (tpr_raw > 1.0 - eps ? 1.0 - eps : tpr_raw);
        double fpr = fpr_raw < eps ? eps : (fpr_raw > 1.0 - eps ? 1.0 - eps : fpr_raw);
        double bns = fabs(tk_probit(tpr) - tk_probit(fpr));
        if (bns > 0) {
          tk_rank_t r = { (int64_t)f, bns };
          #pragma omp critical
          tk_rvec_hmin(per_dim_heaps[b], top_k, r);
        }
      }
    }

    tk_ivec_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    int absent;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++)
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({ ids_union->a[ids_union->n++] = key; }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

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
      if (bit >= 0)
        label_counts->a[(uint64_t)bit % n_hidden]++;
    }

    size_t li = 0;
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden < s)
        li++;
      if (li >= labels->n || labels->a[li] < 0 || (uint64_t)labels->a[li] / n_hidden != s)
        continue;
      size_t li_start = li;
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t byte_idx = f / CHAR_BIT;
        uint8_t bit_idx = f % CHAR_BIT;
        if (!(bitmap_data[sample_offset + byte_idx] & (1u << bit_idx)))
          continue;
        li = li_start;
        while (li < labels->n && labels->a[li] >= 0 && (uint64_t)labels->a[li] / n_hidden == s) {
          uint64_t h = (uint64_t)labels->a[li] % n_hidden;
          tk_iumap_inc(active_counts, (int64_t)(f * n_hidden + h));
          li++;
        }
      }
    }

    double N = (double)n_samples;
    double eps = 0.5 / N;
    tk_rvec_t **per_dim_heaps = (tk_rvec_t **)malloc(n_hidden * sizeof(tk_rvec_t *));
    for (uint64_t h = 0; h < n_hidden; h++)
      per_dim_heaps[h] = tk_rvec_create(0, 0, 0, 0);

    int64_t k, v;
    tk_umap_foreach(active_counts, k, v, ({
      uint64_t f = (uint64_t)k / n_hidden;
      uint64_t b = (uint64_t)k % n_hidden;
      if (f >= n_features || b >= n_hidden) continue;
      double C = (double)feat_counts->a[f];
      double P = (double)label_counts->a[b];
      if (C == 0 || P == 0 || C == N || P == N) continue;
      double a = (double)v;
      double tpr_raw = a / P;
      double fpr_raw = (C - a) / (N - P);
      double tpr = tpr_raw < eps ? eps : (tpr_raw > 1.0 - eps ? 1.0 - eps : tpr_raw);
      double fpr = fpr_raw < eps ? eps : (fpr_raw > 1.0 - eps ? 1.0 - eps : fpr_raw);
      double bns = fabs(tk_probit(tpr) - tk_probit(fpr));
      if (bns > 0) {
        tk_rank_t r = { (int64_t)f, bns };
        tk_rvec_hmin(per_dim_heaps[b], top_k, r);
      }
    }));

    tk_iumap_destroy(active_counts);
    tk_ivec_destroy(feat_counts);
    tk_ivec_destroy(label_counts);

    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_desc(per_dim_heaps[h], 0, per_dim_heaps[h]->n);

    tk_iuset_t *union_set = tk_iuset_create(0, 0);
    int absent;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++)
        tk_iuset_put(union_set, per_dim_heaps[h]->a[i].i, &absent);
    }

    tk_ivec_t *ids_union = tk_ivec_create(L, tk_iuset_size(union_set), 0, 0);
    ids_union->n = 0;
    int64_t key;
    tk_umap_foreach_keys(union_set, key, ({ ids_union->a[ids_union->n++] = key; }));
    tk_ivec_asc(ids_union, 0, ids_union->n);

    tk_iumap_t *id_to_union_idx = tk_iumap_create(0, 0);
    for (uint64_t i = 0; i < ids_union->n; i++) {
      int absent;
      khint_t k = tk_iumap_put(id_to_union_idx, ids_union->a[i], &absent);
      kh_value(id_to_union_idx, k) = (int64_t)i;
    }

    uint64_t total_ids = 0;
    for (uint64_t h = 0; h < n_hidden; h++)
      total_ids += per_dim_heaps[h]->n;

    tk_ivec_t *offsets = tk_ivec_create(L, n_hidden + 1, 0, 0);
    tk_ivec_t *ids = tk_ivec_create(L, total_ids, 0, 0);
    tk_dvec_t *weights = tk_dvec_create(L, total_ids, 0, 0);

    offsets->a[0] = 0;
    ids->n = 0;
    weights->n = 0;
    for (uint64_t h = 0; h < n_hidden; h++) {
      for (uint64_t i = 0; i < per_dim_heaps[h]->n; i++) {
        int64_t feat_id = per_dim_heaps[h]->a[i].i;
        khint_t k = tk_iumap_get(id_to_union_idx, feat_id);
        ids->a[ids->n++] = kh_value(id_to_union_idx, k);
        weights->a[weights->n++] = per_dim_heaps[h]->a[i].d;
      }
      offsets->a[h + 1] = (int64_t)ids->n;
    }
    offsets->n = n_hidden + 1;

    tk_iumap_destroy(id_to_union_idx);
    tk_iuset_destroy(union_set);
    for (uint64_t h = 0; h < n_hidden; h++)
      tk_rvec_destroy(per_dim_heaps[h]);
    free(per_dim_heaps);

  } else {
    tk_ivec_create(L, 0, 0, 0);
    tk_ivec_create(L, 1, 0, 0);
    tk_ivec_create(L, 0, 0, 0);
    tk_dvec_create(L, 0, 0, 0);
  }
}

#endif
