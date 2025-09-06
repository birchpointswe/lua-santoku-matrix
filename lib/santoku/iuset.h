#ifndef TK_IUSET_H
#define TK_IUSET_H

#include <santoku/klib.h>
#include <santoku/ivec.h>
#include <santoku/cvec.h>
#include <santoku/cvec/ext.h>
#include <santoku/rvec.h>
#include <santoku/rvec/ext.h>
#include <santoku/dvec.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

KHASH_INIT(tk_iuset, int64_t, char, 0, kh_int64_hash_func, kh_int64_hash_equal)
typedef khash_t(tk_iuset) tk_iuset_t;

#define tk_iuset_put(...) kh_put(tk_iuset, __VA_ARGS__)
#define tk_iuset_get(...) kh_get(tk_iuset, __VA_ARGS__)
#define tk_iuset_del(...) kh_del(tk_iuset, __VA_ARGS__)
#define tk_iuset_exist(...) kh_exist(__VA_ARGS__)
#define tk_iuset_key(...) kh_key(__VA_ARGS__)
#define tk_iuset_value(...) kh_value(__VA_ARGS__)
#define tk_iuset_begin(...) kh_begin(__VA_ARGS__)
#define tk_iuset_end(...) kh_end(__VA_ARGS__)
#define tk_iuset_size(...) kh_size(__VA_ARGS__)
#define tk_iuset_resize(...) kh_resize(tk_iuset, __VA_ARGS__)
#define tk_iuset_clear(...) kh_clear(tk_iuset, __VA_ARGS__)
#define tk_iuset_destroy(...) kh_destroy(tk_iuset, __VA_ARGS__)
#define tk_iuset_create() kh_init(tk_iuset)
#define tk_iuset_contains(h, v) (tk_iuset_get(h, v) != tk_iuset_end(h))

#define tk_iuset_foreach(h, kvar, code) { khint_t __i; \
	for (__i = kh_begin(h); __i != kh_end(h); __i ++) { \
		if (!kh_exist(h,__i)) continue; \
		(kvar) = kh_key(h,__i); \
		code;	\
	} }

static inline void tk_iuset_union (tk_iuset_t *a, tk_iuset_t *b)
{
  int kha;
  int64_t x;
  tk_iuset_foreach(b, x, ({
    tk_iuset_put(a, x, &kha);
  }));
}

static inline double tk_iuset_jaccard (tk_iuset_t *a, tk_iuset_t *b)
{
  uint64_t intersection = 0;
  uint64_t union_count = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    if (tk_iuset_contains(b, x))
      intersection ++;
  }));
  union_count = tk_iuset_size(a) + tk_iuset_size(b) - intersection;
  if (union_count == 0)
    return 0.0;
  return (double) intersection / (double) union_count;
}

static inline void tk_iuset_intersect (tk_iuset_t *a, tk_iuset_t *b)
{
  khint_t i = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    if (!tk_iuset_contains(b, x))
      kh_del(tk_iuset, a, i);
  }))
}

static inline void tk_iuset_difference (tk_iuset_t *a, tk_iuset_t *b)
{
  khint_t i = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    if (tk_iuset_contains(b, x))
      kh_del(tk_iuset, a, i);
  }))
}


#include <santoku/iumap.h>

static inline void tk_iuset_union_iumap (tk_iuset_t *a, tk_iumap_t *b)
{
  int kha;
  for (khint_t i = tk_iumap_begin(b); i < tk_iumap_end(b); i ++) {
    if (tk_iumap_exist(b, i)) {
      int64_t key = tk_iumap_key(b, i);
      tk_iuset_put(a, key, &kha);
    }
  }
}

static inline void tk_iuset_intersect_iumap (tk_iuset_t *a, tk_iumap_t *b)
{
  khint_t i = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    khint_t khi = tk_iumap_get(b, x);
    if (khi == tk_iumap_end(b))
      kh_del(tk_iuset, a, i);
  }))
}

static inline void tk_iuset_difference_iumap (tk_iuset_t *a, tk_iumap_t *b)
{
  khint_t i = 0;
  int64_t x;
  tk_iuset_foreach(a, x, ({
    khint_t khi = tk_iumap_get(b, x);
    if (khi != tk_iumap_end(b))
      kh_del(tk_iuset, a, i);
  }))
}

static inline tk_iuset_t *tk_iuset_from_ivec (tk_ivec_t *v)
{
  int kha;
  tk_iuset_t *s = tk_iuset_create();
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


    sample_set = tk_iuset_from_ivec(sample_ids);


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
    sample_set = tk_iuset_from_ivec(sample_ids);
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

  if (dest_sample == 0) {
    tk_ivec_clear(dest);
  }


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

    sample_set = tk_iuset_from_ivec(sample_ids);
    


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
  if (sample_ids != NULL && sample_ids->n > 0) {
    sample_set = tk_iuset_from_ivec(sample_ids);
  }

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

    if (dest_stride > 0) {
      total_bits = (dest_sample + 1) * row_stride_bits;
    }
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
    if (total_bytes > old_size) {
      memset(dest_data + old_size, 0, total_bytes - old_size);
    }
  }


  uint64_t write_sample = dest_sample;
  uint64_t write_bit_offset = 0;
  
  if (use_packed) {

    write_bit_offset = dest_sample * row_stride_bits;
  }
  
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
  tk_iuset_foreach(S, k, ({
    out->a[out->n ++] = k;
  }));
  return out;
}

static inline tk_ivec_t *tk_ivec_from_iuset (lua_State *L, tk_iuset_t *s)
{
  tk_ivec_t *v = tk_ivec_create(L, tk_iuset_size(s), 0, 0);
  int64_t x;
  v->n = 0;
  tk_iuset_foreach(s, x, ({
    v->a[v->n ++] = x;
  }))
  return v;
}

static inline tk_ivec_t *tk_ivec_top_select (
  lua_State *L,
  tk_iuset_t *selected
) {

  tk_ivec_t *top_v = tk_ivec_create(L, 0, 0, 0);
  int64_t sel;
  tk_iuset_foreach(selected, sel, ({
    tk_ivec_push(top_v, sel);
  }));
  tk_ivec_shrink(top_v);
  tk_iuset_destroy(selected);
  return top_v;
}

static inline void tk_ivec_select_union (
  lua_State *L,
  tk_iuset_t *selected,
  tk_rvec_t *rankings,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k,
  uint64_t trunc
) {

  int kha;
  uint64_t *offsets = tk_malloc(L, n_hidden * sizeof(uint64_t));
  memset(offsets, 0, n_hidden * sizeof(uint64_t));
  while (tk_iuset_size(selected) < top_k) {
    bool advanced = false;
    for (uint64_t j = 0; j < n_hidden; j ++) {
      tk_rank_t *rankings_h = rankings->a + j * n_visible;
      while (offsets[j] < n_visible) {
        tk_rank_t candidate = rankings_h[offsets[j]];
        offsets[j] ++;
        if (candidate.i >= (int64_t) (n_visible - trunc))
          continue;
        tk_iuset_put(selected, (int64_t) candidate.i, &kha);
        advanced = true;
        break;
      }
      if (tk_iuset_size(selected) >= top_k)
        break;
    }
    if (!advanced)
      break;
  }

  free(offsets);
}

static inline tk_ivec_t *tk_ivec_top_generic (
  lua_State *L,
  tk_dvec_t *scores,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k,
  uint64_t trunc
) {
  tk_iuset_t *selected = tk_iuset_create();
  tk_rvec_t *rankings = tk_rvec_rankings(L, scores, n_visible, n_hidden);
  tk_ivec_select_union(L, selected, rankings, n_visible, n_hidden, top_k, trunc);
  tk_rvec_destroy(rankings);
  return tk_ivec_top_select(L, selected);
}

static inline tk_ivec_t *tk_ivec_bits_top_mi (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k,
  unsigned int n_threads
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);
  tk_dvec_t *scores = tk_ivec_bits_score_mi(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, n_threads);
  int iscores = tk_lua_absindex(L, -1);
  tk_ivec_t *out = tk_ivec_top_generic(L, scores, n_visible, n_hidden, top_k, 0);
  lua_pushvalue(L, iscores);
  return out;
}

static inline tk_ivec_t *tk_ivec_bits_top_chi2 (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k,
  unsigned int n_threads
) {
  tk_ivec_asc(set_bits, 0, set_bits->n);
  tk_dvec_t *scores = tk_ivec_bits_score_chi2(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, n_threads);
  int iscores = tk_lua_absindex(L, -1);
  tk_ivec_t *out = tk_ivec_top_generic(L, scores, n_visible, n_hidden, top_k, 0);
  lua_pushvalue(L, iscores);
  return out;
}

static inline tk_ivec_t *tk_ivec_bits_top_entropy (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_hidden,
  uint64_t top_k,
  unsigned int n_threads
) {
  tk_dvec_t *scores = tk_ivec_bits_score_entropy(L, set_bits, n_samples, n_hidden, n_threads);
  int iscores = tk_lua_absindex(L, -1);
  tk_rvec_t *rankings = tk_rvec_from_dvec(L, scores);
  tk_rvec_kdesc(rankings, top_k, 0, rankings->n);
  tk_ivec_t *out = tk_ivec_from_rvec(L, rankings);
  out->n = top_k < out->n ? top_k : out->n;
  lua_pushvalue(L, iscores);
  lua_remove(L, -3);
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
  uint64_t top_k,
  unsigned int n_threads
) {
  tk_dvec_t *scores = tk_cvec_bits_score_mi(L, bitmap, codes, labels, n_samples, n_features, n_hidden, n_threads);
  int iscores = tk_lua_absindex(L, -1);
  tk_ivec_t *out = tk_ivec_top_generic(L, scores, n_features, n_hidden, top_k, 0);
  lua_pushvalue(L, iscores);
  return out;
}

static inline tk_ivec_t *tk_cvec_bits_top_chi2 (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  uint64_t top_k,
  unsigned int n_threads
) {
  tk_dvec_t *scores = tk_cvec_bits_score_chi2(L, bitmap, codes, labels, n_samples, n_features, n_hidden, n_threads);
  int iscores = tk_lua_absindex(L, -1);
  tk_ivec_t *out = tk_ivec_top_generic(L, scores, n_features, n_hidden, top_k, 0);
  lua_pushvalue(L, iscores);
  return out;
}

static inline tk_ivec_t *tk_cvec_bits_top_entropy (
  lua_State *L,
  tk_cvec_t *codes,
  uint64_t n_samples,
  uint64_t n_hidden,
  uint64_t top_k,
  unsigned int n_threads
) {
  tk_dvec_t *scores = tk_cvec_bits_score_entropy(L, codes, n_samples, n_hidden, n_threads);
  int iscores = tk_lua_absindex(L, -1);
  tk_rvec_t *rankings = tk_rvec_from_dvec(L, scores);
  tk_rvec_kdesc(rankings, top_k, 0, rankings->n);
  tk_ivec_t *out = tk_ivec_from_rvec(L, rankings);
  out->n = top_k < out->n ? top_k : out->n;
  lua_pushvalue(L, iscores);
  lua_remove(L, -3);
  return out;
}

#endif
