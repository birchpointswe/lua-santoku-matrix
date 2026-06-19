#ifndef TK_CVEC_EXT_H
#define TK_CVEC_EXT_H

#if defined(_OPENMP) && !defined(__EMSCRIPTEN__)
#include <omp.h>
#endif
#ifdef __aarch64__
#include <arm_neon.h>
#elif defined(__AVX512F__)
#include <immintrin.h>
#endif
#include <santoku/cvec/base.h>
#include <santoku/iuset.h>
#include <santoku/ivec.h>
#include <santoku/dvec.h>
#include <santoku/rvec.h>
#include <santoku/iumap.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <stdatomic.h>
#include <math.h>

#ifndef TK_CVEC_BITS
#define TK_CVEC_BITS CHAR_BIT
#endif
#ifndef TK_CVEC_BITS_BYTES
#define TK_CVEC_BITS_BYTES(n) (((n) + CHAR_BIT - 1) / CHAR_BIT)
#endif
#ifndef TK_CVEC_BITS_BYTE
#define TK_CVEC_BITS_BYTE(n) ((n) / CHAR_BIT)
#endif
#ifndef TK_CVEC_BITS_BIT
#define TK_CVEC_BITS_BIT(n) ((n) % CHAR_BIT)
#endif

#define TK_CVEC_ZERO_MASK 0x00
#define TK_CVEC_ALL_MASK 0xFF
#define TK_CVEC_POS_MASK 0x55
#define TK_CVEC_NEG_MASK 0xAA

static const uint8_t TK_POPCOUNT_TABLE[256] = {
  0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};

static inline uint8_t tk_popcount8(uint8_t x) {
  return TK_POPCOUNT_TABLE[x];
}


static inline uint64_t tk_cvec_bits_popcount(const uint8_t *data, uint64_t n_bits);
static inline uint64_t tk_cvec_bits_popcount_serial(const uint8_t *data, uint64_t n_bits);
static inline uint64_t tk_cvec_bits_hamming(const uint8_t *a, const uint8_t *b, uint64_t n_bits);
static inline uint64_t tk_cvec_bits_hamming_serial(const uint8_t *a, const uint8_t *b, uint64_t n_bits);
static inline uint64_t tk_cvec_bits_hamming_mask(const uint8_t *a, const uint8_t *b, const uint8_t *mask, uint64_t n_bits);
static inline uint64_t tk_cvec_bits_hamming_mask_serial(const uint8_t *a, const uint8_t *b, const uint8_t *mask, uint64_t n_bits);
static inline void tk_cvec_bits_and(uint8_t *out, const uint8_t *a, const uint8_t *b, uint64_t n_bits);
static inline void tk_cvec_bits_and_serial(uint8_t *out, const uint8_t *a, const uint8_t *b, uint64_t n_bits);
static inline void tk_cvec_bits_or(uint8_t *out, const uint8_t *a, const uint8_t *b, uint64_t n_bits);
static inline void tk_cvec_bits_or_serial(uint8_t *out, const uint8_t *a, const uint8_t *b, uint64_t n_bits);
static inline void tk_cvec_bits_xor(uint8_t *out, const uint8_t *a, const uint8_t *b, uint64_t n_bits);
static inline void tk_cvec_bits_xor_serial(uint8_t *out, const uint8_t *a, const uint8_t *b, uint64_t n_bits);
static inline tk_ivec_t *tk_cvec_bits_to_ivec(lua_State *L, tk_cvec_t *bitmap, uint64_t n_features);
static inline tk_ivec_t *tk_cvec_bits_to_ivec_serial(lua_State *L, tk_cvec_t *bitmap, uint64_t n_features);
static inline tk_cvec_t *tk_cvec_bits_from_ivec(lua_State *L, tk_ivec_t *set_bits, uint64_t n_samples, uint64_t n_features);
static inline tk_cvec_t *tk_cvec_bits_from_ivec_serial(lua_State *L, tk_ivec_t *set_bits, uint64_t n_samples, uint64_t n_features);
static inline int tk_cvec_bits_extend_mapped(tk_cvec_t *base, tk_cvec_t *ext, tk_ivec_t *aids, tk_ivec_t *bids, uint64_t n_base_features, uint64_t n_ext_features, bool project);
static inline int tk_cvec_bits_extend_mapped_serial(tk_cvec_t *base, tk_cvec_t *ext, tk_ivec_t *aids, tk_ivec_t *bids, uint64_t n_base_features, uint64_t n_ext_features, bool project);

static inline void tk_cvec_bits_to_ascii(lua_State *L, const tk_cvec_t *bitmap, uint64_t start_bit, uint64_t end_bit);
static inline void tk_cvec_bits_to_ascii_serial(lua_State *L, const tk_cvec_t *bitmap, uint64_t start_bit, uint64_t end_bit);

static inline tk_cvec_t *tk_cvec_bits_flip_interleave (tk_cvec_t *v, uint64_t n_features) {
  if (n_features == 0)
    return v;

  uint64_t n_samples = v->n / TK_CVEC_BITS_BYTES(n_features);
  if (n_samples == 0)
    return v;

  uint64_t input_bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t output_features = n_features * 2;
  uint64_t output_bytes_per_sample = TK_CVEC_BITS_BYTES(output_features);
  uint64_t output_size = n_samples * output_bytes_per_sample;
  tk_cvec_ensure(v, output_size);
  uint8_t *data = (uint8_t *)v->a;
  uint64_t full_bytes = n_features / CHAR_BIT;
  uint64_t remaining_bits = n_features % CHAR_BIT;

  uint8_t *temp = malloc(input_bytes_per_sample);
  if (!temp) {
    v->n = output_size;
    return NULL;
  }
  for (uint64_t s = n_samples; s > 0; s--) {
    uint64_t idx = s - 1;
    uint64_t input_sample_offset = idx * input_bytes_per_sample;
    uint64_t output_sample_offset = idx * output_bytes_per_sample;

    memcpy(temp, data + input_sample_offset, input_bytes_per_sample);
    memset(data + output_sample_offset, 0, output_bytes_per_sample);
    memcpy(data + output_sample_offset, temp, input_bytes_per_sample);
    uint64_t second_half_bit_offset = n_features;
    uint64_t second_half_byte_offset = second_half_bit_offset / CHAR_BIT;
    uint8_t second_half_bit_shift = second_half_bit_offset & (CHAR_BIT - 1);

    if (second_half_bit_shift == 0) {
      for (uint64_t i = 0; i < full_bytes; i ++) {
        data[output_sample_offset + second_half_byte_offset + i] = ~temp[i];
      }
      if (remaining_bits > 0) {
        uint8_t last_byte = temp[full_bytes];
        uint8_t mask = (1u << remaining_bits) - 1;
        data[output_sample_offset + second_half_byte_offset + full_bytes] = (~last_byte) & mask;
      }
    } else {
      uint8_t carry = 0;
      for (uint64_t i = 0; i < full_bytes; i ++) {
        uint8_t complement = ~temp[i];
        data[output_sample_offset + second_half_byte_offset + i] |= (complement << second_half_bit_shift) | carry;
        carry = complement >> (CHAR_BIT - second_half_bit_shift);
      }
      if (remaining_bits > 0) {
        uint8_t last_byte = temp[full_bytes];
        uint8_t mask = (1u << remaining_bits) - 1;
        uint8_t complement = (~last_byte) & mask;
        data[output_sample_offset + second_half_byte_offset + full_bytes] |= (complement << second_half_bit_shift) | carry;
        if (second_half_bit_shift + remaining_bits > CHAR_BIT) {
          carry = complement >> (CHAR_BIT - second_half_bit_shift);
          data[output_sample_offset + second_half_byte_offset + full_bytes + 1] |= carry;
        }
      } else if (carry != 0) {
        data[output_sample_offset + second_half_byte_offset + full_bytes] |= carry;
      }
    }
  }

  free(temp);
  v->n = output_size;
  return v;
}

static inline int tk_cvec_push_str (tk_cvec_t *v, const char *s)
{
  for (char *p = (char *) s; *p; p ++) {
    if (tk_cvec_push(v, *p) != 0)
      return -1;
  }
  return 0;
}

static inline tk_cvec_t *tk_cvec_bits_extend (
  tk_cvec_t *base,
  tk_cvec_t *ext,
  uint64_t n_base_features,
  uint64_t n_ext_features
) {
  uint64_t n_samples = base->n / TK_CVEC_BITS_BYTES(n_base_features);
  uint64_t n_total_features = n_base_features + n_ext_features;
  uint64_t base_bytes_per_sample = TK_CVEC_BITS_BYTES(n_base_features);
  uint64_t ext_bytes_per_sample = TK_CVEC_BITS_BYTES(n_ext_features);
  uint64_t total_bytes_per_sample = TK_CVEC_BITS_BYTES(n_total_features);
  tk_cvec_ensure(base, n_samples * total_bytes_per_sample);
  uint8_t *base_data = (uint8_t *)base->a;
  uint8_t *ext_data = (uint8_t *)ext->a;
  uint8_t *temp = malloc(total_bytes_per_sample);
  if (!temp)
    return NULL;
  for (int64_t s = (int64_t) n_samples - 1; s >= 0; s--) {
    memset(temp, 0, total_bytes_per_sample);
    uint64_t base_offset = (uint64_t) s * base_bytes_per_sample;
    memcpy(temp, base_data + base_offset, base_bytes_per_sample);
    uint64_t ext_offset = (uint64_t) s * ext_bytes_per_sample;
    uint64_t ext_bit_offset = n_base_features;
    uint64_t ext_byte_offset = ext_bit_offset / CHAR_BIT;
    uint8_t ext_bit_shift = ext_bit_offset % CHAR_BIT;
    if (ext_bit_shift == 0) {
      memcpy(temp + ext_byte_offset, ext_data + ext_offset, ext_bytes_per_sample);
    } else {
      uint8_t carry = 0;
      for (uint64_t i = 0; i < ext_bytes_per_sample; i++) {
        uint8_t byte = ext_data[ext_offset + i];
        temp[ext_byte_offset + i] |= (byte << ext_bit_shift) | carry;
        carry = byte >> (CHAR_BIT - ext_bit_shift);
      }
      if (carry != 0 && ext_byte_offset + ext_bytes_per_sample < total_bytes_per_sample) {
        temp[ext_byte_offset + ext_bytes_per_sample] |= carry;
      }
    }
    uint64_t out_offset = (uint64_t) s * total_bytes_per_sample;
    memcpy(base_data + out_offset, temp, total_bytes_per_sample);
  }
  free(temp);
  base->n = n_samples * total_bytes_per_sample;
  return base;
}

#define TK_GENERATE_SINGLE
#include <santoku/parallel/tpl.h>
#include <santoku/cvec/ext_tpl.h>
#undef TK_GENERATE_SINGLE

#include <santoku/parallel/tpl.h>
#include <santoku/cvec/ext_tpl.h>

static inline void tk_cvec_bits_topk (
  lua_State *L,
  tk_cvec_t *queries,
  tk_cvec_t *corpus,
  uint64_t n_queries,
  uint64_t n_corpus,
  uint64_t n_bits,
  uint64_t k
) {
  uint64_t bytes_per = TK_CVEC_BITS_BYTES(n_bits);
  if (k == 0 || n_queries == 0 || n_corpus == 0 || n_bits == 0) {
    tk_ivec_t *off = tk_ivec_create(L, n_queries + 1);
    off->n = n_queries + 1;
    memset(off->a, 0, (n_queries + 1) * sizeof(int64_t));
    tk_ivec_create(L, 0);
    tk_ivec_create(L, 0);
    return;
  }
  if (k > n_corpus) k = n_corpus;
  uint64_t total = n_queries * k;
  tk_ivec_t *offsets = tk_ivec_create(L, n_queries + 1);
  offsets->n = n_queries + 1;
  for (uint64_t i = 0; i <= n_queries; i++) offsets->a[i] = (int64_t)(i * k);
  tk_ivec_t *indices = tk_ivec_create(L, total);
  indices->n = total;
  tk_ivec_t *distances = tk_ivec_create(L, total);
  distances->n = total;
#if defined(_OPENMP) && !defined(__EMSCRIPTEN__)
  #pragma omp parallel
  {
    tk_rvec_t *heap = tk_rvec_create(NULL, k);
    #pragma omp for schedule(static)
    for (uint64_t q = 0; q < n_queries; q++) {
      tk_rvec_clear(heap);
      const uint8_t *qvec = (const uint8_t *)queries->a + q * bytes_per;
      for (uint64_t c = 0; c < n_corpus; c++) {
        const uint8_t *cvec = (const uint8_t *)corpus->a + c * bytes_per;
        uint64_t dist = tk_cvec_bits_hamming_serial(qvec, cvec, n_bits);
        tk_rvec_hmax(heap, k, tk_rank((int64_t)c, (double)dist));
      }
      tk_rvec_asc(heap, 0, heap->n);
      int64_t *idx_row = indices->a + q * k;
      int64_t *dist_row = distances->a + q * k;
      for (uint64_t h = 0; h < heap->n; h++) {
        idx_row[h] = heap->a[h].i;
        dist_row[h] = (int64_t)(heap->a[h].d);
      }
    }
    tk_rvec_destroy(heap);
  }
#else
  tk_rvec_t *heap = tk_rvec_create(NULL, k);
  for (uint64_t q = 0; q < n_queries; q++) {
    tk_rvec_clear(heap);
    const uint8_t *qvec = (const uint8_t *)queries->a + q * bytes_per;
    for (uint64_t c = 0; c < n_corpus; c++) {
      const uint8_t *cvec = (const uint8_t *)corpus->a + c * bytes_per;
      uint64_t dist = tk_cvec_bits_hamming_serial(qvec, cvec, n_bits);
      tk_rvec_hmax(heap, k, tk_rank((int64_t)c, (double)dist));
    }
    tk_rvec_asc(heap, 0, heap->n);
    int64_t *idx_row = indices->a + q * k;
    int64_t *dist_row = distances->a + q * k;
    for (uint64_t h = 0; h < heap->n; h++) {
      idx_row[h] = heap->a[h].i;
      dist_row[h] = (int64_t)(heap->a[h].d);
    }
  }
  tk_rvec_destroy(heap);
#endif
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

#endif
