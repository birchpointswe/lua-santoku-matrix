#ifndef TK_CVEC_EXT_H
#define TK_CVEC_EXT_H

#include <omp.h>
#include <santoku/cvec/base.h>
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
#ifndef tk_cvec_byte_popcount
#ifdef __GNUC__
#define tk_cvec_byte_popcount(x) __builtin_popcount(x)
#else
static inline int tk_cvec_byte_popcount (unsigned int x) {
  x = x - ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  return (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}
#endif
#endif

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

static inline uint64_t tk_cvec_bits_popcount (
  const uint8_t *data,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  uint64_t count = 0;
  for (uint64_t i = 0; i < full_bytes - (rem_bits > 0); i ++)
    count += (uint64_t) tk_cvec_byte_popcount(data[i]);
  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    count += (uint64_t) tk_cvec_byte_popcount(data[full_bytes - 1] & mask);
  }
  return count;
}

static inline uint64_t tk_cvec_bits_hamming (
  const uint8_t *a,
  const uint8_t *b,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  uint64_t dist = 0;
  for (uint64_t i = 0; i < full_bytes - (rem_bits > 0); i ++)
    dist += (uint64_t) tk_cvec_byte_popcount(a[i] ^ b[i]);
  if (rem_bits > 0) {
    uint8_t x = a[full_bytes - 1] ^ b[full_bytes - 1];
    uint8_t mask = (1U << rem_bits) - 1;
    dist += (uint64_t) tk_cvec_byte_popcount(x & mask);
  }
  return dist;
}

static inline uint64_t tk_cvec_bits_hamming_mask (
  const uint8_t *a,
  const uint8_t *b,
  const uint8_t *mask,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  uint64_t dist = 0;
  for (uint64_t i = 0; i < full_bytes - (rem_bits > 0); i ++)
    dist += (uint64_t) tk_cvec_byte_popcount((a[i] ^ b[i]) & mask[i]);
  if (rem_bits > 0) {
    uint8_t x = (a[full_bytes - 1] ^ b[full_bytes - 1]) & mask[full_bytes - 1];
    uint8_t bit_mask = (1U << rem_bits) - 1;
    dist += (uint64_t) tk_cvec_byte_popcount(x & bit_mask);
  }
  return dist;
}

static inline void tk_cvec_bits_and (
  uint8_t *out,
  const uint8_t *a,
  const uint8_t *b,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  for (uint64_t i = 0; i < full_bytes - (rem_bits > 0); i ++)
    out[i] = a[i] & b[i];
  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    out[full_bytes - 1] = (a[full_bytes - 1] & b[full_bytes - 1]) & mask;
  }
}

static inline void tk_cvec_bits_or (
  uint8_t *out,
  const uint8_t *a,
  const uint8_t *b,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  for (uint64_t i = 0; i < full_bytes - (rem_bits > 0); i ++)
    out[i] = a[i] | b[i];
  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    out[full_bytes - 1] = (a[full_bytes - 1] | b[full_bytes - 1]) & mask;
  }
}

static inline void tk_cvec_bits_xor (
  uint8_t *out,
  const uint8_t *a,
  const uint8_t *b,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  for (uint64_t i = 0; i < full_bytes - (rem_bits > 0); i ++)
    out[i] = a[i] ^ b[i];
  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    out[full_bytes - 1] = (a[full_bytes - 1] ^ b[full_bytes - 1]) & mask;
  }
}

static inline tk_ivec_t *tk_cvec_bits_to_ivec (
  lua_State *L,
  tk_cvec_t *bitmap,
  uint64_t n_features
) {
  uint8_t *data = (uint8_t *)bitmap->a;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t n_samples = bitmap->n / bytes_per_sample;
  uint64_t total_bits = 0;
  #pragma omp parallel for reduction(+:total_bits)
  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t offset = s * bytes_per_sample;
    total_bits += tk_cvec_bits_popcount(data + offset, n_features);
  }
  tk_ivec_t *out = tk_ivec_create(L, total_bits, 0, 0);
  uint64_t *sample_counts = (uint64_t *)calloc(n_samples, sizeof(uint64_t));
  if (!sample_counts) {
    tk_ivec_destroy(out);
    return NULL;
  }

  #pragma omp parallel for
  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t offset = s * bytes_per_sample;
    sample_counts[s] = tk_cvec_bits_popcount(data + offset, n_features);
  }

  uint64_t *offsets = (uint64_t *)malloc((n_samples + 1) * sizeof(uint64_t));
  if (!offsets) {
    free(sample_counts);
    tk_ivec_destroy(out);
    return NULL;
  }
  offsets[0] = 0;
  for (uint64_t s = 0; s < n_samples; s++)
    offsets[s + 1] = offsets[s] + sample_counts[s];

  #pragma omp parallel for
  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t offset = s * bytes_per_sample;
    uint64_t write_pos = offsets[s];
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t byte_idx = f / CHAR_BIT;
      uint8_t bit_idx = f % CHAR_BIT;
      if (data[offset + byte_idx] & (1u << bit_idx)) {
        out->a[write_pos++] = (int64_t) (s * n_features + f);
      }
    }
  }

  out->n = total_bits;
  free(sample_counts);
  free(offsets);
  return out;
}

static inline tk_cvec_t *tk_cvec_bits_from_ivec (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_features
) {
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t total_bytes = n_samples * bytes_per_sample;

  tk_cvec_t *out = tk_cvec_create(L, total_bytes, 0, 0);
  uint8_t *data = (uint8_t *)out->a;
  memset(data, 0, total_bytes);

  #pragma omp parallel for
  for (size_t i = 0; i < set_bits->n; i++) {
    int64_t bit_idx = set_bits->a[i];
    uint64_t sample = (uint64_t) bit_idx / n_features;
    uint64_t feature = (uint64_t) bit_idx % n_features;

    if (bit_idx >= 0 && sample < n_samples && feature < n_features) {
      uint64_t offset = sample * bytes_per_sample;
      uint64_t byte_idx = feature / CHAR_BIT;
      uint8_t bit_pos = feature % CHAR_BIT;
      #pragma omp atomic
      data[offset + byte_idx] |= (1u << bit_pos);
    }
  }

  out->n = total_bytes;
  return out;
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

static inline int tk_cvec_bits_extend_mapped (
  tk_cvec_t *base,
  tk_cvec_t *ext,
  tk_ivec_t *aids,
  tk_ivec_t *bids,
  uint64_t n_base_features,
  uint64_t n_ext_features,
  bool project
) {
  uint64_t n_total_features = n_base_features + n_ext_features;
  uint64_t base_bytes_per_sample = TK_CVEC_BITS_BYTES(n_base_features);
  uint64_t ext_bytes_per_sample = TK_CVEC_BITS_BYTES(n_ext_features);
  uint64_t total_bytes_per_sample = TK_CVEC_BITS_BYTES(n_total_features);
  tk_iumap_t *a_id_to_pos = tk_iumap_from_ivec(0, aids);
  if (!a_id_to_pos)
    return -1;
  uint64_t n_only_b = 0;
  int64_t *b_to_final = (int64_t *)calloc(bids->n, sizeof(int64_t));
  if (!b_to_final) {
    tk_iumap_destroy(a_id_to_pos);
    return -1;
  }
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
  size_t old_aids_n = aids->n;
  if (!project) {
    if (tk_ivec_ensure(aids, final_n_samples) != 0) {
      free(b_to_final);
      tk_iumap_destroy(a_id_to_pos);
      return -1;
    }
    for (size_t i = 0; i < bids->n; i++) {
      if (b_to_final[i] >= (int64_t)old_aids_n) {
        aids->a[aids->n++] = bids->a[i];
      }
    }
  }
  if (tk_cvec_ensure(base, final_n_samples * total_bytes_per_sample) != 0) {
    free(b_to_final);
    tk_iumap_destroy(a_id_to_pos);
    return -1;
  }
  uint8_t *base_data = (uint8_t *)base->a;
  uint8_t *ext_data = (uint8_t *)ext->a;
  tk_iumap_t *b_id_to_pos = tk_iumap_from_ivec(0, bids);
  if (!b_id_to_pos) {
    free(b_to_final);
    tk_iumap_destroy(a_id_to_pos);
    return -1;
  }
  uint8_t *new_data = calloc(final_n_samples, total_bytes_per_sample);
  if (!new_data) {
    free(b_to_final);
    tk_iumap_destroy(a_id_to_pos);
    tk_iumap_destroy(b_id_to_pos);
    return -1;
  }
  #pragma omp parallel for
  for (size_t ai = 0; ai < old_aids_n; ai++) {
    uint64_t dest_offset = ai * total_bytes_per_sample;
    uint64_t src_offset = ai * base_bytes_per_sample;
    memcpy(new_data + dest_offset, base_data + src_offset, base_bytes_per_sample);
    khint_t khi = tk_iumap_get(b_id_to_pos, aids->a[ai]);
    if (khi != tk_iumap_end(b_id_to_pos)) {
      int64_t b_idx = tk_iumap_val(b_id_to_pos, khi);
      uint64_t ext_src_offset = (uint64_t)b_idx * ext_bytes_per_sample;
      uint64_t ext_bit_offset = n_base_features;
      uint64_t ext_byte_offset = ext_bit_offset / CHAR_BIT;
      uint8_t ext_bit_shift = ext_bit_offset % CHAR_BIT;
      if (ext_bit_shift == 0) {
        memcpy(new_data + dest_offset + ext_byte_offset, ext_data + ext_src_offset, ext_bytes_per_sample);
      } else {
        uint8_t carry = 0;
        for (uint64_t i = 0; i < ext_bytes_per_sample; i++) {
          uint8_t byte = ext_data[ext_src_offset + i];
          new_data[dest_offset + ext_byte_offset + i] |= (byte << ext_bit_shift) | carry;
          carry = byte >> (CHAR_BIT - ext_bit_shift);
        }
        if (carry != 0 && ext_byte_offset + ext_bytes_per_sample < total_bytes_per_sample) {
          new_data[dest_offset + ext_byte_offset + ext_bytes_per_sample] |= carry;
        }
      }
    }
  }
  #pragma omp parallel for
  for (size_t bi = 0; bi < bids->n; bi++) {
    int64_t final_pos = b_to_final[bi];
    if (final_pos >= (int64_t)old_aids_n) {
      uint64_t dest_offset = (uint64_t)final_pos * total_bytes_per_sample;
      uint64_t src_offset = bi * ext_bytes_per_sample;
      uint64_t ext_bit_offset = n_base_features;
      uint64_t ext_byte_offset = ext_bit_offset / CHAR_BIT;
      uint8_t ext_bit_shift = ext_bit_offset % CHAR_BIT;
      if (ext_bit_shift == 0) {
        memcpy(new_data + dest_offset + ext_byte_offset, ext_data + src_offset, ext_bytes_per_sample);
      } else {
        uint8_t carry = 0;
        for (uint64_t i = 0; i < ext_bytes_per_sample; i++) {
          uint8_t byte = ext_data[src_offset + i];
          new_data[dest_offset + ext_byte_offset + i] |= (byte << ext_bit_shift) | carry;
          carry = byte >> (CHAR_BIT - ext_bit_shift);
        }
        if (carry != 0 && ext_byte_offset + ext_bytes_per_sample < total_bytes_per_sample) {
          new_data[dest_offset + ext_byte_offset + ext_bytes_per_sample] |= carry;
        }
      }
    }
  }
  memcpy(base_data, new_data, final_n_samples * total_bytes_per_sample);
  free(new_data);
  free(b_to_final);
  tk_iumap_destroy(a_id_to_pos);
  tk_iumap_destroy(b_id_to_pos);
  base->n = final_n_samples * total_bytes_per_sample;
  return 0;
}

static inline void tk_cvec_bits_to_ascii (
  lua_State *L,
  const tk_cvec_t *bitmap,
  uint64_t start_bit,
  uint64_t end_bit
) {
  uint64_t n_bits = end_bit - start_bit;
  tk_cvec_t *ascii = tk_cvec_create(L, n_bits + 1, 0, 0);
  const uint8_t *data = (const uint8_t *)bitmap->a;
  char *out = ascii->a;
  #pragma omp parallel for
  for (uint64_t i = 0; i < n_bits; i++) {
    uint64_t bit_pos = start_bit + i;
    uint64_t byte_idx = bit_pos / 8;
    uint64_t bit_idx = bit_pos % 8;
    out[i] = (data[byte_idx] & (1 << bit_idx)) ? '1' : '0';
  }
  ascii->n = n_bits;
  out[n_bits] = '\0';
  lua_pushstring(L, out);
  lua_remove(L, -2);
}

static inline int tk_cvec_push_str (tk_cvec_t *v, const char *s)
{
  for (char *p = (char *) s; *p; p ++) {
    if (tk_cvec_push(v, *p) != 0)
      return -1;
  }
  return 0;
}

#endif
