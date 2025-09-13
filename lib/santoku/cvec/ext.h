#ifndef TK_CVEC_EXT_H
#define TK_CVEC_EXT_H

#include <santoku/cvec/base.h>
#include <santoku/ivec.h>
#include <santoku/dvec.h>
#include <santoku/rvec.h>
#include <santoku/iumap.h>
#include <santoku/threads.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <stdatomic.h>
#include <math.h>


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


typedef enum {
  TK_CVEC_ENTROPY,
  TK_CVEC_CHI2,
  TK_CVEC_MI,
} tk_cvec_stage_t;

typedef struct {
  uint64_t n_samples, n_features, n_hidden;
  tk_cvec_t *bitmap;
  tk_dvec_t *scores;
  tk_ivec_t *counts, *feat_counts;
  tk_ivec_t *active_counts, *global_counts;
  tk_ivec_t *labels;
  char *codes;
  atomic_ulong *bit_counts;
} tk_cvec_ctx_t;

typedef struct {
  tk_cvec_ctx_t *state;
  uint64_t hfirst, hlast;
  uint64_t sfirst, slast;
  uint64_t ffirst, flast;
} tk_cvec_thread_t;

static inline void tk_cvec_worker (void *dp, int sig)
{
  tk_cvec_thread_t *data = (tk_cvec_thread_t *) dp;
  tk_cvec_ctx_t *state = data->state;
  uint64_t n_samples = state->n_samples;
  uint64_t n_features = state->n_features;
  uint64_t n_hidden = state->n_hidden;
  tk_dvec_t *scores = state->scores;
  tk_ivec_t *counts = state->counts;
  tk_ivec_t *active_counts = state->active_counts;
  tk_ivec_t *global_counts = state->global_counts;
  tk_ivec_t *feat_counts = state->feat_counts;
  atomic_ulong *bit_counts = state->bit_counts;
  char *codes = state->codes;
  uint64_t chunks = TK_CVEC_BITS_BYTES(n_hidden);

  switch ((tk_cvec_stage_t) sig) {

    case TK_CVEC_CHI2:
      for (uint64_t b = data->hfirst; b <= data->hlast; b++) {
        double *scores_b = scores->a + b * n_features;
        for (uint64_t f = 0; f < n_features; f++) {
          int64_t A = active_counts->a[f * n_hidden + b];
          int64_t G = global_counts->a[b];
          int64_t C = feat_counts->a[f];
          if (C == 0 || G == 0 || C == (int64_t) n_samples || G == (int64_t) n_samples) {
            scores_b[f] = 0.0;
            continue;
          }
          int64_t B = G - A;
          int64_t C_ = C - A;
          int64_t D = (int64_t) n_samples - C - B;
          double n = (double) n_samples;
          double E_A = ((double) C * (double) G) / n;
          double E_B = ((double)(n - C) * (double) G) / n;
          double E_C = ((double) C * (double)(n - G)) / n;
          double E_D = ((double)(n - C) * (double)(n - G)) / n;
          double chi2 = 0.0;
          if (E_A > 0)
            chi2 += ((A - E_A)*(A - E_A)) / E_A;
          if (E_B > 0)
            chi2 += ((B - E_B)*(B - E_B)) / E_B;
          if (E_C > 0)
            chi2 += ((C_ - E_C)*(C_ - E_C)) / E_C;
          if (E_D > 0)
            chi2 += ((D - E_D)*(D - E_D)) / E_D;
          scores_b[f] = chi2;
        }
      }
      break;

    case TK_CVEC_MI:
      for (int64_t j = (int64_t) data->hfirst; j <= (int64_t) data->hlast; j++) {
        double *scores_h = scores->a + j * (int64_t) n_features;
        for (int64_t i = 0; i < (int64_t) n_features; i++) {
          int64_t *c = counts->a + i * (int64_t) n_hidden * 4 + j * 4;
          for (int k = 0; k < 4; k++)
            c[k] += 1;
          double total = c[0] + c[1] + c[2] + c[3];
          double mi = 0.0;
          if (total == 0.0)
            continue;
          for (unsigned int o = 0; o < 4; o++) {
            if (c[o] == 0)
              continue;
            double p_fb = c[o] / total;
            unsigned int f = o >> 1;
            unsigned int b = o & 1;
            double pf = (c[2] + c[3]) / total;
            if (f == 0) pf = 1.0 - pf;
            double pb = (c[1] + c[3]) / total;
            if (b == 0) pb = 1.0 - pb;
            double d = pf * pb;
            if (d > 0) mi += p_fb * log2(p_fb / d);
          }
          scores_h[i] = mi;
        }
      }
      break;

    case TK_CVEC_ENTROPY:
      for (uint64_t i = data->sfirst; i <= data->slast; i++) {
        for (uint64_t j = 0; j < n_hidden; j++) {
          uint64_t word = j / CHAR_BIT;
          uint64_t bit = j % CHAR_BIT;
          if (codes[i * chunks + word] & (1 << bit))
            atomic_fetch_add(bit_counts + j, 1);
        }
      }
      break;
  }
}

static inline void tk_cvec_bits_flip_interleave (
  tk_cvec_t *v,
  uint64_t n_features
) {
  if (n_features == 0)
    return;

  uint64_t n_samples = v->n / TK_CVEC_BITS_BYTES(n_features);
  if (n_samples == 0)
    return;

  uint64_t input_bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t output_features = n_features * 2;
  uint64_t output_bytes_per_sample = TK_CVEC_BITS_BYTES(output_features);
  uint64_t output_size = n_samples * output_bytes_per_sample;


  tk_cvec_ensure(v, output_size);
  uint8_t *data = (uint8_t *)v->a;


  uint64_t full_bytes = n_features / CHAR_BIT;
  uint64_t remaining_bits = n_features % CHAR_BIT;



  for (uint64_t s = n_samples; s > 0; s--) {
    uint64_t idx = s - 1;
    uint64_t input_sample_offset = idx * input_bytes_per_sample;
    uint64_t output_sample_offset = idx * output_bytes_per_sample;


    uint8_t *temp = malloc(input_bytes_per_sample);
    if (!temp) {

      v->n = output_size;
      return;
    }
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


    free(temp);
  }


  v->n = output_size;
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
  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t offset = s * bytes_per_sample;
    total_bits += tk_cvec_bits_popcount(data + offset, n_features);
  }

  tk_ivec_t *out = tk_ivec_create(L, total_bits, 0, 0);
  uint64_t idx = 0;

  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t offset = s * bytes_per_sample;
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t byte_idx = f / CHAR_BIT;
      uint8_t bit_idx = f % CHAR_BIT;
      if (data[offset + byte_idx] & (1u << bit_idx)) {
        out->a[idx++] = (int64_t) (s * n_features + f);
      }
    }
  }

  out->n = idx;
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

  for (size_t i = 0; i < set_bits->n; i++) {
    int64_t bit_idx = set_bits->a[i];
    uint64_t sample = (uint64_t) bit_idx / n_features;
    uint64_t feature = (uint64_t) bit_idx % n_features;

    if (bit_idx >= 0 && sample < n_samples && feature < n_features) {
      uint64_t offset = sample * bytes_per_sample;
      uint64_t byte_idx = feature / CHAR_BIT;
      uint8_t bit_pos = feature % CHAR_BIT;
      data[offset + byte_idx] |= (1u << bit_pos);
    }
  }

  out->n = total_bytes;
  return out;
}

static inline void tk_cvec_bits_rearrange (
  tk_cvec_t *bitmap,
  tk_ivec_t *ids,
  uint64_t n_features
) {
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t n_samples = ids->n;
  uint8_t *data = (uint8_t *)bitmap->a;


  uint8_t *temp = malloc(n_samples * bytes_per_sample);
  if (!temp) return;


  for (uint64_t i = 0; i < n_samples; i++) {
    int64_t src_idx = ids->a[i];
    if (src_idx < 0) continue;
    memcpy(temp + i * bytes_per_sample,
           data + (uint64_t) src_idx * bytes_per_sample,
           bytes_per_sample);
  }


  memcpy(data, temp, n_samples * bytes_per_sample);
  free(temp);

  bitmap->n = n_samples * bytes_per_sample;
}

static inline void tk_cvec_bits_extend (
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


  for (int64_t s = (int64_t) n_samples - 1; s >= 0; s--) {
    uint8_t *temp = malloc(total_bytes_per_sample);
    if (!temp) return;
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
    free(temp);
  }

  base->n = n_samples * total_bytes_per_sample;
}

static inline void tk_cvec_bits_extend_mapped (
  tk_cvec_t *base,
  tk_cvec_t *ext,
  tk_ivec_t *aids,
  tk_ivec_t *bids,
  uint64_t n_base_features,
  uint64_t n_ext_features
) {
  uint64_t n_total_features = n_base_features + n_ext_features;
  uint64_t base_bytes_per_sample = TK_CVEC_BITS_BYTES(n_base_features);
  uint64_t ext_bytes_per_sample = TK_CVEC_BITS_BYTES(n_ext_features);
  uint64_t total_bytes_per_sample = TK_CVEC_BITS_BYTES(n_total_features);


  tk_iumap_t *a_id_to_pos = tk_iumap_from_ivec(aids);


  uint64_t n_only_b = 0;
  int64_t *b_to_final = (int64_t *)malloc(bids->n * sizeof(int64_t));

  if (!b_to_final) {
    tk_iumap_destroy(a_id_to_pos);
    return;
  }


  int64_t next_pos = (int64_t)aids->n;
  for (size_t bi = 0; bi < bids->n; bi++) {
    khint_t khi = tk_iumap_get(a_id_to_pos, bids->a[bi]);
    if (khi != tk_iumap_end(a_id_to_pos)) {

      b_to_final[bi] = tk_iumap_value(a_id_to_pos, khi);
    } else {

      b_to_final[bi] = next_pos++;
      n_only_b++;
    }
  }

  uint64_t final_n_samples = aids->n + n_only_b;


  size_t old_aids_n = aids->n;
  tk_ivec_ensure(aids, final_n_samples);

  for (size_t i = 0; i < bids->n; i++) {
    if (b_to_final[i] >= (int64_t)old_aids_n) {
      aids->a[aids->n++] = bids->a[i];
    }
  }


  tk_cvec_ensure(base, final_n_samples * total_bytes_per_sample);
  uint8_t *base_data = (uint8_t *)base->a;
  uint8_t *ext_data = (uint8_t *)ext->a;


  tk_iumap_t *b_id_to_pos = tk_iumap_from_ivec(bids);


  uint8_t *new_data = malloc(final_n_samples * total_bytes_per_sample);
  if (!new_data) {
    free(b_to_final);
    tk_iumap_destroy(a_id_to_pos);
    tk_iumap_destroy(b_id_to_pos);
    return;
  }
  memset(new_data, 0, final_n_samples * total_bytes_per_sample);


  for (size_t ai = 0; ai < old_aids_n; ai++) {
    uint64_t dest_offset = ai * total_bytes_per_sample;
    uint64_t src_offset = ai * base_bytes_per_sample;

    memcpy(new_data + dest_offset, base_data + src_offset, base_bytes_per_sample);


    khint_t khi = tk_iumap_get(b_id_to_pos, aids->a[ai]);
    if (khi != tk_iumap_end(b_id_to_pos)) {
      int64_t b_idx = tk_iumap_value(b_id_to_pos, khi);

      uint64_t ext_src_offset = (uint64_t)b_idx * ext_bytes_per_sample;
      uint64_t ext_bit_offset = n_base_features;
      uint64_t ext_byte_offset = ext_bit_offset / CHAR_BIT;
      uint8_t ext_bit_shift = ext_bit_offset % CHAR_BIT;
      if (ext_bit_shift == 0) {

        memcpy(new_data + dest_offset + ext_byte_offset,
               ext_data + ext_src_offset, ext_bytes_per_sample);
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


  for (size_t bi = 0; bi < bids->n; bi++) {
    int64_t final_pos = b_to_final[bi];
    if (final_pos >= (int64_t)old_aids_n) {

      uint64_t dest_offset = (uint64_t)final_pos * total_bytes_per_sample;
      uint64_t src_offset = bi * ext_bytes_per_sample;


      uint64_t ext_bit_offset = n_base_features;
      uint64_t ext_byte_offset = ext_bit_offset / CHAR_BIT;
      uint8_t ext_bit_shift = ext_bit_offset % CHAR_BIT;
      if (ext_bit_shift == 0) {
        memcpy(new_data + dest_offset + ext_byte_offset,
               ext_data + src_offset, ext_bytes_per_sample);
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
}


static inline tk_dvec_t *tk_cvec_bits_score_chi2 (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  unsigned int n_threads
) {
  tk_cvec_ctx_t ctx;
  tk_cvec_thread_t threads[n_threads];
  tk_threadpool_t *pool = tk_threads_create(L, n_threads, tk_cvec_worker);
  ctx.bitmap = bitmap;
  ctx.codes = codes ? codes->a : NULL;
  ctx.labels = labels;
  ctx.n_features = n_features;
  ctx.n_hidden = n_hidden;
  ctx.n_samples = n_samples;
  ctx.active_counts = tk_ivec_create(L, ctx.n_features * ctx.n_hidden, 0, 0);
  ctx.global_counts = tk_ivec_create(L, ctx.n_hidden, 0, 0);
  ctx.feat_counts = tk_ivec_create(L, ctx.n_features, 0, 0);
  tk_ivec_zero(ctx.active_counts);
  tk_ivec_zero(ctx.global_counts);
  tk_ivec_zero(ctx.feat_counts);

  for (unsigned int i = 0; i < n_threads; i++) {
    tk_cvec_thread_t *data = threads + i;
    pool->threads[i].data = data;
    data->state = &ctx;
    tk_thread_range(i, n_threads, ctx.n_hidden, &data->hfirst, &data->hlast);
  }


  uint8_t *bitmap_data = (uint8_t *)bitmap->a;

  if (ctx.codes != NULL) {

    for (uint64_t s = 0; s < ctx.n_samples; s++) {
      const unsigned char *sample_bitmap =
        (const unsigned char *)(ctx.codes + s * (ctx.n_hidden / CHAR_BIT));
      for (uint64_t chunk = 0; chunk < (ctx.n_hidden / CHAR_BIT); chunk++) {
        unsigned char byte = sample_bitmap[chunk];
        while (byte) {
          int bit = __builtin_ctz(byte);
          uint64_t b = chunk * CHAR_BIT + (unsigned) bit;
          if (b < ctx.n_hidden)
            ctx.global_counts->a[b]++;
          byte &= byte - 1;
        }
      }
    }


    for (uint64_t s = 0; s < n_samples; s++) {
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t bit_idx = s * n_features + f;
        uint64_t byte_idx = bit_idx / CHAR_BIT;
        uint64_t bit_pos = bit_idx % CHAR_BIT;

        if (bitmap_data[byte_idx] & (1 << bit_pos)) {
          ctx.feat_counts->a[f]++;
          const unsigned char *sample_codes =
            (unsigned char *)(ctx.codes + s * (ctx.n_hidden / CHAR_BIT));
          for (uint64_t b = 0; b < ctx.n_hidden; b++) {
            uint64_t chunk = b / CHAR_BIT, bit = b % CHAR_BIT;
            if (sample_codes[chunk] & (1u << bit)) {
              ctx.active_counts->a[f * ctx.n_hidden + b]++;
            }
          }
        }
      }
    }
  } else if (labels != NULL) {

    for (uint64_t s = 0; s < n_samples; s++) {
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t bit_idx = s * n_features + f;
        uint64_t byte_idx = bit_idx / CHAR_BIT;
        uint64_t bit_pos = bit_idx % CHAR_BIT;

        if (bitmap_data[byte_idx] & (1 << bit_pos)) {
          ctx.feat_counts->a[f]++;
          int64_t label = ctx.labels->a[s];
          if (label >= 0 && (size_t) label < ctx.n_hidden)
            ctx.active_counts->a[f * ctx.n_hidden + (uint64_t) label]++;
        }
      }
    }


    tk_ivec_zero(ctx.global_counts);
    for (uint64_t s = 0; s < ctx.n_samples; s++) {
      int64_t label = ctx.labels->a[s];
      if (label >= 0 && (size_t) label < ctx.n_hidden)
        ctx.global_counts->a[label]++;
    }
  }


  ctx.scores = tk_dvec_create(L, ctx.n_hidden * ctx.n_features, 0, 0);
  tk_dvec_zero(ctx.scores);
  tk_threads_signal(pool, TK_CVEC_CHI2, 0);
  tk_threads_destroy(pool);

  tk_ivec_destroy(ctx.active_counts);
  tk_ivec_destroy(ctx.global_counts);
  tk_ivec_destroy(ctx.feat_counts);

  return ctx.scores;
}

static inline tk_dvec_t *tk_cvec_bits_score_mi (
  lua_State *L,
  tk_cvec_t *bitmap,
  tk_cvec_t *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden,
  unsigned int n_threads
) {
  tk_cvec_ctx_t ctx;
  tk_cvec_thread_t threads[n_threads];
  tk_threadpool_t *pool = tk_threads_create(L, n_threads, tk_cvec_worker);
  ctx.bitmap = bitmap;
  ctx.codes = codes ? codes->a : NULL;
  ctx.labels = labels;
  ctx.n_features = n_features;
  ctx.n_hidden = n_hidden;
  ctx.n_samples = n_samples;
  ctx.counts = tk_ivec_create(L, ctx.n_features * ctx.n_hidden * 4, 0, 0);
  tk_ivec_zero(ctx.counts);

  for (unsigned int i = 0; i < n_threads; i++) {
    tk_cvec_thread_t *data = threads + i;
    pool->threads[i].data = data;
    data->state = &ctx;
    tk_thread_range(i, n_threads, ctx.n_hidden, &data->hfirst, &data->hlast);
  }


  uint8_t *bitmap_data = (uint8_t *)bitmap->a;

  if (ctx.codes != NULL) {
    for (uint64_t s = 0; s < n_samples; s++) {
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t bit_idx = s * n_features + f;
        uint64_t byte_idx = bit_idx / CHAR_BIT;
        uint64_t bit_pos = bit_idx % CHAR_BIT;
        bool visible = (bitmap_data[byte_idx] & (1 << bit_pos)) != 0;

        for (uint64_t j = 0; j < ctx.n_hidden; j++) {
          uint64_t chunk = j / CHAR_BIT;
          uint64_t bit = j % CHAR_BIT;
          bool hidden = (ctx.codes[s * (ctx.n_hidden / CHAR_BIT) + chunk] & (1 << bit)) > 0;
          ctx.counts->a[f * ctx.n_hidden * 4 + j * 4 + (visible ? 2 : 0) + (hidden ? 1 : 0)]++;
        }
      }
    }
  } else if (ctx.labels != NULL) {
    for (uint64_t s = 0; s < n_samples; s++) {
      for (uint64_t f = 0; f < n_features; f++) {
        uint64_t bit_idx = s * n_features + f;
        uint64_t byte_idx = bit_idx / CHAR_BIT;
        uint64_t bit_pos = bit_idx % CHAR_BIT;
        bool visible = (bitmap_data[byte_idx] & (1 << bit_pos)) != 0;

        int64_t label = ctx.labels->a[s];
        for (uint64_t j = 0; j < ctx.n_hidden; j++) {
          bool hidden = ((size_t)j == (size_t)label);
          ctx.counts->a[f * ctx.n_hidden * 4 + j * 4 + (visible ? 2 : 0) + (hidden ? 1 : 0)]++;
        }
      }
    }
  }


  ctx.scores = tk_dvec_create(L, ctx.n_hidden * ctx.n_features, 0, 0);
  tk_dvec_zero(ctx.scores);
  tk_threads_signal(pool, TK_CVEC_MI, 0);
  tk_threads_destroy(pool);

  tk_ivec_destroy(ctx.counts);

  return ctx.scores;
}

static inline tk_dvec_t *tk_cvec_bits_score_entropy (
  lua_State *L,
  tk_cvec_t *codes,
  unsigned int n_samples,
  unsigned int n_hidden,
  unsigned int n_threads
) {
  tk_cvec_ctx_t ctx;
  tk_cvec_thread_t threads[n_threads];
  tk_threadpool_t *pool = tk_threads_create(L, n_threads, tk_cvec_worker);
  ctx.codes = codes ? codes->a : NULL;
  ctx.n_samples = n_samples;
  ctx.n_hidden = n_hidden;
  ctx.bit_counts = tk_malloc(L, n_hidden * sizeof(atomic_ulong));
  for (uint64_t i = 0; i < n_hidden; i++)
    atomic_init(ctx.bit_counts + i, 0);

  for (unsigned int i = 0; i < n_threads; i++) {
    tk_cvec_thread_t *data = threads + i;
    pool->threads[i].data = data;
    data->state = &ctx;
    tk_thread_range(i, n_threads, ctx.n_samples, &data->sfirst, &data->slast);
  }


  tk_threads_signal(pool, TK_CVEC_ENTROPY, 0);
  tk_threads_destroy(pool);


  tk_dvec_t *scores = tk_dvec_create(L, ctx.n_hidden, 0, 0);
  for (uint64_t j = 0; j < ctx.n_hidden; j++) {
    double p = (double) ctx.bit_counts[j] / (double) ctx.n_samples;
    double entropy = 0.0;
    if (p > 0.0 && p < 1.0)
      entropy = -(p * log2(p) + (1.0 - p) * log2(1.0 - p));
    scores->a[j] = entropy;
  }

  free(ctx.bit_counts);

  return scores;
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


#endif
