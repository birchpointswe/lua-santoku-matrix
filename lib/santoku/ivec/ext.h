#ifndef TK_IVEC_EXT_H
#define TK_IVEC_EXT_H

#include <santoku/rvec.h>
#include <santoku/threads.h>
#include <santoku/dvec.h>
#include <stdatomic.h>
#include <math.h>
#include <limits.h>


#ifndef TK_CVEC_BITS_BYTES
#define TK_CVEC_BITS_BYTES(n) (((n) + CHAR_BIT - 1) / CHAR_BIT)
#endif


typedef enum {
  TK_IVEC_CHI2,
  TK_IVEC_MI,
  TK_IVEC_ENTROPY,
} tk_ivec_stage_t;

typedef struct {
  uint64_t n_visible, n_hidden, n_samples;
  tk_dvec_t *scores;
  tk_ivec_t *counts, *feat_counts;
  tk_ivec_t *active_counts, *global_counts;
  tk_ivec_t *set_bits;
  tk_ivec_t *labels;
  char *codes;
  atomic_ulong *bit_counts;
} tk_ivec_ctx_t;

typedef struct {
  tk_ivec_ctx_t *state;
  uint64_t hfirst, hlast;
  uint64_t sfirst, slast;
} tk_ivec_thread_t;

static inline void tk_ivec_worker (void *dp, int sig)
{
  tk_ivec_thread_t *data = (tk_ivec_thread_t *) dp;
  tk_ivec_ctx_t *state = data->state;
  uint64_t n_visible = state->n_visible;
  uint64_t n_hidden = state->n_hidden;
  uint64_t n_samples = state->n_samples;
  tk_dvec_t *scores = state->scores;
  tk_ivec_t *counts = state->counts;
  tk_ivec_t *active_counts = state->active_counts;
  tk_ivec_t *global_counts = state->global_counts;
  tk_ivec_t *feat_counts = state->feat_counts;
  switch ((tk_ivec_stage_t) sig) {

    case TK_IVEC_CHI2:
      for (uint64_t b = data->hfirst; b <= data->hlast; b ++) {
        double *scores_b = scores->a + b * n_visible;
        for (uint64_t f = 0; f < n_visible; f ++) {
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

    case TK_IVEC_MI:
      for (int64_t j = (int64_t) data->hfirst; j <= (int64_t) data->hlast; j ++) {
        double *scores_h = scores->a + j * (int64_t) n_visible;
        for (int64_t i = 0; i < (int64_t) n_visible; i ++) {
          int64_t *c = counts->a + i * (int64_t) n_hidden * 4 + j * 4;
          for (int k = 0; k < 4; k ++)
            c[k] += 1;
          double total = c[0] + c[1] + c[2] + c[3];
          double mi = 0.0;
          if (total == 0.0)
            continue;
          for (unsigned int o = 0; o < 4; o ++) {
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
    case TK_IVEC_ENTROPY:

      for (uint64_t i = 0; i < state->set_bits->n; i++) {
        int64_t bit_idx = state->set_bits->a[i];
        if (bit_idx < 0) continue;

        uint64_t sample_idx = (uint64_t)bit_idx / state->n_hidden;
        uint64_t hidden_idx = (uint64_t)bit_idx % state->n_hidden;


        if (hidden_idx >= data->hfirst && hidden_idx <= data->hlast && sample_idx < state->n_samples) {
          atomic_fetch_add(state->bit_counts + hidden_idx, 1);
        }
      }
      break;
  }
}

static inline void tk_ivec_copy_pkeys (
  tk_ivec_t *m0,
  tk_pvec_t *m1,
  int64_t start,
  int64_t end,
  int64_t dest
) {
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

static inline void tk_ivec_copy_pvalues (
  tk_ivec_t *m0,
  tk_pvec_t *m1,
  int64_t start,
  int64_t end,
  int64_t dest
) {
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


static inline void tk_ivec_copy_rkeys (
  tk_ivec_t *m0,
  tk_rvec_t *m1,
  int64_t start,
  int64_t end,
  int64_t dest
) {
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

static inline void tk_ivec_copy_rvalues (
  tk_ivec_t *m0,
  tk_rvec_t *m1,
  int64_t start,
  int64_t end,
  int64_t dest
) {
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

static inline tk_ivec_t *tk_ivec_from_rvec (
  lua_State *L,
  tk_rvec_t *R
) {
  tk_ivec_t *I = tk_ivec_create(L, R->n, 0, 0);
  for (int64_t i = 0; i < (int64_t) R->n; i ++)
    I->a[i] = R->a[i].i;
  return I;
}



static inline void tk_ivec_lookup (
  tk_ivec_t *indices,
  tk_ivec_t *source
) {
  int64_t write_pos = 0;
  for (uint64_t i = 0; i < indices->n; i++) {
    int64_t idx = indices->a[i];

    if (idx >= 0 && idx < (int64_t) source->n) {

      indices->a[write_pos++] = source->a[idx];
    }

  }

  indices->n = (uint64_t) write_pos;
}


static inline tk_rvec_t *tk_rvec_rankings (
  lua_State *L,
  tk_dvec_t *scores,
  uint64_t n_visible,
  uint64_t n_hidden
) {

  tk_rvec_t *rankings = tk_rvec_create(L, n_hidden * n_visible, NULL, NULL);
  for (uint64_t h = 0; h < n_hidden; h ++)
    for (uint64_t v = 0; v < n_visible; v ++)
      rankings->a[h * n_visible + v] = tk_rank((int64_t) v, scores->a[h * n_visible + v]);


  for (uint64_t j = 0; j < n_hidden; j ++)
    tk_rvec_desc(rankings, j * n_visible, (j + 1) * n_visible);


  return rankings;
}


static inline tk_cvec_t *tk_ivec_bits_to_cvec (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_features,
  bool flip_interleave
) {

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
      if (v < 0)
        continue;

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
      if (v < 0)
        continue;
      uint64_t s = (uint64_t) v / n_features;
      uint64_t f = (uint64_t) v % n_features;
      uint64_t byte_off = s * bytes_per_sample + (f / CHAR_BIT);
      uint8_t bit_off = f & (CHAR_BIT - 1);
      out[byte_off] |= (uint8_t) (1u << bit_off);
    }
  }

  return out_cvec;
}


static inline tk_ivec_t *tk_ivec_bits_from_cvec (
  lua_State *L,
  const char *bm,
  uint64_t n_samples,
  uint64_t n_features
) {
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);
  for (uint64_t i = 0; i < n_samples; i ++)
    for (uint64_t j = 0; j < n_features; j ++) {
      uint64_t bit = i * n_features + j;
      uint64_t chunk = bit / CHAR_BIT;
      uint64_t pos = bit % CHAR_BIT;
      if (bm[chunk] & (1 << pos))
        tk_ivec_push(out, (int64_t) bit);
    }
  tk_ivec_shrink(out);
  return out;
}

static inline void tk_ivec_bits_extend (
  tk_ivec_t *base,
  tk_ivec_t *ext,
  uint64_t n_feat,
  uint64_t n_extfeat
) {
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

static inline tk_dvec_t *tk_ivec_bits_score_chi2 (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  unsigned int n_threads
) {
  tk_ivec_ctx_t ctx;
  tk_ivec_thread_t threads[n_threads];
  tk_threadpool_t *pool = tk_threads_create(L, n_threads, tk_ivec_worker);
  ctx.set_bits = set_bits;
  ctx.codes = codes;
  ctx.labels = labels;
  ctx.n_visible = n_visible;
  ctx.n_hidden = n_hidden;
  ctx.n_samples = n_samples;
  ctx.active_counts = tk_ivec_create(L, ctx.n_visible * ctx.n_hidden, 0, 0);
  ctx.global_counts = tk_ivec_create(L, ctx.n_hidden, 0, 0);
  ctx.feat_counts  = tk_ivec_create(L, ctx.n_visible, 0, 0);
  tk_ivec_zero(ctx.active_counts);
  tk_ivec_zero(ctx.global_counts);
  tk_ivec_zero(ctx.feat_counts);
  for (unsigned int i = 0; i < n_threads; i ++) {
    tk_ivec_thread_t *data = threads + i;
    pool->threads[i].data = data;
    data->state = &ctx;
    tk_thread_range(i, n_threads, ctx.n_hidden, &data->hfirst, &data->hlast);
  }



  if (ctx.codes != NULL) {


    for (uint64_t s = 0; s < ctx.n_samples; s++) {
      const unsigned char *sample_bitmap
        = (const unsigned char *) (ctx.codes + s * (ctx.n_hidden / CHAR_BIT));
      for (uint64_t chunk = 0; chunk < (ctx.n_hidden / CHAR_BIT); chunk ++) {
        unsigned char byte = sample_bitmap[chunk];
        while (byte) {
          int bit = __builtin_ctz(byte);
          uint64_t b = chunk * CHAR_BIT + (unsigned) bit;
          if (b < ctx.n_hidden)
            ctx.global_counts->a[b] ++;
          byte &= byte - 1;
        }
      }
    }
    for (uint64_t i = 0; i < ctx.set_bits->n; i ++) {
      if (ctx.set_bits->a[i] < 0)
        continue;
      uint64_t s = (uint64_t) ctx.set_bits->a[i] / ctx.n_visible;
      uint64_t f = (uint64_t) ctx.set_bits->a[i] % ctx.n_visible;
      ctx.feat_counts->a[f] ++;
      const unsigned char *bitmap = (unsigned char *) (ctx.codes + s * (ctx.n_hidden / CHAR_BIT));
      for (uint64_t b = 0; b < ctx.n_hidden; b ++) {
        uint64_t chunk = b / CHAR_BIT, bit = b % CHAR_BIT;
        if (bitmap[chunk] & (1u << bit)) {
          ctx.active_counts->a[f * ctx.n_hidden + b] ++;
        }
      }
    }
  } else if (labels != NULL) {
    for (uint64_t i = 0; i < ctx.set_bits->n; i ++) {
      if (ctx.set_bits->a[i] < 0)
        continue;
      uint64_t s = (uint64_t) ctx.set_bits->a[i] / ctx.n_visible;
      uint64_t f = (uint64_t) ctx.set_bits->a[i] % ctx.n_visible;
      ctx.feat_counts->a[f] ++;
      int64_t label = ctx.labels->a[s];
      if (label < 0)
        continue;
      if ((size_t) label < ctx.n_hidden)
        ctx.active_counts->a[f * ctx.n_hidden + (uint64_t) label] ++;
    }
    tk_ivec_zero(ctx.global_counts);
    for (uint64_t s = 0; s < ctx.n_samples; s ++) {
      int64_t label = ctx.labels->a[s];
      if ((size_t) label < ctx.n_hidden)
        ctx.global_counts->a[label] ++;
    }
  }


  ctx.scores = tk_dvec_create(L, ctx.n_hidden * ctx.n_visible, 0, 0);
  tk_dvec_zero(ctx.scores);
  tk_threads_signal(pool, TK_IVEC_CHI2, 0);
  tk_threads_destroy(pool);

  tk_ivec_destroy(ctx.active_counts);
  tk_ivec_destroy(ctx.global_counts);
  tk_ivec_destroy(ctx.feat_counts);

  return ctx.scores;
}

static inline tk_dvec_t *tk_ivec_bits_score_mi (
  lua_State *L,
  tk_ivec_t *set_bits,
  char *codes,
  tk_ivec_t *labels,
  uint64_t n_samples,
  uint64_t n_visible,
  uint64_t n_hidden,
  unsigned int n_threads
) {
  tk_ivec_ctx_t ctx;
  tk_ivec_thread_t threads[n_threads];
  tk_threadpool_t *pool = tk_threads_create(L, n_threads, tk_ivec_worker);
  ctx.set_bits = set_bits;
  ctx.codes = codes;
  ctx.labels = labels;
  ctx.n_visible = n_visible;
  ctx.n_hidden = n_hidden;
  ctx.n_samples = n_samples;
  ctx.counts = tk_ivec_create(L, ctx.n_visible * ctx.n_hidden * 4, 0, 0);
  tk_ivec_zero(ctx.counts);
  for (unsigned int i = 0; i < n_threads; i ++) {
    tk_ivec_thread_t *data = threads + i;
    pool->threads[i].data = data;
    data->state = &ctx;
    tk_thread_range(i, n_threads, ctx.n_hidden, &data->hfirst, &data->hlast);
  }

  if (ctx.codes != NULL) {
    int64_t last_set_bit = -1;
    for (uint64_t i = 0; i < ctx.set_bits->n; i ++) {
      int64_t set_bit = ctx.set_bits->a[i];
      if (set_bit < 0)
        continue;
      for (int64_t unset_bit = last_set_bit + 1; unset_bit < set_bit; unset_bit ++) {
        uint64_t sample = (uint64_t) unset_bit / ctx.n_visible;
        uint64_t visible = (uint64_t) unset_bit % ctx.n_visible;
        for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
          uint64_t chunk = j / CHAR_BIT;
          uint64_t bit = j % CHAR_BIT;
          bool hidden = (ctx.codes[sample * (ctx.n_hidden / CHAR_BIT) + chunk] & (1 << bit)) > 0;
          ctx.counts->a[visible * ctx.n_hidden * 4 + j * 4 + (hidden ? 1 : 0)] ++;
        }
      }
      uint64_t sample = (uint64_t) set_bit / ctx.n_visible;
      uint64_t visible = (uint64_t) set_bit % ctx.n_visible;
      for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
        uint64_t chunk = j / CHAR_BIT;
        uint64_t bit = j % CHAR_BIT;
        bool hidden = (ctx.codes[sample * (ctx.n_hidden / CHAR_BIT) + chunk] & (1 << bit)) > 0;
        ctx.counts->a[visible * ctx.n_hidden * 4 + j * 4 + 2 + (hidden ? 1 : 0)] ++;
      }
      last_set_bit = set_bit;
    }
    for (uint64_t unset_bit = (uint64_t) last_set_bit + 1; unset_bit < n_samples * n_visible; unset_bit ++) {
      uint64_t sample = unset_bit / ctx.n_visible;
      uint64_t visible = unset_bit % ctx.n_visible;
      for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
        uint64_t chunk = j / CHAR_BIT;
        uint64_t bit = j % CHAR_BIT;
        bool hidden = (ctx.codes[sample * (ctx.n_hidden / CHAR_BIT) + chunk] & (1 << bit)) > 0;
        ctx.counts->a[visible * ctx.n_hidden * 4 + j * 4 + (hidden ? 1 : 0)] ++;
      }
    }
  } else if (ctx.labels != NULL) {
    int64_t last_set_bit = -1;
    for (uint64_t i = 0; i < set_bits->n; i ++) {
      int64_t set_bit = ctx.set_bits->a[i];
      if (set_bit < 0) continue;
      for (uint64_t unset_bit = (uint64_t) last_set_bit + 1; unset_bit < (uint64_t) set_bit; unset_bit ++) {
        uint64_t sample = unset_bit / ctx.n_visible;
        uint64_t visible = unset_bit % ctx.n_visible;
        int64_t label = ctx.labels->a[sample];
        for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
          bool hidden = ((size_t)j == (size_t)label);
          ctx.counts->a[visible * ctx.n_hidden * 4 + j * 4 + (hidden ? 1 : 0)] ++;
        }
      }
      uint64_t sample = (uint64_t) set_bit / ctx.n_visible;
      uint64_t visible = (uint64_t) set_bit % ctx.n_visible;
      int64_t label = ctx.labels->a[sample];
      for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
        bool hidden = ((size_t)j == (size_t)label);
        ctx.counts->a[visible * ctx.n_hidden * 4 + j * 4 + 2 + (hidden ? 1 : 0)] ++;
      }
      last_set_bit = set_bit;
    }
    for (uint64_t unset_bit = (uint64_t) last_set_bit + 1; unset_bit < ctx.n_samples * ctx.n_visible; unset_bit ++) {
      uint64_t sample = unset_bit / ctx.n_visible;
      uint64_t visible = unset_bit % ctx.n_visible;
      int64_t label = ctx.labels->a[sample];
      for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
        bool hidden = ((size_t)j == (size_t)label);
        ctx.counts->a[visible * ctx.n_hidden * 4 + j * 4 + (hidden ? 1 : 0)] ++;
      }
    }
  }


  ctx.scores = tk_dvec_create(L, ctx.n_hidden * ctx.n_visible, 0, 0);
  tk_dvec_zero(ctx.scores);
  tk_threads_signal(pool, TK_IVEC_MI, 0);
  tk_threads_destroy(pool);


  tk_ivec_destroy(ctx.counts);

  return ctx.scores;
}



typedef enum {
  TK_IVEC_JACCARD,
  TK_IVEC_OVERLAP,
  TK_IVEC_DICE,
  TK_IVEC_TVERSKY
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

static inline int64_t tk_ivec_set_find (
  tk_ivec_t *vec,
  int64_t value,
  int64_t *insert_idx
) {
  int64_t left = 0;
  int64_t right = (int64_t)vec->n - 1;
  int64_t mid;
  while (left <= right) {
    mid = left + (right - left) / 2;
    if (vec->a[mid] == value)
      return mid;
    if (vec->a[mid] < value)
      left = mid + 1;
    else
      right = mid - 1;
  }
  if (insert_idx != NULL)
    *insert_idx = left;
  return -1;
}

static inline tk_ivec_t *tk_ivec_set_intersect (
  lua_State *L,
  tk_ivec_t *a,
  tk_ivec_t *b,
  tk_ivec_t *out
) {
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

static inline tk_ivec_t *tk_ivec_set_union (
  lua_State *L,
  tk_ivec_t *a,
  tk_ivec_t *b,
  tk_ivec_t *out
) {

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

static inline tk_dvec_t *tk_ivec_bits_score_entropy (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_hidden,
  unsigned int n_threads
) {
  tk_ivec_ctx_t ctx;
  tk_ivec_thread_t threads[n_threads];
  tk_threadpool_t *pool = tk_threads_create(L, n_threads, tk_ivec_worker);

  ctx.set_bits = set_bits;
  ctx.n_samples = n_samples;
  ctx.n_hidden = n_hidden;
  ctx.bit_counts = tk_malloc(L, n_hidden * sizeof(atomic_ulong));


  for (uint64_t i = 0; i < n_hidden; i++)
    atomic_init(ctx.bit_counts + i, 0);


  for (unsigned int i = 0; i < n_threads; i++) {
    tk_ivec_thread_t *data = threads + i;
    pool->threads[i].data = data;
    data->state = &ctx;
    tk_thread_range(i, n_threads, ctx.n_hidden, &data->hfirst, &data->hlast);
  }


  tk_threads_signal(pool, TK_IVEC_ENTROPY, 0);
  tk_threads_destroy(pool);


  tk_dvec_t *scores = tk_dvec_create(L, n_hidden, 0, 0);
  for (uint64_t j = 0; j < n_hidden; j++) {
    double p = (double) ctx.bit_counts[j] / (double) n_samples;
    double entropy = 0.0;
    if (p > 0.0 && p < 1.0)
      entropy = -(p * log2(p) + (1.0 - p) * log2(1.0 - p));
    scores->a[j] = entropy;
  }

  free(ctx.bit_counts);
  return scores;
}

#endif
