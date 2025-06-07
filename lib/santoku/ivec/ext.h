#ifndef TK_IVEC_EXT_H
#define TK_IVEC_EXT_H

#include <santoku/rvec.h>
#include <santoku/iuset.h>
#include <santoku/threads.h>
#include <stdatomic.h>

typedef enum {
  TK_IVEC_ENTROPY,
  TK_IVEC_CHI2,
  TK_IVEC_MI,
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
  atomic_ulong *bit_counts = state->bit_counts;
  char *codes = state->codes;
  uint64_t chunks = (n_hidden + CHAR_BIT - 1) / CHAR_BIT;
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
      for (uint64_t i = data->sfirst; i <= data->slast; i ++) {
        for (uint64_t j = 0; j < n_hidden; j ++) {
          uint64_t word = j / CHAR_BIT;
          uint64_t bit = j % CHAR_BIT;
          if (codes[i * chunks + word] & (1 << bit))
            atomic_fetch_add(bit_counts + j, 1);
        }
      }
      break;

  }
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

static inline tk_iuset_t *tk_iuset_from_ivec (tk_ivec_t *v)
{
  int kha;
  tk_iuset_t *s = tk_iuset_create();
  for (uint64_t i = 0; i < v->n; i ++)
    tk_iuset_put(s, v->a[i], &kha);
  return s;
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

static inline int tk_ivec_flip_interleave (
  lua_State *L,
  tk_ivec_t *m0,
  uint64_t n_samples,
  uint64_t n_features
) {
  tk_ivec_asc(m0, 0, m0->n);
  size_t total = n_samples * n_features;
  tk_ivec_ensure(L, m0, 2 * total);
  size_t write = m0->n;
  size_t last = 0;
  int64_t s, k;
  size_t n = m0->n;
  for (size_t i_present = 0; i_present < n; i_present ++) {
    int64_t x = m0->a[i_present];
    s = x / (int64_t) n_features;
    k = x % (int64_t) n_features;
    m0->a[i_present] = (s * 2 * (int64_t) n_features) + k;
    for (size_t y = last; y < (size_t) x; y ++) {
      int64_t s = (int64_t) y / (int64_t) n_features;
      int64_t k = (int64_t) y % (int64_t) n_features;
      m0->a[write ++] = (s * 2 * (int64_t) n_features) + (int64_t) n_features + k;
    }
    last = (size_t) x + 1;
  }
  for (size_t y = last; y < total; y ++) {
    int64_t s = (int64_t) y / (int64_t) n_features;
    int64_t k = (int64_t) y % (int64_t) n_features;
    m0->a[write ++] = (s * 2 * (int64_t) n_features) + (int64_t) n_features + k;
  }
  m0->n = write;
  tk_ivec_asc(m0, 0, m0->n);
  return 0;
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
  tk_ivec_shrink(L, top_v);
  tk_iuset_destroy(selected);
  return top_v;
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
    for (uint64_t j = 0; j < n_hidden; j++) {
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

static inline int tk_ivec_filter (
  lua_State *L,
  tk_ivec_t *set_bits,
  tk_ivec_t *top_v,
  uint64_t n_visible
) {
  int64_t *vmap = tk_malloc(L, n_visible * sizeof(int64_t));
  for (unsigned int i = 0; i < n_visible; i ++)
    vmap[i] = -1;
  for (unsigned int i = 0; i < top_v->n; i ++)
    vmap[top_v->a[i]] = i;
  size_t write = 0;
  for (size_t i = 0; i < set_bits->n; i ++) {
    int64_t val = set_bits->a[i];
    if (val < 0)
      continue;
    uint64_t sample = (uint64_t) val / n_visible;
    uint64_t feature = (uint64_t) val % n_visible;
    int64_t new_feature = vmap[feature];
    if (new_feature == -1)
      continue;
    set_bits->a[write ++] = (int64_t) sample * (int64_t) top_v->n + new_feature;
  }
  set_bits->n = write;
  free(vmap);
  return 0;
}

static inline char *tk_ivec_raw_bitmap (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_features,
  size_t *lenp
) {
  size_t len = (n_samples * n_features + CHAR_BIT - 1) / CHAR_BIT;
  char *out = tk_malloc(L, len);
  memset(out, 0, len);
  for (uint64_t i = 0; i < set_bits->n; i ++) {
    int64_t v = set_bits->a[i];
    if (v < 0 || (uint64_t) v >= n_samples * n_features)
      continue;
    out[v / CHAR_BIT] |= (1 << (v % CHAR_BIT));
  }
  *lenp = len;
  return out;
}


static inline tk_ivec_t *tk_ivec_from_bitmap (
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
  tk_ivec_shrink(L, out);
  return out;
}


static inline void tk_ivec_extend_bits (
  lua_State *L,
  tk_ivec_t *base,
  tk_ivec_t *ext,
  uint64_t n_feat,
  uint64_t n_extfeat
) {
  tk_ivec_asc(base, 0, base->n);
  tk_ivec_asc(ext, 0, ext->n);
  size_t total = base->n + ext->n;
  tk_ivec_ensure(L, base, total);
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
    base->a[base->n + i] = (int64_t) (bit - old_off + new_off);
  }
  base->n = total;
  tk_ivec_asc(base, 0, base->n);
}

static inline tk_dvec_t *tk_ivec_score_chi2 (
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



  if (ctx.codes != NULL) {
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
  tk_threads_signal(pool, TK_IVEC_CHI2);
  tk_threads_destroy(pool);

  tk_ivec_destroy(ctx.active_counts);
  tk_ivec_destroy(ctx.global_counts);
  tk_ivec_destroy(ctx.feat_counts);

  return ctx.scores;
}

static inline tk_dvec_t *tk_ivec_score_mi (
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
  tk_threads_signal(pool, TK_IVEC_MI);
  tk_threads_destroy(pool);


  tk_ivec_destroy(ctx.counts);

  return ctx.scores;
}

static inline tk_dvec_t *tk_ivec_score_entropy (
  lua_State *L,
  char *codes,
  unsigned int n_samples,
  unsigned int n_hidden,
  unsigned int n_threads
) {
  tk_ivec_ctx_t ctx;
  tk_ivec_thread_t threads[n_threads];
  tk_threadpool_t *pool = tk_threads_create(L, n_threads, tk_ivec_worker);
  ctx.codes = codes;
  ctx.n_samples = n_samples;
  ctx.n_hidden = n_hidden;
  ctx.bit_counts = tk_malloc(L, n_hidden * sizeof(atomic_ulong));
  for (uint64_t i = 0; i < n_hidden; i ++)
    atomic_init(ctx.bit_counts + i, 0);
  for (unsigned int i = 0; i < n_threads; i ++) {
    tk_ivec_thread_t *data = threads + i;
    pool->threads[i].data = data;
    data->state = &ctx;
    tk_thread_range(i, n_threads, ctx.n_samples, &data->sfirst, &data->slast);
  }


  tk_threads_signal(pool, TK_IVEC_ENTROPY);
  tk_threads_destroy(pool);



  tk_dvec_t *scores = tk_dvec_create(L, ctx.n_hidden, 0, 0);
  for (uint64_t j = 0; j < ctx.n_hidden; j ++) {
    double p = (double) ctx.bit_counts[j] / (double) ctx.n_samples;
    double entropy = 0.0;
    if (p > 0.0 && p < 1.0)
      entropy = -(p * log2(p) + (1.0 - p) * log2(1.0 - p));
    scores->a[j] = entropy;
  }

  free(ctx.bit_counts);

  return scores;
}

static inline tk_ivec_t *tk_ivec_top_mi (
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
  tk_dvec_t *scores = tk_ivec_score_mi(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, n_threads);
  tk_ivec_t *out = tk_ivec_top_generic(L, scores, n_visible, n_hidden, top_k, 0);
  tk_dvec_destroy(scores);
  return out;
}

static inline tk_ivec_t *tk_ivec_top_entropy (
  lua_State *L,
  char *codes,
  uint64_t n_samples,
  uint64_t n_hidden,
  uint64_t top_k,
  unsigned int n_threads
) {
  tk_dvec_t *scores = tk_ivec_score_entropy(L, codes, n_samples, n_hidden, n_threads);
  tk_rvec_t *rankings = tk_rvec_from_dvec(L, scores);
  lua_remove(L, -2);
  tk_rvec_kdesc(rankings, top_k, 0, rankings->n);
  tk_ivec_t *out = tk_ivec_from_rvec(L, rankings);
  lua_remove(L, -2);
  return out;
}

static inline tk_ivec_t *tk_ivec_top_chi2 (
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
  tk_dvec_t *scores = tk_ivec_score_chi2(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, n_threads);
  tk_ivec_t *out = tk_ivec_top_generic(L, scores, n_visible, n_hidden, top_k, 0);
  tk_dvec_destroy(scores);
  return out;
}

#endif
