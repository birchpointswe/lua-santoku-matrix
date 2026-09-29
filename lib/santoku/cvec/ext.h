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
#include <santoku/ivec.h>
#include <santoku/dvec.h>
#include <santoku/rvec.h>
#include <string.h>
#include <limits.h>

#ifndef TK_CVEC_BITS_BYTES
#define TK_CVEC_BITS_BYTES(n) (((n) + CHAR_BIT - 1) / CHAR_BIT)
#endif
#ifndef TK_CVEC_BITS_BIT
#define TK_CVEC_BITS_BIT(n) ((n) % CHAR_BIT)
#endif

#define TK_GENERATE_SINGLE
#include <santoku/parallel/tpl.h>
#include <santoku/cvec/ext_tpl.h>
#undef TK_GENERATE_SINGLE

#include <santoku/parallel/tpl.h>
#include <santoku/cvec/ext_tpl.h>

#define tk_cvec_rank_lt(a, b) ((a).d < (b).d || ((a).d == (b).d && (a).i < (b).i))
KSORT_INIT(tk_cvec_rank, tk_rank_t, tk_cvec_rank_lt)

static inline void tk_cvec_bits_topk_row (
  tk_rvec_t *heap,
  uint64_t k,
  const uint8_t *qvec,
  const uint8_t *corpus,
  uint64_t n_corpus,
  uint64_t bytes_per,
  uint64_t n_bits,
  int64_t *idx_row,
  int64_t *dist_row
) {
  tk_rvec_clear(heap);
  for (uint64_t c = 0; c < n_corpus; c++) {
    uint64_t dist = tk_cvec_bits_hamming_serial(qvec, corpus + c * bytes_per, n_bits);
    tk_rank_t r = tk_rank((int64_t)c, (double)dist);
    if (heap->n < k) {
      tk_rvec_push(heap, r);
      if (heap->n == k)
        ks_heapmake(tk_cvec_rank, k, heap->a);
    } else if (tk_cvec_rank_lt(r, heap->a[0])) {
      heap->a[0] = r;
      ks_heapadjust(tk_cvec_rank, 0, k, heap->a);
    }
  }
  ks_heapsort(tk_cvec_rank, heap->n, heap->a);
  for (uint64_t h = 0; h < heap->n; h++) {
    idx_row[h] = heap->a[h].i;
    dist_row[h] = (int64_t)(heap->a[h].d);
  }
}

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
  const uint8_t *qdata = (const uint8_t *)queries->a;
  const uint8_t *cdata = (const uint8_t *)corpus->a;
#if defined(_OPENMP) && !defined(__EMSCRIPTEN__)
  #pragma omp parallel
  {
    tk_rvec_t *heap = tk_rvec_create(NULL, k);
    #pragma omp for schedule(static)
    for (uint64_t q = 0; q < n_queries; q++)
      tk_cvec_bits_topk_row(heap, k, qdata + q * bytes_per, cdata, n_corpus, bytes_per, n_bits,
        indices->a + q * k, distances->a + q * k);
    tk_rvec_destroy(heap);
  }
#else
  tk_rvec_t *heap = tk_rvec_create(NULL, k);
  for (uint64_t q = 0; q < n_queries; q++)
    tk_cvec_bits_topk_row(heap, k, qdata + q * bytes_per, cdata, n_corpus, bytes_per, n_bits,
      indices->a + q * k, distances->a + q * k);
  tk_rvec_destroy(heap);
#endif
}

static inline void tk_cvec_bits_copy_at (
  uint8_t *dst,
  uint64_t off,
  const uint8_t *src,
  uint64_t n_bits
) {
  uint64_t nb = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem = TK_CVEC_BITS_BIT(n_bits);
  uint64_t o = off / CHAR_BIT;
  unsigned int s = off % CHAR_BIT;
  for (uint64_t j = 0; j < nb; j++) {
    uint8_t x = src[j];
    if (rem > 0 && j == nb - 1)
      x &= (uint8_t)((1u << rem) - 1);
    dst[o + j] |= (uint8_t)(x << s);
    if (s > 0 && (x >> (CHAR_BIT - s)) != 0)
      dst[o + j + 1] |= (uint8_t)(x >> (CHAR_BIT - s));
  }
}

#endif
