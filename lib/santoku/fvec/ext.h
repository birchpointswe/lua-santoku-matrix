#ifndef TK_FVEC_EXT_H
#define TK_FVEC_EXT_H

#if !defined(__EMSCRIPTEN__)
#include <lapacke.h>
#include <cblas.h>
#endif
#if defined(_OPENMP) && !defined(__EMSCRIPTEN__)
#include <omp.h>
#endif
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <santoku/fvec/base.h>
#include <santoku/ivec.h>
#include <santoku/rvec/base.h>

KSORT_INIT_GENERIC(float)

static inline void tk_fvec_mtx_topk (lua_State *L, tk_fvec_t *queries, tk_fvec_t *corpus, uint64_t n_queries, uint64_t n_corpus, uint64_t d, uint64_t k);

#include <santoku/parallel/tpl.h>
#include <santoku/fvec/ext_tpl.h>

#if !defined(__EMSCRIPTEN__)
static inline void tk_fvec_gemv(
  bool transpose, uint64_t rows, uint64_t cols,
  float alpha, float *A, float *x, float beta, float *y
) {
  cblas_sgemv(CblasRowMajor, transpose ? CblasTrans : CblasNoTrans, rows, cols, alpha, A, cols, x, 1, beta, y, 1);
}

static inline void tk_fvec_gemm(
  bool transpose_a, bool transpose_b,
  uint64_t m, uint64_t n, uint64_t k,
  float alpha, float *A, float *B, float beta, float *C
) {
  cblas_sgemm(CblasRowMajor, transpose_a ? CblasTrans : CblasNoTrans, transpose_b ? CblasTrans : CblasNoTrans, m, n, k, alpha, A, transpose_a ? m : k, B, transpose_b ? k : n, beta, C, n);
}

#else

static inline void tk_fvec_gemv(
  bool transpose, uint64_t rows, uint64_t cols,
  float alpha, float *A, float *x, float beta, float *y
) {
  uint64_t out_len = transpose ? cols : rows;
  for (uint64_t i = 0; i < out_len; i++) y[i] = beta == 0.0f ? 0.0f : y[i] * beta;
  if (!transpose) {
    for (uint64_t r = 0; r < rows; r++)
      for (uint64_t c = 0; c < cols; c++)
        y[r] += alpha * A[r * cols + c] * x[c];
  } else {
    for (uint64_t r = 0; r < rows; r++)
      for (uint64_t c = 0; c < cols; c++)
        y[c] += alpha * A[r * cols + c] * x[r];
  }
}

static inline void tk_fvec_gemm(
  bool transpose_a, bool transpose_b,
  uint64_t m, uint64_t n, uint64_t k,
  float alpha, float *A, float *B, float beta, float *C
) {
  for (uint64_t i = 0; i < m * n; i++) C[i] = beta == 0.0f ? 0.0f : C[i] * beta;
  for (uint64_t i = 0; i < m; i++)
    for (uint64_t j = 0; j < n; j++)
      for (uint64_t l = 0; l < k; l++)
        C[i * n + j] += alpha *
          (transpose_a ? A[l * m + i] : A[i * k + l]) *
          (transpose_b ? B[j * k + l] : B[l * n + j]);
}

#endif

#include <santoku/cvec/base.h>
#include <limits.h>
#ifndef TK_CVEC_BITS_BYTES
#define TK_CVEC_BITS_BYTES(n) (((n) + CHAR_BIT - 1) / CHAR_BIT)
#endif

static inline void tk_fvec_mtx_center (
  lua_State *L, tk_fvec_t *data, uint64_t n_cols,
  tk_fvec_t *mean_in, tk_fvec_t **mean_out
) {
  uint64_t N = data->n / n_cols;
  if (mean_in) {
    #pragma omp parallel for
    for (uint64_t d = 0; d < n_cols; d++) {
      float mu = mean_in->a[d];
      for (uint64_t s = 0; s < N; s++)
        data->a[s * n_cols + d] -= mu;
    }
  } else {
    tk_fvec_t *mu = tk_fvec_create(L, n_cols);
    mu->n = n_cols;
    #pragma omp parallel for
    for (uint64_t d = 0; d < n_cols; d++) {
      double sum = 0;
      for (uint64_t s = 0; s < N; s++)
        sum += (double)data->a[s * n_cols + d];
      float m = (float)(sum / (double)N);
      mu->a[d] = m;
      for (uint64_t s = 0; s < N; s++)
        data->a[s * n_cols + d] -= m;
    }
    *mean_out = mu;
  }
}

static inline tk_fvec_t *tk_fvec_mtx_standardize (
  lua_State *L, tk_fvec_t *data, uint64_t N, uint64_t n_cols, bool rms
) {
  tk_fvec_t *w = tk_fvec_create(L, n_cols);
  w->n = n_cols;
  #pragma omp parallel for
  for (uint64_t d = 0; d < n_cols; d++) {
    double sum = 0, sum2 = 0;
    for (uint64_t s = 0; s < N; s++) {
      double v = (double)data->a[s * n_cols + d];
      sum += v; sum2 += v * v;
    }
    double x;
    if (rms) {
      x = sum2 > 0.0 ? sqrt((double)N / sum2) : 0.0;
    } else {
      double m = sum / (double)N;
      double var = sum2 / (double)N - m * m;
      x = var > 1e-24 ? 1.0 / sqrt(var) : 0.0;
    }
    w->a[d] = (float)x;
    for (uint64_t s = 0; s < N; s++)
      data->a[s * n_cols + d] = (float)((double)data->a[s * n_cols + d] * x);
  }
  return w;
}

static inline void tk_fvec_mtx_sign_raw (
  char *out, float *X, uint64_t N, uint64_t stride, uint64_t K
) {
  #pragma omp parallel for
  for (uint64_t i = 0; i < N; i++) {
    float *row = X + i * stride;
    uint8_t *out_row = (uint8_t *)(out + i * TK_CVEC_BITS_BYTES(K));
    uint64_t full_bytes = K / 8;
    for (uint64_t byte_idx = 0; byte_idx < full_bytes; byte_idx++) {
      uint8_t byte_val = 0;
      uint64_t j_base = byte_idx * 8;
      for (uint64_t bit = 0; bit < 8; bit++)
        byte_val |= (uint8_t)((row[j_base + bit] >= 0.0f) << bit);
      out_row[byte_idx] = byte_val;
    }
    uint64_t remaining_start = full_bytes * 8;
    if (remaining_start < K) {
      uint8_t byte_val = 0;
      for (uint64_t j = remaining_start; j < K; j++)
        byte_val |= (uint8_t)((row[j] >= 0.0f) << (j - remaining_start));
      out_row[full_bytes] = byte_val;
    }
  }
}

static inline tk_cvec_t *tk_fvec_mtx_sign (
  lua_State *L, tk_fvec_t *codes, uint64_t n_dims, uint64_t n_trunc
) {
  const size_t N = codes->n / n_dims;
  tk_cvec_t *binary = tk_cvec_create(L, N * TK_CVEC_BITS_BYTES(n_trunc));
  binary->n = N * TK_CVEC_BITS_BYTES(n_trunc);
  tk_cvec_zero(binary);
  tk_fvec_mtx_sign_raw(binary->a, codes->a, N, n_dims, n_trunc);
  return binary;
}

static inline void tk_fvec_mtx_threshold_raw (
  char *out, float *X, float *thresholds, uint64_t N, uint64_t K
) {
  #pragma omp parallel for
  for (uint64_t i = 0; i < N; i++) {
    float *row = X + i * K;
    uint8_t *out_row = (uint8_t *)(out + i * TK_CVEC_BITS_BYTES(K));
    uint64_t full_bytes = K / 8;
    for (uint64_t byte_idx = 0; byte_idx < full_bytes; byte_idx++) {
      uint8_t byte_val = 0;
      uint64_t j_base = byte_idx * 8;
      for (uint64_t bit = 0; bit < 8; bit++)
        byte_val |= (uint8_t)((row[j_base + bit] >= thresholds[j_base + bit]) << bit);
      out_row[byte_idx] = byte_val;
    }
    uint64_t remaining_start = full_bytes * 8;
    if (remaining_start < K) {
      uint8_t byte_val = 0;
      for (uint64_t j = remaining_start; j < K; j++)
        byte_val |= (uint8_t)((row[j] >= thresholds[j]) << (j - remaining_start));
      out_row[full_bytes] = byte_val;
    }
  }
}

static inline tk_cvec_t *tk_fvec_mtx_median (
  lua_State *L, tk_fvec_t *codes, uint64_t n_dims, tk_fvec_t **medians_out
) {
  const uint64_t K = n_dims;
  const size_t N = codes->n / K;
  tk_cvec_t *binary = tk_cvec_create(L, N * TK_CVEC_BITS_BYTES(K));
  binary->n = N * TK_CVEC_BITS_BYTES(K);
  tk_cvec_zero(binary);
  float *medians = malloc(K * sizeof(float));
  float *col_buf = malloc(N * sizeof(float));
  for (uint64_t k = 0; k < K; k++) {
    for (uint64_t i = 0; i < N; i++) col_buf[i] = codes->a[i * K + k];
    ks_introsort(float, N, col_buf);
    medians[k] = (N % 2 == 1) ? col_buf[N / 2] : (col_buf[N / 2 - 1] + col_buf[N / 2]) / 2.0f;
  }
  free(col_buf);
  tk_fvec_mtx_threshold_raw(binary->a, codes->a, medians, N, K);
  if (medians_out) {
    tk_fvec_t *med_vec = tk_fvec_create(L, K);
    med_vec->n = K;
    memcpy(med_vec->a, medians, K * sizeof(float));
    *medians_out = med_vec;
  }
  free(medians);
  return binary;
}

static inline void tk_fvec_round (tk_fvec_t *v, uint64_t start, uint64_t end) {
  if (end > v->n) end = v->n;
  for (uint64_t i = start; i < end; i++) v->a[i] = roundf(v->a[i]);
}

static inline void tk_fvec_trunc (tk_fvec_t *v, uint64_t start, uint64_t end) {
  if (end > v->n) end = v->n;
  for (uint64_t i = start; i < end; i++) v->a[i] = truncf(v->a[i]);
}

static inline void tk_fvec_floor (tk_fvec_t *v, uint64_t start, uint64_t end) {
  if (end > v->n) end = v->n;
  for (uint64_t i = start; i < end; i++) v->a[i] = floorf(v->a[i]);
}

static inline void tk_fvec_ceil (tk_fvec_t *v, uint64_t start, uint64_t end) {
  if (end > v->n) end = v->n;
  for (uint64_t i = start; i < end; i++) v->a[i] = ceilf(v->a[i]);
}

static inline tk_ivec_t *tk_fvec_to_ivec (lua_State *L, tk_fvec_t *v) {
  tk_ivec_t *out = tk_ivec_create(L, v->n);
  for (uint64_t i = 0; i < v->n; i++) out->a[i] = (int64_t)v->a[i];
  return out;
}

static inline tk_dvec_t *tk_fvec_to_dvec (lua_State *L, tk_fvec_t *v, tk_dvec_t *out) {
  if (out == NULL) {
    out = tk_dvec_create(L, v->n);
  } else {
    tk_dvec_ensure(out, v->n);
    out->n = v->n;
  }
  for (uint64_t i = 0; i < v->n; i++) out->a[i] = (double)v->a[i];
  return out;
}

static inline void tk_fvec_mtx_topk (
  lua_State *L, tk_fvec_t *queries, tk_fvec_t *corpus,
  uint64_t n_queries, uint64_t n_corpus, uint64_t d, uint64_t k
) {
  if (k == 0 || n_queries == 0 || n_corpus == 0 || d == 0) {
    tk_ivec_t *off = tk_ivec_create(L, n_queries + 1);
    off->n = n_queries + 1;
    memset(off->a, 0, (n_queries + 1) * sizeof(int64_t));
    tk_ivec_create(L, 0);
    tk_dvec_create(L, 0);
    return;
  }
  if (k > n_corpus) k = n_corpus;
  uint64_t total = n_queries * k;
  tk_ivec_t *offsets = tk_ivec_create(L, n_queries + 1);
  offsets->n = n_queries + 1;
  for (uint64_t i = 0; i <= n_queries; i++) offsets->a[i] = (int64_t)(i * k);
  tk_ivec_t *indices = tk_ivec_create(L, total);
  indices->n = total;
  tk_dvec_t *out_scores = tk_dvec_create(L, total);
  out_scores->n = total;
  uint64_t max_buf = 128ULL * 1024 * 1024 / sizeof(float);
  uint64_t tile = n_corpus > 0 ? max_buf / n_corpus : n_queries;
  if (tile == 0) tile = 1;
  if (tile > n_queries) tile = n_queries;
  float *sbuf = (float *)malloc(tile * n_corpus * sizeof(float));
  if (!sbuf) { luaL_error(L, "mtx_topk: out of memory"); return; }
  for (uint64_t base = 0; base < n_queries; base += tile) {
    uint64_t blk = (base + tile <= n_queries) ? tile : n_queries - base;
#if !defined(__EMSCRIPTEN__)
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
      (int)blk, (int)n_corpus, (int)d,
      1.0f, queries->a + base * d, (int)d,
      corpus->a, (int)d,
      0.0f, sbuf, (int)n_corpus);
#else
    for (uint64_t i = 0; i < blk; i++)
      for (uint64_t j = 0; j < n_corpus; j++) {
        float s = 0.0f;
        const float *q = queries->a + (base + i) * d;
        const float *c = corpus->a + j * d;
        for (uint64_t l = 0; l < d; l++) s += q[l] * c[l];
        sbuf[i * n_corpus + j] = s;
      }
#endif
    #pragma omp parallel
    {
      tk_rvec_t *heap = tk_rvec_create(NULL, k);
      #pragma omp for schedule(static)
      for (uint64_t i = 0; i < blk; i++) {
        tk_rvec_clear(heap);
        float *row = sbuf + i * n_corpus;
        for (uint64_t j = 0; j < n_corpus; j++)
          tk_rvec_hmin(heap, k, tk_rank((int64_t)j, (double)row[j]));
        tk_rvec_desc(heap, 0, heap->n);
        uint64_t qi = base + i;
        int64_t *idx_row = indices->a + qi * k;
        double *sco_row = out_scores->a + qi * k;
        for (uint64_t h = 0; h < heap->n; h++) {
          idx_row[h] = heap->a[h].i;
          sco_row[h] = heap->a[h].d;
        }
      }
      tk_rvec_destroy(heap);
    }
  }
  free(sbuf);
}

#endif
