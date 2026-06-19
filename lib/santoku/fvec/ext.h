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

#define TK_GENERATE_SINGLE
#include <santoku/parallel/tpl.h>
#include <santoku/fvec/ext_tpl.h>
#undef TK_GENERATE_SINGLE

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

static inline float tk_fvec_blas_dot(float *x, float *y, uint64_t n) {
  return cblas_sdot(n, x, 1, y, 1);
}

static inline void tk_fvec_blas_scal(float alpha, float *x, uint64_t n) {
  cblas_sscal(n, alpha, x, 1);
}

static inline void tk_fvec_blas_axpy(float alpha, float *x, float *y, uint64_t n) {
  cblas_saxpy(n, alpha, x, 1, y, 1);
}

static inline void tk_fvec_blas_copy(float *x, float *y, uint64_t n) {
  cblas_scopy(n, x, 1, y, 1);
}

static inline float tk_fvec_blas_nrm2(float *x, uint64_t n) {
  return cblas_snrm2(n, x, 1);
}

#else

static inline void tk_fvec_gemv(
  bool transpose, uint64_t rows, uint64_t cols,
  float alpha, float *A, float *x, float beta, float *y
) {
  uint64_t out_len = transpose ? cols : rows;
  for (uint64_t i = 0; i < out_len; i++) y[i] *= beta;
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
  for (uint64_t i = 0; i < m * n; i++) C[i] *= beta;
  for (uint64_t i = 0; i < m; i++)
    for (uint64_t j = 0; j < n; j++)
      for (uint64_t l = 0; l < k; l++)
        C[i * n + j] += alpha *
          (transpose_a ? A[l * m + i] : A[i * k + l]) *
          (transpose_b ? B[j * k + l] : B[l * n + j]);
}

static inline float tk_fvec_blas_dot(float *x, float *y, uint64_t n) {
  float s = 0; for (uint64_t i = 0; i < n; i++) s += x[i] * y[i]; return s;
}

static inline void tk_fvec_blas_scal(float alpha, float *x, uint64_t n) {
  for (uint64_t i = 0; i < n; i++) x[i] *= alpha;
}

static inline void tk_fvec_blas_axpy(float alpha, float *x, float *y, uint64_t n) {
  for (uint64_t i = 0; i < n; i++) y[i] += alpha * x[i];
}

static inline void tk_fvec_blas_copy(float *x, float *y, uint64_t n) {
  memcpy(y, x, n * sizeof(float));
}

static inline float tk_fvec_blas_nrm2(float *x, uint64_t n) {
  float s = 0; for (uint64_t i = 0; i < n; i++) s += x[i] * x[i]; return sqrtf(s);
}

#endif

static inline float tk_fvec_dot_override(tk_fvec_t *a, tk_fvec_t *b) {
  uint64_t n = a->n < b->n ? a->n : b->n;
  return tk_fvec_blas_dot(a->a, b->a, n);
}

static inline void tk_fvec_scale_override(tk_fvec_t *v, float scale, uint64_t start, uint64_t end) {
  if (end > v->n) { tk_fvec_ensure(v, end); v->n = end; }
  if (end <= start) return;
  tk_fvec_blas_scal(scale, v->a + start, end - start);
}

static inline void tk_fvec_addv_override(tk_fvec_t *a, tk_fvec_t *b, uint64_t start, uint64_t end) {
  if (end > a->n) { tk_fvec_ensure(a, end); a->n = end; }
  if (end > b->n || end <= start) return;
  tk_fvec_blas_axpy(1.0f, b->a + start, a->a + start, end - start);
}

static inline void tk_fvec_multiply_override(tk_fvec_t *a, tk_fvec_t *b, tk_fvec_t *c, uint64_t k, bool transpose_a, bool transpose_b) {
  size_t m = transpose_a ? k : a->n / k;
  size_t n = transpose_b ? k : b->n / k;
  tk_fvec_ensure(c, m * n);
  c->n = m * n;
  tk_fvec_gemm(transpose_a, transpose_b, m, n, k, 1.0f, a->a, b->a, 0.0f, c->a);
}

static inline void tk_fvec_scale_overridev(tk_fvec_t *a, tk_fvec_t *b, uint64_t start, uint64_t end) {
  if (end > a->n) { tk_fvec_ensure(a, end); a->n = end; }
  if (end > b->n || end <= start) return;
  for (size_t i = start; i < end; i++)
    a->a[i] *= b->a[i];
}

static inline tk_fvec_t *tk_fvec_rmags_override(lua_State *L, tk_fvec_t *m0, uint64_t cols) {
  uint64_t rows = m0->n / cols;
  tk_fvec_t *out = tk_fvec_create(L, rows);
  for (uint64_t r = 0; r < rows; r++)
    out->a[r] = tk_fvec_blas_nrm2(m0->a + r * cols, cols);
  out->n = rows;
  return out;
}

static inline tk_fvec_t *tk_fvec_cmags_override(lua_State *L, tk_fvec_t *m0, uint64_t cols) {
  uint64_t rows = m0->n / cols;
  tk_fvec_t *out = tk_fvec_create(L, cols);
  for (uint64_t c = 0; c < cols; c++) {
    float s = 0;
    for (uint64_t r = 0; r < rows; r++) {
      float v = m0->a[r * cols + c];
      s += v * v;
    }
    out->a[c] = sqrtf(s);
  }
  out->n = cols;
  return out;
}

static inline tk_fvec_t *tk_fvec_rsums_override(lua_State *L, tk_fvec_t *m0, uint64_t cols) {
  uint64_t rows = m0->n / cols;
  tk_fvec_t *out = tk_fvec_create(L, rows);
  for (uint64_t r = 0; r < rows; r++) {
    float s = 0;
    for (uint64_t c = 0; c < cols; c++) s += m0->a[r * cols + c];
    out->a[r] = s;
  }
  out->n = rows;
  return out;
}

static inline tk_fvec_t *tk_fvec_csums_override(lua_State *L, tk_fvec_t *m0, uint64_t cols) {
  uint64_t rows = m0->n / cols;
  tk_fvec_t *out = tk_fvec_create(L, cols);
  memset(out->a, 0, cols * sizeof(float));
  for (uint64_t r = 0; r < rows; r++)
    for (uint64_t c = 0; c < cols; c++)
      out->a[c] += m0->a[r * cols + c];
  out->n = cols;
  return out;
}

#include <santoku/iumap.h>
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

static inline void tk_fvec_mtx_zscore (
  lua_State *L, tk_fvec_t *data, uint64_t n_cols,
  tk_fvec_t *istd_in, tk_fvec_t **istd_out
) {
  uint64_t N = data->n / n_cols;
  if (istd_in) {
    #pragma omp parallel for
    for (uint64_t d = 0; d < n_cols; d++) {
      float is = istd_in->a[d];
      for (uint64_t s = 0; s < N; s++)
        data->a[s * n_cols + d] *= is;
    }
  } else {
    tk_fvec_t *is = tk_fvec_create(L, n_cols);
    is->n = n_cols;
    #pragma omp parallel for
    for (uint64_t d = 0; d < n_cols; d++) {
      double sum = 0, sum2 = 0;
      for (uint64_t s = 0; s < N; s++) {
        double v = (double)data->a[s * n_cols + d];
        sum += v; sum2 += v * v;
      }
      double m = sum / (double)N;
      double var = sum2 / (double)N - m * m;
      double istd = var > 1e-24 ? 1.0 / sqrt(var) : 0.0;
      is->a[d] = (float)istd;
      for (uint64_t s = 0; s < N; s++)
        data->a[s * n_cols + d] = (float)((double)data->a[s * n_cols + d] * istd);
    }
    *istd_out = is;
  }
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

static inline tk_fvec_t *tk_fvec_mtx_extend (
  tk_fvec_t *base, tk_fvec_t *ext,
  uint64_t n_base_features, uint64_t n_ext_features
) {
  if (base == NULL || ext == NULL) return NULL;
  uint64_t n_samples = base->n / n_base_features;
  uint64_t n_total_features = n_base_features + n_ext_features;
  tk_fvec_ensure(base, n_samples * n_total_features);
  float *base_data = base->a;
  float *ext_data = ext->a;
  for (int64_t s = (int64_t) n_samples - 1; s >= 0; s--) {
    uint64_t base_offset = (uint64_t) s * n_base_features;
    uint64_t ext_offset = (uint64_t) s * n_ext_features;
    uint64_t out_offset = (uint64_t) s * n_total_features;
    memmove(base_data + out_offset, base_data + base_offset, n_base_features * sizeof(float));
    memcpy(base_data + out_offset + n_base_features, ext_data + ext_offset, n_ext_features * sizeof(float));
  }
  base->n = n_samples * n_total_features;
  return base;
}

static inline int tk_fvec_mtx_extend_mapped (
  tk_fvec_t *base, tk_fvec_t *ext,
  tk_ivec_t *aids, tk_ivec_t *bids,
  uint64_t n_base_features, uint64_t n_ext_features, bool project
) {
  if (base == NULL || ext == NULL || aids == NULL || bids == NULL) return -1;
  uint64_t n_total_features = n_base_features + n_ext_features;
  tk_iumap_t *a_id_to_pos = tk_iumap_from_ivec(0, aids);
  if (!a_id_to_pos) return -1;
  uint64_t n_only_b = 0;
  int64_t *b_to_final = (int64_t *)calloc(bids->n, sizeof(int64_t));
  if (!b_to_final) { tk_iumap_destroy(a_id_to_pos); return -1; }
  int64_t next_pos = (int64_t)aids->n;
  for (size_t bi = 0; bi < bids->n; bi++) {
    khint_t khi = tk_iumap_get(a_id_to_pos, bids->a[bi]);
    if (khi != tk_iumap_end(a_id_to_pos)) {
      b_to_final[bi] = tk_iumap_val(a_id_to_pos, khi);
    } else {
      if (!project) { b_to_final[bi] = next_pos++; n_only_b++; }
      else b_to_final[bi] = -1;
    }
  }
  uint64_t final_n_samples = project ? aids->n : (aids->n + n_only_b);
  size_t old_aids_n = aids->n;
  if (!project) {
    if (tk_ivec_ensure(aids, final_n_samples) != 0) { free(b_to_final); tk_iumap_destroy(a_id_to_pos); return -1; }
    for (size_t i = 0; i < bids->n; i++)
      if (b_to_final[i] >= (int64_t)old_aids_n) aids->a[aids->n++] = bids->a[i];
  }
  if (tk_fvec_ensure(base, final_n_samples * n_total_features) != 0) { free(b_to_final); tk_iumap_destroy(a_id_to_pos); return -1; }
  float *base_data = base->a;
  float *ext_data = ext->a;
  tk_iumap_t *b_id_to_pos = tk_iumap_from_ivec(0, bids);
  if (!b_id_to_pos) { free(b_to_final); tk_iumap_destroy(a_id_to_pos); return -1; }
  float *new_data = calloc(final_n_samples * n_total_features, sizeof(float));
  if (!new_data) { free(b_to_final); tk_iumap_destroy(a_id_to_pos); tk_iumap_destroy(b_id_to_pos); return -1; }
  for (size_t ai = 0; ai < old_aids_n; ai++) {
    uint64_t dest_offset = ai * n_total_features;
    uint64_t src_offset = ai * n_base_features;
    memcpy(new_data + dest_offset, base_data + src_offset, n_base_features * sizeof(float));
    khint_t khi = tk_iumap_get(b_id_to_pos, aids->a[ai]);
    if (khi != tk_iumap_end(b_id_to_pos)) {
      int64_t b_idx = tk_iumap_val(b_id_to_pos, khi);
      uint64_t ext_src_offset = (uint64_t)b_idx * n_ext_features;
      memcpy(new_data + dest_offset + n_base_features, ext_data + ext_src_offset, n_ext_features * sizeof(float));
    }
  }
  for (size_t bi = 0; bi < bids->n; bi++) {
    int64_t final_pos = b_to_final[bi];
    if (final_pos >= (int64_t)old_aids_n) {
      uint64_t dest_offset = (uint64_t)final_pos * n_total_features;
      uint64_t ext_src_offset = bi * n_ext_features;
      memcpy(new_data + dest_offset + n_base_features, ext_data + ext_src_offset, n_ext_features * sizeof(float));
    }
  }
  memcpy(base_data, new_data, final_n_samples * n_total_features * sizeof(float));
  base->n = final_n_samples * n_total_features;
  free(new_data);
  free(b_to_final);
  tk_iumap_destroy(a_id_to_pos);
  tk_iumap_destroy(b_id_to_pos);
  return 0;
}

static inline tk_fvec_t *tk_fvec_mtx_select (
  tk_fvec_t *src_matrix, tk_ivec_t *selected_features,
  tk_ivec_t *sample_ids, uint64_t n_features,
  tk_fvec_t *dest, uint64_t dest_sample, uint64_t dest_stride
) {
  if (src_matrix == NULL) return NULL;
  uint64_t n_samples = src_matrix->n / n_features;
  if (dest == NULL && (selected_features == NULL || selected_features->n == 0) &&
      (sample_ids == NULL || sample_ids->n == 0))
    return src_matrix;
  uint64_t n_output_samples = (sample_ids != NULL && sample_ids->n > 0) ? sample_ids->n : n_samples;
  uint64_t n_selected_features = (selected_features != NULL && selected_features->n > 0) ? selected_features->n : n_features;
  float *src_data = src_matrix->a;
  if (dest != NULL) {
    uint64_t final_stride = (dest_stride > 0) ? dest_stride : n_selected_features;
    if (dest_sample == 0) { tk_fvec_ensure(dest, n_output_samples * final_stride); dest->n = n_output_samples * final_stride; }
    else { tk_fvec_ensure(dest, (dest_sample + n_output_samples) * final_stride); dest->n = (dest_sample + n_output_samples) * final_stride; }
    float *dest_data = dest->a;
    for (uint64_t si = 0; si < n_output_samples; si++) {
      uint64_t s = (sample_ids != NULL && sample_ids->n > 0) ? (uint64_t) sample_ids->a[si] : si;
      if (s >= n_samples) continue;
      uint64_t src_offset = s * n_features;
      uint64_t d_offset = (dest_sample + si) * final_stride;
      if (selected_features != NULL && selected_features->n > 0) {
        for (uint64_t f = 0; f < selected_features->n; f++) {
          int64_t feat_idx = selected_features->a[f];
          if (feat_idx >= 0 && (uint64_t)feat_idx < n_features)
            dest_data[d_offset + f] = src_data[src_offset + (uint64_t)feat_idx];
        }
      } else {
        memcpy(dest_data + d_offset, src_data + src_offset, n_features * sizeof(float));
      }
    }
    return dest;
  } else {
    if (sample_ids != NULL && sample_ids->n > 0) {
      float *tmp = (float *) malloc(n_output_samples * n_selected_features * sizeof(float));
      if (!tmp) return NULL;
      for (uint64_t si = 0; si < n_output_samples; si++) {
        uint64_t s = (uint64_t) sample_ids->a[si];
        if (s >= n_samples) continue;
        uint64_t src_offset = s * n_features;
        uint64_t tmp_offset = si * n_selected_features;
        if (selected_features != NULL && selected_features->n > 0) {
          for (uint64_t f = 0; f < selected_features->n; f++) {
            int64_t feat_idx = selected_features->a[f];
            if (feat_idx >= 0 && (uint64_t)feat_idx < n_features)
              tmp[tmp_offset + f] = src_data[src_offset + (uint64_t)feat_idx];
          }
        } else {
          memcpy(tmp + tmp_offset, src_data + src_offset, n_features * sizeof(float));
        }
      }
      free(src_matrix->a);
      src_matrix->a = tmp;
      src_matrix->n = n_output_samples * n_selected_features;
      src_matrix->m = n_output_samples * n_selected_features;
    } else if (selected_features != NULL && selected_features->n > 0) {
      for (uint64_t si = 0; si < n_samples; si++) {
        uint64_t src_offset = si * n_features;
        uint64_t dst_offset = si * n_selected_features;
        for (uint64_t f = 0; f < selected_features->n; f++) {
          int64_t feat_idx = selected_features->a[f];
          if (feat_idx >= 0 && (uint64_t)feat_idx < n_features)
            src_data[dst_offset + f] = src_data[src_offset + (uint64_t)feat_idx];
        }
      }
      src_matrix->n = n_samples * n_selected_features;
    }
    return src_matrix;
  }
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

static inline void tk_fvec_mtx_standardize (
  lua_State *L, tk_fvec_t *data, uint64_t n_cols,
  tk_fvec_t *mean_in, tk_fvec_t *istd_in,
  tk_fvec_t **mean_out, tk_fvec_t **istd_out
) {
  uint64_t N = data->n / n_cols;
  if (mean_in) {
    #pragma omp parallel for
    for (uint64_t d = 0; d < n_cols; d++) {
      float mu = mean_in->a[d];
      float is = istd_in->a[d];
      for (uint64_t s = 0; s < N; s++)
        data->a[s * n_cols + d] = (data->a[s * n_cols + d] - mu) * is;
    }
  } else {
    tk_fvec_t *mu = tk_fvec_create(L, n_cols);
    mu->n = n_cols;
    tk_fvec_t *is = tk_fvec_create(L, n_cols);
    is->n = n_cols;
    #pragma omp parallel for
    for (uint64_t d = 0; d < n_cols; d++) {
      double sum = 0, sum2 = 0;
      for (uint64_t s = 0; s < N; s++) {
        double v = (double)data->a[s * n_cols + d];
        sum += v; sum2 += v * v;
      }
      double m = sum / (double)N;
      double var = sum2 / (double)N - m * m;
      double istd = var > 1e-24 ? 1.0 / sqrt(var) : 0.0;
      mu->a[d] = (float)m;
      is->a[d] = (float)istd;
      for (uint64_t s = 0; s < N; s++)
        data->a[s * n_cols + d] = (float)(((double)data->a[s * n_cols + d] - m) * istd);
    }
    *mean_out = mu;
    *istd_out = is;
  }
}

#endif
