#ifndef TK_DVEC_EXT_H
#define TK_DVEC_EXT_H

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

#include <santoku/dvec/base.h>
#include <santoku/ivec.h>
#include <santoku/rvec/base.h>

static inline tk_ivec_t *tk_rvec_keys (lua_State *L, tk_rvec_t *P, tk_ivec_t *out);
static inline tk_dvec_t *tk_rvec_values (lua_State *L, tk_rvec_t *P, tk_dvec_t *out);
static inline tk_ivec_t *tk_dvec_mtx_top_variance (lua_State *L, tk_dvec_t *matrix, uint64_t n_samples, uint64_t n_features, uint64_t top_k);
static inline tk_ivec_t *tk_dvec_mtx_top_skewness (lua_State *L, tk_dvec_t *matrix, uint64_t n_samples, uint64_t n_features, uint64_t top_k);
static inline tk_ivec_t *tk_dvec_mtx_top_entropy (lua_State *L, tk_dvec_t *matrix, uint64_t n_samples, uint64_t n_features, uint64_t top_k, uint64_t n_bins);
static inline tk_ivec_t *tk_dvec_mtx_top_bimodality (lua_State *L, tk_dvec_t *matrix, uint64_t n_samples, uint64_t n_features, uint64_t top_k);
static inline tk_ivec_t *tk_dvec_mtx_top_dip (lua_State *L, tk_dvec_t *matrix, uint64_t n_samples, uint64_t n_features, uint64_t top_k);
static inline void tk_dvec_mtx_topk (lua_State *L, tk_dvec_t *queries, tk_dvec_t *corpus, uint64_t n_queries, uint64_t n_corpus, uint64_t d, uint64_t k);

#define TK_GENERATE_SINGLE
#include <santoku/parallel/tpl.h>
#include <santoku/dvec/ext_tpl.h>
#undef TK_GENERATE_SINGLE

#include <santoku/parallel/tpl.h>
#include <santoku/dvec/ext_tpl.h>

static inline size_t tk_dvec_scores_max_curvature (
  double *scores,
  size_t n,
  double *out_val
) {
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? scores[0] : 0.0;
    return n > 0 ? n - 1 : 0;
  }
  double max_curv = 0.0;
  size_t max_idx = 1;
  for (size_t i = 1; i < n - 1; i++) {
    double curv = fabs(scores[i-1] - 2.0 * scores[i] + scores[i+1]);
    if (curv > max_curv) {
      max_curv = curv;
      max_idx = i;
    }
  }

  if (max_curv < 1e-10) {
    if (out_val) *out_val = scores[n - 1];
    return n - 1;
  }
  if (out_val) *out_val = scores[max_idx];
  return max_idx;
}

static inline size_t tk_dvec_scores_lmethod (
  double *scores,
  size_t n,
  double *out_val
) {
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? scores[0] : 0.0;
    return n > 0 ? n - 1 : 0;
  }

  double min_val = scores[0], max_val = scores[0];
  for (size_t i = 1; i < n; i++) {
    if (scores[i] < min_val) min_val = scores[i];
    if (scores[i] > max_val) max_val = scores[i];
  }
  if (max_val - min_val < 1e-10) {
    if (out_val) *out_val = scores[n - 1];
    return n - 1;
  }
  double best_rmse = DBL_MAX;
  size_t best_k = 1;
  for (size_t k = 1; k < n - 1; k++) {
    double sum_x1 = 0.0, sum_y1 = 0.0, sum_xy1 = 0.0, sum_xx1 = 0.0;
    for (size_t i = 0; i <= k; i++) {
      double x = (double)i;
      double y = scores[i];
      sum_x1 += x;
      sum_y1 += y;
      sum_xy1 += x * y;
      sum_xx1 += x * x;
    }
    size_t n1 = k + 1;
    double mean_x1 = sum_x1 / (double)n1;
    double mean_y1 = sum_y1 / (double)n1;
    double slope1 = (sum_xy1 - (double)n1 * mean_x1 * mean_y1) / (sum_xx1 - (double)n1 * mean_x1 * mean_x1 + 1e-10);
    double intercept1 = mean_y1 - slope1 * mean_x1;
    double sum_x2 = 0.0, sum_y2 = 0.0, sum_xy2 = 0.0, sum_xx2 = 0.0;
    for (size_t i = k + 1; i < n; i++) {
      double x = (double)i;
      double y = scores[i];
      sum_x2 += x;
      sum_y2 += y;
      sum_xy2 += x * y;
      sum_xx2 += x * x;
    }
    size_t n2 = n - k - 1;
    double mean_x2 = sum_x2 / (double)n2;
    double mean_y2 = sum_y2 / (double)n2;
    double slope2 = (sum_xy2 - (double)n2 * mean_x2 * mean_y2) / (sum_xx2 - (double)n2 * mean_x2 * mean_x2 + 1e-10);
    double intercept2 = mean_y2 - slope2 * mean_x2;
    double sse = 0.0;
    for (size_t i = 0; i <= k; i++) {
      double pred = slope1 * (double)i + intercept1;
      double err = scores[i] - pred;
      sse += err * err;
    }
    for (size_t i = k + 1; i < n; i++) {
      double pred = slope2 * (double)i + intercept2;
      double err = scores[i] - pred;
      sse += err * err;
    }
    double rmse = sqrt(sse / (double)n);
    if (rmse < best_rmse) {
      best_rmse = rmse;
      best_k = k;
    }
  }
  if (out_val) *out_val = scores[best_k];
  return best_k;
}

static inline size_t tk_dvec_scores_max_gap (
  double *scores,
  size_t n,
  double *out_val
) {
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? scores[0] : 0.0;
    return n > 0 ? n - 1 : 0;
  }
  double max_gap = 0.0;
  size_t max_idx = 0;
  for (size_t i = 0; i < n - 1; i++) {
    double gap = fabs(scores[i + 1] - scores[i]);
    if (gap > max_gap) {
      max_gap = gap;
      max_idx = i;
    }
  }

  if (max_gap < 1e-10) {
    if (out_val) *out_val = scores[n - 1];
    return n - 1;
  }
  if (out_val) *out_val = scores[max_idx];
  return max_idx;
}



static inline size_t tk_dvec_scores_plateau (
  double *scores,
  size_t n,
  double tolerance,
  double *out_val
) {
  if (n == 0) {
    if (out_val) *out_val = 0.0;
    return 0;
  }
  if (n == 1) {
    if (out_val) *out_val = scores[0];
    return 0;
  }
  if (tolerance <= 0.0) tolerance = 0.01;

  double min_score = scores[0], max_score = scores[0];
  for (size_t i = 1; i < n; i++) {
    if (scores[i] < min_score) min_score = scores[i];
    if (scores[i] > max_score) max_score = scores[i];
  }
  double range = max_score - min_score;
  if (range <= 0.0) {
    if (out_val) *out_val = scores[n - 1];
    return n - 1;
  }

  double abs_tolerance = tolerance * range;
  double base = scores[0];
  size_t end_idx = 0;
  for (size_t i = 1; i < n; i++) {
    if (fabs(scores[i] - base) <= abs_tolerance) {
      end_idx = i;
    } else {
      break;
    }
  }
  if (out_val) *out_val = scores[end_idx];
  return end_idx;
}

static inline size_t tk_dvec_scores_kneedle (
  double *scores,
  size_t n,
  double sensitivity,
  double *out_val
) {
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? scores[0] : 0.0;
    return 0;
  }
  if (sensitivity <= 0.0)
    sensitivity = 1.0;
  double min_score = scores[0];
  double max_score = scores[0];
  for (size_t i = 1; i < n; i++) {
    if (scores[i] < min_score) min_score = scores[i];
    if (scores[i] > max_score) max_score = scores[i];
  }
  double score_range = max_score - min_score;
  if (score_range < 1e-10) {
    if (out_val) *out_val = scores[0];
    return 0;
  }
  double *normalized = (double *)malloc(n * sizeof(double));
  if (!normalized) {
    if (out_val) *out_val = scores[0];
    return 0;
  }
  for (size_t i = 0; i < n; i++) {
    double x_norm = (double)i / (double)(n - 1);
    double y_norm = (scores[i] - min_score) / score_range;
    normalized[i] = y_norm - x_norm;
  }
  double *smoothed = (double *)malloc(n * sizeof(double));
  if (!smoothed) {
    free(normalized);
    if (out_val) *out_val = scores[0];
    return 0;
  }
  smoothed[0] = normalized[0];
  smoothed[n - 1] = normalized[n - 1];
  for (size_t i = 1; i < n - 1; i++) {
    smoothed[i] = (normalized[i - 1] + normalized[i] + normalized[i + 1]) / 3.0;
  }
  double max_diff = -DBL_MAX;
  size_t knee_idx = 0;
  for (size_t i = 0; i < n; i++) {
    if (smoothed[i] > max_diff) {
      max_diff = smoothed[i];
      knee_idx = i;
    }
  }
  free(normalized);
  free(smoothed);
  double threshold = max_diff - sensitivity / (double)n;
  size_t final_knee = knee_idx;
  for (size_t i = 0; i < n; i++) {
    double x_norm = (double)i / (double)(n - 1);
    double y_norm = (scores[i] - min_score) / score_range;
    double diff = y_norm - x_norm;
    if (diff >= threshold) {
      final_knee = i;
    }
  }
  if (out_val) *out_val = scores[final_knee];
  return final_knee;
}







static inline size_t tk_dvec_scores_first_gap (
  double *scores,
  size_t n,
  double alpha,
  double *out_val
) {
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? scores[0] : 0.0;
    return n > 0 ? n - 1 : 0;
  }
  if (alpha <= 0.0) alpha = 3.0;

  size_t n_gaps = n - 1;
  double *gaps = (double *)malloc(n_gaps * sizeof(double));
  if (!gaps) {
    if (out_val) *out_val = scores[n - 1];
    return n - 1;
  }

  for (size_t i = 0; i < n_gaps; i++) {
    gaps[i] = fabs(scores[i + 1] - scores[i]);
  }


  for (size_t i = 1; i < n_gaps; i++) {
    double key = gaps[i];
    size_t j = i;
    while (j > 0 && gaps[j - 1] > key) {
      gaps[j] = gaps[j - 1];
      j--;
    }
    gaps[j] = key;
  }

  double median_gap;
  if (n_gaps % 2 == 1) {
    median_gap = gaps[n_gaps / 2];
  } else {
    median_gap = (gaps[n_gaps / 2 - 1] + gaps[n_gaps / 2]) / 2.0;
  }
  free(gaps);

  if (median_gap <= 0.0) {
    double max_gap = 0.0;
    size_t max_idx = n - 1;
    for (size_t i = 0; i < n - 1; i++) {
      double gap = fabs(scores[i + 1] - scores[i]);
      if (gap > max_gap) {
        max_gap = gap;
        max_idx = i;
      }
    }
    if (out_val) *out_val = scores[max_idx];
    return max_idx;
  }

  double threshold = alpha * median_gap;

  for (size_t i = 0; i < n - 1; i++) {
    double gap = fabs(scores[i + 1] - scores[i]);
    if (gap > threshold) {
      if (out_val) *out_val = scores[i];
      return i;
    }
  }

  if (out_val) *out_val = scores[n - 1];
  return n - 1;
}




static inline size_t tk_dvec_scores_otsu (
  double *scores,
  size_t n,
  double *out_val
) {
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? scores[0] : 0.0;
    return n > 0 ? n - 1 : 0;
  }


  double min_val = scores[0], max_val = scores[0];
  for (size_t i = 1; i < n; i++) {
    if (scores[i] < min_val) min_val = scores[i];
    if (scores[i] > max_val) max_val = scores[i];
  }
  if (max_val - min_val < 1e-10) {
    if (out_val) *out_val = scores[n - 1];
    return n - 1;
  }


  double total_sum = 0.0;
  for (size_t i = 0; i < n; i++) {
    total_sum += scores[i];
  }




  double best_variance = -1.0;
  size_t best_k = 0;
  double sum0 = 0.0;

  for (size_t k = 0; k < n - 1; k++) {
    sum0 += scores[k];
    double sum1 = total_sum - sum0;

    size_t n0 = k + 1;
    size_t n1 = n - n0;

    double w0 = (double)n0 / (double)n;
    double w1 = (double)n1 / (double)n;

    double mean0 = sum0 / (double)n0;
    double mean1 = sum1 / (double)n1;


    double variance = w0 * w1 * (mean0 - mean1) * (mean0 - mean1);

    if (variance > best_variance) {
      best_variance = variance;
      best_k = k;
    }
  }

  if (out_val) *out_val = scores[best_k];
  return best_k;
}

static inline void tk_dvec_scores_tolerance (
  double *scores,
  size_t n,
  double tolerance,
  size_t *out_start,
  size_t *out_end
) {
  if (n == 0) {
    *out_start = 0;
    *out_end = 0;
    return;
  }
  if (n == 1) {
    *out_start = 0;
    *out_end = 0;
    return;
  }
  size_t best_start = 0;
  size_t best_end = 0;
  size_t best_len = 1;
  for (size_t i = 0; i < n; i++) {
    double base = scores[i];
    size_t span_end = i;
    for (size_t j = i + 1; j < n; j++) {
      if (fabs(scores[j] - base) <= tolerance) {
        span_end = j;
      } else {
        break;
      }
    }
    size_t span_len = span_end - i + 1;
    if (span_len > best_len) {
      best_len = span_len;
      best_start = i;
      best_end = span_end;
    }
  }
  *out_start = best_start;
  *out_end = best_end;
}

#if !defined(__EMSCRIPTEN__)

static inline void tk_dvec_gemv(
  bool transpose, uint64_t rows, uint64_t cols,
  double alpha, double *A, double *x, double beta, double *y
) {
  cblas_dgemv(CblasRowMajor, transpose ? CblasTrans : CblasNoTrans, rows, cols, alpha, A, cols, x, 1, beta, y, 1);
}

static inline void tk_dvec_gemm(
  bool transpose_a, bool transpose_b,
  uint64_t m, uint64_t n, uint64_t k,
  double alpha, double *A, double *B, double beta, double *C
) {
  cblas_dgemm(CblasRowMajor, transpose_a ? CblasTrans : CblasNoTrans, transpose_b ? CblasTrans : CblasNoTrans, m, n, k, alpha, A, transpose_a ? m : k, B, transpose_b ? k : n, beta, C, n);
}

static inline double tk_dvec_blas_dot(double *x, double *y, uint64_t n) {
  return cblas_ddot(n, x, 1, y, 1);
}

static inline void tk_dvec_blas_scal(double alpha, double *x, uint64_t n) {
  cblas_dscal(n, alpha, x, 1);
}

static inline void tk_dvec_blas_axpy(double alpha, double *x, double *y, uint64_t n) {
  cblas_daxpy(n, alpha, x, 1, y, 1);
}

static inline void tk_dvec_blas_copy(double *x, double *y, uint64_t n) {
  cblas_dcopy(n, x, 1, y, 1);
}

static inline double tk_dvec_blas_nrm2(double *x, uint64_t n) {
  return cblas_dnrm2(n, x, 1);
}

#else

static inline void tk_dvec_gemv(
  bool transpose, uint64_t rows, uint64_t cols,
  double alpha, double *A, double *x, double beta, double *y
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

static inline void tk_dvec_gemm(
  bool transpose_a, bool transpose_b,
  uint64_t m, uint64_t n, uint64_t k,
  double alpha, double *A, double *B, double beta, double *C
) {
  for (uint64_t i = 0; i < m * n; i++) C[i] *= beta;
  for (uint64_t i = 0; i < m; i++)
    for (uint64_t j = 0; j < n; j++)
      for (uint64_t l = 0; l < k; l++)
        C[i * n + j] += alpha *
          (transpose_a ? A[l * m + i] : A[i * k + l]) *
          (transpose_b ? B[j * k + l] : B[l * n + j]);
}

static inline double tk_dvec_blas_dot(double *x, double *y, uint64_t n) {
  double s = 0; for (uint64_t i = 0; i < n; i++) s += x[i] * y[i]; return s;
}

static inline void tk_dvec_blas_scal(double alpha, double *x, uint64_t n) {
  for (uint64_t i = 0; i < n; i++) x[i] *= alpha;
}

static inline void tk_dvec_blas_axpy(double alpha, double *x, double *y, uint64_t n) {
  for (uint64_t i = 0; i < n; i++) y[i] += alpha * x[i];
}

static inline void tk_dvec_blas_copy(double *x, double *y, uint64_t n) {
  memcpy(y, x, n * sizeof(double));
}

static inline double tk_dvec_blas_nrm2(double *x, uint64_t n) {
  double s = 0; for (uint64_t i = 0; i < n; i++) s += x[i] * x[i]; return sqrt(s);
}

#endif

static inline double tk_dvec_dot_override(tk_dvec_t *a, tk_dvec_t *b) {
  uint64_t n = a->n < b->n ? a->n : b->n;
  return tk_dvec_blas_dot(a->a, b->a, n);
}

static inline void tk_dvec_scale_override(tk_dvec_t *v, double scale, uint64_t start, uint64_t end) {
  if (end > v->n) {
    tk_dvec_ensure(v, end);
    v->n = end;
  }
  if (end <= start) return;
  tk_dvec_blas_scal(scale, v->a + start, end - start);
}

static inline void tk_dvec_addv_override(tk_dvec_t *a, tk_dvec_t *b, uint64_t start, uint64_t end) {
  if (end > a->n) {
    tk_dvec_ensure(a, end);
    a->n = end;
  }
  if (end > b->n || end <= start) return;
  tk_dvec_blas_axpy(1.0, b->a + start, a->a + start, end - start);
}

static inline void tk_dvec_multiply_override(tk_dvec_t *a, tk_dvec_t *b, tk_dvec_t *c, uint64_t k, bool transpose_a, bool transpose_b) {
  size_t m = transpose_a ? k : a->n / k;
  size_t n = transpose_b ? k : b->n / k;
  tk_dvec_ensure(c, m * n);
  c->n = m * n;
  tk_dvec_gemm(transpose_a, transpose_b, m, n, k, 1.0, a->a, b->a, 0.0, c->a);
}

static inline void tk_dvec_scale_overridev(tk_dvec_t *a, tk_dvec_t *b, uint64_t start, uint64_t end) {
  if (end > a->n) {
    tk_dvec_ensure(a, end);
    a->n = end;
  }
  if (end > b->n || end <= start) return;
  for (size_t i = start; i < end; i++)
    a->a[i] *= b->a[i];
}

static inline tk_dvec_t *tk_dvec_rmags_override(lua_State *L, tk_dvec_t *m0, uint64_t cols) {
  uint64_t rows = m0->n / cols;
  tk_dvec_t *out = tk_dvec_create(L, rows);
  for (uint64_t r = 0; r < rows; r++)
    out->a[r] = tk_dvec_blas_nrm2(m0->a + r * cols, cols);
  out->n = rows;
  return out;
}

static inline tk_dvec_t *tk_dvec_cmags_override(lua_State *L, tk_dvec_t *m0, uint64_t cols) {
  uint64_t rows = m0->n / cols;
  tk_dvec_t *out = tk_dvec_create(L, cols);
  for (uint64_t c = 0; c < cols; c++) {
    double s = 0;
    for (uint64_t r = 0; r < rows; r++) {
      double v = m0->a[r * cols + c];
      s += v * v;
    }
    out->a[c] = sqrt(s);
  }
  out->n = cols;
  return out;
}

static inline tk_dvec_t *tk_dvec_rsums_override(lua_State *L, tk_dvec_t *m0, uint64_t cols) {
  uint64_t rows = m0->n / cols;
  tk_dvec_t *out = tk_dvec_create(L, rows);
  for (uint64_t r = 0; r < rows; r++) {
    double s = 0;
    for (uint64_t c = 0; c < cols; c++) s += m0->a[r * cols + c];
    out->a[r] = s;
  }
  out->n = rows;
  return out;
}

static inline tk_dvec_t *tk_dvec_csums_override(lua_State *L, tk_dvec_t *m0, uint64_t cols) {
  uint64_t rows = m0->n / cols;
  tk_dvec_t *out = tk_dvec_create(L, cols);
  memset(out->a, 0, cols * sizeof(double));
  for (uint64_t r = 0; r < rows; r++)
    for (uint64_t c = 0; c < cols; c++)
      out->a[c] += m0->a[r * cols + c];
  out->n = cols;
  return out;
}

#include <santoku/iumap.h>
#include <santoku/cvec/ext.h>

static inline void tk_dvec_mtx_center (
  lua_State *L,
  tk_dvec_t *data,
  uint64_t n_cols,
  tk_dvec_t *mean_in,
  tk_dvec_t **mean_out
) {
  uint64_t N = data->n / n_cols;
  if (mean_in) {
    #pragma omp parallel for
    for (uint64_t d = 0; d < n_cols; d++) {
      double mu = mean_in->a[d];
      for (uint64_t s = 0; s < N; s++)
        data->a[s * n_cols + d] -= mu;
    }
  } else {
    tk_dvec_t *mu = tk_dvec_create(L, n_cols);
    mu->n = n_cols;
    #pragma omp parallel for
    for (uint64_t d = 0; d < n_cols; d++) {
      double sum = 0;
      for (uint64_t s = 0; s < N; s++)
        sum += data->a[s * n_cols + d];
      double m = sum / (double)N;
      mu->a[d] = m;
      for (uint64_t s = 0; s < N; s++)
        data->a[s * n_cols + d] -= m;
    }
    *mean_out = mu;
  }
}

static inline void tk_dvec_mtx_zscore (
  lua_State *L,
  tk_dvec_t *data,
  uint64_t n_cols,
  tk_dvec_t *istd_in,
  tk_dvec_t **istd_out
) {
  uint64_t N = data->n / n_cols;
  if (istd_in) {
    #pragma omp parallel for
    for (uint64_t d = 0; d < n_cols; d++) {
      double is = istd_in->a[d];
      for (uint64_t s = 0; s < N; s++)
        data->a[s * n_cols + d] *= is;
    }
  } else {
    tk_dvec_t *is = tk_dvec_create(L, n_cols);
    is->n = n_cols;
    #pragma omp parallel for
    for (uint64_t d = 0; d < n_cols; d++) {
      double sum = 0, sum2 = 0;
      for (uint64_t s = 0; s < N; s++) {
        double v = data->a[s * n_cols + d];
        sum += v;
        sum2 += v * v;
      }
      double m = sum / (double)N;
      double var = sum2 / (double)N - m * m;
      double istd = var > 1e-24 ? 1.0 / sqrt(var) : 0.0;
      is->a[d] = istd;
      for (uint64_t s = 0; s < N; s++)
        data->a[s * n_cols + d] *= istd;
    }
    *istd_out = is;
  }
}

static inline void tk_dvec_mtx_sign_raw (
  char *out,
  double *X,
  uint64_t N,
  uint64_t stride,
  uint64_t K
) {
  #pragma omp parallel for
  for (uint64_t i = 0; i < N; i++) {
    double *row = X + i * stride;
    uint8_t *out_row = (uint8_t *)(out + i * TK_CVEC_BITS_BYTES(K));
    uint64_t full_bytes = K / 8;
    for (uint64_t byte_idx = 0; byte_idx < full_bytes; byte_idx++) {
      uint8_t byte_val = 0;
      uint64_t j_base = byte_idx * 8;
      for (uint64_t bit = 0; bit < 8; bit++)
        byte_val |= (uint8_t)((row[j_base + bit] >= 0.0) << bit);
      out_row[byte_idx] = byte_val;
    }
    uint64_t remaining_start = full_bytes * 8;
    if (remaining_start < K) {
      uint8_t byte_val = 0;
      for (uint64_t j = remaining_start; j < K; j++)
        byte_val |= (uint8_t)((row[j] >= 0.0) << (j - remaining_start));
      out_row[full_bytes] = byte_val;
    }
  }
}

static inline tk_cvec_t *tk_dvec_mtx_sign (
  lua_State *L,
  tk_dvec_t *codes,
  uint64_t n_dims,
  uint64_t n_trunc
) {
  const size_t N = codes->n / n_dims;
  tk_cvec_t *binary = tk_cvec_create(L, N * TK_CVEC_BITS_BYTES(n_trunc));
  binary->n = N * TK_CVEC_BITS_BYTES(n_trunc);
  tk_cvec_zero(binary);
  tk_dvec_mtx_sign_raw(binary->a, codes->a, N, n_dims, n_trunc);
  return binary;
}

static inline void tk_dvec_mtx_threshold_raw (
  char *out,
  double *X,
  double *thresholds,
  uint64_t N,
  uint64_t K
) {
  #pragma omp parallel for
  for (uint64_t i = 0; i < N; i++) {
    double *row = X + i * K;
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

static int tk_dvec_mtx_cmp_double (const void *a, const void *b) {
  double da = *(const double *)a;
  double db = *(const double *)b;
  return (da > db) - (da < db);
}

static inline tk_cvec_t *tk_dvec_mtx_median (
  lua_State *L,
  tk_dvec_t *codes,
  uint64_t n_dims,
  tk_dvec_t **medians_out
) {
  const uint64_t K = n_dims;
  const size_t N = codes->n / K;
  tk_cvec_t *binary = tk_cvec_create(L, N * TK_CVEC_BITS_BYTES(K));
  binary->n = N * TK_CVEC_BITS_BYTES(K);
  tk_cvec_zero(binary);
  double *medians = malloc(K * sizeof(double));
  double *col_buf = malloc(N * sizeof(double));
  for (uint64_t k = 0; k < K; k++) {
    for (uint64_t i = 0; i < N; i++)
      col_buf[i] = codes->a[i * K + k];
    qsort(col_buf, N, sizeof(double), tk_dvec_mtx_cmp_double);
    medians[k] = (N % 2 == 1) ? col_buf[N / 2] : (col_buf[N / 2 - 1] + col_buf[N / 2]) / 2.0;
  }
  free(col_buf);
  tk_dvec_mtx_threshold_raw(binary->a, codes->a, medians, N, K);
  if (medians_out) {
    tk_dvec_t *med_vec = tk_dvec_create(L, K);
    med_vec->n = K;
    memcpy(med_vec->a, medians, K * sizeof(double));
    *medians_out = med_vec;
  }
  free(medians);
  return binary;
}

static inline void tk_dvec_round (tk_dvec_t *v, uint64_t start, uint64_t end) {
  if (end > v->n) end = v->n;
  for (uint64_t i = start; i < end; i++)
    v->a[i] = round(v->a[i]);
}

static inline void tk_dvec_trunc (tk_dvec_t *v, uint64_t start, uint64_t end) {
  if (end > v->n) end = v->n;
  for (uint64_t i = start; i < end; i++)
    v->a[i] = trunc(v->a[i]);
}

static inline void tk_dvec_floor (tk_dvec_t *v, uint64_t start, uint64_t end) {
  if (end > v->n) end = v->n;
  for (uint64_t i = start; i < end; i++)
    v->a[i] = floor(v->a[i]);
}

static inline void tk_dvec_ceil (tk_dvec_t *v, uint64_t start, uint64_t end) {
  if (end > v->n) end = v->n;
  for (uint64_t i = start; i < end; i++)
    v->a[i] = ceil(v->a[i]);
}

static inline tk_ivec_t *tk_dvec_to_ivec (lua_State *L, tk_dvec_t *v) {
  tk_ivec_t *out = tk_ivec_create(L, v->n);
  for (uint64_t i = 0; i < v->n; i++)
    out->a[i] = (int64_t)v->a[i];
  return out;
}

static inline tk_fvec_t *tk_dvec_to_fvec (lua_State *L, tk_dvec_t *v, tk_fvec_t *out) {
  if (out == NULL) {
    out = tk_fvec_create(L, v->n);
  } else {
    tk_fvec_ensure(out, v->n);
    out->n = v->n;
  }
  for (uint64_t i = 0; i < v->n; i++)
    out->a[i] = (float)v->a[i];
  return out;
}

static inline tk_dvec_t *tk_dvec_mtx_extend (
  tk_dvec_t *base,
  tk_dvec_t *ext,
  uint64_t n_base_features,
  uint64_t n_ext_features
) {
  if (base == NULL || ext == NULL)
    return NULL;
  uint64_t n_samples = base->n / n_base_features;
  uint64_t n_total_features = n_base_features + n_ext_features;
  tk_dvec_ensure(base, n_samples * n_total_features);
  double *base_data = base->a;
  double *ext_data = ext->a;
  for (int64_t s = (int64_t) n_samples - 1; s >= 0; s--) {
    uint64_t base_offset = (uint64_t) s * n_base_features;
    uint64_t ext_offset = (uint64_t) s * n_ext_features;
    uint64_t out_offset = (uint64_t) s * n_total_features;
    memmove(base_data + out_offset, base_data + base_offset, n_base_features * sizeof(double));
    memcpy(base_data + out_offset + n_base_features, ext_data + ext_offset, n_ext_features * sizeof(double));
  }
  base->n = n_samples * n_total_features;
  return base;
}

static inline int tk_dvec_mtx_extend_mapped (
  tk_dvec_t *base,
  tk_dvec_t *ext,
  tk_ivec_t *aids,
  tk_ivec_t *bids,
  uint64_t n_base_features,
  uint64_t n_ext_features,
  bool project
) {
  if (base == NULL || ext == NULL || aids == NULL || bids == NULL)
    return -1;
  uint64_t n_total_features = n_base_features + n_ext_features;
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
  if (tk_dvec_ensure(base, final_n_samples * n_total_features) != 0) {
    free(b_to_final);
    tk_iumap_destroy(a_id_to_pos);
    return -1;
  }
  double *base_data = base->a;
  double *ext_data = ext->a;
  tk_iumap_t *b_id_to_pos = tk_iumap_from_ivec(0, bids);
  if (!b_id_to_pos) {
    free(b_to_final);
    tk_iumap_destroy(a_id_to_pos);
    return -1;
  }
  double *new_data = calloc(final_n_samples * n_total_features, sizeof(double));
  if (!new_data) {
    free(b_to_final);
    tk_iumap_destroy(a_id_to_pos);
    tk_iumap_destroy(b_id_to_pos);
    return -1;
  }
  for (size_t ai = 0; ai < old_aids_n; ai++) {
    uint64_t dest_offset = ai * n_total_features;
    uint64_t src_offset = ai * n_base_features;
    memcpy(new_data + dest_offset, base_data + src_offset, n_base_features * sizeof(double));
    khint_t khi = tk_iumap_get(b_id_to_pos, aids->a[ai]);
    if (khi != tk_iumap_end(b_id_to_pos)) {
      int64_t b_idx = tk_iumap_val(b_id_to_pos, khi);
      uint64_t ext_src_offset = (uint64_t)b_idx * n_ext_features;
      memcpy(new_data + dest_offset + n_base_features, ext_data + ext_src_offset, n_ext_features * sizeof(double));
    }
  }
  for (size_t bi = 0; bi < bids->n; bi++) {
    int64_t final_pos = b_to_final[bi];
    if (final_pos >= (int64_t)old_aids_n) {
      uint64_t dest_offset = (uint64_t)final_pos * n_total_features;
      uint64_t ext_src_offset = bi * n_ext_features;
      memcpy(new_data + dest_offset + n_base_features, ext_data + ext_src_offset, n_ext_features * sizeof(double));
    }
  }
  memcpy(base_data, new_data, final_n_samples * n_total_features * sizeof(double));
  base->n = final_n_samples * n_total_features;
  free(new_data);
  free(b_to_final);
  tk_iumap_destroy(a_id_to_pos);
  tk_iumap_destroy(b_id_to_pos);
  return 0;
}

static inline tk_ivec_t *tk_dvec_mtx_top_variance (
  lua_State *L,
  tk_dvec_t *matrix,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t top_k
) {
  if (matrix == NULL || n_samples == 0 || n_features == 0)
    return NULL;
  double *data = matrix->a;
  tk_rvec_t *top_heap = tk_rvec_create(0, 0);
  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    double sum = 0.0;
    for (uint64_t s = 0; s < n_samples; s++) {
      sum += data[s * n_features + f];
    }
    double mean = sum / (double)n_samples;
    double var_sum = 0.0;
    for (uint64_t s = 0; s < n_samples; s++) {
      double diff = data[s * n_features + f] - mean;
      var_sum += diff * diff;
    }
    double variance = var_sum / (double)n_samples;
    tk_rank_t r = { (int64_t)f, variance };
    #pragma omp critical
    tk_rvec_hmin(top_heap, top_k, r);
  }
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}

static inline tk_ivec_t *tk_dvec_mtx_top_skewness (
  lua_State *L,
  tk_dvec_t *matrix,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t top_k
) {
  if (matrix == NULL || n_samples == 0 || n_features == 0)
    return NULL;
  double *data = matrix->a;
  tk_rvec_t *top_heap = tk_rvec_create(0, 0);
  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    double sum = 0.0;
    for (uint64_t s = 0; s < n_samples; s++) {
      sum += data[s * n_features + f];
    }
    double mean = sum / (double)n_samples;
    double m2_sum = 0.0;
    double m3_sum = 0.0;
    for (uint64_t s = 0; s < n_samples; s++) {
      double diff = data[s * n_features + f] - mean;
      double diff2 = diff * diff;
      m2_sum += diff2;
      m3_sum += diff2 * diff;
    }
    double variance = m2_sum / (double)n_samples;
    double m3 = m3_sum / (double)n_samples;
    double skewness = 0.0;
    if (variance > 1e-10) {
      double std_dev = sqrt(variance);
      skewness = m3 / (std_dev * std_dev * std_dev);
    }
    tk_rank_t r = { (int64_t)f, -fabs(skewness) };
    #pragma omp critical
    tk_rvec_hmin(top_heap, top_k, r);
  }
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}

static inline tk_ivec_t *tk_dvec_mtx_top_entropy (
  lua_State *L,
  tk_dvec_t *matrix,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t top_k,
  uint64_t n_bins
) {
  if (matrix == NULL || n_samples == 0 || n_features == 0)
    return NULL;
  if (n_bins == 0)
    n_bins = 32;
  double *data = matrix->a;
  tk_rvec_t *top_heap = tk_rvec_create(0, 0);
  #pragma omp parallel
  {
    uint64_t *bin_counts = (uint64_t *)calloc(n_bins, sizeof(uint64_t));
    if (!bin_counts) {
      #pragma omp critical
      {
        tk_rvec_destroy(top_heap);
        luaL_error(L, "ESBER: failed to allocate bin buffer");
      }
    }
    #pragma omp for
    for (uint64_t f = 0; f < n_features; f++) {
      double min_val = data[0 * n_features + f];
      double max_val = min_val;
      for (uint64_t s = 1; s < n_samples; s++) {
        double val = data[s * n_features + f];
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
      }
      double range = max_val - min_val;
      if (range < 1e-10) {
        tk_rank_t r = { (int64_t)f, 0.0 };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
        continue;
      }
      memset(bin_counts, 0, n_bins * sizeof(uint64_t));
      for (uint64_t s = 0; s < n_samples; s++) {
        double val = data[s * n_features + f];
        double normalized = (val - min_val) / range;
        uint64_t bin_idx = (uint64_t)(normalized * (double)(n_bins - 1));
        if (bin_idx >= n_bins) bin_idx = n_bins - 1;
        bin_counts[bin_idx]++;
      }
      double entropy = 0.0;
      for (uint64_t b = 0; b < n_bins; b++) {
        if (bin_counts[b] > 0) {
          double p = (double)bin_counts[b] / (double)n_samples;
          entropy -= p * log2(p);
        }
      }
      tk_rank_t r = { (int64_t)f, entropy };
      #pragma omp critical
      tk_rvec_hmin(top_heap, top_k, r);
    }
    free(bin_counts);
  }
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}

static inline tk_ivec_t *tk_dvec_mtx_top_bimodality (
  lua_State *L,
  tk_dvec_t *matrix,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t top_k
) {
  if (matrix == NULL || n_samples == 0 || n_features == 0)
    return NULL;
  double *data = matrix->a;
  tk_rvec_t *top_heap = tk_rvec_create(0, 0);
  #pragma omp parallel for
  for (uint64_t f = 0; f < n_features; f++) {
    double sum = 0.0;
    for (uint64_t s = 0; s < n_samples; s++) {
      sum += data[s * n_features + f];
    }
    double mean = sum / (double)n_samples;
    double m2_sum = 0.0;
    double m3_sum = 0.0;
    double m4_sum = 0.0;
    for (uint64_t s = 0; s < n_samples; s++) {
      double diff = data[s * n_features + f] - mean;
      double diff2 = diff * diff;
      double diff3 = diff2 * diff;
      double diff4 = diff2 * diff2;
      m2_sum += diff2;
      m3_sum += diff3;
      m4_sum += diff4;
    }
    double variance = m2_sum / (double)n_samples;
    double m3 = m3_sum / (double)n_samples;
    double m4 = m4_sum / (double)n_samples;
    double bimodality = 0.0;
    if (variance > 1e-10 && n_samples > 3) {
      double std_dev = sqrt(variance);
      double skewness = m3 / (std_dev * std_dev * std_dev);
      double kurtosis = m4 / (variance * variance);
      double excess_kurtosis = kurtosis - 3.0;
      double n = (double)n_samples;
      double sample_size_correction = 3.0 * (n - 1.0) * (n - 1.0) / ((n - 2.0) * (n - 3.0));
      bimodality = (skewness * skewness + 1.0) / (excess_kurtosis + sample_size_correction);
    }
    tk_rank_t r = { (int64_t)f, -bimodality };
    #pragma omp critical
    tk_rvec_hmin(top_heap, top_k, r);
  }
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}

static inline tk_ivec_t *tk_dvec_mtx_top_dip (
  lua_State *L,
  tk_dvec_t *matrix,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t top_k
) {
  if (matrix == NULL || n_samples == 0 || n_features == 0 || n_samples < 5)
    return NULL;
  double *data = matrix->a;
  tk_rvec_t *top_heap = tk_rvec_create(0, 0);
  uint64_t n_bins = n_samples < 100 ? 16 : 32;
  #pragma omp parallel
  {
    uint64_t *bins = (uint64_t *)calloc(n_bins, sizeof(uint64_t));
    double *smoothed = (double *)calloc(n_bins, sizeof(double));
    if (!bins || !smoothed) {
      #pragma omp critical
      {
        tk_rvec_destroy(top_heap);
        if (bins) free(bins);
        if (smoothed) free(smoothed);
        luaL_error(L, "Multimodality test: failed to allocate buffers");
      }
    }
    #pragma omp for
    for (uint64_t f = 0; f < n_features; f++) {
      double min_val = data[0 * n_features + f];
      double max_val = min_val;
      for (uint64_t s = 1; s < n_samples; s++) {
        double val = data[s * n_features + f];
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
      }
      double range = max_val - min_val;
      if (range < 1e-10) {
        tk_rank_t r = { (int64_t)f, 0.0 };
        #pragma omp critical
        tk_rvec_hmin(top_heap, top_k, r);
        continue;
      }
      memset(bins, 0, n_bins * sizeof(uint64_t));
      for (uint64_t s = 0; s < n_samples; s++) {
        double val = data[s * n_features + f];
        double normalized = (val - min_val) / range;
        uint64_t bin_idx = (uint64_t)(normalized * (double)(n_bins - 1));
        if (bin_idx >= n_bins) bin_idx = n_bins - 1;
        bins[bin_idx]++;
      }
      for (uint64_t b = 0; b < n_bins; b++) {
        smoothed[b] = (double)bins[b] / (double)n_samples;
      }
      for (uint64_t pass = 0; pass < 2; pass++) {
        for (uint64_t b = 1; b < n_bins - 1; b++) {
          double val = (smoothed[b-1] + 2.0 * smoothed[b] + smoothed[b+1]) / 4.0;
          smoothed[b] = val;
        }
      }
      uint64_t n_local_min = 0;
      double deepest_valley = 0.0;
      for (uint64_t b = 1; b < n_bins - 1; b++) {
        if (smoothed[b] < smoothed[b-1] && smoothed[b] < smoothed[b+1]) {
          n_local_min++;
          double left_peak = smoothed[b-1];
          double right_peak = smoothed[b+1];
          for (int64_t i = (int64_t)b - 2; i >= 0; i--) {
            if (smoothed[i] > left_peak) left_peak = smoothed[i];
          }
          for (uint64_t i = b + 2; i < n_bins; i++) {
            if (smoothed[i] > right_peak) right_peak = smoothed[i];
          }
          double valley_depth = fmin(left_peak - smoothed[b], right_peak - smoothed[b]);
          if (valley_depth > deepest_valley) {
            deepest_valley = valley_depth;
          }
        }
      }
      double multimodality_score = deepest_valley * (double)n_local_min;
      tk_rank_t r = { (int64_t)f, multimodality_score };
      #pragma omp critical
      tk_rvec_hmin(top_heap, top_k, r);
    }
    free(bins);
    free(smoothed);
  }
  tk_rvec_desc(top_heap, 0, top_heap->n);
  tk_ivec_t *out = tk_ivec_create(L, 0);
  tk_dvec_t *weights = tk_dvec_create(L, 0);
  tk_rvec_keys(L, top_heap, out);
  tk_rvec_values(L, top_heap, weights);
  tk_rvec_destroy(top_heap);
  return out;
}

static inline void tk_dvec_mtx_topk (
  lua_State *L,
  tk_dvec_t *queries,
  tk_dvec_t *corpus,
  uint64_t n_queries,
  uint64_t n_corpus,
  uint64_t d,
  uint64_t k
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
  for (uint64_t i = 0; i <= n_queries; i++)
    offsets->a[i] = (int64_t)(i * k);
  tk_ivec_t *indices = tk_ivec_create(L, total);
  indices->n = total;
  tk_dvec_t *out_scores = tk_dvec_create(L, total);
  out_scores->n = total;
  uint64_t max_buf = 128ULL * 1024 * 1024 / sizeof(double);
  uint64_t tile = n_corpus > 0 ? max_buf / n_corpus : n_queries;
  if (tile == 0) tile = 1;
  if (tile > n_queries) tile = n_queries;
  double *sbuf = (double *)malloc(tile * n_corpus * sizeof(double));
  if (!sbuf) {
    luaL_error(L, "mtx_topk: out of memory");
    return;
  }
  for (uint64_t base = 0; base < n_queries; base += tile) {
    uint64_t blk = (base + tile <= n_queries) ? tile : n_queries - base;
#if !defined(__EMSCRIPTEN__)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
      (int)blk, (int)n_corpus, (int)d,
      1.0, queries->a + base * d, (int)d,
      corpus->a, (int)d,
      0.0, sbuf, (int)n_corpus);
#else
    for (uint64_t i = 0; i < blk; i++)
      for (uint64_t j = 0; j < n_corpus; j++) {
        double s = 0.0;
        const double *q = queries->a + (base + i) * d;
        const double *c = corpus->a + j * d;
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
        double *row = sbuf + i * n_corpus;
        for (uint64_t j = 0; j < n_corpus; j++)
          tk_rvec_hmin(heap, k, tk_rank((int64_t)j, row[j]));
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

static inline void tk_dvec_mtx_standardize (
  lua_State *L,
  tk_dvec_t *data,
  uint64_t n_cols,
  tk_dvec_t *mean_in,
  tk_dvec_t *istd_in,
  tk_dvec_t **mean_out,
  tk_dvec_t **istd_out
) {
  uint64_t N = data->n / n_cols;
  if (mean_in) {
    #pragma omp parallel for
    for (uint64_t d = 0; d < n_cols; d++) {
      double mu = mean_in->a[d];
      double is = istd_in->a[d];
      for (uint64_t s = 0; s < N; s++)
        data->a[s * n_cols + d] = (data->a[s * n_cols + d] - mu) * is;
    }
  } else {
    tk_dvec_t *mu = tk_dvec_create(L, n_cols);
    mu->n = n_cols;
    tk_dvec_t *is = tk_dvec_create(L, n_cols);
    is->n = n_cols;
    #pragma omp parallel for
    for (uint64_t d = 0; d < n_cols; d++) {
      double sum = 0, sum2 = 0;
      for (uint64_t s = 0; s < N; s++) {
        double v = data->a[s * n_cols + d];
        sum += v;
        sum2 += v * v;
      }
      double m = sum / (double)N;
      double var = sum2 / (double)N - m * m;
      double istd = var > 1e-24 ? 1.0 / sqrt(var) : 0.0;
      mu->a[d] = m;
      is->a[d] = istd;
      for (uint64_t s = 0; s < N; s++)
        data->a[s * n_cols + d] = (data->a[s * n_cols + d] - m) * istd;
    }
    *mean_out = mu;
    *istd_out = is;
  }
}

#endif
