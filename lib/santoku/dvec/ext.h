#ifndef TK_DVEC_EXT_H
#define TK_DVEC_EXT_H

#include <lapacke.h>
#include <cblas.h>
#include <omp.h>

#include <santoku/dvec/base.h>
#include <santoku/ivec.h>

#define TK_GENERATE_SINGLE
#include <santoku/parallel/tpl.h>
#include <santoku/dvec/ext_tpl.h>
#undef TK_GENERATE_SINGLE

#include <santoku/parallel/tpl.h>
#include <santoku/dvec/ext_tpl.h>

static inline size_t tk_dvec_scores_kaiser (
  double *scores,
  size_t n,
  double *out_val
) {
  if (n == 0) {
    if (out_val) *out_val = 0.0;
    return 0;
  }
  if (n == 1) {
    if (out_val) *out_val = scores[0];
    return 1;
  }
  double sum = 0.0;
  for (size_t i = 0; i < n; i++)
    sum += scores[i];
  double mean = sum / (double)n;
  for (size_t i = 0; i < n; i++) {
    if (scores[i] < mean) {
      if (out_val) *out_val = scores[i];
      return i;
    }
  }
  if (out_val) *out_val = (n > 0) ? scores[n-1] : 0.0;
  return n;
}

static inline size_t tk_dvec_scores_max_curvature (
  double *scores,
  size_t n,
  double *out_val
) {
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? scores[0] : 0.0;
    return 0;
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
    return 0;
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
    return 0;
  }
  double max_gap = 0.0;
  size_t max_idx = 0;
  for (size_t i = 0; i < n - 1; i++) {
    double gap = scores[i + 1] - scores[i];
    if (gap > max_gap) {
      max_gap = gap;
      max_idx = i;
    }
  }
  if (out_val) *out_val = scores[max_idx];
  return max_idx;
}

static inline size_t tk_dvec_scores_max_drop (
  double *scores,
  size_t n,
  double *out_val
) {
  if (n < 2) {
    if (out_val) *out_val = (n > 0) ? scores[0] : 0.0;
    return 0;
  }
  double max_drop = 0.0;
  size_t max_idx = 0;
  for (size_t i = 0; i < n - 1; i++) {
    double drop = scores[i] - scores[i + 1];
    if (drop > max_drop) {
      max_drop = drop;
      max_idx = i;
    }
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
  double base = scores[0];
  size_t end_idx = 0;
  for (size_t i = 1; i < n; i++) {
    if (fabs(scores[i] - base) <= tolerance) {
      end_idx = i;
    } else {
      break;
    }
  }
  if (out_val) *out_val = scores[end_idx];
  return end_idx;
}

static inline size_t tk_dvec_scores_max_acceleration (
  double *scores,
  size_t n,
  double *out_val
) {
  if (n < 3) {
    if (out_val) *out_val = (n > 0) ? scores[0] : 0.0;
    return 0;
  }
  double max_accel = -DBL_MAX;
  size_t max_idx = 0;
  for (size_t i = 0; i < n - 2; i++) {
    double gap1 = scores[i + 1] - scores[i];
    double gap2 = scores[i + 2] - scores[i + 1];
    double accel = gap2 - gap1;
    if (accel > max_accel) {
      max_accel = accel;
      max_idx = i;
    }
  }
  if (out_val) *out_val = scores[max_idx];
  return max_idx;
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

static inline void tk_dvec_gemv(
  bool transpose,
  uint64_t rows,
  uint64_t cols,
  double alpha,
  double *A,
  double *x,
  double beta,
  double *y
) {
  cblas_dgemv(CblasRowMajor, transpose ? CblasTrans : CblasNoTrans, rows, cols, alpha, A, cols, x, 1, beta, y, 1);
}

static inline void tk_dvec_gemm(
  bool transpose_a,
  bool transpose_b,
  uint64_t m,
  uint64_t n,
  uint64_t k,
  double alpha,
  double *A,
  double *B,
  double beta,
  double *C
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

static inline double tk_dvec_dot_override(tk_dvec_t *a, tk_dvec_t *b) {
  uint64_t n = a->n < b->n ? a->n : b->n;
  return cblas_ddot(n, a->a, 1, b->a, 1);
}

static inline void tk_dvec_scale_override(tk_dvec_t *v, double scale, uint64_t start, uint64_t end) {
  if (end > v->n) {
    tk_dvec_ensure(v, end);
    v->n = end;
  }
  if (end <= start) return;
  cblas_dscal(end - start, scale, v->a + start, 1);
}

static inline void tk_dvec_addv_override(tk_dvec_t *a, tk_dvec_t *b, uint64_t start, uint64_t end) {
  if (end > a->n) {
    tk_dvec_ensure(a, end);
    a->n = end;
  }
  if (end > b->n || end <= start) return;
  cblas_daxpy(end - start, 1.0, b->a + start, 1, a->a + start, 1);
}

static inline void tk_dvec_multiply_override(tk_dvec_t *a, tk_dvec_t *b, tk_dvec_t *c, uint64_t k, bool transpose_a, bool transpose_b) {
  size_t m = transpose_a ? k : a->n / k;
  size_t n = transpose_b ? k : b->n / k;
  tk_dvec_ensure(c, m * n);
  c->n = m * n;
  cblas_dgemm(CblasRowMajor,
              transpose_a ? CblasTrans : CblasNoTrans,
              transpose_b ? CblasTrans : CblasNoTrans,
              m, n, k,
              1.0, a->a, transpose_a ? m : k,
              b->a, transpose_b ? k : n,
              0.0, c->a, n);
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
  tk_dvec_t *out = tk_dvec_create(L, rows, 0, 0);
  for (uint64_t r = 0; r < rows; r++) {
    out->a[r] = cblas_dnrm2(cols, m0->a + r * cols, 1);
  }
  out->n = rows;
  return out;
}

static inline tk_dvec_t *tk_dvec_cmags_override(lua_State *L, tk_dvec_t *m0, uint64_t cols) {
  uint64_t rows = m0->n / cols;
  tk_dvec_t *out = tk_dvec_create(L, cols, 0, 0);
  for (uint64_t c = 0; c < cols; c++) {
    out->a[c] = cblas_dnrm2(rows, m0->a + c, cols);
  }
  out->n = cols;
  return out;
}

static inline tk_dvec_t *tk_dvec_rsums_override(lua_State *L, tk_dvec_t *m0, uint64_t cols) {
  uint64_t rows = m0->n / cols;
  tk_dvec_t *out = tk_dvec_create(L, rows, 0, 0);
  tk_dvec_t *ones = tk_dvec_create(NULL, cols, 0, 0);
  for (uint64_t i = 0; i < cols; i++) ones->a[i] = 1.0;
  ones->n = cols;
  cblas_dgemv(CblasRowMajor, CblasNoTrans, rows, cols, 1.0, m0->a, cols, ones->a, 1, 0.0, out->a, 1);
  out->n = rows;
  tk_dvec_destroy(ones);
  return out;
}

static inline tk_dvec_t *tk_dvec_csums_override(lua_State *L, tk_dvec_t *m0, uint64_t cols) {
  uint64_t rows = m0->n / cols;
  tk_dvec_t *out = tk_dvec_create(L, cols, 0, 0);
  tk_dvec_t *ones = tk_dvec_create(NULL, rows, 0, 0);
  for (uint64_t i = 0; i < rows; i++) ones->a[i] = 1.0;
  ones->n = rows;
  cblas_dgemv(CblasRowMajor, CblasTrans, rows, cols, 1.0, m0->a, cols, ones->a, 1, 0.0, out->a, 1);
  out->n = cols;
  tk_dvec_destroy(ones);
  return out;
}

#include <santoku/iumap.h>

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

#endif
