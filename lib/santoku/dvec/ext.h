#ifndef TK_DVEC_EXT_H
#define TK_DVEC_EXT_H

#include <santoku/dvec/base.h>
#include <santoku/ivec.h>

static inline void tk_dvec_center (
  double *M,
  size_t N,
  size_t K
) {
  for (size_t j = 0; j < K; j ++) {
    double mu = 0.0;
    for (size_t i = 0; i < N; i ++) mu += M[i * K + j];
    mu /= (double)N;
    for (size_t i = 0; i < N; i ++) M[i * K + j] -= mu;
  }
}

static inline void tk_dvec_rnorml2 (
  double *M,
  size_t N,
  size_t K
) {
  for (size_t i = 0; i < N; i++) {
    double sum = 0.0;
    for (size_t j = 0; j < K; j++) {
      double val = M[i * K + j];
      sum += val * val;
    }
    double norm = sqrt(sum);
    if (norm > 0.0) {
      for (size_t j = 0; j < K; j++) {
        M[i * K + j] /= norm;
      }
    }
  }
}

static inline tk_dvec_t *tk_dvec_multiply_bits (
  lua_State *L,
  tk_dvec_t *P,
  tk_ivec_t *raw_features,
  uint64_t n_samples,
  uint64_t n_features,
  uint64_t n_hidden
) {
  tk_dvec_t *V = tk_dvec_create(L, n_samples * n_hidden, 0, 0);
  tk_dvec_zero(V);
  for (size_t i = 0; i < raw_features->n; i ++) {
    int64_t b = raw_features->a[i];
    int64_t s = b / (int64_t) n_features;
    int64_t f = b % (int64_t) n_features;
    double *pa = P->a + f * (int64_t) n_hidden;
    double *va = V->a + s * (int64_t) n_hidden;
    for (uint64_t k = 0; k < n_hidden; k ++)
      va[k] += pa[k];
  }
  return V;
}



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
    double slope1 = (sum_xy1 - (double)n1 * mean_x1 * mean_y1) /
                    (sum_xx1 - (double)n1 * mean_x1 * mean_x1 + 1e-10);
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
    double slope2 = (sum_xy2 - (double)n2 * mean_x2 * mean_y2) /
                    (sum_xx2 - (double)n2 * mean_x2 * mean_x2 + 1e-10);
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

#endif
