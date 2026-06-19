#ifndef TK_IVEC_EXT_H
#define TK_IVEC_EXT_H

#if defined(_OPENMP) && !defined(__EMSCRIPTEN__)
#include <omp.h>
#endif
#include <santoku/rvec.h>
#include <santoku/dvec.h>
#include <santoku/iumap.h>
#include <stdatomic.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifndef TK_CVEC_BITS_BYTES
#define TK_CVEC_BITS_BYTES(n) (((n) + CHAR_BIT - 1) / CHAR_BIT)
#endif

static inline tk_dvec_t *tk_ivec_to_dvec (lua_State *L, tk_ivec_t *v) {
  tk_dvec_t *out = tk_dvec_create(L, v->n);
  for (uint64_t i = 0; i < v->n; i++)
    out->a[i] = (double)v->a[i];
  return out;
}

static inline tk_svec_t *tk_ivec_to_svec (lua_State *L, tk_ivec_t *v) {
  tk_svec_t *out = tk_svec_create(L, v->n);
  for (uint64_t i = 0; i < v->n; i++)
    out->a[i] = (int32_t)v->a[i];
  return out;
}

static inline void tk_ivec_lookup (tk_ivec_t *indices, tk_ivec_t *source) {
  int64_t write_pos = 0;
  for (uint64_t i = 0; i < indices->n; i++) {
    int64_t idx = indices->a[i];
    if (idx >= 0 && idx < (int64_t) source->n) {
      indices->a[write_pos++] = source->a[idx];
    }
  }
  indices->n = (uint64_t) write_pos;
}

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

static inline int64_t tk_ivec_set_find (int64_t *arr, int64_t start, int64_t end, int64_t value) {
  int64_t left = start;
  int64_t right = end - 1;
  while (left <= right) {
    int64_t mid = left + (right - left) / 2;
    if (arr[mid] == value)
      return mid;
    if (arr[mid] < value)
      left = mid + 1;
    else
      right = mid - 1;
  }
  return -(left + 1);
}
static inline void tk_ivec_set_insert (tk_ivec_t *vec, int64_t pos, int64_t value) {
  if (pos < 0) pos = -(pos + 1);
  if (pos < 0) pos = 0;
  if (pos > (int64_t)vec->n) pos = (int64_t)vec->n;
  tk_ivec_ensure(vec, vec->n + 1);
  if (pos < (int64_t)vec->n) {
    memmove(vec->a + pos + 1, vec->a + pos, (vec->n - (size_t) pos) * sizeof(int64_t));
  }
  vec->a[pos] = value;
  vec->n++;
}

static inline tk_ivec_t *tk_ivec_set_intersect (lua_State *L, tk_ivec_t *a, tk_ivec_t *b, tk_ivec_t *out) {
  if (out == NULL) {
    size_t min_size = a->n < b->n ? a->n : b->n;
    out = tk_ivec_create(L, min_size);
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

static inline tk_ivec_t *tk_ivec_set_union (lua_State *L, tk_ivec_t *a, tk_ivec_t *b, tk_ivec_t *out) {
  if (out == NULL) {
    out = tk_ivec_create(L, a->n + b->n);
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


#endif
