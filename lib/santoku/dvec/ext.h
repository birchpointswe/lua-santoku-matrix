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

#endif
