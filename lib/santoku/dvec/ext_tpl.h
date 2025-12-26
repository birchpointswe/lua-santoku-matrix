


#ifndef tk_parallel_sfx
#error "Must include santoku/parallel/tpl.h before this template"
#endif

static inline void tk_parallel_sfx(tk_dvec_center) (
  double *M,
  size_t N,
  size_t K
) {

  TK_PARALLEL_FOR(schedule(static))
  for (size_t j = 0; j < K; j ++) {
    double mu = 0.0;
    for (size_t i = 0; i < N; i ++) mu += M[i * K + j];
    mu /= (double)N;
    for (size_t i = 0; i < N; i ++) M[i * K + j] -= mu;
  }
}

static inline void tk_parallel_sfx(tk_dvec_rnorml2) (
  double *M,
  size_t N,
  size_t K
) {
  TK_PARALLEL_FOR(schedule(static))
  for (size_t i = 0; i < N; i++) {
    double *row = M + i * K;
#if !defined(__EMSCRIPTEN__)
    double norm = cblas_dnrm2(K, row, 1);
    if (norm > 0.0) {
      cblas_dscal(K, 1.0 / norm, row, 1);
    }
#else
    double norm = 0.0;
    for (size_t j = 0; j < K; j++) norm += row[j] * row[j];
    norm = sqrt(norm);
    if (norm > 0.0) {
      double inv = 1.0 / norm;
      for (size_t j = 0; j < K; j++) row[j] *= inv;
    }
#endif
  }
}
