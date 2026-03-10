#ifndef tk_parallel_sfx
#error "Must include santoku/parallel/tpl.h before this template"
#endif

static inline void tk_parallel_sfx(tk_fvec_center) (
  float *M,
  size_t N,
  size_t K
) {
  TK_PARALLEL_FOR(schedule(static))
  for (size_t j = 0; j < K; j ++) {
    float mu = 0.0f;
    for (size_t i = 0; i < N; i ++) mu += M[i * K + j];
    mu /= (float)N;
    for (size_t i = 0; i < N; i ++) M[i * K + j] -= mu;
  }
}

static inline void tk_parallel_sfx(tk_fvec_rnorml2) (
  float *M,
  size_t N,
  size_t K
) {
  TK_PARALLEL_FOR(schedule(static))
  for (size_t i = 0; i < N; i++) {
    float *row = M + i * K;
#if !defined(__EMSCRIPTEN__)
    float norm = cblas_snrm2(K, row, 1);
    if (norm > 0.0f) {
      cblas_sscal(K, 1.0f / norm, row, 1);
    }
#else
    float norm = 0.0f;
    for (size_t j = 0; j < K; j++) norm += row[j] * row[j];
    norm = sqrtf(norm);
    if (norm > 0.0f) {
      float inv = 1.0f / norm;
      for (size_t j = 0; j < K; j++) row[j] *= inv;
    }
#endif
  }
}
