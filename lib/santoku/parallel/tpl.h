#undef TK_PARALLEL
#undef TK_PARALLEL_FOR
#undef TK_FOR
#undef TK_BARRIER
#undef TK_SINGLE_BEGIN
#undef TK_SINGLE_END
#undef TK_CRITICAL_BEGIN
#undef TK_CRITICAL_END
#undef TK_ATOMIC
#undef tk_parallel_sfx
#undef tk_parallel_sfx_helper

#ifndef TK_STR_HELPER
#define TK_STR_HELPER(...) #__VA_ARGS__
#define TK_STR(...) TK_STR_HELPER(__VA_ARGS__)
#endif

#ifdef TK_GENERATE_SINGLE
  #define TK_PARALLEL
  #define TK_PARALLEL_FOR(...)
  #define TK_FOR(...)
  #define TK_BARRIER
  #define TK_SINGLE_BEGIN
  #define TK_SINGLE_END
  #define TK_CRITICAL_BEGIN
  #define TK_CRITICAL_END
  #define TK_ATOMIC
  #define tk_parallel_sfx_helper(base) base##_serial
  #define tk_parallel_sfx(base) tk_parallel_sfx_helper(base)
#else
  #define TK_PARALLEL _Pragma("omp parallel")
  #define TK_PARALLEL_FOR(...) _Pragma(TK_STR(omp parallel for __VA_ARGS__))
  #define TK_FOR(...) _Pragma(TK_STR(omp for __VA_ARGS__))
  #define TK_BARRIER _Pragma("omp barrier")
  #define TK_SINGLE_BEGIN _Pragma("omp single") {
  #define TK_SINGLE_END }
  #define TK_CRITICAL_BEGIN _Pragma("omp critical") {
  #define TK_CRITICAL_END }
  #define TK_ATOMIC _Pragma("omp atomic")
  #define tk_parallel_sfx(base) base
#endif
