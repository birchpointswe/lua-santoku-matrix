#ifndef tk_parallel_sfx
#error "Must include santoku/parallel/tpl.h before this template"
#endif

#if defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__)
  #include <immintrin.h>
#elif defined(__AVX2__)
  #include <immintrin.h>
#elif defined(__ARM_NEON)
  #include <arm_neon.h>
#endif

static inline uint64_t tk_cvec_bits_popcount_serial(const uint8_t *data, uint64_t n_bits);

static inline uint64_t tk_parallel_sfx(tk_cvec_bits_popcount) (
  const uint8_t * __restrict__ data,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  uint64_t main_bytes = full_bytes - (rem_bits > 0 ? 1 : 0);
  uint64_t count = 0;

  uint64_t n_chunks = main_bytes / 8;
  uint64_t remaining = main_bytes % 8;
  (void) remaining;
  const uint64_t *data64 = (const uint64_t *)data;

#if defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__)
  uint64_t vec_chunks = n_chunks / 8;

  TK_PARALLEL_FOR(reduction(+:count))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m512i vdata = _mm512_loadu_si512(&data64[v * 8]);
    __m512i vpop = _mm512_popcnt_epi64(vdata);
    count += (uint64_t) _mm512_reduce_add_epi64(vpop);
  }

  for (uint64_t i = vec_chunks * 8; i < n_chunks; i++)
    count += (uint64_t) __builtin_popcountll(data64[i]);

#elif defined(__AVX2__)
  uint64_t vec_chunks = n_chunks / 4;

  const __m256i lookup = _mm256_set_epi8(
    4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0,
    4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0
  );
  const __m256i low_mask = _mm256_set1_epi8(0x0F);
  const __m256i zero = _mm256_setzero_si256();

  TK_PARALLEL_FOR(reduction(+:count))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m256i vdata = _mm256_loadu_si256((__m256i*)&data64[v * 4]);

    __m256i v_low = _mm256_and_si256(vdata, low_mask);
    __m256i v_pop_low = _mm256_shuffle_epi8(lookup, v_low);
    __m256i v_high = _mm256_and_si256(_mm256_srli_epi16(vdata, 4), low_mask);
    __m256i v_pop_high = _mm256_shuffle_epi8(lookup, v_high);

    __m256i v_pop_bytes = _mm256_add_epi8(v_pop_low, v_pop_high);
    __m256i v_pop_sum = _mm256_sad_epu8(v_pop_bytes, zero);

    count += (uint64_t) _mm256_extract_epi64(v_pop_sum, 0);
    count += (uint64_t) _mm256_extract_epi64(v_pop_sum, 1);
    count += (uint64_t) _mm256_extract_epi64(v_pop_sum, 2);
    count += (uint64_t) _mm256_extract_epi64(v_pop_sum, 3);
  }

  for (uint64_t i = vec_chunks * 4; i < n_chunks; i++)
    count += (uint64_t) __builtin_popcountll(data64[i]);

#elif defined(__ARM_NEON)
  uint64_t vec_chunks = n_chunks / 2;

  TK_PARALLEL_FOR(reduction(+:count))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    uint64x2_t vdata = vld1q_u64(&data64[v * 2]);
    uint8x16_t vcnt = vcntq_u8(vreinterpretq_u8_u64(vdata));
    count += (uint64_t) vaddlvq_u8(vcnt);
  }

  for (uint64_t i = vec_chunks * 2; i < n_chunks; i++)
    count += (uint64_t) __builtin_popcountll(data64[i]);

#else
  TK_PARALLEL_FOR(reduction(+:count))
  for (uint64_t i = 0; i < n_chunks; i ++)
    count += (uint64_t) __builtin_popcountll(data64[i]);
#endif

  for (uint64_t i = n_chunks * 8; i < main_bytes; i ++)
    count += (uint64_t) __builtin_popcount(data[i]);

  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    count += (uint64_t) __builtin_popcount(data[full_bytes - 1] & mask);
  }
  return count;
}

static inline uint64_t tk_parallel_sfx(tk_cvec_bits_hamming) (
  const uint8_t * __restrict__ a,
  const uint8_t * __restrict__ b,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  uint64_t main_bytes = full_bytes - (rem_bits > 0 ? 1 : 0);
  uint64_t dist = 0;

  uint64_t n_chunks = main_bytes / 8;
  const uint64_t *a64 = (const uint64_t *)a;
  const uint64_t *b64 = (const uint64_t *)b;

#if defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__)
  uint64_t i = 0;
  uint64_t vec_chunks = n_chunks / 8;

  TK_PARALLEL_FOR(reduction(+:dist))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m512i va = _mm512_loadu_si512(&a64[v * 8]);
    __m512i vb = _mm512_loadu_si512(&b64[v * 8]);
    __m512i vxor = _mm512_xor_si512(va, vb);
    __m512i vpop = _mm512_popcnt_epi64(vxor);
    dist += (uint64_t) _mm512_reduce_add_epi64(vpop);
  }

  for (i = vec_chunks * 8; i < n_chunks; i++) {
    uint64_t xor_val = a64[i] ^ b64[i];
    dist += (uint64_t) __builtin_popcountll(xor_val);
  }

#elif defined(__AVX2__)
  uint64_t i = 0;
  uint64_t vec_chunks = n_chunks / 4;
  const __m256i lookup = _mm256_set_epi8(
    4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0, 4, 3, 3,
    2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0
  );
  const __m256i low_mask = _mm256_set1_epi8(0x0F);
  const __m256i zero = _mm256_setzero_si256();

  TK_PARALLEL_FOR(reduction(+:dist))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m256i va = _mm256_loadu_si256((__m256i*)&a64[v * 4]);
    __m256i vb = _mm256_loadu_si256((__m256i*)&b64[v * 4]);
    __m256i vxor = _mm256_xor_si256(va, vb);

    __m256i v_low = _mm256_and_si256(vxor, low_mask);
    __m256i v_pop_low = _mm256_shuffle_epi8(lookup, v_low);

    __m256i v_high = _mm256_and_si256(_mm256_srli_epi16(vxor, 4), low_mask);
    __m256i v_pop_high = _mm256_shuffle_epi8(lookup, v_high);

    __m256i v_pop_bytes = _mm256_add_epi8(v_pop_low, v_pop_high);

    __m256i v_pop_sum = _mm256_sad_epu8(v_pop_bytes, zero);

    dist += (uint64_t) _mm256_extract_epi64(v_pop_sum, 0);
    dist += (uint64_t) _mm256_extract_epi64(v_pop_sum, 1);
    dist += (uint64_t) _mm256_extract_epi64(v_pop_sum, 2);
    dist += (uint64_t) _mm256_extract_epi64(v_pop_sum, 3);
  }

  for (i = vec_chunks * 4; i < n_chunks; i++) {
    uint64_t xor_val = a64[i] ^ b64[i];
    dist += (uint64_t) __builtin_popcountll(xor_val);
  }

  for (uint64_t bi = n_chunks * 8; bi < main_bytes; bi++) {
    dist += (uint64_t) __builtin_popcount(a[bi] ^ b[bi]);
  }

#elif defined(__ARM_NEON)
  uint64_t i = 0;
  uint64_t vec_chunks = n_chunks / 2;

  TK_PARALLEL_FOR(reduction(+:dist))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    uint64x2_t va = vld1q_u64(&a64[v * 2]);
    uint64x2_t vb = vld1q_u64(&b64[v * 2]);
    uint64x2_t vxor = veorq_u64(va, vb);

    uint8x16_t vcnt = vcntq_u8(vreinterpretq_u8_u64(vxor));
    dist += (uint64_t) vaddlvq_u8(vcnt);
  }

  for (i = vec_chunks * 2; i < n_chunks; i++) {
    uint64_t xor_val = a64[i] ^ b64[i];
    dist += (uint64_t) __builtin_popcountll(xor_val);
  }

  for (uint64_t bi = n_chunks * 8; bi < main_bytes; bi++) {
    dist += (uint64_t) __builtin_popcount(a[bi] ^ b[bi]);
  }

#else
  TK_PARALLEL_FOR(reduction(+:dist))
  for (uint64_t i = 0; i < main_bytes; i++) {
    dist += (uint64_t) __builtin_popcount(a[i] ^ b[i]);
  }
#endif

  if (rem_bits > 0) {
    uint8_t x = a[full_bytes - 1] ^ b[full_bytes - 1];
    uint8_t mask = (1U << rem_bits) - 1;
    dist += (uint64_t) __builtin_popcount(x & mask);
  }
  return dist;
}

static inline uint64_t tk_parallel_sfx(tk_cvec_bits_hamming_mask) (
  const uint8_t * __restrict__ a,
  const uint8_t * __restrict__ b,
  const uint8_t * __restrict__ mask,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  uint64_t main_bytes = full_bytes - (rem_bits > 0 ? 1 : 0);
  uint64_t dist = 0;

  uint64_t n_chunks = main_bytes / 8;
  uint64_t remaining = main_bytes % 8;
  const uint64_t *a64 = (const uint64_t *)a;
  const uint64_t *b64 = (const uint64_t *)b;
  const uint64_t *mask64 = (const uint64_t *)mask;

#if defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__)
  uint64_t vec_chunks = n_chunks / 8;

  TK_PARALLEL_FOR(reduction(+:dist))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m512i va = _mm512_loadu_si512(&a64[v * 8]);
    __m512i vb = _mm512_loadu_si512(&b64[v * 8]);
    __m512i vm = _mm512_loadu_si512(&mask64[v * 8]);
    __m512i vxor = _mm512_xor_si512(va, vb);
    __m512i vmasked = _mm512_and_si512(vxor, vm);
    __m512i vpop = _mm512_popcnt_epi64(vmasked);
    dist += (uint64_t) _mm512_reduce_add_epi64(vpop);
  }

  for (uint64_t i = vec_chunks * 8; i < n_chunks; i++) {
    uint64_t masked_xor = (a64[i] ^ b64[i]) & mask64[i];
    dist += (uint64_t) __builtin_popcountll(masked_xor);
  }

  for (uint64_t i = n_chunks * 8; i < main_bytes; i++) {
    dist += (uint64_t) __builtin_popcount((a[i] ^ b[i]) & mask[i]);
  }

#elif defined(__AVX2__)
  uint64_t vec_chunks = n_chunks / 4;

  const __m256i lookup = _mm256_set_epi8(
    4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0,
    4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0
  );
  const __m256i low_mask = _mm256_set1_epi8(0x0F);
  const __m256i zero = _mm256_setzero_si256();

  TK_PARALLEL_FOR(reduction(+:dist))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m256i va = _mm256_loadu_si256((__m256i*)&a64[v * 4]);
    __m256i vb = _mm256_loadu_si256((__m256i*)&b64[v * 4]);
    __m256i vm = _mm256_loadu_si256((__m256i*)&mask64[v * 4]);

    __m256i vxor = _mm256_xor_si256(va, vb);
    __m256i v_to_pop = _mm256_and_si256(vxor, vm);

    __m256i v_low = _mm256_and_si256(v_to_pop, low_mask);
    __m256i v_pop_low = _mm256_shuffle_epi8(lookup, v_low);
    __m256i v_high = _mm256_and_si256(_mm256_srli_epi16(v_to_pop, 4), low_mask);
    __m256i v_pop_high = _mm256_shuffle_epi8(lookup, v_high);

    __m256i v_pop_bytes = _mm256_add_epi8(v_pop_low, v_pop_high);
    __m256i v_pop_sum = _mm256_sad_epu8(v_pop_bytes, zero);

    dist += (uint64_t) _mm256_extract_epi64(v_pop_sum, 0);
    dist += (uint64_t) _mm256_extract_epi64(v_pop_sum, 1);
    dist += (uint64_t) _mm256_extract_epi64(v_pop_sum, 2);
    dist += (uint64_t) _mm256_extract_epi64(v_pop_sum, 3);
  }

  for (uint64_t i = vec_chunks * 4; i < n_chunks; i++) {
    uint64_t masked_xor = (a64[i] ^ b64[i]) & mask64[i];
    dist += (uint64_t) __builtin_popcountll(masked_xor);
  }

#elif defined(__ARM_NEON)
  uint64_t vec_chunks = n_chunks / 2;

  TK_PARALLEL_FOR(reduction(+:dist))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    uint64x2_t va = vld1q_u64(&a64[v * 2]);
    uint64x2_t vb = vld1q_u64(&b64[v * 2]);
    uint64x2_t vm = vld1q_u64(&mask64[v * 2]);
    uint64x2_t vxor = veorq_u64(va, vb);
    uint64x2_t vmasked = vandq_u64(vxor, vm);

    uint8x16_t vcnt = vcntq_u8(vreinterpretq_u8_u64(vmasked));
    dist += (uint64_t) vaddlvq_u8(vcnt);
  }

  for (uint64_t i = vec_chunks * 2; i < n_chunks; i++) {
    uint64_t masked_xor = (a64[i] ^ b64[i]) & mask64[i];
    dist += (uint64_t) __builtin_popcountll(masked_xor);
  }

#else
  TK_PARALLEL_FOR(reduction(+:dist))
  for (uint64_t i = 0; i < n_chunks; i ++) {
    uint64_t masked_xor = (a64[i] ^ b64[i]) & mask64[i];
    dist += (uint64_t) __builtin_popcountll(masked_xor);
  }
#endif

  for (uint64_t i = n_chunks * 8; i < n_chunks * 8 + remaining; i ++)
    dist += (uint64_t) __builtin_popcount((a[i] ^ b[i]) & mask[i]);

  if (rem_bits > 0) {
    uint8_t x = a[full_bytes - 1] ^ b[full_bytes - 1];
    uint8_t m = mask[full_bytes - 1] & ((1U << rem_bits) - 1);
    dist += (uint64_t) __builtin_popcount(x & m);
  }
  return dist;
}

static inline void tk_parallel_sfx(tk_cvec_bits_popcount_andnot) (
  const uint8_t * __restrict__ a,
  const uint8_t * __restrict__ b,
  uint64_t n_bits,
  uint64_t * __restrict__ pop_a_out,
  uint64_t * __restrict__ pop_andnot_out
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  uint64_t main_bytes = full_bytes - (rem_bits > 0 ? 1 : 0);

  uint64_t pop_a = 0;
  uint64_t pop_andnot = 0;

  uint64_t n_chunks = main_bytes / 8;
  uint64_t remaining = main_bytes % 8;
  const uint64_t *a64 = (const uint64_t *)a;
  const uint64_t *b64 = (const uint64_t *)b;

#if defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__)
  uint64_t vec_chunks = n_chunks / 8;

  TK_PARALLEL_FOR(reduction(+:pop_a, pop_andnot))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m512i va = _mm512_loadu_si512(&a64[v * 8]);
    __m512i vb = _mm512_loadu_si512(&b64[v * 8]);

    __m512i vandnot = _mm512_andnot_si512(vb, va);

    __m512i vpop_a = _mm512_popcnt_epi64(va);
    __m512i vpop_andnot = _mm512_popcnt_epi64(vandnot);

    pop_a += (uint64_t) _mm512_reduce_add_epi64(vpop_a);
    pop_andnot += (uint64_t) _mm512_reduce_add_epi64(vpop_andnot);
  }

  for (uint64_t i = vec_chunks * 8; i < n_chunks; i++) {
    uint64_t va = a64[i];
    uint64_t vandnot = va & ~b64[i];
    pop_a += (uint64_t) __builtin_popcountll(va);
    pop_andnot += (uint64_t) __builtin_popcountll(vandnot);
  }

#elif defined(__AVX2__)
  uint64_t vec_chunks = n_chunks / 4;

  const __m256i lookup = _mm256_set_epi8(
      4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0,
      4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0
  );
  const __m256i low_mask = _mm256_set1_epi8(0x0F);
  const __m256i zero = _mm256_setzero_si256();

  TK_PARALLEL_FOR(reduction(+:pop_a, pop_andnot))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m256i va = _mm256_loadu_si256((__m256i*)&a64[v * 4]);
    __m256i vb = _mm256_loadu_si256((__m256i*)&b64[v * 4]);

    __m256i vandnot = _mm256_andnot_si256(vb, va);

    __m256i va_low = _mm256_and_si256(va, low_mask);
    __m256i va_pop_low = _mm256_shuffle_epi8(lookup, va_low);
    __m256i va_high = _mm256_and_si256(_mm256_srli_epi16(va, 4), low_mask);
    __m256i va_pop_high = _mm256_shuffle_epi8(lookup, va_high);
    __m256i va_pop_bytes = _mm256_add_epi8(va_pop_low, va_pop_high);
    __m256i va_pop_sum = _mm256_sad_epu8(va_pop_bytes, zero);

    __m256i vn_low = _mm256_and_si256(vandnot, low_mask);
    __m256i vn_pop_low = _mm256_shuffle_epi8(lookup, vn_low);
    __m256i vn_high = _mm256_and_si256(_mm256_srli_epi16(vandnot, 4), low_mask);
    __m256i vn_pop_high = _mm256_shuffle_epi8(lookup, vn_high);
    __m256i vn_pop_bytes = _mm256_add_epi8(vn_pop_low, vn_pop_high);
    __m256i vn_pop_sum = _mm256_sad_epu8(vn_pop_bytes, zero);

    pop_a += (uint64_t) _mm256_extract_epi64(va_pop_sum, 0) +
             (uint64_t) _mm256_extract_epi64(va_pop_sum, 1) +
             (uint64_t) _mm256_extract_epi64(va_pop_sum, 2) +
             (uint64_t) _mm256_extract_epi64(va_pop_sum, 3);

    pop_andnot += (uint64_t) _mm256_extract_epi64(vn_pop_sum, 0) +
                  (uint64_t) _mm256_extract_epi64(vn_pop_sum, 1) +
                  (uint64_t) _mm256_extract_epi64(vn_pop_sum, 2) +
                  (uint64_t) _mm256_extract_epi64(vn_pop_sum, 3);
  }

  for (uint64_t i = vec_chunks * 4; i < n_chunks; i++) {
    uint64_t va = a64[i];
    uint64_t vandnot = va & ~b64[i];
    pop_a += (uint64_t) __builtin_popcountll(va);
    pop_andnot += (uint64_t) __builtin_popcountll(vandnot);
  }

#elif defined(__ARM_NEON)
  uint64_t vec_chunks = n_chunks / 2;

  TK_PARALLEL_FOR(reduction(+:pop_a, pop_andnot))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    uint64x2_t va = vld1q_u64(&a64[v * 2]);
    uint64x2_t vb = vld1q_u64(&b64[v * 2]);

    uint64x2_t vandnot = vbicq_u64(va, vb);

    uint8x16_t vcnt_a = vcntq_u8(vreinterpretq_u8_u64(va));
    uint8x16_t vcnt_andnot = vcntq_u8(vreinterpretq_u8_u64(vandnot));

    pop_a += (uint64_t) vaddlvq_u8(vcnt_a);
    pop_andnot += (uint64_t) vaddlvq_u8(vcnt_andnot);
  }

  for (uint64_t i = vec_chunks * 2; i < n_chunks; i++) {
    uint64_t va = a64[i];
    uint64_t vandnot = va & ~b64[i];
    pop_a += (uint64_t) __builtin_popcountll(va);
    pop_andnot += (uint64_t) __builtin_popcountll(vandnot);
  }

#else
  TK_PARALLEL_FOR(reduction(+:pop_a, pop_andnot))
  for (uint64_t i = 0; i < n_chunks; i ++) {
    uint64_t va = a64[i];
    uint64_t vandnot = va & ~b64[i];
    pop_a += (uint64_t) __builtin_popcountll(va);
    pop_andnot += (uint64_t) __builtin_popcountll(vandnot);
  }
#endif

  for (uint64_t i = n_chunks * 8; i < n_chunks * 8 + remaining; i ++) {
    uint8_t va = a[i];
    uint8_t vandnot = va & ~b[i];
    pop_a += (uint64_t) __builtin_popcount(va);
    pop_andnot += (uint64_t) __builtin_popcount(vandnot);
  }

  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    uint8_t va = a[full_bytes - 1] & mask;
    uint8_t vandnot = va & ~b[full_bytes - 1];

    pop_a += (uint64_t) __builtin_popcount(va);
    pop_andnot += (uint64_t) __builtin_popcount(vandnot);
  }

  *pop_a_out = pop_a;
  *pop_andnot_out = pop_andnot;
}

static inline void tk_parallel_sfx(tk_cvec_bits_andnot) (
  uint8_t * __restrict__ out,
  const uint8_t * __restrict__ a,
  const uint8_t * __restrict__ b,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  uint64_t main_bytes = full_bytes - (rem_bits > 0 ? 1 : 0);

  uint64_t n_chunks = main_bytes / 8;
  uint64_t remaining = main_bytes % 8;
  uint64_t *out64 = (uint64_t *)out;
  const uint64_t *a64 = (const uint64_t *)a;
  const uint64_t *b64 = (const uint64_t *)b;

#if defined(__AVX512F__)
  uint64_t vec_chunks = n_chunks / 8;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m512i va = _mm512_loadu_si512(&a64[v * 8]);
    __m512i vb = _mm512_loadu_si512(&b64[v * 8]);
    __m512i vout = _mm512_andnot_si512(vb, va);
    _mm512_storeu_si512(&out64[v * 8], vout);
  }

  for (uint64_t i = vec_chunks * 8; i < n_chunks; i++)
    out64[i] = a64[i] & ~b64[i];

#elif defined(__AVX2__)
  uint64_t vec_chunks = n_chunks / 4;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m256i va = _mm256_loadu_si256((__m256i*)&a64[v * 4]);
    __m256i vb = _mm256_loadu_si256((__m256i*)&b64[v * 4]);
    __m256i vout = _mm256_andnot_si256(vb, va);
    _mm256_storeu_si256((__m256i*)&out64[v * 4], vout);
  }

  for (uint64_t i = vec_chunks * 4; i < n_chunks; i++)
    out64[i] = a64[i] & ~b64[i];

#elif defined(__ARM_NEON)
  uint64_t vec_chunks = n_chunks / 2;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    uint64x2_t va = vld1q_u64(&a64[v * 2]);
    uint64x2_t vb = vld1q_u64(&b64[v * 2]);
    uint64x2_t vout = vbicq_u64(va, vb);
    vst1q_u64(&out64[v * 2], vout);
  }

  for (uint64_t i = vec_chunks * 2; i < n_chunks; i++)
    out64[i] = a64[i] & ~b64[i];

#else
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n_chunks; i ++)
    out64[i] = a64[i] & ~b64[i];
#endif

  for (uint64_t i = n_chunks * 8; i < n_chunks * 8 + remaining; i ++)
    out[i] = a[i] & ~b[i];

  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    out[full_bytes - 1] = (a[full_bytes - 1] & ~b[full_bytes - 1]) & mask;
  }
}

static inline void tk_parallel_sfx(tk_cvec_bits_and) (
  uint8_t * __restrict__ out,
  const uint8_t * __restrict__ a,
  const uint8_t * __restrict__ b,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  uint64_t main_bytes = full_bytes - (rem_bits > 0 ? 1 : 0);

  uint64_t n_chunks = main_bytes / 8;
  uint64_t remaining = main_bytes % 8;
  uint64_t *out64 = (uint64_t *)out;
  const uint64_t *a64 = (const uint64_t *)a;
  const uint64_t *b64 = (const uint64_t *)b;

#if defined(__AVX512F__)
  uint64_t vec_chunks = n_chunks / 8;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m512i va = _mm512_loadu_si512(&a64[v * 8]);
    __m512i vb = _mm512_loadu_si512(&b64[v * 8]);
    __m512i vout = _mm512_and_si512(va, vb);
    _mm512_storeu_si512(&out64[v * 8], vout);
  }

  for (uint64_t i = vec_chunks * 8; i < n_chunks; i++)
    out64[i] = a64[i] & b64[i];

#elif defined(__AVX2__)
  uint64_t vec_chunks = n_chunks / 4;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m256i va = _mm256_loadu_si256((__m256i*)&a64[v * 4]);
    __m256i vb = _mm256_loadu_si256((__m256i*)&b64[v * 4]);
    __m256i vout = _mm256_and_si256(va, vb);
    _mm256_storeu_si256((__m256i*)&out64[v * 4], vout);
  }

  for (uint64_t i = vec_chunks * 4; i < n_chunks; i++)
    out64[i] = a64[i] & b64[i];

#elif defined(__ARM_NEON)
  uint64_t vec_chunks = n_chunks / 2;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    uint64x2_t va = vld1q_u64(&a64[v * 2]);
    uint64x2_t vb = vld1q_u64(&b64[v * 2]);
    uint64x2_t vout = vandq_u64(va, vb);
    vst1q_u64(&out64[v * 2], vout);
  }

  for (uint64_t i = vec_chunks * 2; i < n_chunks; i++)
    out64[i] = a64[i] & b64[i];

#else
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n_chunks; i ++)
    out64[i] = a64[i] & b64[i];
#endif

  for (uint64_t i = n_chunks * 8; i < n_chunks * 8 + remaining; i ++)
    out[i] = a[i] & b[i];

  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    out[full_bytes - 1] = (a[full_bytes - 1] & b[full_bytes - 1]) & mask;
  }
}

static inline void tk_parallel_sfx(tk_cvec_bits_or) (
  uint8_t * __restrict__ out,
  const uint8_t * __restrict__ a,
  const uint8_t * __restrict__ b,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  uint64_t main_bytes = full_bytes - (rem_bits > 0 ? 1 : 0);

  uint64_t n_chunks = main_bytes / 8;
  uint64_t remaining = main_bytes % 8;
  uint64_t *out64 = (uint64_t *)out;
  const uint64_t *a64 = (const uint64_t *)a;
  const uint64_t *b64 = (const uint64_t *)b;

#if defined(__AVX512F__)
  uint64_t vec_chunks = n_chunks / 8;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m512i va = _mm512_loadu_si512(&a64[v * 8]);
    __m512i vb = _mm512_loadu_si512(&b64[v * 8]);
    __m512i vout = _mm512_or_si512(va, vb);
    _mm512_storeu_si512(&out64[v * 8], vout);
  }

  for (uint64_t i = vec_chunks * 8; i < n_chunks; i++)
    out64[i] = a64[i] | b64[i];

#elif defined(__AVX2__)
  uint64_t vec_chunks = n_chunks / 4;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m256i va = _mm256_loadu_si256((__m256i*)&a64[v * 4]);
    __m256i vb = _mm256_loadu_si256((__m256i*)&b64[v * 4]);
    __m256i vout = _mm256_or_si256(va, vb);
    _mm256_storeu_si256((__m256i*)&out64[v * 4], vout);
  }

  for (uint64_t i = vec_chunks * 4; i < n_chunks; i++)
    out64[i] = a64[i] | b64[i];

#elif defined(__ARM_NEON)
  uint64_t vec_chunks = n_chunks / 2;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    uint64x2_t va = vld1q_u64(&a64[v * 2]);
    uint64x2_t vb = vld1q_u64(&b64[v * 2]);
    uint64x2_t vout = vorrq_u64(va, vb);
    vst1q_u64(&out64[v * 2], vout);
  }

  for (uint64_t i = vec_chunks * 2; i < n_chunks; i++)
    out64[i] = a64[i] | b64[i];

#else
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n_chunks; i ++)
    out64[i] = a64[i] | b64[i];
#endif

  for (uint64_t i = n_chunks * 8; i < n_chunks * 8 + remaining; i ++)
    out[i] = a[i] | b[i];

  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    out[full_bytes - 1] = (a[full_bytes - 1] | b[full_bytes - 1]) & mask;
  }
}

static inline void tk_parallel_sfx(tk_cvec_bits_xor) (
  uint8_t * __restrict__ out,
  const uint8_t * __restrict__ a,
  const uint8_t * __restrict__ b,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  uint64_t main_bytes = full_bytes - (rem_bits > 0 ? 1 : 0);

  uint64_t n_chunks = main_bytes / 8;
  uint64_t remaining = main_bytes % 8;
  uint64_t *out64 = (uint64_t *)out;
  const uint64_t *a64 = (const uint64_t *)a;
  const uint64_t *b64 = (const uint64_t *)b;

#if defined(__AVX512F__)
  uint64_t vec_chunks = n_chunks / 8;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m512i va = _mm512_loadu_si512(&a64[v * 8]);
    __m512i vb = _mm512_loadu_si512(&b64[v * 8]);
    __m512i vout = _mm512_xor_si512(va, vb);
    _mm512_storeu_si512(&out64[v * 8], vout);
  }

  for (uint64_t i = vec_chunks * 8; i < n_chunks; i++)
    out64[i] = a64[i] ^ b64[i];

#elif defined(__AVX2__)
  uint64_t vec_chunks = n_chunks / 4;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    __m256i va = _mm256_loadu_si256((__m256i*)&a64[v * 4]);
    __m256i vb = _mm256_loadu_si256((__m256i*)&b64[v * 4]);
    __m256i vout = _mm256_xor_si256(va, vb);
    _mm256_storeu_si256((__m256i*)&out64[v * 4], vout);
  }

  for (uint64_t i = vec_chunks * 4; i < n_chunks; i++)
    out64[i] = a64[i] ^ b64[i];

#elif defined(__ARM_NEON)
  uint64_t vec_chunks = n_chunks / 2;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t v = 0; v < vec_chunks; v++) {
    uint64x2_t va = vld1q_u64(&a64[v * 2]);
    uint64x2_t vb = vld1q_u64(&b64[v * 2]);
    uint64x2_t vout = veorq_u64(va, vb);
    vst1q_u64(&out64[v * 2], vout);
  }

  for (uint64_t i = vec_chunks * 2; i < n_chunks; i++)
    out64[i] = a64[i] ^ b64[i];

#else
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n_chunks; i ++)
    out64[i] = a64[i] ^ b64[i];
#endif

  for (uint64_t i = n_chunks * 8; i < n_chunks * 8 + remaining; i ++)
    out[i] = a[i] ^ b[i];

  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    out[full_bytes - 1] = (a[full_bytes - 1] ^ b[full_bytes - 1]) & mask;
  }
}

static inline tk_ivec_t *tk_parallel_sfx(tk_cvec_bits_to_ivec) (
  lua_State *L,
  tk_cvec_t *bitmap,
  uint64_t n_features
) {
  uint8_t *data = (uint8_t *)bitmap->a;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t n_samples = bitmap->n / bytes_per_sample;
  uint64_t total_bits = 0;
  TK_PARALLEL_FOR(reduction(+:total_bits))
  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t offset = s * bytes_per_sample;
    total_bits += tk_cvec_bits_popcount_serial(data + offset, n_features);
  }
  tk_ivec_t *out = tk_ivec_create(L, total_bits, 0, 0);
  uint64_t *sample_counts = (uint64_t *)calloc(n_samples, sizeof(uint64_t));
  if (!sample_counts) {
    tk_ivec_destroy(out);
    return NULL;
  }

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t offset = s * bytes_per_sample;
    sample_counts[s] = tk_cvec_bits_popcount_serial(data + offset, n_features);
  }

  uint64_t *offsets = (uint64_t *)malloc((n_samples + 1) * sizeof(uint64_t));
  if (!offsets) {
    free(sample_counts);
    tk_ivec_destroy(out);
    return NULL;
  }
  offsets[0] = 0;
  for (uint64_t s = 0; s < n_samples; s++)
    offsets[s + 1] = offsets[s] + sample_counts[s];

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t s = 0; s < n_samples; s++) {
    uint64_t offset = s * bytes_per_sample;
    uint64_t write_pos = offsets[s];
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t byte_idx = f / CHAR_BIT;
      uint8_t bit_idx = f % CHAR_BIT;
      if (data[offset + byte_idx] & (1u << bit_idx)) {
        out->a[write_pos++] = (int64_t) (s * n_features + f);
      }
    }
  }

  out->n = total_bits;
  free(sample_counts);
  free(offsets);
  return out;
}

static inline tk_cvec_t *tk_parallel_sfx(tk_cvec_bits_from_ivec) (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_features
) {
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  uint64_t total_bytes = n_samples * bytes_per_sample;

  tk_cvec_t *out = tk_cvec_create(L, total_bytes, 0, 0);
  uint8_t *data = (uint8_t *)out->a;
  memset(data, 0, total_bytes);

  TK_PARALLEL_FOR(schedule(static))
  for (size_t i = 0; i < set_bits->n; i++) {
    int64_t bit_idx = set_bits->a[i];
    uint64_t sample = (uint64_t) bit_idx / n_features;
    uint64_t feature = (uint64_t) bit_idx % n_features;

    if (bit_idx >= 0 && sample < n_samples && feature < n_features) {
      uint64_t offset = sample * bytes_per_sample;
      uint64_t byte_idx = feature / CHAR_BIT;
      uint8_t bit_pos = feature % CHAR_BIT;
      TK_ATOMIC
      data[offset + byte_idx] |= (1u << bit_pos);
    }
  }

  out->n = total_bytes;
  return out;
}

static inline int tk_parallel_sfx(tk_cvec_bits_extend_mapped) (
  tk_cvec_t *base,
  tk_cvec_t *ext,
  tk_ivec_t *aids,
  tk_ivec_t *bids,
  uint64_t n_base_features,
  uint64_t n_ext_features,
  bool project
) {
  uint64_t n_total_features = n_base_features + n_ext_features;
  uint64_t base_bytes_per_sample = TK_CVEC_BITS_BYTES(n_base_features);
  uint64_t ext_bytes_per_sample = TK_CVEC_BITS_BYTES(n_ext_features);
  uint64_t total_bytes_per_sample = TK_CVEC_BITS_BYTES(n_total_features);
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
  if (tk_cvec_ensure(base, final_n_samples * total_bytes_per_sample) != 0) {
    free(b_to_final);
    tk_iumap_destroy(a_id_to_pos);
    return -1;
  }
  uint8_t *base_data = (uint8_t *)base->a;
  uint8_t *ext_data = (uint8_t *)ext->a;
  tk_iumap_t *b_id_to_pos = tk_iumap_from_ivec(0, bids);
  if (!b_id_to_pos) {
    free(b_to_final);
    tk_iumap_destroy(a_id_to_pos);
    return -1;
  }
  uint8_t *new_data = calloc(final_n_samples, total_bytes_per_sample);
  if (!new_data) {
    free(b_to_final);
    tk_iumap_destroy(a_id_to_pos);
    tk_iumap_destroy(b_id_to_pos);
    return -1;
  }
  TK_PARALLEL_FOR(schedule(static))
  for (size_t ai = 0; ai < old_aids_n; ai++) {
    uint64_t dest_offset = ai * total_bytes_per_sample;
    uint64_t src_offset = ai * base_bytes_per_sample;
    memcpy(new_data + dest_offset, base_data + src_offset, base_bytes_per_sample);
    khint_t khi = tk_iumap_get(b_id_to_pos, aids->a[ai]);
    if (khi != tk_iumap_end(b_id_to_pos)) {
      int64_t b_idx = tk_iumap_val(b_id_to_pos, khi);
      uint64_t ext_src_offset = (uint64_t)b_idx * ext_bytes_per_sample;
      uint64_t ext_bit_offset = n_base_features;
      uint64_t ext_byte_offset = ext_bit_offset / CHAR_BIT;
      uint8_t ext_bit_shift = ext_bit_offset % CHAR_BIT;
      if (ext_bit_shift == 0) {
        memcpy(new_data + dest_offset + ext_byte_offset, ext_data + ext_src_offset, ext_bytes_per_sample);
      } else {
        uint8_t carry = 0;
        for (uint64_t i = 0; i < ext_bytes_per_sample; i++) {
          uint8_t byte = ext_data[ext_src_offset + i];
          new_data[dest_offset + ext_byte_offset + i] |= (byte << ext_bit_shift) | carry;
          carry = byte >> (CHAR_BIT - ext_bit_shift);
        }
        if (carry != 0 && ext_byte_offset + ext_bytes_per_sample < total_bytes_per_sample) {
          new_data[dest_offset + ext_byte_offset + ext_bytes_per_sample] |= carry;
        }
      }
    }
  }
  TK_PARALLEL_FOR(schedule(static))
  for (size_t bi = 0; bi < bids->n; bi++) {
    int64_t final_pos = b_to_final[bi];
    if (final_pos >= (int64_t)old_aids_n) {
      uint64_t dest_offset = (uint64_t)final_pos * total_bytes_per_sample;
      uint64_t src_offset = bi * ext_bytes_per_sample;
      uint64_t ext_bit_offset = n_base_features;
      uint64_t ext_byte_offset = ext_bit_offset / CHAR_BIT;
      uint8_t ext_bit_shift = ext_bit_offset % CHAR_BIT;
      if (ext_bit_shift == 0) {
        memcpy(new_data + dest_offset + ext_byte_offset, ext_data + src_offset, ext_bytes_per_sample);
      } else {
        uint8_t carry = 0;
        for (uint64_t i = 0; i < ext_bytes_per_sample; i++) {
          uint8_t byte = ext_data[src_offset + i];
          new_data[dest_offset + ext_byte_offset + i] |= (byte << ext_bit_shift) | carry;
          carry = byte >> (CHAR_BIT - ext_bit_shift);
        }
        if (carry != 0 && ext_byte_offset + ext_bytes_per_sample < total_bytes_per_sample) {
          new_data[dest_offset + ext_byte_offset + ext_bytes_per_sample] |= carry;
        }
      }
    }
  }
  memcpy(base_data, new_data, final_n_samples * total_bytes_per_sample);
  free(new_data);
  free(b_to_final);
  tk_iumap_destroy(a_id_to_pos);
  tk_iumap_destroy(b_id_to_pos);
  base->n = final_n_samples * total_bytes_per_sample;
  return 0;
}

static inline void tk_parallel_sfx(tk_cvec_bits_to_ascii) (
  lua_State *L,
  const tk_cvec_t *bitmap,
  uint64_t start_bit,
  uint64_t end_bit
) {
  uint64_t n_bits = end_bit - start_bit;
  tk_cvec_t *ascii = tk_cvec_create(L, n_bits + 1, 0, 0);
  const uint8_t *data = (const uint8_t *)bitmap->a;
  char *out = ascii->a;
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n_bits; i++) {
    uint64_t bit_pos = start_bit + i;
    uint64_t byte_idx = bit_pos / 8;
    uint64_t bit_idx = bit_pos % 8;
    out[i] = (data[byte_idx] & (1 << bit_idx)) ? '1' : '0';
  }
  ascii->n = n_bits;
  out[n_bits] = '\0';
  lua_pushstring(L, out);
  lua_remove(L, -2);
}
