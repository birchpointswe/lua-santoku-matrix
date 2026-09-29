#ifndef tk_parallel_sfx
#error "Must include santoku/parallel/tpl.h before this template"
#endif

static inline uint64_t tk_parallel_sfx(tk_cvec_bits_popcount) (
  const uint8_t * __restrict__ data,
  uint64_t n_bits
) {
  uint64_t full_bytes = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t rem_bits = TK_CVEC_BITS_BIT(n_bits);
  uint64_t main_bytes = full_bytes - (rem_bits > 0 ? 1 : 0);

  uint64_t count = 0;

#if defined(__AVX512VPOPCNTDQ__)
  uint64_t n512 = main_bytes / 64;
  __m512i vacc = _mm512_setzero_si512();
  for (uint64_t i = 0; i < n512; i++)
    vacc = _mm512_add_epi64(vacc, _mm512_popcnt_epi64(_mm512_loadu_si512(data + i * 64)));
  count = (uint64_t)_mm512_reduce_add_epi64(vacc);
  for (uint64_t i = n512 * 64; i < main_bytes; i++)
    count += (uint64_t)__builtin_popcount(data[i]);
#elif defined(__aarch64__)
  uint64_t n16 = main_bytes / 16, ni = 0;
  while (ni < n16) {
    uint16x8_t vacc = vdupq_n_u16(0);
    uint64_t end = ni + 4095; if (end > n16) end = n16;
    for (; ni < end; ni++)
      vacc = vpadalq_u8(vacc, vcntq_u8(vld1q_u8(data + ni * 16)));
    count += vaddlvq_u16(vacc);
  }
  for (uint64_t i = n16 * 16; i < main_bytes; i++)
    count += (uint64_t)__builtin_popcount(data[i]);
#else
  uint64_t c0 = 0, c1 = 0, c2 = 0, c3 = 0;
  uint64_t n64 = main_bytes / 8;
  uint64_t n64_4 = (n64 / 4) * 4;

  TK_PARALLEL_FOR(reduction(+:c0,c1,c2,c3))
  for (uint64_t i = 0; i < n64_4; i += 4) {
    uint64_t d0, d1, d2, d3;
    memcpy(&d0, &data[(i + 0) * 8], 8);
    memcpy(&d1, &data[(i + 1) * 8], 8);
    memcpy(&d2, &data[(i + 2) * 8], 8);
    memcpy(&d3, &data[(i + 3) * 8], 8);
    c0 += (uint64_t)__builtin_popcountll(d0);
    c1 += (uint64_t)__builtin_popcountll(d1);
    c2 += (uint64_t)__builtin_popcountll(d2);
    c3 += (uint64_t)__builtin_popcountll(d3);
  }
  count = c0 + c1 + c2 + c3;

  for (uint64_t i = n64_4; i < n64; i++) {
    uint64_t chunk;
    memcpy(&chunk, &data[i * 8], 8);
    count += (uint64_t)__builtin_popcountll(chunk);
  }

  uint64_t offset = n64 * 8;
  for (uint64_t i = offset; i < main_bytes; i++)
    count += (uint64_t)__builtin_popcount(data[i]);
#endif

  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    count += (uint64_t)__builtin_popcount(data[full_bytes - 1] & mask);
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

#if defined(__AVX512VPOPCNTDQ__)
  uint64_t n512 = main_bytes / 64;
  __m512i vacc = _mm512_setzero_si512();
  for (uint64_t i = 0; i < n512; i++)
    vacc = _mm512_add_epi64(vacc, _mm512_popcnt_epi64(_mm512_xor_si512(
      _mm512_loadu_si512(a + i * 64), _mm512_loadu_si512(b + i * 64))));
  dist = (uint64_t)_mm512_reduce_add_epi64(vacc);
  for (uint64_t i = n512 * 64; i < main_bytes; i++)
    dist += (uint64_t)__builtin_popcount(a[i] ^ b[i]);
#elif defined(__aarch64__)
  uint64_t n16 = main_bytes / 16, ni = 0;
  while (ni < n16) {
    uint16x8_t vacc = vdupq_n_u16(0);
    uint64_t end = ni + 4095; if (end > n16) end = n16;
    for (; ni < end; ni++)
      vacc = vpadalq_u8(vacc, vcntq_u8(veorq_u8(vld1q_u8(a + ni * 16), vld1q_u8(b + ni * 16))));
    dist += vaddlvq_u16(vacc);
  }
  for (uint64_t i = n16 * 16; i < main_bytes; i++)
    dist += (uint64_t)__builtin_popcount(a[i] ^ b[i]);
#else
  uint64_t d0 = 0, d1 = 0, d2 = 0, d3 = 0;
  uint64_t n64 = main_bytes / 8;
  uint64_t n64_4 = (n64 / 4) * 4;

  TK_PARALLEL_FOR(reduction(+:d0,d1,d2,d3))
  for (uint64_t i = 0; i < n64_4; i += 4) {
    uint64_t a0, a1, a2, a3, b0, b1, b2, b3;
    memcpy(&a0, &a[(i + 0) * 8], 8);
    memcpy(&a1, &a[(i + 1) * 8], 8);
    memcpy(&a2, &a[(i + 2) * 8], 8);
    memcpy(&a3, &a[(i + 3) * 8], 8);
    memcpy(&b0, &b[(i + 0) * 8], 8);
    memcpy(&b1, &b[(i + 1) * 8], 8);
    memcpy(&b2, &b[(i + 2) * 8], 8);
    memcpy(&b3, &b[(i + 3) * 8], 8);
    d0 += (uint64_t)__builtin_popcountll(a0 ^ b0);
    d1 += (uint64_t)__builtin_popcountll(a1 ^ b1);
    d2 += (uint64_t)__builtin_popcountll(a2 ^ b2);
    d3 += (uint64_t)__builtin_popcountll(a3 ^ b3);
  }
  dist = d0 + d1 + d2 + d3;

  for (uint64_t i = n64_4; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 8], 8);
    memcpy(&b_chunk, &b[i * 8], 8);
    dist += (uint64_t)__builtin_popcountll(a_chunk ^ b_chunk);
  }

  uint64_t offset = n64 * 8;
  for (uint64_t i = offset; i < main_bytes; i++)
    dist += (uint64_t)__builtin_popcount(a[i] ^ b[i]);
#endif

  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    dist += (uint64_t)__builtin_popcount((a[full_bytes - 1] ^ b[full_bytes - 1]) & mask);
  }

  return dist;
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

#if defined(__AVX512F__)
  uint64_t n512 = main_bytes / 64;
  for (uint64_t i = 0; i < n512; i++)
    _mm512_storeu_si512(out + i * 64, _mm512_andnot_si512(
      _mm512_loadu_si512(b + i * 64), _mm512_loadu_si512(a + i * 64)));
  uint64_t offset = n512 * 64;
#elif defined(__SIZEOF_INT128__)
  uint64_t n128 = main_bytes / 16;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n128; i++) {
    __uint128_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 16], sizeof(__uint128_t));
    memcpy(&b_chunk, &b[i * 16], sizeof(__uint128_t));
    __uint128_t result = a_chunk & ~b_chunk;
    memcpy(&out[i * 16], &result, sizeof(__uint128_t));
  }

  uint64_t offset = n128 * 16;
  uint64_t n64 = (main_bytes - offset) / 8;

  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[offset + i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[offset + i * 8], sizeof(uint64_t));
    uint64_t result = a_chunk & ~b_chunk;
    memcpy(&out[offset + i * 8], &result, sizeof(uint64_t));
  }

  offset += n64 * 8;
#else
  uint64_t n64 = main_bytes / 8;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[i * 8], sizeof(uint64_t));
    uint64_t result = a_chunk & ~b_chunk;
    memcpy(&out[i * 8], &result, sizeof(uint64_t));
  }

  uint64_t offset = n64 * 8;
#endif

  for (uint64_t i = offset; i < main_bytes; i++)
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

#ifdef __SIZEOF_INT128__
  uint64_t n128 = main_bytes / 16;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n128; i++) {
    __uint128_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 16], sizeof(__uint128_t));
    memcpy(&b_chunk, &b[i * 16], sizeof(__uint128_t));
    __uint128_t result = a_chunk & b_chunk;
    memcpy(&out[i * 16], &result, sizeof(__uint128_t));
  }

  uint64_t offset = n128 * 16;
  uint64_t n64 = (main_bytes - offset) / 8;

  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[offset + i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[offset + i * 8], sizeof(uint64_t));
    uint64_t result = a_chunk & b_chunk;
    memcpy(&out[offset + i * 8], &result, sizeof(uint64_t));
  }

  offset += n64 * 8;
#else
  uint64_t n64 = main_bytes / 8;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[i * 8], sizeof(uint64_t));
    uint64_t result = a_chunk & b_chunk;
    memcpy(&out[i * 8], &result, sizeof(uint64_t));
  }

  uint64_t offset = n64 * 8;
#endif

  for (uint64_t i = offset; i < main_bytes; i++)
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

#ifdef __SIZEOF_INT128__
  uint64_t n128 = main_bytes / 16;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n128; i++) {
    __uint128_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 16], sizeof(__uint128_t));
    memcpy(&b_chunk, &b[i * 16], sizeof(__uint128_t));
    __uint128_t result = a_chunk | b_chunk;
    memcpy(&out[i * 16], &result, sizeof(__uint128_t));
  }

  uint64_t offset = n128 * 16;
  uint64_t n64 = (main_bytes - offset) / 8;

  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[offset + i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[offset + i * 8], sizeof(uint64_t));
    uint64_t result = a_chunk | b_chunk;
    memcpy(&out[offset + i * 8], &result, sizeof(uint64_t));
  }

  offset += n64 * 8;
#else
  uint64_t n64 = main_bytes / 8;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[i * 8], sizeof(uint64_t));
    uint64_t result = a_chunk | b_chunk;
    memcpy(&out[i * 8], &result, sizeof(uint64_t));
  }

  uint64_t offset = n64 * 8;
#endif

  for (uint64_t i = offset; i < main_bytes; i++)
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

#ifdef __SIZEOF_INT128__
  uint64_t n128 = main_bytes / 16;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n128; i++) {
    __uint128_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 16], sizeof(__uint128_t));
    memcpy(&b_chunk, &b[i * 16], sizeof(__uint128_t));
    __uint128_t result = a_chunk ^ b_chunk;
    memcpy(&out[i * 16], &result, sizeof(__uint128_t));
  }

  uint64_t offset = n128 * 16;
  uint64_t n64 = (main_bytes - offset) / 8;

  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[offset + i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[offset + i * 8], sizeof(uint64_t));
    uint64_t result = a_chunk ^ b_chunk;
    memcpy(&out[offset + i * 8], &result, sizeof(uint64_t));
  }

  offset += n64 * 8;
#else
  uint64_t n64 = main_bytes / 8;

  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[i * 8], sizeof(uint64_t));
    uint64_t result = a_chunk ^ b_chunk;
    memcpy(&out[i * 8], &result, sizeof(uint64_t));
  }

  uint64_t offset = n64 * 8;
#endif

  for (uint64_t i = offset; i < main_bytes; i++)
    out[i] = a[i] ^ b[i];

  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    out[full_bytes - 1] = (a[full_bytes - 1] ^ b[full_bytes - 1]) & mask;
  }
}
