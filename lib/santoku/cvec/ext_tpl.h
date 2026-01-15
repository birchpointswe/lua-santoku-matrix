#ifndef tk_parallel_sfx
#error "Must include santoku/parallel/tpl.h before this template"
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

#ifdef __SIZEOF_INT128__
  uint64_t n128 = main_bytes / 16;

  TK_PARALLEL_FOR(reduction(+:count))
  for (uint64_t i = 0; i < n128; i++) {
    __uint128_t chunk;
    memcpy(&chunk, &data[i * 16], sizeof(__uint128_t));
    uint64_t low = (uint64_t)chunk;
    uint64_t high = (uint64_t)(chunk >> 64);
    count += (uint64_t)__builtin_popcountll(low) + (uint64_t)__builtin_popcountll(high);
  }

  uint64_t offset = n128 * 16;
  uint64_t n64 = (main_bytes - offset) / 8;

  for (uint64_t i = 0; i < n64; i++) {
    uint64_t chunk;
    memcpy(&chunk, &data[offset + i * 8], sizeof(uint64_t));
    count += (uint64_t)__builtin_popcountll(chunk);
  }
  offset += n64 * 8;
#else
  uint64_t n64 = main_bytes / 8;

  TK_PARALLEL_FOR(reduction(+:count))
  for (uint64_t i = 0; i < n64; i++) {
    uint64_t chunk;
    memcpy(&chunk, &data[i * 8], sizeof(uint64_t));
    count += (uint64_t)__builtin_popcountll(chunk);
  }
  uint64_t offset = n64 * 8;
#endif

  for (uint64_t i = offset; i < main_bytes; i++)
    count += tk_popcount8(data[i]);

  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    count += tk_popcount8(data[full_bytes - 1] & mask);
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

#ifdef __SIZEOF_INT128__
  uint64_t n128 = main_bytes / 16;

  TK_PARALLEL_FOR(reduction(+:dist))
  for (uint64_t i = 0; i < n128; i++) {
    __uint128_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 16], sizeof(__uint128_t));
    memcpy(&b_chunk, &b[i * 16], sizeof(__uint128_t));
    __uint128_t xor_chunk = a_chunk ^ b_chunk;
    uint64_t low = (uint64_t)xor_chunk;
    uint64_t high = (uint64_t)(xor_chunk >> 64);
    dist += (uint64_t)__builtin_popcountll(low) + (uint64_t)__builtin_popcountll(high);
  }

  uint64_t offset = n128 * 16;
  uint64_t n64 = (main_bytes - offset) / 8;

  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[offset + i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[offset + i * 8], sizeof(uint64_t));
    uint64_t xor_chunk = a_chunk ^ b_chunk;
    dist += (uint64_t)__builtin_popcountll(xor_chunk);
  }
  offset += n64 * 8;
#else
  uint64_t n64 = main_bytes / 8;

  TK_PARALLEL_FOR(reduction(+:dist))
  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[i * 8], sizeof(uint64_t));
    uint64_t xor_chunk = a_chunk ^ b_chunk;
    dist += (uint64_t)__builtin_popcountll(xor_chunk);
  }
  uint64_t offset = n64 * 8;
#endif

  for (uint64_t i = offset; i < main_bytes; i++)
    dist += tk_popcount8(a[i] ^ b[i]);

  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    dist += tk_popcount8((a[full_bytes - 1] ^ b[full_bytes - 1]) & mask);
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

#ifdef __SIZEOF_INT128__
  uint64_t n128 = main_bytes / 16;

  TK_PARALLEL_FOR(reduction(+:dist))
  for (uint64_t i = 0; i < n128; i++) {
    __uint128_t a_chunk, b_chunk, m_chunk;
    memcpy(&a_chunk, &a[i * 16], sizeof(__uint128_t));
    memcpy(&b_chunk, &b[i * 16], sizeof(__uint128_t));
    memcpy(&m_chunk, &mask[i * 16], sizeof(__uint128_t));
    __uint128_t masked = (a_chunk ^ b_chunk) & m_chunk;
    uint64_t low = (uint64_t)masked;
    uint64_t high = (uint64_t)(masked >> 64);
    dist += (uint64_t)__builtin_popcountll(low) + (uint64_t)__builtin_popcountll(high);
  }

  uint64_t offset = n128 * 16;
  uint64_t n64 = (main_bytes - offset) / 8;

  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk, m_chunk;
    memcpy(&a_chunk, &a[offset + i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[offset + i * 8], sizeof(uint64_t));
    memcpy(&m_chunk, &mask[offset + i * 8], sizeof(uint64_t));
    uint64_t masked = (a_chunk ^ b_chunk) & m_chunk;
    dist += (uint64_t)__builtin_popcountll(masked);
  }
  offset += n64 * 8;
#else
  uint64_t n64 = main_bytes / 8;

  TK_PARALLEL_FOR(reduction(+:dist))
  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk, m_chunk;
    memcpy(&a_chunk, &a[i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[i * 8], sizeof(uint64_t));
    memcpy(&m_chunk, &mask[i * 8], sizeof(uint64_t));
    uint64_t masked = (a_chunk ^ b_chunk) & m_chunk;
    dist += (uint64_t)__builtin_popcountll(masked);
  }
  uint64_t offset = n64 * 8;
#endif

  for (uint64_t i = offset; i < main_bytes; i++)
    dist += tk_popcount8((a[i] ^ b[i]) & mask[i]);

  if (rem_bits > 0) {
    uint8_t m = mask[full_bytes - 1] & ((1U << rem_bits) - 1);
    dist += tk_popcount8((a[full_bytes - 1] ^ b[full_bytes - 1]) & m);
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

#ifdef __SIZEOF_INT128__
  uint64_t n128 = main_bytes / 16;

  TK_PARALLEL_FOR(reduction(+:pop_a,pop_andnot))
  for (uint64_t i = 0; i < n128; i++) {
    __uint128_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 16], sizeof(__uint128_t));
    memcpy(&b_chunk, &b[i * 16], sizeof(__uint128_t));
    __uint128_t va = a_chunk;
    __uint128_t vandnot = va & ~b_chunk;
    uint64_t a_low = (uint64_t)va;
    uint64_t a_high = (uint64_t)(va >> 64);
    uint64_t n_low = (uint64_t)vandnot;
    uint64_t n_high = (uint64_t)(vandnot >> 64);
    pop_a += (uint64_t)__builtin_popcountll(a_low) + (uint64_t)__builtin_popcountll(a_high);
    pop_andnot += (uint64_t)__builtin_popcountll(n_low) + (uint64_t)__builtin_popcountll(n_high);
  }

  uint64_t offset = n128 * 16;
  uint64_t n64 = (main_bytes - offset) / 8;

  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[offset + i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[offset + i * 8], sizeof(uint64_t));
    uint64_t va = a_chunk;
    uint64_t vandnot = va & ~b_chunk;
    pop_a += (uint64_t)__builtin_popcountll(va);
    pop_andnot += (uint64_t)__builtin_popcountll(vandnot);
  }
  offset += n64 * 8;
#else
  uint64_t n64 = main_bytes / 8;

  TK_PARALLEL_FOR(reduction(+:pop_a,pop_andnot))
  for (uint64_t i = 0; i < n64; i++) {
    uint64_t a_chunk, b_chunk;
    memcpy(&a_chunk, &a[i * 8], sizeof(uint64_t));
    memcpy(&b_chunk, &b[i * 8], sizeof(uint64_t));
    uint64_t va = a_chunk;
    uint64_t vandnot = va & ~b_chunk;
    pop_a += (uint64_t)__builtin_popcountll(va);
    pop_andnot += (uint64_t)__builtin_popcountll(vandnot);
  }
  uint64_t offset = n64 * 8;
#endif

  for (uint64_t i = offset; i < main_bytes; i++) {
    uint8_t va = a[i];
    uint8_t vandnot = va & ~b[i];
    pop_a += tk_popcount8(va);
    pop_andnot += tk_popcount8(vandnot);
  }

  if (rem_bits > 0) {
    uint8_t mask = (1U << rem_bits) - 1;
    uint8_t va = a[full_bytes - 1] & mask;
    uint8_t vandnot = va & ~b[full_bytes - 1];

    pop_a += tk_popcount8(va);
    pop_andnot += tk_popcount8(vandnot);
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

#ifdef __SIZEOF_INT128__
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
  if (base == NULL || ext == NULL || aids == NULL || bids == NULL)
    return -1;

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
