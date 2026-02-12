#ifndef tk_parallel_sfx
#error "Must include santoku/parallel/tpl.h before this template"
#endif

static inline uint64_t tk_cvec_bits_popcount(const uint8_t *data, uint64_t n_bits);
static inline uint64_t tk_cvec_bits_popcount_serial(const uint8_t *data, uint64_t n_bits);

static inline tk_ivec_t *tk_parallel_sfx(tk_ivec_bits_from_cvec) (lua_State *L, const char *bm, uint64_t n_samples, uint64_t n_features) {
  tk_ivec_t *out = tk_ivec_create(L, 0, 0, 0);

  uint64_t total_bits_in_bitmap = n_samples * n_features;
  uint64_t total_bits = tk_parallel_sfx(tk_cvec_bits_popcount)((const uint8_t *) bm, total_bits_in_bitmap);
  if (tk_ivec_ensure(out, total_bits) != 0) {
    tk_ivec_destroy(out);
    return NULL;
  }
  uint64_t *sample_counts = (uint64_t *)calloc(n_samples, sizeof(uint64_t));
  if (!sample_counts) {
    tk_ivec_destroy(out);
    return NULL;
  }

  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n_samples; i ++) {
    const uint8_t *sample_start = (const uint8_t *)(bm + i * bytes_per_sample);
    sample_counts[i] = tk_cvec_bits_popcount_serial(sample_start, n_features);
  }
  uint64_t *offsets = (uint64_t *)malloc((n_samples + 1) * sizeof(uint64_t));
  if (!offsets) {
    free(sample_counts);
    tk_ivec_destroy(out);
    return NULL;
  }
  offsets[0] = 0;
  for (uint64_t i = 0; i < n_samples; i++)
    offsets[i + 1] = offsets[i] + sample_counts[i];
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n_samples; i ++) {
    uint64_t write_pos = offsets[i];
    for (uint64_t j = 0; j < n_features; j ++) {
      uint64_t bit = i * n_features + j;
      uint64_t chunk = bit / CHAR_BIT;
      uint64_t pos = bit % CHAR_BIT;
      if (bm[chunk] & (1 << pos)) {
        out->a[write_pos++] = (int64_t) bit;
      }
    }
  }
  out->n = total_bits;
  free(sample_counts);
  free(offsets);
  tk_ivec_shrink(out);
  return out;
}

static inline tk_cvec_t *tk_parallel_sfx(tk_ivec_bits_to_cvec) (lua_State *L, tk_ivec_t *set_bits, uint64_t n_samples, uint64_t n_features, bool flip_interleave) {
  uint64_t output_features = flip_interleave ? (n_features * 2) : n_features;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(output_features);
  size_t len = n_samples * bytes_per_sample;
  tk_cvec_t *out_cvec = tk_cvec_create(L, len, NULL, NULL);
  uint8_t *out = (uint8_t *)out_cvec->a;

  uint64_t *counts = (uint64_t *)calloc(n_samples, sizeof(uint64_t));
  if (!counts) return out_cvec;
  for (uint64_t idx = 0; idx < set_bits->n; idx++) {
    int64_t v = set_bits->a[idx];
    if (v < 0) continue;
    uint64_t s = (uint64_t)v / n_features;
    if (s < n_samples) counts[s]++;
  }
  uint64_t *sample_offsets = (uint64_t *)malloc((n_samples + 1) * sizeof(uint64_t));
  if (!sample_offsets) { free(counts); return out_cvec; }
  sample_offsets[0] = 0;
  for (uint64_t s = 0; s < n_samples; s++)
    sample_offsets[s + 1] = sample_offsets[s] + counts[s];
  uint64_t *binned = (uint64_t *)malloc(sample_offsets[n_samples] * sizeof(uint64_t));
  if (!binned) { free(counts); free(sample_offsets); return out_cvec; }
  memset(counts, 0, n_samples * sizeof(uint64_t));
  for (uint64_t idx = 0; idx < set_bits->n; idx++) {
    int64_t v = set_bits->a[idx];
    if (v < 0) continue;
    uint64_t s = (uint64_t)v / n_features;
    if (s >= n_samples) continue;
    binned[sample_offsets[s] + counts[s]] = (uint64_t)v % n_features;
    counts[s]++;
  }
  free(counts);

  if (flip_interleave) {
    TK_PARALLEL_FOR(schedule(static))
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      memset(out + sample_offset, 0, bytes_per_sample);
      uint64_t neg_byte_offset = n_features / CHAR_BIT;
      uint8_t neg_bit_shift = n_features & (CHAR_BIT - 1);
      if (neg_bit_shift == 0) {
        uint64_t full_bytes = n_features / CHAR_BIT;
        memset(out + sample_offset + neg_byte_offset, 0xFF, full_bytes);
        uint64_t remaining_bits = n_features % CHAR_BIT;
        if (remaining_bits > 0)
          out[sample_offset + neg_byte_offset + full_bytes] = (1u << remaining_bits) - 1;
      } else {
        for (uint64_t k = 0; k < n_features; k++) {
          uint64_t bit_pos = n_features + k;
          uint64_t byte_off = sample_offset + (bit_pos / CHAR_BIT);
          uint8_t bit_off = bit_pos & (CHAR_BIT - 1);
          out[byte_off] |= (1u << bit_off);
        }
      }
      for (uint64_t bi = sample_offsets[s]; bi < sample_offsets[s + 1]; bi++) {
        uint64_t k = binned[bi];
        uint64_t pos_byte = sample_offset + (k / CHAR_BIT);
        uint8_t pos_bit = k & (CHAR_BIT - 1);
        out[pos_byte] |= (1u << pos_bit);
        uint64_t neg_bit_pos = n_features + k;
        uint64_t neg_byte = sample_offset + (neg_bit_pos / CHAR_BIT);
        uint8_t neg_bit = neg_bit_pos & (CHAR_BIT - 1);
        out[neg_byte] &= ~(1u << neg_bit);
      }
    }
  } else {
    memset(out, 0, len);
    TK_PARALLEL_FOR(schedule(static))
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t bi = sample_offsets[s]; bi < sample_offsets[s + 1]; bi++) {
        uint64_t f = binned[bi];
        uint64_t byte_off = sample_offset + (f / CHAR_BIT);
        uint8_t bit_off = f & (CHAR_BIT - 1);
        out[byte_off] |= (1u << bit_off);
      }
    }
  }

  free(binned);
  free(sample_offsets);
  return out_cvec;
}

static inline uint64_t tk_parallel_sfx(tk_ivec_bits_to_cvec_grouped) (
  lua_State *L,
  tk_ivec_t *set_bits,
  uint64_t n_samples,
  uint64_t n_visible,
  tk_ivec_t *offsets,
  tk_ivec_t *features,
  bool flip_interleave
) {
  uint64_t n_classes = offsets->n > 0 ? offsets->n - 1 : 0;
  if (n_classes == 0) {
    tk_cvec_create(L, 0, NULL, NULL);
    return 0;
  }
  uint64_t max_k = 0;
  for (uint64_t c = 0; c < n_classes; c++) {
    uint64_t class_k = (uint64_t)(offsets->a[c + 1] - offsets->a[c]);
    if (class_k > max_k) max_k = class_k;
  }
  if (max_k == 0) {
    tk_cvec_create(L, 0, NULL, NULL);
    return 0;
  }
  uint64_t bytes_per_class = flip_interleave ? TK_CVEC_BITS_BYTES(max_k * 2) : TK_CVEC_BITS_BYTES(max_k);
  uint64_t bytes_per_sample = n_classes * bytes_per_class;
  uint64_t bits_per_class = bytes_per_class * CHAR_BIT;
  size_t len = n_samples * bytes_per_sample;
  tk_cvec_t *out_cvec = tk_cvec_create(L, len, NULL, NULL);
  uint8_t *out = (uint8_t *)out_cvec->a;

  tk_iumap_t *feat_map = tk_iumap_create(0, 0);
  if (!feat_map) return max_k;
  for (uint64_t c = 0; c < n_classes; c++) {
    uint64_t start = (uint64_t)offsets->a[c];
    uint64_t end = (uint64_t)offsets->a[c + 1];
    for (uint64_t i = start; i < end; i++) {
      int64_t fid = features->a[i];
      if (fid < 0 || (uint64_t)fid >= n_visible) continue;
      uint64_t local = i - start;
      int64_t key = (int64_t)((uint64_t)fid * n_classes + c);
      int absent;
      khint_t kit = tk_iumap_put(feat_map, key, &absent);
      kh_value(feat_map, kit) = (int64_t)local;
    }
  }

  uint64_t *counts = (uint64_t *)calloc(n_samples, sizeof(uint64_t));
  if (!counts) { tk_iumap_destroy(feat_map); return max_k; }
  for (uint64_t idx = 0; idx < set_bits->n; idx++) {
    int64_t v = set_bits->a[idx];
    if (v < 0) continue;
    uint64_t s = (uint64_t)v / n_visible;
    if (s < n_samples) counts[s]++;
  }
  uint64_t *sample_offsets = (uint64_t *)malloc((n_samples + 1) * sizeof(uint64_t));
  if (!sample_offsets) { free(counts); tk_iumap_destroy(feat_map); return max_k; }
  sample_offsets[0] = 0;
  for (uint64_t s = 0; s < n_samples; s++)
    sample_offsets[s + 1] = sample_offsets[s] + counts[s];
  uint64_t *binned = (uint64_t *)malloc(sample_offsets[n_samples] * sizeof(uint64_t));
  if (!binned) { free(counts); free(sample_offsets); tk_iumap_destroy(feat_map); return max_k; }
  memset(counts, 0, n_samples * sizeof(uint64_t));
  for (uint64_t idx = 0; idx < set_bits->n; idx++) {
    int64_t v = set_bits->a[idx];
    if (v < 0) continue;
    uint64_t s = (uint64_t)v / n_visible;
    if (s >= n_samples) continue;
    binned[sample_offsets[s] + counts[s]] = (uint64_t)v % n_visible;
    counts[s]++;
  }
  free(counts);

  if (flip_interleave) {
    TK_PARALLEL_FOR(schedule(static))
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      memset(out + sample_offset, 0, bytes_per_sample);
      for (uint64_t c = 0; c < n_classes; c++) {
        for (uint64_t local = 0; local < max_k; local++) {
          uint64_t neg_bit = c * bits_per_class + max_k + local;
          uint64_t byte_off = sample_offset + (neg_bit / CHAR_BIT);
          uint8_t bit_off = neg_bit & (CHAR_BIT - 1);
          out[byte_off] |= (1u << bit_off);
        }
      }
      for (uint64_t bi = sample_offsets[s]; bi < sample_offsets[s + 1]; bi++) {
        uint64_t f = binned[bi];
        for (uint64_t c = 0; c < n_classes; c++) {
          int64_t key = (int64_t)(f * n_classes + c);
          khint_t kit = tk_iumap_get(feat_map, key);
          if (kit == kh_end(feat_map)) continue;
          int64_t local = kh_value(feat_map, kit);
          uint64_t pos_bit = c * bits_per_class + (uint64_t)local;
          uint64_t pos_byte = sample_offset + (pos_bit / CHAR_BIT);
          uint8_t pos_bit_off = pos_bit & (CHAR_BIT - 1);
          out[pos_byte] |= (1u << pos_bit_off);
          uint64_t neg_bit = c * bits_per_class + max_k + (uint64_t)local;
          uint64_t neg_byte = sample_offset + (neg_bit / CHAR_BIT);
          uint8_t neg_bit_off = neg_bit & (CHAR_BIT - 1);
          out[neg_byte] &= ~(1u << neg_bit_off);
        }
      }
    }
  } else {
    memset(out, 0, len);
    TK_PARALLEL_FOR(schedule(static))
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t sample_offset = s * bytes_per_sample;
      for (uint64_t bi = sample_offsets[s]; bi < sample_offsets[s + 1]; bi++) {
        uint64_t f = binned[bi];
        for (uint64_t c = 0; c < n_classes; c++) {
          int64_t key = (int64_t)(f * n_classes + c);
          khint_t kit = tk_iumap_get(feat_map, key);
          if (kit == kh_end(feat_map)) continue;
          int64_t local = kh_value(feat_map, kit);
          uint64_t out_idx = c * max_k + (uint64_t)local;
          uint64_t byte_off = sample_offset + (out_idx / CHAR_BIT);
          uint8_t bit_off = out_idx & (CHAR_BIT - 1);
          out[byte_off] |= (1u << bit_off);
        }
      }
    }
  }

  free(binned);
  free(sample_offsets);
  tk_iumap_destroy(feat_map);
  return max_k;
}

static inline tk_ivec_t *tk_parallel_sfx(tk_ivec_bits_extend) (tk_ivec_t *base, tk_ivec_t *ext, uint64_t n_base_features, uint64_t n_ext_features) {
  if (base == NULL || ext == NULL)
    return NULL;
  size_t total = base->n + ext->n;
  if (tk_ivec_ensure(base, total) != 0)
    return NULL;
  TK_PARALLEL_FOR(schedule(static))
  for (size_t i = 0; i < base->n; i ++) {
    uint64_t bit = (uint64_t) base->a[i];
    uint64_t sample = bit / n_base_features;
    uint64_t old_off = sample * n_base_features;
    uint64_t new_off = sample * (n_base_features + n_ext_features);
    base->a[i] = (int64_t) (bit - old_off + new_off);
  }
  TK_PARALLEL_FOR(schedule(static))
  for (size_t i = 0; i < ext->n; i ++) {
    uint64_t bit = (uint64_t) ext->a[i];
    uint64_t sample = bit / n_ext_features;
    uint64_t old_off = sample * n_ext_features;
    uint64_t new_off = sample * (n_base_features + n_ext_features);
    base->a[base->n + i] = (int64_t) (bit - old_off + new_off + n_base_features);
  }
  base->n = total;
  tk_ivec_asc(base, 0, base->n);
  return base;
}

static inline int tk_parallel_sfx(tk_ivec_bits_extend_mapped) (tk_ivec_t *base, tk_ivec_t *ext, tk_ivec_t *aids, tk_ivec_t *bids, uint64_t n_base_features, uint64_t n_ext_features, bool project) {
  if (base == NULL || ext == NULL || aids == NULL || bids == NULL)
    return -1;

  tk_iumap_t *a_id_to_pos = tk_iumap_from_ivec(0, aids);
  if (!a_id_to_pos)
    return -1;
  uint64_t n_only_b = 0;
  int64_t *b_to_final = (int64_t *)calloc(bids->n, sizeof(int64_t));
  int64_t *a_to_final = (int64_t *)calloc(aids->n, sizeof(int64_t));

  if (!b_to_final || !a_to_final) {
    free(b_to_final);
    free(a_to_final);
    tk_iumap_destroy(a_id_to_pos);
    return -1;
  }
  TK_PARALLEL_FOR(schedule(static))
  for (size_t i = 0; i < aids->n; i++)
    a_to_final[i] = (int64_t)i;
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
  uint64_t n_total_features = n_base_features + n_ext_features;
  size_t old_aids_n = aids->n;
  if (!project) {
    if (tk_ivec_ensure(aids, final_n_samples) != 0) {
      free(b_to_final);
      free(a_to_final);
      tk_iumap_destroy(a_id_to_pos);
      return -1;
    }
    for (size_t i = 0; i < bids->n; i++) {
      if (b_to_final[i] >= (int64_t)old_aids_n)
        aids->a[aids->n++] = bids->a[i];
    }
  }
  size_t max_bits = 0;
  TK_PARALLEL_FOR(reduction(+:max_bits))
  for (size_t i = 0; i < base->n; i++) {
    uint64_t bit = (uint64_t)base->a[i];
    uint64_t sample = bit / n_base_features;
    if (sample < old_aids_n)
      max_bits++;
  }
  TK_PARALLEL_FOR(reduction(+:max_bits))
  for (size_t i = 0; i < ext->n; i++) {
    uint64_t bit = (uint64_t)ext->a[i];
    uint64_t sample = bit / n_ext_features;
    if (sample < bids->n && b_to_final[sample] >= 0)
      max_bits++;
  }
  if (tk_ivec_ensure(base, max_bits) != 0) {
    free(b_to_final);
    free(a_to_final);
    tk_iumap_destroy(a_id_to_pos);
    return -1;
  }
  TK_PARALLEL_FOR(schedule(static))
  for (int64_t i = (int64_t)base->n - 1; i >= 0; i--) {
    uint64_t bit = (uint64_t)base->a[i];
    uint64_t sample = bit / n_base_features;
    uint64_t feature = bit % n_base_features;
    if (sample < old_aids_n) {
      uint64_t new_sample_pos = (uint64_t)a_to_final[sample];
      base->a[i] = (int64_t)(new_sample_pos * n_total_features + feature);
    }
  }
  size_t insert_pos = base->n;
  for (size_t i = 0; i < ext->n; i++) {
    uint64_t bit = (uint64_t)ext->a[i];
    uint64_t b_sample = bit / n_ext_features;
    uint64_t feature = bit % n_ext_features;
    if (b_sample < bids->n && b_to_final[b_sample] >= 0) {
      uint64_t final_sample = (uint64_t)b_to_final[b_sample];
      base->a[insert_pos++] = (int64_t)(final_sample * n_total_features + n_base_features + feature);
    }
  }
  base->n = insert_pos;
  tk_ivec_asc(base, 0, base->n);
  free(b_to_final);
  free(a_to_final);
  tk_iumap_destroy(a_id_to_pos);
  return 0;
}

static inline tk_ivec_t *tk_parallel_sfx(tk_ivec_bits_bipartite) (
  lua_State *L,
  tk_ivec_t *src,
  uint64_t n_docs,
  uint64_t n_labels,
  const char *mode,
  tk_ivec_t *tokens,
  uint64_t n_tokens,
  uint64_t *out_n_features
) {
  if (src == NULL) return NULL;
  uint64_t n_edges = src->n;

  if (strcmp(mode, "project") == 0) {
    *out_n_features = n_labels;
    uint64_t out_size = n_edges + n_labels;
    tk_ivec_t *out = tk_ivec_create(L, out_size, 0, 0);
    if (!out || tk_ivec_ensure(out, out_size) != 0) {
      if (out) { tk_ivec_destroy(out); lua_pop(L, 1); }
      return NULL;
    }
    uint64_t wp = 0;
    for (uint64_t e = 0; e < n_edges; e++)
      out->a[wp++] = src->a[e];
    for (uint64_t l = 0; l < n_labels; l++)
      out->a[wp++] = (int64_t)((n_docs + l) * n_labels + l);
    out->n = wp;
    tk_ivec_asc(out, 0, out->n);
    return out;
  }

  if (strcmp(mode, "adjacency") == 0) {
    uint64_t n_features = n_docs + n_labels;
    *out_n_features = n_features;
    uint64_t out_size = n_docs + n_labels + 2 * n_edges;
    tk_ivec_t *out = tk_ivec_create(L, out_size, 0, 0);
    if (!out || tk_ivec_ensure(out, out_size) != 0) {
      if (out) { tk_ivec_destroy(out); lua_pop(L, 1); }
      return NULL;
    }
    uint64_t wp = 0;
    for (uint64_t d = 0; d < n_docs; d++)
      out->a[wp++] = (int64_t)(d * n_features + d);
    for (uint64_t l = 0; l < n_labels; l++)
      out->a[wp++] = (int64_t)((n_docs + l) * n_features + (n_docs + l));
    for (uint64_t e = 0; e < n_edges; e++) {
      uint64_t bit = (uint64_t)src->a[e];
      uint64_t doc = bit / n_labels;
      uint64_t label = bit % n_labels;
      out->a[wp++] = (int64_t)(doc * n_features + (n_docs + label));
      out->a[wp++] = (int64_t)((n_docs + label) * n_features + doc);
    }
    out->n = wp;
    tk_ivec_asc(out, 0, out->n);
    return out;
  }

  if (strcmp(mode, "inherit") == 0) {

  if (tokens == NULL)
    return NULL;

  uint64_t n_features = n_tokens;
  *out_n_features = n_features;

  uint64_t *r1_offsets = (uint64_t *)calloc(n_docs + 1, sizeof(uint64_t));
  if (!r1_offsets) return NULL;
  for (uint64_t i = 0; i < tokens->n; i++) {
    uint64_t doc = (uint64_t)tokens->a[i] / n_tokens;
    if (doc < n_docs)
      r1_offsets[doc + 1]++;
  }
  for (uint64_t d = 0; d < n_docs; d++)
    r1_offsets[d + 1] += r1_offsets[d];

  uint64_t *lab_offsets = (uint64_t *)calloc(n_labels + 1, sizeof(uint64_t));
  uint64_t *lab_docs = (uint64_t *)malloc(n_edges * sizeof(uint64_t));
  if (!lab_offsets || !lab_docs) {
    free(r1_offsets); free(lab_offsets); free(lab_docs);
    return NULL;
  }
  for (uint64_t e = 0; e < n_edges; e++) {
    uint64_t label = (uint64_t)src->a[e] % n_labels;
    lab_offsets[label + 1]++;
  }
  for (uint64_t l = 0; l < n_labels; l++)
    lab_offsets[l + 1] += lab_offsets[l];
  uint64_t *cursors = (uint64_t *)malloc(n_labels * sizeof(uint64_t));
  if (!cursors) {
    free(r1_offsets); free(lab_offsets); free(lab_docs);
    return NULL;
  }
  memcpy(cursors, lab_offsets, n_labels * sizeof(uint64_t));
  for (uint64_t e = 0; e < n_edges; e++) {
    uint64_t bit = (uint64_t)src->a[e];
    uint64_t doc = bit / n_labels;
    uint64_t label = bit % n_labels;
    lab_docs[cursors[label]++] = doc;
  }
  free(cursors);

  uint64_t bm_bytes = TK_CVEC_BITS_BYTES(n_tokens);
  uint8_t *bm = (uint8_t *)calloc(1, bm_bytes);
  if (!bm) {
    free(r1_offsets); free(lab_offsets); free(lab_docs);
    return NULL;
  }
  uint64_t total_inherited = 0;
  for (uint64_t l = 0; l < n_labels; l++) {
    memset(bm, 0, bm_bytes);
    for (uint64_t j = lab_offsets[l]; j < lab_offsets[l + 1]; j++) {
      uint64_t doc = lab_docs[j];
      for (uint64_t k = r1_offsets[doc]; k < r1_offsets[doc + 1]; k++) {
        uint64_t tok = (uint64_t)tokens->a[k] % n_tokens;
        bm[tok / CHAR_BIT] |= (1 << (tok % CHAR_BIT));
      }
    }
    for (uint64_t b = 0; b < bm_bytes; b++)
      total_inherited += (uint64_t)__builtin_popcount((unsigned int)bm[b]);
  }

  uint64_t out_size = tokens->n + total_inherited;
  tk_ivec_t *out = tk_ivec_create(L, out_size, 0, 0);
  if (!out || tk_ivec_ensure(out, out_size) != 0) {
    if (out) { tk_ivec_destroy(out); lua_pop(L, 1); }
    free(bm); free(r1_offsets); free(lab_offsets); free(lab_docs);
    return NULL;
  }

  uint64_t wp = 0;
  for (uint64_t i = 0; i < tokens->n; i++) {
    uint64_t bit = (uint64_t)tokens->a[i];
    uint64_t doc = bit / n_tokens;
    uint64_t tok = bit % n_tokens;
    out->a[wp++] = (int64_t)(doc * n_features + tok);
  }

  for (uint64_t l = 0; l < n_labels; l++) {
    memset(bm, 0, bm_bytes);
    for (uint64_t j = lab_offsets[l]; j < lab_offsets[l + 1]; j++) {
      uint64_t doc = lab_docs[j];
      for (uint64_t k = r1_offsets[doc]; k < r1_offsets[doc + 1]; k++) {
        uint64_t tok = (uint64_t)tokens->a[k] % n_tokens;
        uint64_t byte = tok / CHAR_BIT;
        uint64_t bit_pos = tok % CHAR_BIT;
        if (!(bm[byte] & (1 << bit_pos))) {
          bm[byte] |= (1 << bit_pos);
          out->a[wp++] = (int64_t)((n_docs + l) * n_features + tok);
        }
      }
    }
  }

  out->n = wp;
  free(bm); free(r1_offsets); free(lab_offsets); free(lab_docs);
  tk_ivec_asc(out, 0, out->n);
  return out;

  }

  return NULL;
}

static inline void tk_parallel_sfx(tk_ivec_bits_to_hv) (
  lua_State *L,
  tk_ivec_t *offsets,
  tk_ivec_t *features,
  tk_cvec_t *tokens,
  uint64_t n_tokens,
  uint64_t hv_size,
  uint64_t shift_stride
) {
  uint64_t n_samples = offsets->n > 0 ? offsets->n - 1 : 0;
  uint64_t hv_bytes = TK_CVEC_BITS_BYTES(hv_size);
  uint64_t total = n_samples * hv_bytes;
  tk_cvec_t *out_cvec = tk_cvec_create(L, total, NULL, NULL);
  uint8_t *out = (uint8_t *)out_cvec->a;
  memset(out, 0, total);
  uint8_t *tok = (uint8_t *)tokens->a;
  TK_PARALLEL_FOR(schedule(static))
  for (uint64_t i = 0; i < n_samples; i++) {
    uint64_t start = (uint64_t)offsets->a[i];
    uint64_t end = (uint64_t)offsets->a[i + 1];
    uint8_t *row = out + i * hv_bytes;
    for (uint64_t j = start; j < end; j++) {
      int64_t fid = features->a[j];
      if (fid < 0 || (uint64_t)fid >= n_tokens) continue;
      uint8_t *token = tok + (uint64_t)fid * hv_bytes;
      uint64_t shift = ((j - start) * shift_stride) % hv_size;
      if (shift == 0) {
        for (uint64_t b = 0; b < hv_bytes; b++)
          row[b] |= token[b];
      } else {
        tk_hv_shift_or(row, token, hv_size, hv_bytes, shift);
      }
    }
  }
}

