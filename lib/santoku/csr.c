#include <santoku/ivec.h>
#include <santoku/fvec.h>
#include <santoku/dvec.h>
#include <santoku/cvec.h>
#include <santoku/cvec/ext.h>
#include <santoku/iumap.h>
#include <santoku/iumap/ext.h>
#include <santoku/iuset.h>
#include <santoku/iuset/ext.h>
#include <santoku/rvec.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#if defined(_OPENMP) && !defined(__EMSCRIPTEN__)
#include <omp.h>
#endif

#ifndef TK_CVEC_BITS_BYTES
#define TK_CVEC_BITS_BYTES(n) (((n) + CHAR_BIT - 1) / CHAR_BIT)
#endif

static int tk_csr_from_cvec_lua (lua_State *L)
{
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(n_features);
  const uint8_t *data = (const uint8_t *)bitmap->a;
  uint64_t total = 0;
  for (uint64_t s = 0; s < n_samples; s++) {
    const uint8_t *row = data + s * bytes_per_sample;
    for (uint64_t b = 0; b < bytes_per_sample; b++)
      total += (uint64_t)__builtin_popcount((unsigned int)row[b]);
  }
  tk_ivec_t *offsets = tk_ivec_create(L, n_samples + 1);
  tk_ivec_t *neighbors = tk_ivec_create(L, total);
  neighbors->n = total;
  offsets->a[0] = 0;
  uint64_t pos = 0;
  for (uint64_t s = 0; s < n_samples; s++) {
    const uint8_t *row = data + s * bytes_per_sample;
    for (uint64_t f = 0; f < n_features; f++) {
      uint64_t byte_idx = f / CHAR_BIT;
      uint64_t bit_pos = f % CHAR_BIT;
      if (row[byte_idx] & (1u << bit_pos))
        neighbors->a[pos++] = (int64_t)f;
    }
    offsets->a[s + 1] = (int64_t)pos;
  }
  return 2;
}

static int tk_csr_to_cvec_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *offsets = tk_ivec_peek(L, 1, "offsets");
  tk_ivec_t *neighbors = tk_ivec_peek(L, 2, "neighbors");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  bool flip_interleave = (nargs >= 5) ? lua_toboolean(L, 5) : false;
  uint64_t out_features = flip_interleave ? (n_features * 2) : n_features;
  uint64_t bytes_per_sample = TK_CVEC_BITS_BYTES(out_features);
  size_t len = n_samples * bytes_per_sample;
  tk_cvec_t *out = tk_cvec_create(L, len);
  uint8_t *buf = (uint8_t *)out->a;
  if (flip_interleave) {
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t base = s * bytes_per_sample;
      memset(buf + base, 0, bytes_per_sample);
      for (uint64_t k = 0; k < n_features; k++) {
        uint64_t bit_pos = n_features + k;
        buf[base + bit_pos / CHAR_BIT] |= (1u << (bit_pos % CHAR_BIT));
      }
      int64_t lo = offsets->a[s], hi = offsets->a[s + 1];
      for (int64_t j = lo; j < hi; j++) {
        uint64_t f = (uint64_t)neighbors->a[j];
        buf[base + f / CHAR_BIT] |= (1u << (f % CHAR_BIT));
        uint64_t neg = n_features + f;
        buf[base + neg / CHAR_BIT] &= ~(1u << (neg % CHAR_BIT));
      }
    }
  } else {
    memset(buf, 0, len);
    for (uint64_t s = 0; s < n_samples; s++) {
      uint64_t base = s * bytes_per_sample;
      int64_t lo = offsets->a[s], hi = offsets->a[s + 1];
      for (int64_t j = lo; j < hi; j++) {
        uint64_t f = (uint64_t)neighbors->a[j];
        buf[base + f / CHAR_BIT] |= (1u << (f % CHAR_BIT));
      }
    }
  }
  return 1;
}

static int tk_csr_from_dvec_lua (lua_State *L)
{
  tk_dvec_t *dense = tk_dvec_peek(L, 1, "dense");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t total = 0;
  for (uint64_t i = 0; i < n_samples * n_features; i++) {
    if (dense->a[i] != 0.0)
      total++;
  }
  tk_ivec_t *offsets = tk_ivec_create(L, n_samples + 1);
  tk_ivec_t *neighbors = tk_ivec_create(L, total);
  tk_dvec_t *values = tk_dvec_create(L, total);
  neighbors->n = total;
  values->n = total;
  offsets->a[0] = 0;
  uint64_t pos = 0;
  for (uint64_t s = 0; s < n_samples; s++) {
    for (uint64_t f = 0; f < n_features; f++) {
      double v = dense->a[s * n_features + f];
      if (v != 0.0) {
        neighbors->a[pos] = (int64_t)f;
        values->a[pos] = v;
        pos++;
      }
    }
    offsets->a[s + 1] = (int64_t)pos;
  }
  return 3;
}

static int tk_csr_to_dvec_lua (lua_State *L)
{
  tk_ivec_t *offsets = tk_ivec_peek(L, 1, "offsets");
  tk_ivec_t *neighbors = tk_ivec_peek(L, 2, "neighbors");
  tk_dvec_t *values = lua_isnil(L, 3) ? NULL : tk_dvec_peek(L, 3, "values");
  uint64_t n_samples = tk_lua_checkunsigned(L, 4, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 5, "n_features");
  uint64_t total = n_samples * n_features;
  tk_dvec_t *out = tk_dvec_create(L, total);
  memset(out->a, 0, total * sizeof(double));
  for (uint64_t s = 0; s < n_samples; s++) {
    int64_t lo = offsets->a[s], hi = offsets->a[s + 1];
    for (int64_t j = lo; j < hi; j++) {
      uint64_t f = (uint64_t)neighbors->a[j];
      out->a[s * n_features + f] = values ? values->a[j] : 1.0;
    }
  }
  return 1;
}

static int tk_csr_from_fvec_lua (lua_State *L)
{
  tk_fvec_t *dense = tk_fvec_peek(L, 1, "dense");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t total = 0;
  for (uint64_t i = 0; i < n_samples * n_features; i++) {
    if (dense->a[i] != 0.0f)
      total++;
  }
  tk_ivec_t *offsets = tk_ivec_create(L, n_samples + 1);
  tk_ivec_t *neighbors = tk_ivec_create(L, total);
  tk_fvec_t *values = tk_fvec_create(L, total);
  neighbors->n = total;
  values->n = total;
  offsets->a[0] = 0;
  uint64_t pos = 0;
  for (uint64_t s = 0; s < n_samples; s++) {
    for (uint64_t f = 0; f < n_features; f++) {
      float v = dense->a[s * n_features + f];
      if (v != 0.0f) {
        neighbors->a[pos] = (int64_t)f;
        values->a[pos] = v;
        pos++;
      }
    }
    offsets->a[s + 1] = (int64_t)pos;
  }
  return 3;
}

static int tk_csr_to_fvec_lua (lua_State *L)
{
  tk_ivec_t *offsets = tk_ivec_peek(L, 1, "offsets");
  tk_ivec_t *neighbors = tk_ivec_peek(L, 2, "neighbors");
  tk_fvec_t *values = lua_isnil(L, 3) ? NULL : tk_fvec_peek(L, 3, "values");
  uint64_t n_samples = tk_lua_checkunsigned(L, 4, "n_samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 5, "n_features");
  uint64_t total = n_samples * n_features;
  tk_fvec_t *out = tk_fvec_create(L, total);
  memset(out->a, 0, total * sizeof(float));
  for (uint64_t s = 0; s < n_samples; s++) {
    int64_t lo = offsets->a[s], hi = offsets->a[s + 1];
    for (int64_t j = lo; j < hi; j++) {
      uint64_t f = (uint64_t)neighbors->a[j];
      out->a[s * n_features + f] = values ? values->a[j] : 1.0f;
    }
  }
  return 1;
}

static int tk_csr_extend_lua (lua_State *L)
{
  tk_ivec_t *off1 = tk_ivec_peek(L, 1, "off1");
  tk_ivec_t *nbr1 = tk_ivec_peek(L, 2, "nbr1");
  tk_ivec_t *off2 = tk_ivec_peek(L, 3, "off2");
  tk_ivec_t *nbr2 = tk_ivec_peek(L, 4, "nbr2");
  uint64_t n_feat1 = tk_lua_checkunsigned(L, 5, "n_feat1");
  uint64_t n_samples = off1->n - 1;
  uint64_t total = nbr1->n + nbr2->n;
  tk_ivec_t *out_off = tk_ivec_create(L, n_samples + 1);
  tk_ivec_t *out_nbr = tk_ivec_create(L, total);
  out_nbr->n = total;
  out_off->a[0] = 0;
  uint64_t pos = 0;
  for (uint64_t s = 0; s < n_samples; s++) {
    int64_t lo1 = off1->a[s], hi1 = off1->a[s + 1];
    for (int64_t j = lo1; j < hi1; j++)
      out_nbr->a[pos++] = nbr1->a[j];
    int64_t lo2 = off2->a[s], hi2 = off2->a[s + 1];
    for (int64_t j = lo2; j < hi2; j++)
      out_nbr->a[pos++] = (int64_t)((uint64_t)nbr2->a[j] + n_feat1);
    out_off->a[s + 1] = (int64_t)pos;
  }
  return 2;
}

static int tk_csr_select_lua (lua_State *L)
{
  tk_ivec_t *offsets = tk_ivec_peek(L, 1, "offsets");
  tk_ivec_t *neighbors = tk_ivec_peek(L, 2, "neighbors");
  tk_ivec_t *remap_ids = tk_ivec_peek(L, 3, "remap_ids");
  tk_iumap_t *inverse = tk_iumap_from_ivec(L, remap_ids);
  if (!inverse)
    return luaL_error(L, "select: allocation failed");
  uint64_t n_rows = offsets->n - 1;
  tk_ivec_t *new_off = tk_ivec_create(L, n_rows + 1);
  tk_ivec_t *new_nbr = tk_ivec_create(L, neighbors->n);
  new_off->a[0] = 0;
  uint64_t pos = 0;
  for (uint64_t r = 0; r < n_rows; r++) {
    int64_t lo = offsets->a[r], hi = offsets->a[r + 1];
    for (int64_t j = lo; j < hi; j++) {
      int64_t new_id = tk_iumap_get_or(inverse, neighbors->a[j], -1);
      if (new_id >= 0)
        new_nbr->a[pos++] = new_id;
    }
    new_off->a[r + 1] = (int64_t)pos;
  }
  new_nbr->n = pos;
  tk_iumap_destroy(inverse);
  return 2;
}

static int tk_csr_seq_select_lua (lua_State *L)
{
  tk_ivec_t *tokens = tk_ivec_peek(L, 1, "tokens");
  tk_ivec_t *offsets = tk_ivec_peek(L, 2, "offsets");
  tk_ivec_t *keep_ids = tk_ivec_peek(L, 3, "keep_ids");
  tk_iumap_t *inverse = tk_iumap_from_ivec(L, keep_ids);
  if (!inverse)
    return luaL_error(L, "seq_select: allocation failed");
  uint64_t n_docs = offsets->n - 1;
  tk_ivec_t *new_tok = tk_ivec_create(L, tokens->n);
  tk_ivec_t *new_off = tk_ivec_create(L, n_docs + 1);
  bool has_values = false;
  bool is_fvec = false;
  tk_dvec_t *dvalues = NULL;
  tk_fvec_t *fvalues = NULL;
  tk_dvec_t *new_dval = NULL;
  tk_fvec_t *new_fval = NULL;
  if (tk_lua_testuserdata(L, 4, "tk_fvec_t")) {
    fvalues = tk_fvec_peek(L, 4, "values");
    new_fval = tk_fvec_create(L, tokens->n);
    has_values = true;
    is_fvec = true;
  } else {
    dvalues = tk_dvec_peekopt(L, 4);
    if (dvalues) {
      new_dval = tk_dvec_create(L, tokens->n);
      has_values = true;
    }
  }
  new_off->a[0] = 0;
  uint64_t pos = 0;
  for (uint64_t d = 0; d < n_docs; d++) {
    int64_t lo = offsets->a[d], hi = offsets->a[d + 1];
    for (int64_t j = lo; j < hi; j++) {
      int64_t new_id = tk_iumap_get_or(inverse, tokens->a[j], -1);
      if (new_id < 0) continue;
      new_tok->a[pos] = new_id;
      if (is_fvec) new_fval->a[pos] = fvalues->a[j];
      else if (new_dval) new_dval->a[pos] = dvalues->a[j];
      pos++;
    }
    new_off->a[d + 1] = (int64_t)pos;
  }
  new_tok->n = pos;
  if (new_fval) new_fval->n = pos;
  if (new_dval) new_dval->n = pos;
  tk_iumap_destroy(inverse);
  return has_values ? 3 : 2;
}

static int tk_csr_merge_lua (lua_State *L)
{
  tk_ivec_t *off1 = tk_ivec_peek(L, 1, "off1");
  tk_ivec_t *nbr1 = tk_ivec_peek(L, 2, "nbr1");
  tk_ivec_t *off2 = tk_ivec_peek(L, 3, "off2");
  tk_ivec_t *nbr2 = tk_ivec_peek(L, 4, "nbr2");
  uint64_t n = off1->n - 1;
  uint64_t total = nbr1->n + nbr2->n;
  tk_ivec_t *out_off = tk_ivec_create(L, n + 1);
  tk_ivec_t *out_nbr = tk_ivec_create(L, total);
  out_nbr->n = total;
  out_off->a[0] = 0;
  uint64_t pos = 0;
  for (uint64_t i = 0; i < n; i++) {
    int64_t s1 = off1->a[i], e1 = off1->a[i + 1];
    int64_t s2 = off2->a[i], e2 = off2->a[i + 1];
    int64_t i1 = s1, i2 = s2;
    while (i1 < e1 && i2 < e2) {
      if (nbr1->a[i1] <= nbr2->a[i2])
        out_nbr->a[pos++] = nbr1->a[i1++];
      else
        out_nbr->a[pos++] = nbr2->a[i2++];
    }
    while (i1 < e1) out_nbr->a[pos++] = nbr1->a[i1++];
    while (i2 < e2) out_nbr->a[pos++] = nbr2->a[i2++];
    out_off->a[i + 1] = (int64_t)pos;
  }
  return 2;
}

static int tk_csr_subsample_lua (lua_State *L)
{
  tk_ivec_t *offsets = tk_ivec_peek(L, 1, "offsets");
  tk_ivec_t *neighbors = tk_ivec_peek(L, 2, "neighbors");
  tk_ivec_t *sample_ids = tk_ivec_peek(L, 3, "sample_ids");
  int64_t n = (int64_t)sample_ids->n;
  int64_t total = 0;
  for (int64_t i = 0; i < n; i++) {
    int64_t sid = sample_ids->a[i];
    total += offsets->a[sid + 1] - offsets->a[sid];
  }
  tk_ivec_t *new_off = tk_ivec_create(L, (uint64_t)(n + 1));
  tk_ivec_t *new_nbr = tk_ivec_create(L, (uint64_t)total);
  new_nbr->n = (uint64_t)total;
  int64_t pos = 0;
  for (int64_t i = 0; i < n; i++) {
    int64_t sid = sample_ids->a[i];
    int64_t lo = offsets->a[sid], hi = offsets->a[sid + 1];
    new_off->a[i] = pos;
    for (int64_t j = lo; j < hi; j++)
      new_nbr->a[pos++] = neighbors->a[j];
  }
  new_off->a[n] = pos;
  return 2;
}

static int tk_csr_transpose_lua (lua_State *L)
{
  tk_ivec_t *offsets = tk_ivec_peek(L, 1, "offsets");
  tk_ivec_t *tokens = tk_ivec_peek(L, 2, "tokens");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "n_samples");
  uint64_t n_tokens = tk_lua_checkunsigned(L, 4, "n_tokens");
  bool has_values = false;
  bool is_fvec = false;
  tk_dvec_t *dvalues = NULL;
  tk_fvec_t *fvalues = NULL;
  if (tk_lua_testuserdata(L, 5, "tk_fvec_t")) {
    fvalues = tk_fvec_peek(L, 5, "values");
    has_values = true;
    is_fvec = true;
  } else {
    dvalues = tk_dvec_peekopt(L, 5);
    if (dvalues) has_values = true;
  }
  uint64_t nnz = tokens->n;
  int64_t *counts = (int64_t *)calloc(n_tokens + 1, sizeof(int64_t));
  if (!counts)
    return luaL_error(L, "transpose: allocation failed");
  for (uint64_t i = 0; i < nnz; i++)
    counts[tokens->a[i] + 1]++;
  for (uint64_t t = 0; t < n_tokens; t++)
    counts[t + 1] += counts[t];
  tk_ivec_t *csc_off = tk_ivec_create(L, n_tokens + 1);
  csc_off->n = n_tokens + 1;
  memcpy(csc_off->a, counts, (n_tokens + 1) * sizeof(int64_t));
  tk_ivec_t *csc_rows = tk_ivec_create(L, nnz);
  csc_rows->n = nnz;
  tk_fvec_t *csc_fvals = NULL;
  tk_dvec_t *csc_dvals = NULL;
  if (is_fvec) {
    csc_fvals = tk_fvec_create(L, nnz);
    csc_fvals->n = nnz;
  } else if (dvalues) {
    csc_dvals = tk_dvec_create(L, nnz);
    csc_dvals->n = nnz;
  }
  for (uint64_t s = 0; s < n_samples; s++) {
    for (int64_t j = offsets->a[s]; j < offsets->a[s + 1]; j++) {
      int64_t tok = tokens->a[j];
      int64_t pos = counts[tok]++;
      csc_rows->a[pos] = (int64_t)s;
      if (csc_fvals) csc_fvals->a[pos] = fvalues->a[j];
      else if (csc_dvals) csc_dvals->a[pos] = dvalues->a[j];
    }
  }
  free(counts);
  return has_values ? 3 : 2;
}

static int tk_csr_complement_lua (lua_State *L)
{
  tk_ivec_t *subset = tk_ivec_peek(L, 1, "subset");
  uint64_t n = tk_lua_checkunsigned(L, 2, "n");
  uint64_t cn = n - subset->n;
  tk_ivec_t *out = tk_ivec_create(L, cn);
  out->n = cn;
  uint64_t si = 0, oi = 0;
  for (uint64_t i = 0; i < n; i++) {
    if (si < subset->n && (uint64_t)subset->a[si] == i)
      si++;
    else
      out->a[oi++] = (int64_t)i;
  }
  return 1;
}

static int tk_csr_sort_desc_lua (lua_State *L)
{
  tk_ivec_t *off = tk_ivec_peek(L, 1, "offsets");
  tk_ivec_t *nbr = tk_ivec_peek(L, 2, "neighbors");
  uint64_t n = off->n - 1;
  tk_ivec_t *out_n = tk_ivec_create(L, nbr->n);
  memcpy(out_n->a, nbr->a, nbr->n * sizeof(int64_t));
  if (tk_lua_testuserdata(L, 3, "tk_fvec_t")) {
    tk_fvec_t *scores = tk_fvec_peek(L, 3, "scores");
    tk_fvec_t *out_s = tk_fvec_create(L, scores->n);
    memcpy(out_s->a, scores->a, scores->n * sizeof(float));
    #pragma omp parallel for schedule(dynamic, 64)
    for (uint64_t i = 0; i < n; i++) {
      int64_t s = off->a[i], e = off->a[i + 1];
      for (int64_t j = s + 1; j < e; j++) {
        float ks = out_s->a[j];
        int64_t kn = out_n->a[j];
        int64_t p = j - 1;
        while (p >= s && out_s->a[p] < ks) {
          out_s->a[p + 1] = out_s->a[p];
          out_n->a[p + 1] = out_n->a[p];
          p--;
        }
        out_s->a[p + 1] = ks;
        out_n->a[p + 1] = kn;
      }
    }
  } else {
    tk_dvec_t *scores = tk_dvec_peek(L, 3, "scores");
    tk_dvec_t *out_s = tk_dvec_create(L, scores->n);
    memcpy(out_s->a, scores->a, scores->n * sizeof(double));
    #pragma omp parallel for schedule(dynamic, 64)
    for (uint64_t i = 0; i < n; i++) {
      int64_t s = off->a[i], e = off->a[i + 1];
      for (int64_t j = s + 1; j < e; j++) {
        double ks = out_s->a[j];
        int64_t kn = out_n->a[j];
        int64_t p = j - 1;
        while (p >= s && out_s->a[p] < ks) {
          out_s->a[p + 1] = out_s->a[p];
          out_n->a[p + 1] = out_n->a[p];
          p--;
        }
        out_s->a[p + 1] = ks;
        out_n->a[p + 1] = kn;
      }
    }
  }
  return 2;
}

static int tk_csr_scatter_fixed_k_lua (lua_State *L)
{
  tk_ivec_t *dst_nbr = tk_ivec_peek(L, 1, "dst_neighbors");
  tk_ivec_t *src_off = tk_ivec_peek(L, 3, "src_offsets");
  tk_ivec_t *src_nbr = tk_ivec_peek(L, 4, "src_neighbors");
  tk_ivec_t *val_ids = tk_ivec_peek(L, 6, "val_ids");
  uint64_t k = tk_lua_checkunsigned(L, 7, "k");
  uint64_t val_n = val_ids->n;
  if (tk_lua_testuserdata(L, 2, "tk_fvec_t")) {
    tk_fvec_t *dst_sco = tk_fvec_peek(L, 2, "dst_scores");
    tk_fvec_t *src_sco = tk_fvec_peek(L, 5, "src_scores");
    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < val_n; i++) {
      int64_t orig = val_ids->a[i];
      int64_t dst = orig * (int64_t)k;
      int64_t src = src_off->a[i];
      for (uint64_t j = 0; j < k; j++) {
        dst_nbr->a[dst + (int64_t)j] = src_nbr->a[src + (int64_t)j];
        dst_sco->a[dst + (int64_t)j] = src_sco->a[src + (int64_t)j];
      }
    }
  } else {
    tk_dvec_t *dst_sco = tk_dvec_peek(L, 2, "dst_scores");
    tk_dvec_t *src_sco = tk_dvec_peek(L, 5, "src_scores");
    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < val_n; i++) {
      int64_t orig = val_ids->a[i];
      int64_t dst = orig * (int64_t)k;
      int64_t src = src_off->a[i];
      for (uint64_t j = 0; j < k; j++) {
        dst_nbr->a[dst + (int64_t)j] = src_nbr->a[src + (int64_t)j];
        dst_sco->a[dst + (int64_t)j] = src_sco->a[src + (int64_t)j];
      }
    }
  }
  return 0;
}

static int tk_csr_scatter_rows_lua (lua_State *L)
{
  tk_ivec_t *val_ids = tk_ivec_peek(L, 3, "val_ids");
  uint64_t stride = tk_lua_checkunsigned(L, 4, "stride");
  uint64_t val_n = val_ids->n;
  if (tk_lua_testuserdata(L, 1, "tk_fvec_t")) {
    tk_fvec_t *dst = tk_fvec_peek(L, 1, "dst");
    tk_fvec_t *src = tk_fvec_peek(L, 2, "src");
    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < val_n; i++) {
      int64_t orig = val_ids->a[i];
      memcpy(dst->a + orig * (int64_t)stride,
             src->a + (int64_t)(i * stride),
             stride * sizeof(float));
    }
  } else {
    tk_dvec_t *dst = tk_dvec_peek(L, 1, "dst");
    tk_dvec_t *src = tk_dvec_peek(L, 2, "src");
    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < val_n; i++) {
      int64_t orig = val_ids->a[i];
      memcpy(dst->a + orig * (int64_t)stride,
             src->a + (int64_t)(i * stride),
             stride * sizeof(double));
    }
  }
  return 0;
}

static int tk_csr_sample_select_lua (lua_State *L)
{
  tk_ivec_t *offsets = tk_ivec_peek(L, 1, "offsets");
  tk_ivec_t *neighbors = tk_ivec_peek(L, 2, "neighbors");
  tk_ivec_t *sample_ids = tk_ivec_peek(L, 3, "sample_ids");
  tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n = sample_ids->n;
  int64_t total = 0;
  for (uint64_t i = 0; i < n; i++) {
    int64_t sid = sample_ids->a[i];
    total += offsets->a[sid + 1] - offsets->a[sid];
  }
  tk_ivec_t *new_off = tk_ivec_create(L, n + 1);
  tk_ivec_t *new_nbr = tk_ivec_create(L, (uint64_t)total);
  tk_iuset_t *feat_set = tk_iuset_create(0, 0);
  int64_t pos = 0;
  new_off->a[0] = 0;
  for (uint64_t i = 0; i < n; i++) {
    int64_t sid = sample_ids->a[i];
    int64_t lo = offsets->a[sid], hi = offsets->a[sid + 1];
    for (int64_t j = lo; j < hi; j++) {
      new_nbr->a[pos] = neighbors->a[j];
      int absent;
      tk_iuset_put(feat_set, neighbors->a[j], &absent);
      pos++;
    }
    new_off->a[i + 1] = pos;
  }
  new_nbr->n = (uint64_t)pos;
  uint64_t n_local = tk_iuset_size(feat_set);
  tk_ivec_t *local_feats = tk_ivec_create(L, n_local);
  local_feats->n = n_local;
  uint64_t fi = 0;
  int64_t fkey;
  tk_umap_foreach_keys(feat_set, fkey, ({
    local_feats->a[fi++] = fkey;
  }));
  tk_iuset_destroy(feat_set);
  tk_ivec_asc(local_feats, 0, local_feats->n);
  tk_iumap_t *remap = tk_iumap_create(0, 0);
  for (uint64_t i = 0; i < local_feats->n; i++) {
    int absent;
    uint32_t k = tk_iumap_put(remap, local_feats->a[i], &absent);
    tk_iumap_setval(remap, k, (int64_t)i);
  }
  for (uint64_t i = 0; i < new_nbr->n; i++)
    new_nbr->a[i] = tk_iumap_get_or(remap, new_nbr->a[i], new_nbr->a[i]);
  tk_iumap_destroy(remap);
  lua_pushinteger(L, (lua_Integer)n_local);
  return 4;
}

static tk_ivec_t *tk_csr_dedup_to_flat (
  lua_State *L,
  tk_ivec_t *offsets,
  tk_ivec_t *neighbors,
  uint64_t n_samples,
  uint64_t n_features
) {
  tk_ivec_t *out = tk_ivec_create(L, neighbors->n);
  char *seen = (char *)calloc(n_features, 1);
  uint64_t pos = 0;
  for (uint64_t s = 0; s < n_samples; s++) {
    int64_t lo = offsets->a[s], hi = offsets->a[s + 1];
    for (int64_t j = lo; j < hi; j++) {
      uint64_t f = (uint64_t)neighbors->a[j];
      if (!seen[f]) {
        seen[f] = 1;
        out->a[pos++] = (int64_t)(s * n_features + f);
      }
    }
    for (int64_t j = lo; j < hi; j++)
      seen[(uint64_t)neighbors->a[j]] = 0;
  }
  out->n = pos;
  free(seen);
  return out;
}

static int tk_csr_top_chi2_lua (lua_State *L)
{
  lua_settop(L, 9);
  tk_ivec_t *feat_offsets = tk_ivec_peek(L, 1, "feat_offsets");
  tk_ivec_t *feat_neighbors = tk_ivec_peek(L, 2, "feat_neighbors");
  uint64_t n_features = tk_lua_checkunsigned(L, 5, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 6, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 7) ? n_features : tk_lua_checkunsigned(L, 7, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 8) ? 0 : tk_lua_checkunsigned(L, 8, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 9));
  uint64_t n_samples = feat_offsets->n - 1;
  tk_ivec_t *labels = NULL;
  char *codes = NULL;
  if (tk_lua_testuserdata(L, 3, "tk_cvec_t")) {
    tk_cvec_t *cvec = tk_cvec_peek(L, 3, "codes");
    codes = cvec->a;
  } else if (!lua_isnil(L, 3)) {
    tk_ivec_t *label_off = tk_ivec_peek(L, 3, "label_offsets");
    tk_ivec_t *label_nbr = tk_ivec_peek(L, 4, "label_neighbors");
    labels = tk_ivec_create(L, label_nbr->n);
    labels->n = label_nbr->n;
    for (uint64_t s = 0; s < n_samples; s++) {
      int64_t lo = label_off->a[s], hi = label_off->a[s + 1];
      for (int64_t j = lo; j < hi; j++)
        labels->a[j] = (int64_t)(s * n_hidden + (uint64_t)label_nbr->a[j]);
    }
  }
  tk_ivec_t *set_bits = tk_csr_dedup_to_flat(L, feat_offsets, feat_neighbors, n_samples, n_features);
  tk_ivec_bits_top_chi2(L, set_bits, codes, labels, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static int tk_csr_top_mi_lua (lua_State *L)
{
  lua_settop(L, 9);
  tk_ivec_t *feat_offsets = tk_ivec_peek(L, 1, "feat_offsets");
  tk_ivec_t *feat_neighbors = tk_ivec_peek(L, 2, "feat_neighbors");
  uint64_t n_features = tk_lua_checkunsigned(L, 5, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 6, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 7) ? n_features : tk_lua_checkunsigned(L, 7, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 8) ? 0 : tk_lua_checkunsigned(L, 8, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 9));
  uint64_t n_samples = feat_offsets->n - 1;
  tk_ivec_t *labels = NULL;
  char *codes = NULL;
  if (tk_lua_testuserdata(L, 3, "tk_cvec_t")) {
    tk_cvec_t *cvec = tk_cvec_peek(L, 3, "codes");
    codes = cvec->a;
  } else if (!lua_isnil(L, 3)) {
    tk_ivec_t *label_off = tk_ivec_peek(L, 3, "label_offsets");
    tk_ivec_t *label_nbr = tk_ivec_peek(L, 4, "label_neighbors");
    labels = tk_ivec_create(L, label_nbr->n);
    labels->n = label_nbr->n;
    for (uint64_t s = 0; s < n_samples; s++) {
      int64_t lo = label_off->a[s], hi = label_off->a[s + 1];
      for (int64_t j = lo; j < hi; j++)
        labels->a[j] = (int64_t)(s * n_hidden + (uint64_t)label_nbr->a[j]);
    }
  }
  tk_ivec_t *set_bits = tk_csr_dedup_to_flat(L, feat_offsets, feat_neighbors, n_samples, n_features);
  tk_ivec_bits_top_mi(L, set_bits, codes, labels, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static int tk_csr_top_bns_lua (lua_State *L)
{
  lua_settop(L, 9);
  tk_ivec_t *feat_offsets = tk_ivec_peek(L, 1, "feat_offsets");
  tk_ivec_t *feat_neighbors = tk_ivec_peek(L, 2, "feat_neighbors");
  uint64_t n_features = tk_lua_checkunsigned(L, 5, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 6, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 7) ? n_features : tk_lua_checkunsigned(L, 7, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 8) ? 0 : tk_lua_checkunsigned(L, 8, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 9));
  uint64_t n_samples = feat_offsets->n - 1;
  tk_ivec_t *labels = NULL;
  char *codes = NULL;
  if (tk_lua_testuserdata(L, 3, "tk_cvec_t")) {
    tk_cvec_t *cvec = tk_cvec_peek(L, 3, "codes");
    codes = cvec->a;
  } else if (!lua_isnil(L, 3)) {
    tk_ivec_t *label_off = tk_ivec_peek(L, 3, "label_offsets");
    tk_ivec_t *label_nbr = tk_ivec_peek(L, 4, "label_neighbors");
    labels = tk_ivec_create(L, label_nbr->n);
    labels->n = label_nbr->n;
    for (uint64_t s = 0; s < n_samples; s++) {
      int64_t lo = label_off->a[s], hi = label_off->a[s + 1];
      for (int64_t j = lo; j < hi; j++)
        labels->a[j] = (int64_t)(s * n_hidden + (uint64_t)label_nbr->a[j]);
    }
  }
  tk_ivec_t *set_bits = tk_csr_dedup_to_flat(L, feat_offsets, feat_neighbors, n_samples, n_features);
  tk_ivec_bits_top_bns(L, set_bits, codes, labels, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static int tk_csr_top_entropy_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_ivec_t *feat_offsets = tk_ivec_peek(L, 1, "feat_offsets");
  tk_ivec_t *feat_neighbors = tk_ivec_peek(L, 2, "feat_neighbors");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t top_k = lua_isnil(L, 4) ? n_features : tk_lua_checkunsigned(L, 4, "top_k");
  uint64_t n_samples = feat_offsets->n - 1;
  tk_ivec_t *set_bits = tk_csr_dedup_to_flat(L, feat_offsets, feat_neighbors, n_samples, n_features);
  tk_ivec_bits_top_entropy(L, set_bits, n_samples, n_features, top_k);
  return 2;
}

static int tk_csr_top_idf_lua (lua_State *L)
{
  lua_settop(L, 6);
  tk_ivec_t *feat_offsets = tk_ivec_peek(L, 1, "feat_offsets");
  tk_ivec_t *feat_neighbors = tk_ivec_peek(L, 2, "feat_neighbors");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t top_k = lua_isnil(L, 4) ? n_features : tk_lua_checkunsigned(L, 4, "top_k");
  double min_df = tk_lua_optnumber(L, 5, "min_df", 0.0);
  double max_df = tk_lua_optnumber(L, 6, "max_df", 1.0);
  uint64_t n_samples = feat_offsets->n - 1;
  tk_ivec_t *set_bits = tk_csr_dedup_to_flat(L, feat_offsets, feat_neighbors, n_samples, n_features);
  if (!tk_ivec_bits_top_idf(L, set_bits, n_samples, n_features, min_df, max_df, top_k))
    return tk_lua_verror(L, 2, "top_idf", "allocation failed");
  return 2;
}

static int tk_csr_top_df_lua (lua_State *L)
{
  lua_settop(L, 6);
  tk_ivec_t *feat_offsets = tk_ivec_peek(L, 1, "feat_offsets");
  tk_ivec_t *feat_neighbors = tk_ivec_peek(L, 2, "feat_neighbors");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "n_features");
  uint64_t top_k = lua_isnil(L, 4) ? n_features : tk_lua_checkunsigned(L, 4, "top_k");
  double min_df = tk_lua_optnumber(L, 5, "min_df", 0.0);
  double max_df = tk_lua_optnumber(L, 6, "max_df", 1.0);
  uint64_t n_samples = feat_offsets->n - 1;
  tk_ivec_t *set_bits = tk_csr_dedup_to_flat(L, feat_offsets, feat_neighbors, n_samples, n_features);
  if (!tk_ivec_bits_top_df(L, set_bits, n_samples, n_features, min_df, max_df, top_k))
    return tk_lua_verror(L, 2, "top_df", "allocation failed");
  return 2;
}

static int tk_csr_top_reg_f_lua (lua_State *L)
{
  lua_settop(L, 8);
  tk_ivec_t *feat_offsets = tk_ivec_peek(L, 1, "feat_offsets");
  tk_ivec_t *feat_neighbors = tk_ivec_peek(L, 2, "feat_neighbors");
  tk_dvec_t *targets = tk_dvec_peek(L, 3, "targets");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  uint64_t n_samples = feat_offsets->n - 1;
  tk_ivec_t *set_bits = tk_csr_dedup_to_flat(L, feat_offsets, feat_neighbors, n_samples, n_features);
  tk_ivec_bits_top_reg_f(L, set_bits, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static int tk_csr_top_reg_pearson_lua (lua_State *L)
{
  lua_settop(L, 8);
  tk_ivec_t *feat_offsets = tk_ivec_peek(L, 1, "feat_offsets");
  tk_ivec_t *feat_neighbors = tk_ivec_peek(L, 2, "feat_neighbors");
  tk_dvec_t *targets = tk_dvec_peek(L, 3, "targets");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  uint64_t n_samples = feat_offsets->n - 1;
  tk_ivec_t *set_bits = tk_csr_dedup_to_flat(L, feat_offsets, feat_neighbors, n_samples, n_features);
  tk_ivec_bits_top_reg_pearson(L, set_bits, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static int tk_csr_top_reg_r2_lua (lua_State *L)
{
  lua_settop(L, 8);
  tk_ivec_t *feat_offsets = tk_ivec_peek(L, 1, "feat_offsets");
  tk_ivec_t *feat_neighbors = tk_ivec_peek(L, 2, "feat_neighbors");
  tk_dvec_t *targets = tk_dvec_peek(L, 3, "targets");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  uint64_t n_samples = feat_offsets->n - 1;
  tk_ivec_t *set_bits = tk_csr_dedup_to_flat(L, feat_offsets, feat_neighbors, n_samples, n_features);
  tk_ivec_bits_top_reg_r2(L, set_bits, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static int tk_csr_top_reg_cohen_lua (lua_State *L)
{
  lua_settop(L, 8);
  tk_ivec_t *feat_offsets = tk_ivec_peek(L, 1, "feat_offsets");
  tk_ivec_t *feat_neighbors = tk_ivec_peek(L, 2, "feat_neighbors");
  tk_dvec_t *targets = tk_dvec_peek(L, 3, "targets");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  uint64_t n_samples = feat_offsets->n - 1;
  tk_ivec_t *set_bits = tk_csr_dedup_to_flat(L, feat_offsets, feat_neighbors, n_samples, n_features);
  tk_ivec_bits_top_reg_cohen(L, set_bits, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static int tk_csr_top_reg_auc_lua (lua_State *L)
{
  lua_settop(L, 8);
  tk_ivec_t *feat_offsets = tk_ivec_peek(L, 1, "feat_offsets");
  tk_ivec_t *feat_neighbors = tk_ivec_peek(L, 2, "feat_neighbors");
  tk_dvec_t *targets = tk_dvec_peek(L, 3, "targets");
  uint64_t n_features = tk_lua_checkunsigned(L, 4, "n_features");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "n_hidden");
  uint64_t per_class_k = lua_isnil(L, 6) ? n_features : tk_lua_checkunsigned(L, 6, "per_class_k");
  uint64_t union_top_k = lua_isnil(L, 7) ? 0 : tk_lua_checkunsigned(L, 7, "union_top_k");
  tk_pool_t pool = tk_pool_from_string(lua_tostring(L, 8));
  uint64_t n_samples = feat_offsets->n - 1;
  tk_ivec_t *set_bits = tk_csr_dedup_to_flat(L, feat_offsets, feat_neighbors, n_samples, n_features);
  tk_ivec_bits_top_reg_auc(L, set_bits, targets, n_samples, n_features, n_hidden, per_class_k, union_top_k, pool);
  return 5;
}

static luaL_Reg tk_csr_fns[] = {
  { "from_cvec", tk_csr_from_cvec_lua },
  { "to_cvec", tk_csr_to_cvec_lua },
  { "from_dvec", tk_csr_from_dvec_lua },
  { "to_dvec", tk_csr_to_dvec_lua },
  { "from_fvec", tk_csr_from_fvec_lua },
  { "to_fvec", tk_csr_to_fvec_lua },
  { "extend", tk_csr_extend_lua },
  { "select", tk_csr_select_lua },
  { "seq_select", tk_csr_seq_select_lua },
  { "merge", tk_csr_merge_lua },
  { "subsample", tk_csr_subsample_lua },
  { "transpose", tk_csr_transpose_lua },
  { "complement", tk_csr_complement_lua },
  { "sort_desc", tk_csr_sort_desc_lua },
  { "scatter_fixed_k", tk_csr_scatter_fixed_k_lua },
  { "scatter_rows", tk_csr_scatter_rows_lua },
  { "sample_select", tk_csr_sample_select_lua },
  { "top_chi2", tk_csr_top_chi2_lua },
  { "top_mi", tk_csr_top_mi_lua },
  { "top_bns", tk_csr_top_bns_lua },
  { "top_entropy", tk_csr_top_entropy_lua },
  { "top_idf", tk_csr_top_idf_lua },
  { "top_df", tk_csr_top_df_lua },
  { "top_reg_f", tk_csr_top_reg_f_lua },
  { "top_reg_pearson", tk_csr_top_reg_pearson_lua },
  { "top_reg_r2", tk_csr_top_reg_r2_lua },
  { "top_reg_cohen", tk_csr_top_reg_cohen_lua },
  { "top_reg_auc", tk_csr_top_reg_auc_lua },
  { NULL, NULL }
};

int luaopen_santoku_csr (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_csr_fns);
  return 1;
}
