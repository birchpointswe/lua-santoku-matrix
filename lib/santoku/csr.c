#include <santoku/iuset.h>
#include <santoku/csr.h>
#include <santoku/mtx.h>
#include <santoku/iumap.h>
#include <santoku/iumap/ext.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#ifndef TK_CVEC_BITS_BYTES
#define TK_CVEC_BITS_BYTES(n) (((n) + CHAR_BIT - 1) / CHAR_BIT)
#endif

static inline void *tk_csr_new_values (lua_State *L, tk_tag_t tag, uint64_t n)
{
  switch (tag) {
    case TK_TAG_I32: return tk_svec_create(L, n);
    case TK_TAG_I64: return tk_ivec_create(L, n);
    case TK_TAG_F32: return tk_fvec_create(L, n);
    case TK_TAG_F64: return tk_dvec_create(L, n);
    case TK_TAG_U8: return tk_cvec_create(L, n);
    default:
      tk_lua_verror(L, 2, "csr", "unsupported value type");
      return NULL;
  }
}

static inline void tk_csr_vals_grow (lua_State *L, tk_csr_t *X, uint64_t total)
{
  int rc;
  switch (X->tag) {
    case TK_TAG_I32: rc = tk_svec_ensure((tk_svec_t *) X->values, total); if (rc == 0) ((tk_svec_t *) X->values)->n = total; break;
    case TK_TAG_I64: rc = tk_ivec_ensure((tk_ivec_t *) X->values, total); if (rc == 0) ((tk_ivec_t *) X->values)->n = total; break;
    case TK_TAG_F32: rc = tk_fvec_ensure((tk_fvec_t *) X->values, total); if (rc == 0) ((tk_fvec_t *) X->values)->n = total; break;
    case TK_TAG_F64: rc = tk_dvec_ensure((tk_dvec_t *) X->values, total); if (rc == 0) ((tk_dvec_t *) X->values)->n = total; break;
    case TK_TAG_U8: rc = tk_cvec_ensure((tk_cvec_t *) X->values, total); if (rc == 0) ((tk_cvec_t *) X->values)->n = total; break;
    default: rc = 0; break;
  }
  if (rc != 0)
    tk_lua_verror(L, 2, "csr", "allocation failed");
}

static inline void tk_csr_vals_push (lua_State *L, tk_csr_t *X, double x)
{
  switch (X->tag) {
    case TK_TAG_I32: tk_svec_push((tk_svec_t *) X->values, (int32_t) x); break;
    case TK_TAG_I64: tk_ivec_push((tk_ivec_t *) X->values, (int64_t) x); break;
    case TK_TAG_F32: tk_fvec_push((tk_fvec_t *) X->values, (float) x); break;
    case TK_TAG_F64: tk_dvec_push((tk_dvec_t *) X->values, x); break;
    case TK_TAG_U8: tk_cvec_push((tk_cvec_t *) X->values, (char) (unsigned char) x); break;
    default:
      tk_lua_verror(L, 2, "csr", "push with value on a binary csr");
      break;
  }
}

static inline tk_tag_t tk_csr_tag_of_vals (lua_State *L, int i, void **child)
{
  void *p;
  if ((p = tk_fvec_peekopt(L, i)) != NULL) { *child = p; return TK_TAG_F32; }
  if ((p = tk_dvec_peekopt(L, i)) != NULL) { *child = p; return TK_TAG_F64; }
  if ((p = tk_ivec_peekopt(L, i)) != NULL) { *child = p; return TK_TAG_I64; }
  if ((p = tk_svec_peekopt(L, i)) != NULL) { *child = p; return TK_TAG_I32; }
  if ((p = tk_cvec_peekopt(L, i)) != NULL) { *child = p; return TK_TAG_U8; }
  return TK_TAG_NONE;
}

static int tk_csr_create_lua (lua_State *L)
{
  lua_settop(L, 1);
  luaL_checktype(L, 1, LUA_TTABLE);
  lua_getfield(L, 1, "n_cols");
  uint64_t n_cols = lua_isnil(L, -1) ? 0 : tk_lua_checkunsigned(L, -1, "n_cols");
  lua_pop(L, 1);
  lua_getfield(L, 1, "offsets");
  if (!lua_isnil(L, -1)) {

    tk_ivec_t *off = tk_ivec_peek(L, -1, "offsets");
    int io = lua_gettop(L);
    lua_getfield(L, 1, "neighbors");
    tk_ivec_t *nbr = tk_ivec_peek(L, -1, "neighbors");
    int in_ = lua_gettop(L);
    lua_getfield(L, 1, "values");
    tk_tag_t tag = TK_TAG_NONE;
    void *vals = NULL;
    int iv = 0;
    if (!lua_isnil(L, -1)) {
      tag = tk_csr_tag_of_vals(L, -1, &vals);
      if (tag == TK_TAG_NONE)
        return tk_lua_verror(L, 3, "csr", "values", "expected a santoku vec");
      iv = lua_gettop(L);
    }
    if (off->n < 1 || off->a[0] != 0)
      return tk_lua_verror(L, 3, "csr", "offsets", "must start at 0");
    tk_csr_push(L, tag, n_cols, io, off, in_, nbr, iv, vals);
    return 1;
  }
  lua_pop(L, 1);

  lua_getfield(L, 1, "values");
  tk_tag_t tag = TK_TAG_NONE;
  if (lua_type(L, -1) == LUA_TSTRING)
    tag = tk_tag_from_string(lua_tostring(L, -1));
  lua_pop(L, 1);
  tk_ivec_t *off = tk_ivec_create(L, 1);
  off->a[0] = 0;
  int io = lua_gettop(L);
  tk_ivec_t *nbr = tk_ivec_create(L, 0);
  int in_ = lua_gettop(L);
  void *vals = NULL;
  int iv = 0;
  if (tag != TK_TAG_NONE) {
    vals = tk_csr_new_values(L, tag, 0);
    iv = lua_gettop(L);
  }
  tk_csr_push(L, tag, n_cols, io, off, in_, nbr, iv, vals);
  return 1;
}


static int tk_csr_from_classes_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_ivec_t *cls = tk_ivec_peek(L, 1, "classes");
  uint64_t n_cols = 0;
  if (t >= 2 && !lua_isnil(L, 2))
    n_cols = tk_lua_checkunsigned(L, 2, "n_cols");
  else {
    for (uint64_t i = 0; i < cls->n; i ++)
      if (cls->a[i] >= 0 && (uint64_t) cls->a[i] + 1 > n_cols)
        n_cols = (uint64_t) cls->a[i] + 1;
  }
  tk_ivec_t *off = tk_ivec_create(L, cls->n + 1);
  int io = lua_gettop(L);
  tk_ivec_t *nbr = tk_ivec_create(L, cls->n);
  int in_ = lua_gettop(L);
  for (uint64_t i = 0; i < cls->n; i ++) {
    off->a[i] = (int64_t) i;
    nbr->a[i] = cls->a[i];
  }
  off->a[cls->n] = (int64_t) cls->n;
  tk_csr_push(L, TK_TAG_NONE, n_cols, io, off, in_, nbr, 0, NULL);
  return 1;
}


static int tk_csr_from_mask_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_ivec_t *mask = tk_ivec_peek(L, 1, "mask");
  tk_ivec_t *off = tk_ivec_create(L, mask->n + 1);
  int io = lua_gettop(L);
  tk_ivec_t *nbr = tk_ivec_create(L, 0);
  int in_ = lua_gettop(L);
  off->a[0] = 0;
  for (uint64_t i = 0; i < mask->n; i ++) {
    if (mask->a[i] != 0)
      tk_ivec_push(nbr, 0);
    off->a[i + 1] = (int64_t) nbr->n;
  }
  tk_csr_push(L, TK_TAG_NONE, 1, io, off, in_, nbr, 0, NULL);
  return 1;
}


static int tk_csr_from_bits_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_cvec_t *bitmap = tk_cvec_peek(L, 1, "bitmap");
  uint64_t n_rows = tk_lua_checkunsigned(L, 2, "n_rows");
  uint64_t n_cols = tk_lua_checkunsigned(L, 3, "n_cols");
  uint64_t bps = TK_CVEC_BITS_BYTES(n_cols);
  const uint8_t *data = (const uint8_t *) bitmap->a;
  uint64_t total = 0;
  for (uint64_t s = 0; s < n_rows; s ++) {
    const uint8_t *row = data + s * bps;
    for (uint64_t b = 0; b < bps; b ++)
      total += (uint64_t) __builtin_popcount((unsigned int) row[b]);
  }
  tk_ivec_t *off = tk_ivec_create(L, n_rows + 1);
  int io = lua_gettop(L);
  tk_ivec_t *nbr = tk_ivec_create(L, total);
  int in_ = lua_gettop(L);
  nbr->n = total;
  off->a[0] = 0;
  uint64_t pos = 0;
  for (uint64_t s = 0; s < n_rows; s ++) {
    const uint8_t *row = data + s * bps;
    for (uint64_t f = 0; f < n_cols; f ++)
      if (row[f / CHAR_BIT] & (1u << (f % CHAR_BIT)))
        nbr->a[pos ++] = (int64_t) f;
    off->a[s + 1] = (int64_t) pos;
  }
  tk_csr_push(L, TK_TAG_NONE, n_cols, io, off, in_, nbr, 0, NULL);
  return 1;
}

static int tk_csr_to_bits_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  bool flip = lua_toboolean(L, 2);
  uint64_t n_rows = tk_csr_rows(X);
  uint64_t out_cols = flip ? X->n_cols * 2 : X->n_cols;
  uint64_t bps = TK_CVEC_BITS_BYTES(out_cols);
  size_t len = n_rows * bps;
  tk_cvec_t *out = tk_cvec_create(L, len);
  uint8_t *buf = (uint8_t *) out->a;
  memset(buf, 0, len);
  for (uint64_t s = 0; s < n_rows; s ++) {
    uint64_t base = s * bps;
    if (flip)
      for (uint64_t k = 0; k < X->n_cols; k ++) {
        uint64_t bp = X->n_cols + k;
        buf[base + bp / CHAR_BIT] |= (1u << (bp % CHAR_BIT));
      }
    for (int64_t j = X->offsets->a[s]; j < X->offsets->a[s + 1]; j ++) {
      uint64_t f = (uint64_t) X->neighbors->a[j];
      buf[base + f / CHAR_BIT] |= (1u << (f % CHAR_BIT));
      if (flip) {
        uint64_t neg = X->n_cols + f;
        buf[base + neg / CHAR_BIT] &= ~(1u << (neg % CHAR_BIT));
      }
    }
  }
  return 1;
}

static int tk_csr_to_dense_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_mtx_t *out = NULL;
  tk_tag_t tag;
  if (lua_type(L, 2) == LUA_TUSERDATA) {
    out = tk_mtx_peek(L, 2, "out");
    tag = out->tag;
  } else {
    tag = lua_isnil(L, 2)
      ? (X->tag == TK_TAG_NONE ? TK_TAG_F64 : X->tag)
      : tk_tag_from_string(luaL_checkstring(L, 2));
    if (!lua_isnil(L, 3)) {
      out = tk_mtx_peek(L, 3, "out");
      if (out->tag != tag)
        return tk_lua_verror(L, 2, "csr", "to_dense: out type mismatch");
    }
  }
  if (tag == TK_TAG_NONE || tag == TK_TAG_BITS)
    return tk_lua_verror(L, 2, "csr", "to_dense: bad element type");
  uint64_t n_rows = tk_csr_rows(X);
  tk_lua_require_mod(L, "santoku.mtx");
  tk_mtx_t *M;
  if (out != NULL) {
    M = out;
    tk_mtx_reshape(L, M, n_rows, X->n_cols);
    lua_settop(L, lua_type(L, 2) == LUA_TUSERDATA ? 2 : 3);
  } else {
    M = tk_mtx_push_new(L, tag, n_rows, X->n_cols);
  }
  memset(tk_mtx_ptr(M), 0, tk_tag_size(tag) * n_rows * X->n_cols);
  for (uint64_t s = 0; s < n_rows; s ++)
    for (int64_t j = X->offsets->a[s]; j < X->offsets->a[s + 1]; j ++)
      tk_mtx_set1(M, s * X->n_cols + (uint64_t) X->neighbors->a[j],
        X->tag == TK_TAG_NONE ? 1.0 : tk_csr_val1(X, (uint64_t) j));
  return 1;
}

static int tk_csr_shape_lua (lua_State *L)
{
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  lua_pushinteger(L, (lua_Integer) tk_csr_rows(X));
  lua_pushinteger(L, (lua_Integer) X->n_cols);
  return 2;
}

static int tk_csr_nnz_lua (lua_State *L)
{
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  lua_pushinteger(L, (lua_Integer) tk_csr_nnz(X));
  return 1;
}

static int tk_csr_type_lua (lua_State *L)
{
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  lua_pushstring(L, tk_tag_name(X->tag));
  return 1;
}

static int tk_csr_offsets_lua (lua_State *L)
{
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_eph_get(L, 1, X->offsets);
  return 1;
}

static int tk_csr_neighbors_lua (lua_State *L)
{
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_eph_get(L, 1, X->neighbors);
  return 1;
}

static int tk_csr_values_lua (lua_State *L)
{
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  if (X->values == NULL)
    return 0;
  tk_eph_get(L, 1, X->values);
  return 1;
}

static int tk_csr_push_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  int64_t col = tk_lua_checkinteger(L, 2, "col");
  tk_ivec_push(X->neighbors, col);
  if (X->tag != TK_TAG_NONE)
    tk_csr_vals_push(L, X, t >= 3 ? luaL_checknumber(L, 3) : 1.0);
  else if (t >= 3)
    return tk_lua_verror(L, 2, "csr", "push with value on a binary csr");
  lua_settop(L, 1);
  return 1;
}

static int tk_csr_row_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_ivec_push(X->offsets, (int64_t) X->neighbors->n);
  lua_settop(L, 1);
  return 1;
}

static int tk_csr_rows_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_ivec_t *ids = tk_ivec_peek(L, 2, "ids");
  uint64_t n_rows = tk_csr_rows(X);
  uint64_t total = 0;
  for (uint64_t i = 0; i < ids->n; i ++) {
    int64_t s = ids->a[i];
    if (s < 0 || (uint64_t) s >= n_rows)
      return tk_lua_verror(L, 2, "csr", "row id out of range");
    total += (uint64_t) (X->offsets->a[s + 1] - X->offsets->a[s]);
  }
  tk_ivec_t *off = tk_ivec_create(L, ids->n + 1);
  int io = lua_gettop(L);
  tk_ivec_t *nbr = tk_ivec_create(L, total);
  int in_ = lua_gettop(L);
  nbr->n = total;
  void *vals = NULL;
  int iv = 0;
  if (X->tag != TK_TAG_NONE) {
    vals = tk_csr_new_values(L, X->tag, total);
    iv = lua_gettop(L);
  }
  size_t esz = tk_tag_size(X->tag);
  const char *vsrc = X->tag != TK_TAG_NONE ? (const char *) tk_csr_val_ptr(X) : NULL;
  uint64_t pos = 0;
  for (uint64_t i = 0; i < ids->n; i ++) {
    int64_t s = ids->a[i];
    int64_t lo = X->offsets->a[s], hi = X->offsets->a[s + 1];
    off->a[i] = (int64_t) pos;
    memcpy(nbr->a + pos, X->neighbors->a + lo, (size_t) (hi - lo) * sizeof(int64_t));
    if (vsrc) {
      char *vdst;
      switch (X->tag) {
        case TK_TAG_I32: vdst = (char *) ((tk_svec_t *) vals)->a; break;
        case TK_TAG_I64: vdst = (char *) ((tk_ivec_t *) vals)->a; break;
        case TK_TAG_F32: vdst = (char *) ((tk_fvec_t *) vals)->a; break;
        case TK_TAG_F64: vdst = (char *) ((tk_dvec_t *) vals)->a; break;
        default: vdst = ((tk_cvec_t *) vals)->a; break;
      }
      memcpy(vdst + pos * esz, vsrc + (uint64_t) lo * esz, (size_t) (hi - lo) * esz);
    }
    pos += (uint64_t) (hi - lo);
  }
  off->a[ids->n] = (int64_t) pos;
  tk_csr_push(L, X->tag, X->n_cols, io, off, in_, nbr, iv, vals);
  return 1;
}

static int tk_csr_select_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_ivec_t *remap = tk_ivec_peek(L, 2, "cols");
  tk_iumap_t *inverse = tk_iumap_from_ivec(L, remap);
  if (!inverse)
    return tk_lua_verror(L, 2, "csr", "allocation failed");
  int imap = lua_gettop(L);
  uint64_t n_rows = tk_csr_rows(X);
  tk_ivec_t *off = tk_ivec_create(L, n_rows + 1);
  int io = lua_gettop(L);
  tk_ivec_t *nbr = tk_ivec_create(L, X->neighbors->n);
  int in_ = lua_gettop(L);
  void *vals = NULL;
  int iv = 0;
  if (X->tag != TK_TAG_NONE) {
    vals = tk_csr_new_values(L, X->tag, X->neighbors->n);
    iv = lua_gettop(L);
  }
  off->a[0] = 0;
  uint64_t pos = 0;
  for (uint64_t r = 0; r < n_rows; r ++) {
    for (int64_t j = X->offsets->a[r]; j < X->offsets->a[r + 1]; j ++) {
      int64_t new_id = tk_iumap_get_or(inverse, X->neighbors->a[j], -1);
      if (new_id >= 0) {
        if (X->tag != TK_TAG_NONE)
          switch (X->tag) {
            case TK_TAG_I32: ((tk_svec_t *) vals)->a[pos] = ((tk_svec_t *) X->values)->a[j]; break;
            case TK_TAG_I64: ((tk_ivec_t *) vals)->a[pos] = ((tk_ivec_t *) X->values)->a[j]; break;
            case TK_TAG_F32: ((tk_fvec_t *) vals)->a[pos] = ((tk_fvec_t *) X->values)->a[j]; break;
            case TK_TAG_F64: ((tk_dvec_t *) vals)->a[pos] = ((tk_dvec_t *) X->values)->a[j]; break;
            default: ((tk_cvec_t *) vals)->a[pos] = ((tk_cvec_t *) X->values)->a[j]; break;
          }
        nbr->a[pos ++] = new_id;
      }
    }
    off->a[r + 1] = (int64_t) pos;
  }
  nbr->n = pos;
  if (X->tag != TK_TAG_NONE)
    switch (X->tag) {
      case TK_TAG_I32: ((tk_svec_t *) vals)->n = pos; break;
      case TK_TAG_I64: ((tk_ivec_t *) vals)->n = pos; break;
      case TK_TAG_F32: ((tk_fvec_t *) vals)->n = pos; break;
      case TK_TAG_F64: ((tk_dvec_t *) vals)->n = pos; break;
      default: ((tk_cvec_t *) vals)->n = pos; break;
    }
  tk_iumap_destroy(inverse);
  tk_csr_push(L, X->tag, remap->n, io, off, in_, nbr, iv, vals);
  lua_remove(L, imap);
  return 1;
}


static int tk_csr_hcat_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_csr_t *Y = tk_csr_peek(L, 2, "other");
  if (X->tag != Y->tag)
    return tk_lua_verror(L, 2, "csr", "hcat requires matching value types");
  uint64_t n_rows = tk_csr_rows(X);
  if (n_rows != tk_csr_rows(Y))
    return tk_lua_verror(L, 2, "csr", "hcat requires matching row counts");
  uint64_t na = X->neighbors->n, nb = Y->neighbors->n, total = na + nb;
  if (tk_ivec_ensure(X->neighbors, total) != 0)
    return tk_lua_verror(L, 2, "csr", "allocation failed");
  X->neighbors->n = total;
  if (X->tag != TK_TAG_NONE)
    tk_csr_vals_grow(L, X, total);
  size_t esz = tk_tag_size(X->tag);
  char *vdst = X->tag != TK_TAG_NONE ? (char *) tk_csr_val_ptr(X) : NULL;
  const char *vsrc = Y->tag != TK_TAG_NONE ? (const char *) tk_csr_val_ptr(Y) : NULL;
  int64_t shift = (int64_t) X->n_cols;
  for (uint64_t r = n_rows; r > 0; r --) {
    uint64_t i = r - 1;
    int64_t alo = X->offsets->a[i], ahi = X->offsets->a[i + 1];
    int64_t blo = Y->offsets->a[i], bhi = Y->offsets->a[i + 1];
    int64_t dst = alo + blo;
    memmove(X->neighbors->a + dst, X->neighbors->a + alo, (size_t) (ahi - alo) * sizeof(int64_t));
    if (vdst)
      memmove(vdst + (uint64_t) dst * esz, vdst + (uint64_t) alo * esz, (size_t) (ahi - alo) * esz);
    int64_t bdst = dst + (ahi - alo);
    for (int64_t j = blo; j < bhi; j ++)
      X->neighbors->a[bdst + (j - blo)] = Y->neighbors->a[j] + shift;
    if (vdst && vsrc)
      memcpy(vdst + (uint64_t) bdst * esz, vsrc + (uint64_t) blo * esz, (size_t) (bhi - blo) * esz);
  }
  for (uint64_t r = 0; r <= n_rows; r ++)
    X->offsets->a[r] += Y->offsets->a[r];
  X->n_cols += Y->n_cols;
  lua_settop(L, 1);
  return 1;
}

static int tk_csr_transpose_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  uint64_t n_rows = tk_csr_rows(X);
  uint64_t nnz = X->neighbors->n;
  int64_t *counts = (int64_t *) calloc(X->n_cols + 1, sizeof(int64_t));
  if (!counts)
    return tk_lua_verror(L, 2, "csr", "allocation failed");
  for (uint64_t i = 0; i < nnz; i ++)
    counts[X->neighbors->a[i] + 1] ++;
  for (uint64_t t = 0; t < X->n_cols; t ++)
    counts[t + 1] += counts[t];
  tk_ivec_t *off = tk_ivec_create(L, X->n_cols + 1);
  int io = lua_gettop(L);
  off->n = X->n_cols + 1;
  memcpy(off->a, counts, (X->n_cols + 1) * sizeof(int64_t));
  tk_ivec_t *nbr = tk_ivec_create(L, nnz);
  int in_ = lua_gettop(L);
  nbr->n = nnz;
  void *vals = NULL;
  int iv = 0;
  if (X->tag != TK_TAG_NONE) {
    vals = tk_csr_new_values(L, X->tag, nnz);
    iv = lua_gettop(L);
  }
  size_t esz = tk_tag_size(X->tag);
  const char *vsrc = X->tag != TK_TAG_NONE ? (const char *) tk_csr_val_ptr(X) : NULL;
  char *vdst = NULL;
  if (vals != NULL)
    switch (X->tag) {
      case TK_TAG_I32: vdst = (char *) ((tk_svec_t *) vals)->a; break;
      case TK_TAG_I64: vdst = (char *) ((tk_ivec_t *) vals)->a; break;
      case TK_TAG_F32: vdst = (char *) ((tk_fvec_t *) vals)->a; break;
      case TK_TAG_F64: vdst = (char *) ((tk_dvec_t *) vals)->a; break;
      default: vdst = ((tk_cvec_t *) vals)->a; break;
    }
  for (uint64_t s = 0; s < n_rows; s ++)
    for (int64_t j = X->offsets->a[s]; j < X->offsets->a[s + 1]; j ++) {
      int64_t tok = X->neighbors->a[j];
      int64_t pos = counts[tok] ++;
      nbr->a[pos] = (int64_t) s;
      if (vdst)
        memcpy(vdst + (uint64_t) pos * esz, vsrc + (uint64_t) j * esz, esz);
    }
  free(counts);
  tk_csr_push(L, X->tag, n_rows, io, off, in_, nbr, iv, vals);
  return 1;
}


static inline void tk_csr_materialize (lua_State *L, tk_csr_t *X, int ix)
{
  if (X->tag != TK_TAG_NONE)
    return;
  tk_fvec_t *vals = tk_fvec_create(L, X->neighbors->n);
  for (uint64_t i = 0; i < X->neighbors->n; i ++)
    vals->a[i] = 1.0f;
  X->tag = TK_TAG_F32;
  X->values = vals;
  tk_eph_anchor(L, ix, lua_gettop(L), vals);
  lua_pop(L, 1);
}

static int tk_csr_normalize_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_csr_materialize(L, X, 1);
  uint64_t n_rows = tk_csr_rows(X);
  for (uint64_t r = 0; r < n_rows; r ++) {
    int64_t lo = X->offsets->a[r], hi = X->offsets->a[r + 1];
    double ss = 0.0;
    for (int64_t j = lo; j < hi; j ++) {
      double v = tk_csr_val1(X, (uint64_t) j);
      ss += v * v;
    }
    if (ss > 0.0) {
      double inv = 1.0 / sqrt(ss);
      for (int64_t j = lo; j < hi; j ++)
        tk_csr_setval1(X, (uint64_t) j, tk_csr_val1(X, (uint64_t) j) * inv);
    }
  }
  lua_settop(L, 1);
  return 1;
}

static int tk_csr_scale_cols_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_fvec_t *wf = tk_fvec_peekopt(L, 2);
  tk_dvec_t *wd = wf == NULL ? tk_dvec_peek(L, 2, "weights") : NULL;
  uint64_t wn = wf != NULL ? wf->n : wd->n;
  if (wn < X->n_cols)
    return tk_lua_verror(L, 2, "csr", "scale_cols: weights shorter than n_cols");
  tk_csr_materialize(L, X, 1);
  for (uint64_t i = 0; i < X->neighbors->n; i ++) {
    double w = wf != NULL ? (double) wf->a[X->neighbors->a[i]] : wd->a[X->neighbors->a[i]];
    tk_csr_setval1(X, i, tk_csr_val1(X, i) * w);
  }
  lua_settop(L, 1);
  return 1;
}

static int tk_csr_sumsq_cols_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_ivec_t *bounds = lua_isnil(L, 2) ? NULL : tk_ivec_peek(L, 2, "bounds");
  if (bounds == NULL) {
    tk_dvec_t *out = tk_dvec_create(L, X->n_cols);
    memset(out->a, 0, X->n_cols * sizeof(double));
    for (uint64_t i = 0; i < X->neighbors->n; i ++) {
      double v = X->tag == TK_TAG_NONE ? 1.0 : tk_csr_val1(X, i);
      out->a[X->neighbors->a[i]] += v * v;
    }
    return 1;
  }
  uint64_t nb = bounds->n > 0 ? bounds->n - 1 : 0;
  tk_dvec_t *out = tk_dvec_create(L, nb);
  memset(out->a, 0, nb * sizeof(double));
  for (uint64_t i = 0; i < X->neighbors->n; i ++) {
    int64_t c = X->neighbors->a[i];

    uint64_t lo = 0, hi = nb;
    while (lo + 1 < hi) {
      uint64_t mid = (lo + hi) / 2;
      if (bounds->a[mid] <= c) lo = mid; else hi = mid;
    }
    if (c >= bounds->a[lo] && c < bounds->a[lo + 1]) {
      double v = X->tag == TK_TAG_NONE ? 1.0 : tk_csr_val1(X, i);
      out->a[lo] += v * v;
    }
  }
  return 1;
}

static int tk_csr_eq_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *a = tk_csr_peek(L, 1, "csr");
  tk_csr_t *b = tk_csr_peek(L, 2, "other");
  bool r = a->tag == b->tag && a->n_cols == b->n_cols
    && a->offsets->n == b->offsets->n && a->neighbors->n == b->neighbors->n
    && tk_ivec_eq(a->offsets, b->offsets, 0, a->offsets->n)
    && tk_ivec_eq(a->neighbors, b->neighbors, 0, a->neighbors->n);
  if (r && a->tag != TK_TAG_NONE)
    r = memcmp(tk_csr_val_ptr(a), tk_csr_val_ptr(b), tk_tag_size(a->tag) * a->neighbors->n) == 0;
  lua_pushboolean(L, r);
  return 1;
}

static int tk_csr_persist_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  FILE *fh = tk_lua_fopen(L, luaL_checkstring(L, 2), "w");
  char magic[4] = { 'T', 'K', 'c', 's' };
  uint8_t version = 1;
  uint8_t tag = (uint8_t) X->tag;
  uint64_t no = X->offsets->n, nn = X->neighbors->n;
  tk_lua_fwrite(L, magic, 4, 1, fh);
  tk_lua_fwrite(L, (char *) &version, 1, 1, fh);
  tk_lua_fwrite(L, (char *) &tag, 1, 1, fh);
  tk_lua_fwrite(L, (char *) &X->n_cols, sizeof(uint64_t), 1, fh);
  tk_lua_fwrite(L, (char *) &no, sizeof(uint64_t), 1, fh);
  tk_lua_fwrite(L, (char *) &nn, sizeof(uint64_t), 1, fh);
  tk_lua_fwrite(L, (char *) X->offsets->a, sizeof(int64_t) * no, 1, fh);
  tk_lua_fwrite(L, (char *) X->neighbors->a, sizeof(int64_t) * nn, 1, fh);
  if (X->tag != TK_TAG_NONE)
    tk_lua_fwrite(L, (char *) tk_csr_val_ptr(X), tk_tag_size(X->tag) * nn, 1, fh);
  tk_lua_fclose(L, fh);
  return 0;
}

static int tk_csr_load_lua (lua_State *L)
{
  lua_settop(L, 1);
  FILE *fh = tk_lua_fopen(L, luaL_checkstring(L, 1), "r");
  char magic[4];
  uint8_t version, tag8;
  uint64_t n_cols, no, nn;
  tk_lua_fread(L, magic, 4, 1, fh);
  if (memcmp(magic, "TKcs", 4) != 0)
    return tk_lua_verror(L, 2, "csr", "load: bad magic");
  tk_lua_fread(L, (char *) &version, 1, 1, fh);
  if (version != 1)
    return tk_lua_verror(L, 2, "csr", "load: unsupported version");
  tk_lua_fread(L, (char *) &tag8, 1, 1, fh);
  tk_lua_fread(L, (char *) &n_cols, sizeof(uint64_t), 1, fh);
  tk_lua_fread(L, (char *) &no, sizeof(uint64_t), 1, fh);
  tk_lua_fread(L, (char *) &nn, sizeof(uint64_t), 1, fh);
  tk_ivec_t *off = tk_ivec_create(L, no);
  int io = lua_gettop(L);
  tk_lua_fread(L, (char *) off->a, sizeof(int64_t) * no, 1, fh);
  tk_ivec_t *nbr = tk_ivec_create(L, nn);
  int in_ = lua_gettop(L);
  tk_lua_fread(L, (char *) nbr->a, sizeof(int64_t) * nn, 1, fh);
  void *vals = NULL;
  int iv = 0;
  if ((tk_tag_t) tag8 != TK_TAG_NONE) {
    vals = tk_csr_new_values(L, (tk_tag_t) tag8, nn);
    iv = lua_gettop(L);
    char *vdst;
    switch ((tk_tag_t) tag8) {
      case TK_TAG_I32: vdst = (char *) ((tk_svec_t *) vals)->a; break;
      case TK_TAG_I64: vdst = (char *) ((tk_ivec_t *) vals)->a; break;
      case TK_TAG_F32: vdst = (char *) ((tk_fvec_t *) vals)->a; break;
      case TK_TAG_F64: vdst = (char *) ((tk_dvec_t *) vals)->a; break;
      default: vdst = ((tk_cvec_t *) vals)->a; break;
    }
    tk_lua_fread(L, vdst, tk_tag_size((tk_tag_t) tag8) * nn, 1, fh);
  }
  tk_lua_fclose(L, fh);
  tk_csr_push(L, (tk_tag_t) tag8, n_cols, io, off, in_, nbr, iv, vals);
  return 1;
}

static luaL_Reg tk_csr_mt_fns[] = {
  { "shape", tk_csr_shape_lua },
  { "nnz", tk_csr_nnz_lua },
  { "type", tk_csr_type_lua },
  { "offsets", tk_csr_offsets_lua },
  { "neighbors", tk_csr_neighbors_lua },
  { "values", tk_csr_values_lua },
  { "push", tk_csr_push_lua },
  { "row", tk_csr_row_lua },
  { "rows", tk_csr_rows_lua },
  { "select", tk_csr_select_lua },
  { "hcat", tk_csr_hcat_lua },
  { "transpose", tk_csr_transpose_lua },
  { "normalize", tk_csr_normalize_lua },
  { "scale_cols", tk_csr_scale_cols_lua },
  { "sumsq_cols", tk_csr_sumsq_cols_lua },
  { "to_bits", tk_csr_to_bits_lua },
  { "to_dense", tk_csr_to_dense_lua },
  { "eq", tk_csr_eq_lua },
  { "persist", tk_csr_persist_lua },
  { NULL, NULL }
};

static luaL_Reg tk_csr_module_fns[] = {
  { "create", tk_csr_create_lua },
  { "from_classes", tk_csr_from_classes_lua },
  { "from_mask", tk_csr_from_mask_lua },
  { "from_bits", tk_csr_from_bits_lua },
  { "load", tk_csr_load_lua },
  { NULL, NULL }
};

int luaopen_santoku_csr (lua_State *L)
{
  luaL_newmetatable(L, TK_CSR_MT);
  lua_newtable(L);
  luaL_register(L, NULL, tk_csr_mt_fns);
  lua_setfield(L, -2, "__index");
  lua_pop(L, 1);
  lua_newtable(L);
  luaL_register(L, NULL, tk_csr_module_fns);
  return 1;
}
