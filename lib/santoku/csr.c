#include <santoku/iuset.h>
#include <santoku/csr.h>
#include <santoku/mtx.h>
#include <santoku/iumap.h>
#include <santoku/iumap/ext.h>
#include <santoku/rvec.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

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
    void *nbr;
    tk_tag_t ntag;
    tk_svec_t *nbr_s = tk_svec_peekopt(L, -1);
    if (nbr_s != NULL) { nbr = nbr_s; ntag = TK_TAG_I32; }
    else { nbr = tk_ivec_peek(L, -1, "neighbors"); ntag = TK_TAG_I64; }
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
    tk_csr_push(L, tag, ntag, n_cols, io, off, in_, nbr, iv, vals);
    return 1;
  }
  lua_pop(L, 1);

  lua_getfield(L, 1, "values");
  tk_tag_t tag = TK_TAG_NONE;
  if (lua_type(L, -1) == LUA_TSTRING)
    tag = tk_tag_from_string(lua_tostring(L, -1));
  lua_pop(L, 1);
  lua_getfield(L, 1, "neighbors");
  tk_tag_t ntag = TK_TAG_I64;
  if (lua_type(L, -1) == LUA_TSTRING) {
    tk_tag_t nt = tk_tag_from_string(lua_tostring(L, -1));
    if (nt == TK_TAG_I32 || nt == TK_TAG_I64) ntag = nt;
  }
  lua_pop(L, 1);
  tk_ivec_t *off = tk_ivec_create(L, 1);
  off->a[0] = 0;
  int io = lua_gettop(L);
  void *nbr = tk_csr_new_nbr(L, ntag, 0);
  int in_ = lua_gettop(L);
  void *vals = NULL;
  int iv = 0;
  if (tag != TK_TAG_NONE) {
    vals = tk_csr_new_values(L, tag, 0);
    iv = lua_gettop(L);
  }
  tk_csr_push(L, tag, ntag, n_cols, io, off, in_, nbr, iv, vals);
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
  tk_csr_push(L, TK_TAG_NONE, TK_TAG_I64, n_cols, io, off, in_, nbr, 0, NULL);
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
  tk_csr_push(L, TK_TAG_NONE, TK_TAG_I64, 1, io, off, in_, nbr, 0, NULL);
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
  tk_csr_push(L, TK_TAG_NONE, TK_TAG_I64, n_cols, io, off, in_, nbr, 0, NULL);
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
      uint64_t f = (uint64_t) tk_csr_nbr(X, (uint64_t) j);
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
      tk_mtx_set1(M, s * X->n_cols + (uint64_t) tk_csr_nbr(X, (uint64_t) j),
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
  tk_csr_nbr_push(X, col);
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
  tk_ivec_push(X->offsets, (int64_t) tk_csr_nbr_n(X));
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
  void *nbr = tk_csr_new_nbr(L, X->ntag, total);
  int in_ = lua_gettop(L);
  tk_nbr_setn(nbr, X->ntag, total);
  void *vals = NULL;
  int iv = 0;
  if (X->tag != TK_TAG_NONE) {
    vals = tk_csr_new_values(L, X->tag, total);
    iv = lua_gettop(L);
  }
  size_t esz = tk_tag_size(X->tag);
  size_t nesz = tk_nbr_esz(X->ntag);
  char *ndst = (char *) tk_nbr_aptr(nbr, X->ntag);
  const char *nsrc = (const char *) tk_csr_nbr_ptr(X);
  const char *vsrc = X->tag != TK_TAG_NONE ? (const char *) tk_csr_val_ptr(X) : NULL;
  uint64_t pos = 0;
  for (uint64_t i = 0; i < ids->n; i ++) {
    int64_t s = ids->a[i];
    int64_t lo = X->offsets->a[s], hi = X->offsets->a[s + 1];
    off->a[i] = (int64_t) pos;
    memcpy(ndst + pos * nesz, nsrc + (uint64_t) lo * nesz, (size_t) (hi - lo) * nesz);
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
  tk_csr_push(L, X->tag, X->ntag, X->n_cols, io, off, in_, nbr, iv, vals);
  return 1;
}


static int tk_csr_append_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_csr_t *B = tk_csr_peek(L, 2, "other");
  if (X->n_cols != B->n_cols)
    return tk_lua_verror(L, 2, "csr", "append requires matching n_cols");
  if (X->ntag != B->ntag)
    return tk_lua_verror(L, 2, "csr", "append requires matching neighbor type");
  if (X->tag != B->tag)
    return tk_lua_verror(L, 2, "csr", "append requires matching value type");
  uint64_t base = tk_csr_nnz(X);
  uint64_t nb = tk_csr_nnz(B);
  size_t nesz = tk_nbr_esz(X->ntag);
  if (tk_csr_nbr_ensure(X, base + nb) != 0)
    return tk_lua_verror(L, 2, "csr", "allocation failed");
  memcpy((char *) tk_csr_nbr_ptr(X) + base * nesz, tk_csr_nbr_ptr(B), nb * nesz);
  tk_csr_nbr_setn(X, base + nb);
  if (X->tag != TK_TAG_NONE) {
    size_t esz = tk_tag_size(X->tag);
    tk_csr_vals_grow(L, X, base + nb);
    memcpy((char *) tk_csr_val_ptr(X) + base * esz, tk_csr_val_ptr(B), nb * esz);
  }
  uint64_t bd = tk_csr_rows(B);
  for (uint64_t d = 1; d <= bd; d ++)
    tk_ivec_push(X->offsets, (int64_t) (base + (uint64_t) B->offsets->a[d]));
  lua_settop(L, 1);
  return 1;
}


static int tk_csr_clone_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  uint64_t no = X->offsets->n;
  uint64_t nn = tk_csr_nnz(X);
  tk_ivec_t *off = tk_ivec_create(L, no);
  int io = lua_gettop(L);
  memcpy(off->a, X->offsets->a, no * sizeof(int64_t));
  void *nbr = tk_csr_new_nbr(L, X->ntag, nn);
  int in_ = lua_gettop(L);
  memcpy(tk_nbr_aptr(nbr, X->ntag), tk_csr_nbr_ptr(X), nn * tk_nbr_esz(X->ntag));
  void *vals = NULL;
  int iv = 0;
  if (X->tag != TK_TAG_NONE) {
    vals = tk_csr_new_values(L, X->tag, nn);
    iv = lua_gettop(L);
    char *vdst;
    switch (X->tag) {
      case TK_TAG_I32: vdst = (char *) ((tk_svec_t *) vals)->a; break;
      case TK_TAG_I64: vdst = (char *) ((tk_ivec_t *) vals)->a; break;
      case TK_TAG_F32: vdst = (char *) ((tk_fvec_t *) vals)->a; break;
      case TK_TAG_F64: vdst = (char *) ((tk_dvec_t *) vals)->a; break;
      default: vdst = ((tk_cvec_t *) vals)->a; break;
    }
    memcpy(vdst, tk_csr_val_ptr(X), nn * tk_tag_size(X->tag));
  }
  tk_csr_push(L, X->tag, X->ntag, X->n_cols, io, off, in_, nbr, iv, vals);
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
  void *nbr = tk_csr_new_nbr(L, X->ntag, tk_csr_nbr_n(X));
  int in_ = lua_gettop(L);
  void *vals = NULL;
  int iv = 0;
  if (X->tag != TK_TAG_NONE) {
    vals = tk_csr_new_values(L, X->tag, tk_csr_nbr_n(X));
    iv = lua_gettop(L);
  }
  off->a[0] = 0;
  uint64_t pos = 0;
  for (uint64_t r = 0; r < n_rows; r ++) {
    for (int64_t j = X->offsets->a[r]; j < X->offsets->a[r + 1]; j ++) {
      int64_t new_id = tk_iumap_get_or(inverse, tk_csr_nbr(X, (uint64_t) j), -1);
      if (new_id >= 0) {
        if (X->tag != TK_TAG_NONE)
          switch (X->tag) {
            case TK_TAG_I32: ((tk_svec_t *) vals)->a[pos] = ((tk_svec_t *) X->values)->a[j]; break;
            case TK_TAG_I64: ((tk_ivec_t *) vals)->a[pos] = ((tk_ivec_t *) X->values)->a[j]; break;
            case TK_TAG_F32: ((tk_fvec_t *) vals)->a[pos] = ((tk_fvec_t *) X->values)->a[j]; break;
            case TK_TAG_F64: ((tk_dvec_t *) vals)->a[pos] = ((tk_dvec_t *) X->values)->a[j]; break;
            default: ((tk_cvec_t *) vals)->a[pos] = ((tk_cvec_t *) X->values)->a[j]; break;
          }
        tk_nbr_set(nbr, X->ntag, pos ++, new_id);
      }
    }
    off->a[r + 1] = (int64_t) pos;
  }
  tk_nbr_setn(nbr, X->ntag, pos);
  if (X->tag != TK_TAG_NONE)
    switch (X->tag) {
      case TK_TAG_I32: ((tk_svec_t *) vals)->n = pos; break;
      case TK_TAG_I64: ((tk_ivec_t *) vals)->n = pos; break;
      case TK_TAG_F32: ((tk_fvec_t *) vals)->n = pos; break;
      case TK_TAG_F64: ((tk_dvec_t *) vals)->n = pos; break;
      default: ((tk_cvec_t *) vals)->n = pos; break;
    }
  tk_iumap_destroy(inverse);
  tk_csr_push(L, X->tag, X->ntag, remap->n, io, off, in_, nbr, iv, vals);
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
  if (X->ntag != Y->ntag)
    return tk_lua_verror(L, 2, "csr", "hcat requires matching neighbor types");
  uint64_t n_rows = tk_csr_rows(X);
  if (n_rows != tk_csr_rows(Y))
    return tk_lua_verror(L, 2, "csr", "hcat requires matching row counts");
  uint64_t na = tk_csr_nbr_n(X), nb = tk_csr_nbr_n(Y), total = na + nb;
  if (tk_csr_nbr_ensure(X, total) != 0)
    return tk_lua_verror(L, 2, "csr", "allocation failed");
  tk_csr_nbr_setn(X, total);
  if (X->tag != TK_TAG_NONE)
    tk_csr_vals_grow(L, X, total);
  size_t esz = tk_tag_size(X->tag);
  size_t nesz = tk_nbr_esz(X->ntag);
  char *ndst = (char *) tk_csr_nbr_ptr(X);
  char *vdst = X->tag != TK_TAG_NONE ? (char *) tk_csr_val_ptr(X) : NULL;
  const char *vsrc = Y->tag != TK_TAG_NONE ? (const char *) tk_csr_val_ptr(Y) : NULL;
  int64_t shift = (int64_t) X->n_cols;
  for (uint64_t r = n_rows; r > 0; r --) {
    uint64_t i = r - 1;
    int64_t alo = X->offsets->a[i], ahi = X->offsets->a[i + 1];
    int64_t blo = Y->offsets->a[i], bhi = Y->offsets->a[i + 1];
    int64_t dst = alo + blo;
    memmove(ndst + (uint64_t) dst * nesz, ndst + (uint64_t) alo * nesz, (size_t) (ahi - alo) * nesz);
    if (vdst)
      memmove(vdst + (uint64_t) dst * esz, vdst + (uint64_t) alo * esz, (size_t) (ahi - alo) * esz);
    int64_t bdst = dst + (ahi - alo);
    for (int64_t j = blo; j < bhi; j ++)
      tk_csr_setnbr(X, (uint64_t) (bdst + (j - blo)), tk_csr_nbr(Y, (uint64_t) j) + shift);
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
  uint64_t nnz = tk_csr_nbr_n(X);
  int64_t *counts = (int64_t *) calloc(X->n_cols + 1, sizeof(int64_t));
  if (!counts)
    return tk_lua_verror(L, 2, "csr", "allocation failed");
  for (uint64_t i = 0; i < nnz; i ++)
    counts[tk_csr_nbr(X, i) + 1] ++;
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
      int64_t tok = tk_csr_nbr(X, (uint64_t) j);
      int64_t pos = counts[tok] ++;
      nbr->a[pos] = (int64_t) s;
      if (vdst)
        memcpy(vdst + (uint64_t) pos * esz, vsrc + (uint64_t) j * esz, esz);
    }
  free(counts);
  tk_csr_push(L, X->tag, TK_TAG_I64, n_rows, io, off, in_, nbr, iv, vals);
  return 1;
}


static inline void tk_csr_materialize (lua_State *L, tk_csr_t *X, int ix)
{
  if (X->tag != TK_TAG_NONE)
    return;
  uint64_t nn = tk_csr_nbr_n(X);
  tk_fvec_t *vals = tk_fvec_create(L, nn);
  for (uint64_t i = 0; i < nn; i ++)
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
  uint64_t nn = tk_csr_nbr_n(X);
  for (uint64_t i = 0; i < nn; i ++) {
    int64_t c = tk_csr_nbr(X, i);
    double w = wf != NULL ? (double) wf->a[c] : wd->a[c];
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
  uint64_t nn = tk_csr_nbr_n(X);
  if (bounds == NULL) {
    tk_dvec_t *out = tk_dvec_create(L, X->n_cols);
    memset(out->a, 0, X->n_cols * sizeof(double));
    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < nn; i ++) {
      double v = X->tag == TK_TAG_NONE ? 1.0 : tk_csr_val1(X, i);
      double vv = v * v;
      int64_t c = tk_csr_nbr(X, i);
      #pragma omp atomic
      out->a[c] += vv;
    }
    return 1;
  }
  uint64_t nb = bounds->n > 0 ? bounds->n - 1 : 0;
  tk_dvec_t *out = tk_dvec_create(L, nb);
  memset(out->a, 0, nb * sizeof(double));
  for (uint64_t i = 0; i < nn; i ++) {
    int64_t c = tk_csr_nbr(X, i);

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

static int tk_csr_nnz_cols_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_ivec_t *bounds = lua_isnil(L, 2) ? NULL : tk_ivec_peek(L, 2, "bounds");
  uint64_t nn = tk_csr_nbr_n(X);
  if (bounds == NULL) {
    tk_ivec_t *out = tk_ivec_create(L, X->n_cols);
    memset(out->a, 0, X->n_cols * sizeof(int64_t));
    #pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < nn; i ++) {
      int64_t c = tk_csr_nbr(X, i);
      #pragma omp atomic
      out->a[c] += 1;
    }
    return 1;
  }
  uint64_t nb = bounds->n > 0 ? bounds->n - 1 : 0;
  tk_ivec_t *out = tk_ivec_create(L, nb);
  memset(out->a, 0, nb * sizeof(int64_t));
  for (uint64_t i = 0; i < nn; i ++) {
    int64_t c = tk_csr_nbr(X, i);
    uint64_t lo = 0, hi = nb;
    while (lo + 1 < hi) {
      uint64_t mid = (lo + hi) / 2;
      if (bounds->a[mid] <= c) lo = mid; else hi = mid;
    }
    if (c >= bounds->a[lo] && c < bounds->a[lo + 1])
      out->a[lo] += 1;
  }
  return 1;
}

static inline void tk_csr_scale_by_cols (tk_csr_t *X, tk_fvec_t *wf, tk_dvec_t *wd)
{
  uint64_t nn = tk_csr_nbr_n(X);
  for (uint64_t i = 0; i < nn; i ++) {
    int64_t c = tk_csr_nbr(X, i);
    double w = wf != NULL ? (double) wf->a[c] : wd->a[c];
    tk_csr_setval1(X, i, tk_csr_val1(X, i) * w);
  }
}

static inline double tk_csr_probit (double p)
{
  if (p <= 0.0) return -1e10;
  if (p >= 1.0) return 1e10;
  static const double a[] = {
    -3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
     1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00 };
  static const double b[] = {
    -5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
     6.680131188771972e+01, -1.328068155288572e+01 };
  static const double c[] = {
    -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
    -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
  static const double d[] = {
    7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
    3.754408661907416e+00 };
  double plow = 0.02425, phigh = 1.0 - plow;
  double q, r;
  if (p < plow) {
    q = sqrt(-2.0 * log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
  } else if (p <= phigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
  } else {
    q = sqrt(-2.0 * log(1.0 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
  }
}

#define TK_CSR_BNS_EPS 0.5

static inline double tk_csr_bns_score (double N, double C, double P, double A)
{
  if (C <= 0 || C >= N || P <= 0 || P >= N) return 0.0;
  double tpr = (A + TK_CSR_BNS_EPS) / (P + 2.0 * TK_CSR_BNS_EPS);
  double fpr = (C - A + TK_CSR_BNS_EPS) / (N - P + 2.0 * TK_CSR_BNS_EPS);
  return fabs(tk_csr_probit(tpr) - tk_csr_probit(fpr));
}





static int tk_csr_bns_lua (lua_State *L)
{
  lua_settop(L, 3);
  int noscale = lua_toboolean(L, 3);
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_csr_t *Y = tk_csr_peekopt(L, 2);
  tk_csr_materialize(L, X, 1);
  if (Y == NULL) {
    tk_fvec_t *wf = tk_fvec_peekopt(L, 2);
    tk_dvec_t *wd = wf == NULL ? tk_dvec_peek(L, 2, "weights") : NULL;
    uint64_t wn = wf != NULL ? wf->n : wd->n;
    if (wn < X->n_cols)
      return tk_lua_verror(L, 2, "csr", "bns: weights shorter than n_cols");
    tk_csr_scale_by_cols(X, wf, wd);
    lua_pushvalue(L, 2);
    return 1;
  }
  uint64_t nc = X->n_cols, n_labels = Y->n_cols, n_rows = tk_csr_rows(X);
  if (Y->offsets->n != X->offsets->n)
    return tk_lua_verror(L, 2, "csr", "bns: labels row count mismatch");
  double N = (double) n_rows;
  uint32_t *doc_freq = (uint32_t *) calloc(nc, sizeof(uint32_t));
  uint32_t *label_freq = (uint32_t *) calloc(n_labels, sizeof(uint32_t));

  int64_t *lbl_off = (int64_t *) malloc((n_labels + 1) * sizeof(int64_t));
  if (!doc_freq || !label_freq || !lbl_off) {
    free(doc_freq); free(label_freq); free(lbl_off);
    return tk_lua_verror(L, 2, "csr", "bns: alloc failed");
  }
  uint64_t xnn = tk_csr_nbr_n(X);
  for (uint64_t j = 0; j < xnn; j ++)
    doc_freq[tk_csr_nbr(X, j)] ++;
  for (uint64_t d = 0; d < n_rows; d ++)
    for (int64_t j = Y->offsets->a[d]; j < Y->offsets->a[d + 1]; j ++) {
      uint64_t b = (uint64_t) tk_csr_nbr(Y, (uint64_t) j);
      if (b < n_labels) label_freq[b] ++;
    }
  lbl_off[0] = 0;
  for (uint64_t b = 0; b < n_labels; b ++)
    lbl_off[b + 1] = lbl_off[b] + label_freq[b];
  uint32_t *lbl_docs = (uint32_t *) malloc((uint64_t) lbl_off[n_labels] * sizeof(uint32_t));
  uint32_t *lbl_pos = (uint32_t *) calloc(n_labels, sizeof(uint32_t));
  if (!lbl_docs || !lbl_pos) {
    free(doc_freq); free(label_freq); free(lbl_off); free(lbl_docs); free(lbl_pos);
    return tk_lua_verror(L, 2, "csr", "bns: alloc failed");
  }
  for (uint64_t d = 0; d < n_rows; d ++)
    for (int64_t j = Y->offsets->a[d]; j < Y->offsets->a[d + 1]; j ++) {
      uint64_t b = (uint64_t) tk_csr_nbr(Y, (uint64_t) j);
      if (b < n_labels) lbl_docs[lbl_off[b] + lbl_pos[b] ++] = (uint32_t) d;
    }
  free(lbl_pos);
  tk_fvec_t *w = tk_fvec_create(L, nc);
  w->n = nc;
  memset(w->a, 0, nc * sizeof(float));
  float *cooc = (float *) calloc(nc, sizeof(float));
  int32_t *touched = (int32_t *) malloc(nc * sizeof(int32_t));
  if (!cooc || !touched) {
    free(doc_freq); free(label_freq); free(lbl_off); free(lbl_docs); free(cooc); free(touched);
    return tk_lua_verror(L, 2, "csr", "bns: alloc failed");
  }
  for (uint64_t b = 0; b < n_labels; b ++) {
    double P = (double) label_freq[b];
    if (P <= 0.0 || P >= N) continue;
    uint32_t n_touched = 0;
    for (int64_t di = lbl_off[b]; di < lbl_off[b + 1]; di ++) {
      uint32_t dd = lbl_docs[di];
      for (int64_t j = X->offsets->a[dd]; j < X->offsets->a[dd + 1]; j ++) {
        int32_t f = (int32_t) tk_csr_nbr(X, (uint64_t) j);
        if (cooc[f] == 0.0f) touched[n_touched ++] = f;
        cooc[f] += 1.0f;
      }
    }
    for (uint32_t i = 0; i < n_touched; i ++) {
      int32_t f = touched[i];
      float sc = (float) tk_csr_bns_score(N, (double) doc_freq[f], P, (double) cooc[f]);
      if (sc > w->a[f]) w->a[f] = sc;
      cooc[f] = 0.0f;
    }
  }
  free(cooc); free(touched);
  free(doc_freq); free(label_freq); free(lbl_off); free(lbl_docs);
  if (!noscale)
    tk_csr_scale_by_cols(X, w, NULL);
  return 1;
}

#define TK_CSR_AUC_EPS 1e-4






static int tk_csr_auc_spearman (lua_State *L, tk_csr_t *X, tk_dvec_t *Yt)
{
  tk_csr_materialize(L, X, 1);
  uint64_t nc = X->n_cols, n_rows = tk_csr_rows(X);
  if (Yt->n < n_rows)
    return tk_lua_verror(L, 2, "csr", "auc: targets shorter than rows");
  double N = (double) n_rows;
  uint64_t nn = tk_csr_nbr_n(X);
  int has_vals = X->tag != TK_TAG_NONE;
  uint32_t *doc_freq = (uint32_t *) calloc(nc ? nc : 1, sizeof(uint32_t));
  uint64_t *col_off = (uint64_t *) malloc((nc + 1) * sizeof(uint64_t));
  uint32_t *col_cur = (uint32_t *) calloc(nc ? nc : 1, sizeof(uint32_t));
  uint32_t *row_of = (uint32_t *) malloc((nn ? nn : 1) * sizeof(uint32_t));
  float *val_of = (float *) malloc((nn ? nn : 1) * sizeof(float));
  double *yrk = (double *) malloc((n_rows ? n_rows : 1) * sizeof(double));
  if (!doc_freq || !col_off || !col_cur || !row_of || !val_of || !yrk) {
    free(doc_freq); free(col_off); free(col_cur);
    free(row_of); free(val_of); free(yrk);
    return tk_lua_verror(L, 2, "csr", "auc: alloc failed");
  }
  for (uint64_t j = 0; j < nn; j ++)
    doc_freq[tk_csr_nbr(X, j)] ++;
  col_off[0] = 0;
  for (uint64_t c = 0; c < nc; c ++)
    col_off[c + 1] = col_off[c] + doc_freq[c];
  for (uint64_t d = 0; d < n_rows; d ++)
    for (int64_t j = X->offsets->a[d]; j < X->offsets->a[d + 1]; j ++) {
      int64_t c = tk_csr_nbr(X, (uint64_t) j);
      uint64_t p = col_off[c] + col_cur[c] ++;
      row_of[p] = (uint32_t) d;
      val_of[p] = has_vals ? (float) tk_csr_val1(X, (uint64_t) j) : 1.0f;
    }
  tk_rvec_t *Sy = tk_rvec_create(L, n_rows ? n_rows : 1);
  for (uint64_t d = 0; d < n_rows; d ++)
    tk_rvec_push(Sy, tk_rank((int64_t) d, Yt->a[d]));
  tk_rvec_asc(Sy, 0, Sy->n);
  {
    uint64_t p = 0;
    while (p < n_rows) {
      double v = Sy->a[p].d;
      uint64_t q = p;
      while (q + 1 < n_rows && Sy->a[q + 1].d == v) q ++;
      double mr = ((double) (p + 1) + (double) (q + 1)) / 2.0;
      for (uint64_t t = p; t <= q; t ++) yrk[Sy->a[t].i] = mr;
      p = q + 1;
    }
  }
  double mean = (N + 1.0) / 2.0;
  double Sy_sum = N * (N + 1.0) / 2.0;
  double Syy = 0.0;
  for (uint64_t d = 0; d < n_rows; d ++) Syy += yrk[d] * yrk[d];
  double var_y = Syy / N - mean * mean;
  tk_fvec_t *w = tk_fvec_create(L, nc);
  int iw = lua_gettop(L);
  w->n = nc;
  memset(w->a, 0, nc * sizeof(float));
  uint64_t maxcol = 0;
  for (uint64_t c = 0; c < nc; c ++) if (doc_freq[c] > maxcol) maxcol = doc_freq[c];
  tk_rvec_t *S = tk_rvec_create(L, maxcol ? maxcol : 1);
  double *rk = (double *) malloc((maxcol ? maxcol : 1) * sizeof(double));
  if (!rk) {
    free(doc_freq); free(col_off); free(col_cur);
    free(row_of); free(val_of); free(yrk);
    return tk_lua_verror(L, 2, "csr", "auc: alloc failed");
  }
  for (uint64_t c = 0; c < nc; c ++) {
    uint64_t cn = doc_freq[c];
    tk_rvec_clear(S);
    for (uint64_t p = 0; p < cn; p ++)
      tk_rvec_push(S, tk_rank((int64_t) row_of[col_off[c] + p], (double) val_of[col_off[c] + p]));
    tk_rvec_asc(S, 0, S->n);
    uint64_t neg = 0, ez = 0;
    for (uint64_t p = 0; p < cn; p ++) {
      if (S->a[p].d < 0.0) neg ++;
      else if (S->a[p].d == 0.0) ez ++;
      else break;
    }
    uint64_t z = (n_rows - cn) + ez;
    double midrank_zero = (double) neg + ((double) z + 1.0) / 2.0;
    uint64_t p = 0;
    while (p < cn) {
      double v = S->a[p].d;
      uint64_t q = p;
      while (q + 1 < cn && S->a[q + 1].d == v) q ++;
      double mr;
      if (v == 0.0) mr = midrank_zero;
      else if (v < 0.0) mr = ((double) (p + 1) + (double) (q + 1)) / 2.0;
      else mr = ((double) (p + 1 - ez) + (double) z + (double) (q + 1 - ez) + (double) z) / 2.0;
      for (uint64_t t = p; t <= q; t ++) rk[t] = mr;
      p = q + 1;
    }
    double sxy = 0.0, sxx = 0.0, sy_expl = 0.0;
    for (uint64_t t = 0; t < cn; t ++) {
      double yr = yrk[S->a[t].i];
      sxy += rk[t] * yr;
      sxx += rk[t] * rk[t];
      sy_expl += yr;
    }
    sxy += midrank_zero * (Sy_sum - sy_expl);
    sxx += (double) (n_rows - cn) * midrank_zero * midrank_zero;
    double var_x = sxx / N - mean * mean;
    if (var_x > 0.0 && var_y > 0.0) {
      double rho = (sxy / N - mean * mean) / sqrt(var_x * var_y);
      double auc = (rho + 1.0) / 2.0;
      auc = (auc + TK_CSR_AUC_EPS) / (1.0 + 2.0 * TK_CSR_AUC_EPS);
      w->a[c] = (float) (fabs(tk_csr_probit(auc)) * M_SQRT2);
    }
  }
  free(rk);
  free(doc_freq); free(col_off); free(col_cur);
  free(row_of); free(val_of); free(yrk);
  lua_pushvalue(L, iw);
  return 1;
}









static int tk_csr_auc_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_dvec_t *Yt = tk_dvec_peekopt(L, 2);
  if (Yt != NULL)
    return tk_csr_auc_spearman(L, X, Yt);
  tk_csr_t *Y = tk_csr_peek(L, 2, "labels");
  tk_csr_materialize(L, X, 1);
  uint64_t nc = X->n_cols, n_labels = Y->n_cols, n_rows = tk_csr_rows(X);
  if (Y->offsets->n != X->offsets->n)
    return tk_lua_verror(L, 2, "csr", "auc: labels row count mismatch");
  double N = (double) n_rows;
  uint64_t nn = tk_csr_nbr_n(X);
  int has_vals = X->tag != TK_TAG_NONE;
  uint32_t *doc_freq = (uint32_t *) calloc(nc ? nc : 1, sizeof(uint32_t));
  uint32_t *label_freq = (uint32_t *) calloc(n_labels ? n_labels : 1, sizeof(uint32_t));
  uint64_t *col_off = (uint64_t *) malloc((nc + 1) * sizeof(uint64_t));
  uint32_t *col_cur = (uint32_t *) calloc(nc ? nc : 1, sizeof(uint32_t));
  uint32_t *row_of = (uint32_t *) malloc((nn ? nn : 1) * sizeof(uint32_t));
  float *val_of = (float *) malloc((nn ? nn : 1) * sizeof(float));
  double *rank_sum = (double *) calloc(n_labels ? n_labels : 1, sizeof(double));
  uint32_t *cnt = (uint32_t *) calloc(n_labels ? n_labels : 1, sizeof(uint32_t));
  if (!doc_freq || !label_freq || !col_off || !col_cur || !row_of || !val_of || !rank_sum || !cnt) {
    free(doc_freq); free(label_freq); free(col_off); free(col_cur);
    free(row_of); free(val_of); free(rank_sum); free(cnt);
    return tk_lua_verror(L, 2, "csr", "auc: alloc failed");
  }
  for (uint64_t j = 0; j < nn; j ++)
    doc_freq[tk_csr_nbr(X, j)] ++;
  col_off[0] = 0;
  for (uint64_t c = 0; c < nc; c ++)
    col_off[c + 1] = col_off[c] + doc_freq[c];
  for (uint64_t d = 0; d < n_rows; d ++)
    for (int64_t j = X->offsets->a[d]; j < X->offsets->a[d + 1]; j ++) {
      int64_t c = tk_csr_nbr(X, (uint64_t) j);
      uint64_t p = col_off[c] + col_cur[c] ++;
      row_of[p] = (uint32_t) d;
      val_of[p] = has_vals ? (float) tk_csr_val1(X, (uint64_t) j) : 1.0f;
    }
  for (uint64_t d = 0; d < n_rows; d ++)
    for (int64_t j = Y->offsets->a[d]; j < Y->offsets->a[d + 1]; j ++) {
      uint64_t b = (uint64_t) tk_csr_nbr(Y, (uint64_t) j);
      if (b < n_labels) label_freq[b] ++;
    }
  tk_fvec_t *w = tk_fvec_create(L, nc);
  int iw = lua_gettop(L);
  w->n = nc;
  memset(w->a, 0, nc * sizeof(float));
  uint64_t maxcol = 0;
  for (uint64_t c = 0; c < nc; c ++) if (doc_freq[c] > maxcol) maxcol = doc_freq[c];
  tk_rvec_t *S = tk_rvec_create(L, maxcol ? maxcol : 1);
  double *rk = (double *) malloc((maxcol ? maxcol : 1) * sizeof(double));
  if (!rk) {
    free(doc_freq); free(label_freq); free(col_off); free(col_cur);
    free(row_of); free(val_of); free(rank_sum); free(cnt);
    return tk_lua_verror(L, 2, "csr", "auc: alloc failed");
  }
  for (uint64_t c = 0; c < nc; c ++) {
    uint64_t cn = doc_freq[c];
    tk_rvec_clear(S);
    for (uint64_t p = 0; p < cn; p ++)
      tk_rvec_push(S, tk_rank((int64_t) row_of[col_off[c] + p], (double) val_of[col_off[c] + p]));
    tk_rvec_asc(S, 0, S->n);
    uint64_t neg = 0, ez = 0;
    for (uint64_t p = 0; p < cn; p ++) {
      if (S->a[p].d < 0.0) neg ++;
      else if (S->a[p].d == 0.0) ez ++;
      else break;
    }
    uint64_t z = (n_rows - cn) + ez;
    double midrank_zero = (double) neg + ((double) z + 1.0) / 2.0;
    uint64_t p = 0;
    while (p < cn) {
      double v = S->a[p].d;
      uint64_t q = p;
      while (q + 1 < cn && S->a[q + 1].d == v) q ++;
      double mr;
      if (v == 0.0) mr = midrank_zero;
      else if (v < 0.0) mr = ((double) (p + 1) + (double) (q + 1)) / 2.0;
      else mr = ((double) (p + 1 - ez) + (double) z + (double) (q + 1 - ez) + (double) z) / 2.0;
      for (uint64_t t = p; t <= q; t ++) rk[t] = mr;
      p = q + 1;
    }
    for (uint64_t t = 0; t < cn; t ++) {
      uint32_t d = (uint32_t) S->a[t].i;
      for (int64_t j = Y->offsets->a[d]; j < Y->offsets->a[d + 1]; j ++) {
        uint64_t b = (uint64_t) tk_csr_nbr(Y, (uint64_t) j);
        if (b < n_labels) { rank_sum[b] += rk[t]; cnt[b] ++; }
      }
    }
    float best = 0.0f;
    for (uint64_t b = 0; b < n_labels; b ++) {
      double P = (double) label_freq[b];
      if (P > 0.0 && P < N) {
        double rs = rank_sum[b] + (P - (double) cnt[b]) * midrank_zero;
        double u = rs - P * (P + 1.0) / 2.0;
        double auc = u / (P * (N - P));
        auc = (auc + TK_CSR_AUC_EPS) / (1.0 + 2.0 * TK_CSR_AUC_EPS);
        float sc = (float) (fabs(tk_csr_probit(auc)) * M_SQRT2);
        if (sc > best) best = sc;
      }
      rank_sum[b] = 0.0; cnt[b] = 0;
    }
    w->a[c] = best;
  }
  free(rk);
  free(doc_freq); free(label_freq); free(col_off); free(col_cur);
  free(row_of); free(val_of); free(rank_sum); free(cnt);
  lua_pushvalue(L, iw);
  return 1;
}




static int tk_csr_standardize_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_fvec_t *wf = tk_fvec_peekopt(L, 2);
  tk_dvec_t *wd = wf == NULL ? tk_dvec_peekopt(L, 2) : NULL;
  tk_csr_materialize(L, X, 1);
  if (wf != NULL || wd != NULL) {
    uint64_t wn = wf != NULL ? wf->n : wd->n;
    if (wn < X->n_cols)
      return tk_lua_verror(L, 2, "csr", "standardize: weights shorter than n_cols");
    tk_csr_scale_by_cols(X, wf, wd);
    lua_pushvalue(L, 2);
    return 1;
  }
  uint64_t nc = X->n_cols, n_rows = tk_csr_rows(X);
  double *sum = (double *) calloc(nc, sizeof(double));
  double *ssq = (double *) calloc(nc, sizeof(double));
  if (!sum || !ssq) { free(sum); free(ssq); return tk_lua_verror(L, 2, "csr", "standardize: alloc failed"); }
  uint64_t nn = tk_csr_nbr_n(X);
  for (uint64_t i = 0; i < nn; i ++) {
    double v = tk_csr_val1(X, i);
    int64_t c = tk_csr_nbr(X, i);
    sum[c] += v; ssq[c] += v * v;
  }
  tk_fvec_t *w = tk_fvec_create(L, nc);
  w->n = nc;
  double n = (double) n_rows;
  for (uint64_t c = 0; c < nc; c ++) {
    double mean = n > 0 ? sum[c] / n : 0.0;
    double var = n > 0 ? ssq[c] / n - mean * mean : 0.0;
    double sd = sqrt(var > 0.0 ? var : 0.0);
    w->a[c] = sd > 1e-10 ? (float) (1.0 / sd) : 0.0f;
  }
  free(sum); free(ssq);
  tk_csr_scale_by_cols(X, w, NULL);
  return 1;
}




static int tk_csr_idf_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  tk_fvec_t *wf = tk_fvec_peekopt(L, 2);
  tk_dvec_t *wd = wf == NULL ? tk_dvec_peekopt(L, 2) : NULL;
  tk_csr_materialize(L, X, 1);
  if (wf != NULL || wd != NULL) {
    uint64_t wn = wf != NULL ? wf->n : wd->n;
    if (wn < X->n_cols)
      return tk_lua_verror(L, 2, "csr", "idf: weights shorter than n_cols");
    tk_csr_scale_by_cols(X, wf, wd);
    lua_pushvalue(L, 2);
    return 1;
  }
  uint64_t nc = X->n_cols, n_rows = tk_csr_rows(X);
  uint32_t *df = (uint32_t *) calloc(nc, sizeof(uint32_t));
  if (!df) return tk_lua_verror(L, 2, "csr", "idf: alloc failed");
  uint64_t nn = tk_csr_nbr_n(X);
  for (uint64_t i = 0; i < nn; i ++)
    df[tk_csr_nbr(X, i)] ++;
  tk_fvec_t *w = tk_fvec_create(L, nc);
  w->n = nc;
  double N = (double) n_rows;
  for (uint64_t c = 0; c < nc; c ++) {
    double d = (double) df[c];
    w->a[c] = (float) log((N - d + 0.5) / (d + 0.5));
  }
  free(df);
  tk_csr_scale_by_cols(X, w, NULL);
  return 1;
}

static int tk_csr_eq_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *a = tk_csr_peek(L, 1, "csr");
  tk_csr_t *b = tk_csr_peek(L, 2, "other");
  bool r = a->tag == b->tag && a->ntag == b->ntag && a->n_cols == b->n_cols
    && a->offsets->n == b->offsets->n && tk_csr_nbr_n(a) == tk_csr_nbr_n(b)
    && tk_ivec_eq(a->offsets, b->offsets, 0, a->offsets->n)
    && memcmp(tk_csr_nbr_ptr(a), tk_csr_nbr_ptr(b), tk_nbr_esz(a->ntag) * tk_csr_nbr_n(a)) == 0;
  if (r && a->tag != TK_TAG_NONE)
    r = memcmp(tk_csr_val_ptr(a), tk_csr_val_ptr(b), tk_tag_size(a->tag) * tk_csr_nbr_n(a)) == 0;
  lua_pushboolean(L, r);
  return 1;
}

static int tk_csr_persist_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  FILE *fh = tk_lua_fopen(L, luaL_checkstring(L, 2), "w");
  char magic[4] = { 'T', 'K', 'c', 's' };
  uint8_t version = 2;
  uint8_t tag = (uint8_t) X->tag;
  uint8_t ntag = (uint8_t) X->ntag;
  uint64_t no = X->offsets->n, nn = tk_csr_nbr_n(X);
  tk_lua_fwrite(L, magic, 4, 1, fh);
  tk_lua_fwrite(L, (char *) &version, 1, 1, fh);
  tk_lua_fwrite(L, (char *) &tag, 1, 1, fh);
  tk_lua_fwrite(L, (char *) &ntag, 1, 1, fh);
  tk_lua_fwrite(L, (char *) &X->n_cols, sizeof(uint64_t), 1, fh);
  tk_lua_fwrite(L, (char *) &no, sizeof(uint64_t), 1, fh);
  tk_lua_fwrite(L, (char *) &nn, sizeof(uint64_t), 1, fh);
  tk_lua_fwrite(L, (char *) X->offsets->a, sizeof(int64_t) * no, 1, fh);
  tk_lua_fwrite(L, (char *) tk_csr_nbr_ptr(X), tk_nbr_esz(X->ntag) * nn, 1, fh);
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
  uint8_t version, tag8, ntag8 = (uint8_t) TK_TAG_I64;
  uint64_t n_cols, no, nn;
  tk_lua_fread(L, magic, 4, 1, fh);
  if (memcmp(magic, "TKcs", 4) != 0)
    return tk_lua_verror(L, 2, "csr", "load: bad magic");
  tk_lua_fread(L, (char *) &version, 1, 1, fh);
  if (version != 1 && version != 2)
    return tk_lua_verror(L, 2, "csr", "load: unsupported version");
  tk_lua_fread(L, (char *) &tag8, 1, 1, fh);
  if (version >= 2)
    tk_lua_fread(L, (char *) &ntag8, 1, 1, fh);
  tk_lua_fread(L, (char *) &n_cols, sizeof(uint64_t), 1, fh);
  tk_lua_fread(L, (char *) &no, sizeof(uint64_t), 1, fh);
  tk_lua_fread(L, (char *) &nn, sizeof(uint64_t), 1, fh);
  tk_ivec_t *off = tk_ivec_create(L, no);
  int io = lua_gettop(L);
  tk_lua_fread(L, (char *) off->a, sizeof(int64_t) * no, 1, fh);
  void *nbr = tk_csr_new_nbr(L, (tk_tag_t) ntag8, nn);
  int in_ = lua_gettop(L);
  tk_lua_fread(L, (char *) tk_nbr_aptr(nbr, (tk_tag_t) ntag8), tk_nbr_esz((tk_tag_t) ntag8) * nn, 1, fh);
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
  tk_csr_push(L, (tk_tag_t) tag8, (tk_tag_t) ntag8, n_cols, io, off, in_, nbr, iv, vals);
  return 1;
}

static int tk_csr_i32_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");
  if (X->ntag == TK_TAG_I32) { lua_pushvalue(L, 1); return 1; }
  uint64_t n_off = X->offsets->n;
  uint64_t total = tk_csr_nbr_n(X);
  tk_ivec_t *off = tk_ivec_create(L, n_off);
  int io = lua_gettop(L);
  for (uint64_t i = 0; i < n_off; i ++) off->a[i] = X->offsets->a[i];
  off->n = n_off;
  void *nbr = tk_csr_new_nbr(L, TK_TAG_I32, total);
  int in_ = lua_gettop(L);
  tk_nbr_setn(nbr, TK_TAG_I32, total);
  const int64_t *nsrc = (const int64_t *) tk_csr_nbr_ptr(X);
  int32_t *ndst = (int32_t *) tk_nbr_aptr(nbr, TK_TAG_I32);
  for (uint64_t i = 0; i < total; i ++) {
    int64_t v = nsrc[i];
    if (v < 0 || v > (int64_t) INT32_MAX)
      return tk_lua_verror(L, 2, "csr", "i32: neighbor index exceeds int32 range");
    ndst[i] = (int32_t) v;
  }
  void *vals = NULL;
  int iv = 0;
  if (X->tag != TK_TAG_NONE) {
    vals = tk_csr_new_values(L, X->tag, total);
    iv = lua_gettop(L);
    char *vdst;
    switch (X->tag) {
      case TK_TAG_I32: vdst = (char *) ((tk_svec_t *) vals)->a; break;
      case TK_TAG_I64: vdst = (char *) ((tk_ivec_t *) vals)->a; break;
      case TK_TAG_F32: vdst = (char *) ((tk_fvec_t *) vals)->a; break;
      case TK_TAG_F64: vdst = (char *) ((tk_dvec_t *) vals)->a; break;
      default: vdst = ((tk_cvec_t *) vals)->a; break;
    }
    memcpy(vdst, tk_csr_val_ptr(X), tk_tag_size(X->tag) * total);
  }
  tk_csr_push(L, X->tag, TK_TAG_I32, X->n_cols, io, off, in_, nbr, iv, vals);
  return 1;
}

static int tk_csr_clear_lua (lua_State *L)
{
  tk_csr_t *X = tk_csr_peek(L, 1, "csr");


  if (X->offsets->m >= 1) { X->offsets->a[0] = 0; X->offsets->n = 1; }
  else X->offsets->n = 0;
  tk_csr_nbr_setn(X, 0);
  if (X->values) switch (X->tag) {
    case TK_TAG_I32: ((tk_svec_t *) X->values)->n = 0; break;
    case TK_TAG_I64: ((tk_ivec_t *) X->values)->n = 0; break;
    case TK_TAG_F32: ((tk_fvec_t *) X->values)->n = 0; break;
    case TK_TAG_F64: ((tk_dvec_t *) X->values)->n = 0; break;
    case TK_TAG_U8:  ((tk_cvec_t *) X->values)->n = 0; break;
    default: break;
  }
  return 0;
}

static luaL_Reg tk_csr_mt_fns[] = {
  { "shape", tk_csr_shape_lua },
  { "clear", tk_csr_clear_lua },
  { "nnz", tk_csr_nnz_lua },
  { "type", tk_csr_type_lua },
  { "offsets", tk_csr_offsets_lua },
  { "neighbors", tk_csr_neighbors_lua },
  { "values", tk_csr_values_lua },
  { "push", tk_csr_push_lua },
  { "row", tk_csr_row_lua },
  { "rows", tk_csr_rows_lua },
  { "append", tk_csr_append_lua },
  { "clone", tk_csr_clone_lua },
  { "select", tk_csr_select_lua },
  { "hcat", tk_csr_hcat_lua },
  { "transpose", tk_csr_transpose_lua },
  { "normalize", tk_csr_normalize_lua },
  { "scale_cols", tk_csr_scale_cols_lua },
  { "sumsq_cols", tk_csr_sumsq_cols_lua },
  { "nnz_cols", tk_csr_nnz_cols_lua },
  { "standardize", tk_csr_standardize_lua },
  { "idf", tk_csr_idf_lua },
  { "bns", tk_csr_bns_lua },
  { "auc", tk_csr_auc_lua },
  { "to_bits", tk_csr_to_bits_lua },
  { "to_dense", tk_csr_to_dense_lua },
  { "i32", tk_csr_i32_lua },
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
