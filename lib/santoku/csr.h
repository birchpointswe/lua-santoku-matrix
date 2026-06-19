#ifndef TK_CSR_H
#define TK_CSR_H

#include <santoku/shaped.h>
#include <santoku/ivec.h>
#include <santoku/svec.h>
#include <santoku/dvec.h>
#include <santoku/fvec.h>
#include <santoku/cvec.h>

#define TK_CSR_MT "tk_csr_t"

typedef struct {
  tk_tag_t tag;
  tk_tag_t ntag;
  uint64_t n_cols;
  tk_ivec_t *offsets;
  void *neighbors;
  void *values;
} tk_csr_t;

static inline tk_csr_t *tk_csr_peek (lua_State *L, int i, const char *name)
{
  tk_csr_t *X = (tk_csr_t *) luaL_checkudata(L, i, TK_CSR_MT);
  if (X == NULL)
    tk_lua_verror(L, 2, name, "expected a csr");
  return X;
}

static inline tk_csr_t *tk_csr_peekopt (lua_State *L, int i)
{
  if (lua_type(L, i) != LUA_TUSERDATA)
    return NULL;
  void *p = lua_touserdata(L, i);
  if (!lua_getmetatable(L, i))
    return NULL;
  luaL_getmetatable(L, TK_CSR_MT);
  bool ok = lua_rawequal(L, -1, -2);
  lua_pop(L, 2);
  return ok ? (tk_csr_t *) p : NULL;
}

static inline uint64_t tk_csr_rows (tk_csr_t *X)
{
  return X->offsets->n > 0 ? X->offsets->n - 1 : 0;
}




static inline int64_t tk_nbr_get (void *nv, tk_tag_t ntag, uint64_t i) {
  return ntag == TK_TAG_I32 ? (int64_t) ((tk_svec_t *) nv)->a[i] : ((tk_ivec_t *) nv)->a[i];
}
static inline void tk_nbr_set (void *nv, tk_tag_t ntag, uint64_t i, int64_t v) {
  if (ntag == TK_TAG_I32) ((tk_svec_t *) nv)->a[i] = (int32_t) v;
  else ((tk_ivec_t *) nv)->a[i] = v;
}
static inline uint64_t tk_nbr_n (void *nv, tk_tag_t ntag) {
  return ntag == TK_TAG_I32 ? ((tk_svec_t *) nv)->n : ((tk_ivec_t *) nv)->n;
}
static inline void tk_nbr_setn (void *nv, tk_tag_t ntag, uint64_t n) {
  if (ntag == TK_TAG_I32) ((tk_svec_t *) nv)->n = n; else ((tk_ivec_t *) nv)->n = n;
}
static inline void *tk_nbr_aptr (void *nv, tk_tag_t ntag) {
  return ntag == TK_TAG_I32 ? (void *) ((tk_svec_t *) nv)->a : (void *) ((tk_ivec_t *) nv)->a;
}
static inline size_t tk_nbr_esz (tk_tag_t ntag) {
  return ntag == TK_TAG_I32 ? sizeof(int32_t) : sizeof(int64_t);
}
static inline int tk_nbr_ensure (void *nv, tk_tag_t ntag, uint64_t n) {
  return ntag == TK_TAG_I32 ? tk_svec_ensure((tk_svec_t *) nv, n) : tk_ivec_ensure((tk_ivec_t *) nv, n);
}
static inline void *tk_csr_new_nbr (lua_State *L, tk_tag_t ntag, uint64_t n) {
  return ntag == TK_TAG_I32 ? (void *) tk_svec_create(L, n) : (void *) tk_ivec_create(L, n);
}
static inline int64_t tk_csr_nbr (tk_csr_t *X, uint64_t i) { return tk_nbr_get(X->neighbors, X->ntag, i); }
static inline void tk_csr_setnbr (tk_csr_t *X, uint64_t i, int64_t v) { tk_nbr_set(X->neighbors, X->ntag, i, v); }
static inline uint64_t tk_csr_nbr_n (tk_csr_t *X) { return tk_nbr_n(X->neighbors, X->ntag); }
static inline void tk_csr_nbr_setn (tk_csr_t *X, uint64_t n) { tk_nbr_setn(X->neighbors, X->ntag, n); }
static inline void *tk_csr_nbr_ptr (tk_csr_t *X) { return tk_nbr_aptr(X->neighbors, X->ntag); }
static inline int tk_csr_nbr_ensure (tk_csr_t *X, uint64_t n) { return tk_nbr_ensure(X->neighbors, X->ntag, n); }
static inline void tk_csr_nbr_push (tk_csr_t *X, int64_t v) {
  if (X->ntag == TK_TAG_I32) tk_svec_push((tk_svec_t *) X->neighbors, (int32_t) v);
  else tk_ivec_push((tk_ivec_t *) X->neighbors, v);
}

static inline uint64_t tk_csr_nnz (tk_csr_t *X)
{
  return tk_csr_nbr_n(X);
}

static inline void *tk_csr_val_ptr (tk_csr_t *X)
{
  switch (X->tag) {
    case TK_TAG_I32: return ((tk_svec_t *) X->values)->a;
    case TK_TAG_I64: return ((tk_ivec_t *) X->values)->a;
    case TK_TAG_F32: return ((tk_fvec_t *) X->values)->a;
    case TK_TAG_F64: return ((tk_dvec_t *) X->values)->a;
    case TK_TAG_U8: return ((tk_cvec_t *) X->values)->a;
    default: return NULL;
  }
}

static inline double tk_csr_val1 (tk_csr_t *X, uint64_t i)
{
  switch (X->tag) {
    case TK_TAG_I32: return (double) ((tk_svec_t *) X->values)->a[i];
    case TK_TAG_I64: return (double) ((tk_ivec_t *) X->values)->a[i];
    case TK_TAG_F32: return (double) ((tk_fvec_t *) X->values)->a[i];
    case TK_TAG_F64: return ((tk_dvec_t *) X->values)->a[i];
    case TK_TAG_U8: return (double) (unsigned char) ((tk_cvec_t *) X->values)->a[i];
    default: return 1.0;
  }
}

static inline void tk_csr_setval1 (tk_csr_t *X, uint64_t i, double x)
{
  switch (X->tag) {
    case TK_TAG_I32: ((tk_svec_t *) X->values)->a[i] = (int32_t) x; break;
    case TK_TAG_I64: ((tk_ivec_t *) X->values)->a[i] = (int64_t) x; break;
    case TK_TAG_F32: ((tk_fvec_t *) X->values)->a[i] = (float) x; break;
    case TK_TAG_F64: ((tk_dvec_t *) X->values)->a[i] = x; break;
    case TK_TAG_U8: ((tk_cvec_t *) X->values)->a[i] = (char) (unsigned char) x; break;
    default: break;
  }
}

static inline uint64_t tk_csr_val_len (tk_csr_t *X)
{
  switch (X->tag) {
    case TK_TAG_I32: return ((tk_svec_t *) X->values)->n;
    case TK_TAG_I64: return ((tk_ivec_t *) X->values)->n;
    case TK_TAG_F32: return ((tk_fvec_t *) X->values)->n;
    case TK_TAG_F64: return ((tk_dvec_t *) X->values)->n;
    case TK_TAG_U8: return ((tk_cvec_t *) X->values)->n;
    default: return 0;
  }
}




static inline tk_csr_t *tk_csr_push (lua_State *L, tk_tag_t tag, tk_tag_t ntag, uint64_t n_cols,
  int io, tk_ivec_t *offsets, int in_, void *neighbors, int iv, void *values)
{
  tk_csr_t *X = (tk_csr_t *) lua_newuserdata(L, sizeof(tk_csr_t));
  X->tag = tag;
  X->ntag = ntag;
  X->n_cols = n_cols;
  X->offsets = offsets;
  X->neighbors = neighbors;
  X->values = values;
  luaL_getmetatable(L, TK_CSR_MT);
  lua_setmetatable(L, -2);
  int ix = lua_gettop(L);
  tk_eph_init(L, ix);
  tk_eph_anchor(L, ix, io, offsets);
  tk_eph_anchor(L, ix, in_, neighbors);
  if (iv != 0 && values != NULL)
    tk_eph_anchor(L, ix, iv, values);
  return X;
}

#endif
