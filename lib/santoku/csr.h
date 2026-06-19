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
  uint64_t n_cols;
  tk_ivec_t *offsets;
  tk_ivec_t *neighbors;
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

static inline uint64_t tk_csr_nnz (tk_csr_t *X)
{
  return X->neighbors->n;
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



static inline tk_csr_t *tk_csr_push (lua_State *L, tk_tag_t tag, uint64_t n_cols,
  int io, tk_ivec_t *offsets, int in_, tk_ivec_t *neighbors, int iv, void *values)
{
  tk_csr_t *X = (tk_csr_t *) lua_newuserdata(L, sizeof(tk_csr_t));
  X->tag = tag;
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
