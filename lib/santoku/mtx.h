#ifndef TK_MTX_H
#define TK_MTX_H

#include <santoku/shaped.h>
#include <santoku/ivec.h>
#include <santoku/svec.h>
#include <santoku/dvec.h>
#include <santoku/fvec.h>
#include <santoku/cvec.h>
#include <limits.h>

#define TK_MTX_MT "tk_mtx_t"

typedef struct {
  tk_tag_t tag;
  uint64_t n_rows;
  uint64_t n_cols;
  void *v;
} tk_mtx_t;

static inline tk_mtx_t *tk_mtx_peek (lua_State *L, int i, const char *name)
{
  tk_mtx_t *M = (tk_mtx_t *) luaL_checkudata(L, i, TK_MTX_MT);
  if (M == NULL)
    tk_lua_verror(L, 2, name, "expected a mtx");
  return M;
}

static inline tk_mtx_t *tk_mtx_peekopt (lua_State *L, int i)
{
  if (lua_type(L, i) != LUA_TUSERDATA)
    return NULL;
  void *p = lua_touserdata(L, i);
  if (!lua_getmetatable(L, i))
    return NULL;
  luaL_getmetatable(L, TK_MTX_MT);
  bool ok = lua_rawequal(L, -1, -2);
  lua_pop(L, 2);
  return ok ? (tk_mtx_t *) p : NULL;
}

static inline void *tk_mtx_ptr (tk_mtx_t *M)
{
  switch (M->tag) {
    case TK_TAG_I32: return ((tk_svec_t *) M->v)->a;
    case TK_TAG_I64: return ((tk_ivec_t *) M->v)->a;
    case TK_TAG_F32: return ((tk_fvec_t *) M->v)->a;
    case TK_TAG_F64: return ((tk_dvec_t *) M->v)->a;
    default: return ((tk_cvec_t *) M->v)->a;
  }
}

#ifndef TK_CVEC_BITS_BYTES
#define TK_CVEC_BITS_BYTES(n) (((n) + CHAR_BIT - 1) / CHAR_BIT)
#endif

static inline uint64_t tk_mtx_rowbytes (tk_mtx_t *M)
{
  return M->tag == TK_TAG_BITS
    ? TK_CVEC_BITS_BYTES(M->n_cols)
    : tk_tag_size(M->tag) * M->n_cols;
}

static inline uint64_t tk_mtx_len (tk_mtx_t *M)
{
  switch (M->tag) {
    case TK_TAG_I32: return ((tk_svec_t *) M->v)->n;
    case TK_TAG_I64: return ((tk_ivec_t *) M->v)->n;
    case TK_TAG_F32: return ((tk_fvec_t *) M->v)->n;
    case TK_TAG_F64: return ((tk_dvec_t *) M->v)->n;
    default: return ((tk_cvec_t *) M->v)->n;
  }
}

static inline double tk_mtx_get1 (tk_mtx_t *M, uint64_t i)
{
  switch (M->tag) {
    case TK_TAG_I32: return (double) ((tk_svec_t *) M->v)->a[i];
    case TK_TAG_I64: return (double) ((tk_ivec_t *) M->v)->a[i];
    case TK_TAG_F32: return (double) ((tk_fvec_t *) M->v)->a[i];
    case TK_TAG_F64: return ((tk_dvec_t *) M->v)->a[i];
    default: return (double) (unsigned char) ((tk_cvec_t *) M->v)->a[i];
  }
}

static inline void tk_mtx_set1 (tk_mtx_t *M, uint64_t i, double x)
{
  switch (M->tag) {
    case TK_TAG_I32: ((tk_svec_t *) M->v)->a[i] = (int32_t) x; break;
    case TK_TAG_I64: ((tk_ivec_t *) M->v)->a[i] = (int64_t) x; break;
    case TK_TAG_F32: ((tk_fvec_t *) M->v)->a[i] = (float) x; break;
    case TK_TAG_F64: ((tk_dvec_t *) M->v)->a[i] = x; break;
    default: ((tk_cvec_t *) M->v)->a[i] = (char) (unsigned char) x; break;
  }
}


static inline void *tk_mtx_new_child (lua_State *L, tk_tag_t tag, uint64_t n)
{
  switch (tag) {
    case TK_TAG_I32: return tk_svec_create(L, n);
    case TK_TAG_I64: return tk_ivec_create(L, n);
    case TK_TAG_F32: return tk_fvec_create(L, n);
    case TK_TAG_F64: return tk_dvec_create(L, n);
    case TK_TAG_U8: return tk_cvec_create(L, n);
    case TK_TAG_BITS: return tk_cvec_create(L, n);
    default:
      tk_lua_verror(L, 2, "mtx", "unsupported element type");
      return NULL;
  }
}



static inline tk_mtx_t *tk_mtx_push (lua_State *L, tk_tag_t tag, uint64_t rows, uint64_t cols, int ic, void *child)
{
  tk_mtx_t *M = (tk_mtx_t *) lua_newuserdata(L, sizeof(tk_mtx_t));
  M->tag = tag;
  M->n_rows = rows;
  M->n_cols = cols;
  M->v = child;
  luaL_getmetatable(L, TK_MTX_MT);
  lua_setmetatable(L, -2);
  tk_eph_init(L, lua_gettop(L));
  tk_eph_anchor(L, lua_gettop(L), ic, child);
  return M;
}

static inline void tk_mtx_grow (lua_State *L, tk_mtx_t *M, uint64_t total)
{
  int rc;
  switch (M->tag) {
    case TK_TAG_I32: rc = tk_svec_ensure((tk_svec_t *) M->v, total); if (rc == 0) ((tk_svec_t *) M->v)->n = total; break;
    case TK_TAG_I64: rc = tk_ivec_ensure((tk_ivec_t *) M->v, total); if (rc == 0) ((tk_ivec_t *) M->v)->n = total; break;
    case TK_TAG_F32: rc = tk_fvec_ensure((tk_fvec_t *) M->v, total); if (rc == 0) ((tk_fvec_t *) M->v)->n = total; break;
    case TK_TAG_F64: rc = tk_dvec_ensure((tk_dvec_t *) M->v, total); if (rc == 0) ((tk_dvec_t *) M->v)->n = total; break;
    default: rc = tk_cvec_ensure((tk_cvec_t *) M->v, total); if (rc == 0) ((tk_cvec_t *) M->v)->n = total; break;
  }
  if (rc != 0)
    tk_lua_verror(L, 2, "mtx", "allocation failed");
}


static inline void tk_mtx_reshape (lua_State *L, tk_mtx_t *M, uint64_t rows, uint64_t cols)
{
  M->n_rows = rows;
  M->n_cols = cols;
  tk_mtx_grow(L, M, M->tag == TK_TAG_BITS ? rows * TK_CVEC_BITS_BYTES(cols) : rows * cols);
}


static inline tk_mtx_t *tk_mtx_push_new (lua_State *L, tk_tag_t tag, uint64_t rows, uint64_t cols)
{
  uint64_t n = tag == TK_TAG_BITS ? rows * TK_CVEC_BITS_BYTES(cols) : rows * cols;
  void *child = tk_mtx_new_child(L, tag, n);
  int ic = lua_gettop(L);
  tk_mtx_t *M = tk_mtx_push(L, tag, rows, cols, ic, child);
  lua_remove(L, ic);
  return M;
}

#endif
