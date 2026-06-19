#include <santoku/iuset.h>
#include <santoku/mtx.h>
#include <santoku/csr.h>
#include <string.h>

static inline uint64_t tk_mtx_checkidx (lua_State *L, tk_mtx_t *M, int ir, int ic)
{
  uint64_t r = tk_lua_checkunsigned(L, ir, "row");
  uint64_t c = tk_lua_checkunsigned(L, ic, "col");
  if (M->tag == TK_TAG_BITS)
    tk_lua_verror(L, 2, "mtx", "get/set not supported for bits layout");
  if (r >= M->n_rows || c >= M->n_cols)
    tk_lua_verror(L, 2, "mtx", "index out of range");
  return r * M->n_cols + c;
}

static inline bool tk_mtx_axis_row (lua_State *L, int i)
{
  const char *a = luaL_checkstring(L, i);
  if (strcmp(a, "row") == 0) return true;
  if (strcmp(a, "col") == 0) return false;
  tk_lua_verror(L, 2, "mtx", "axis must be \"row\" or \"col\"");
  return true;
}






static inline tk_tag_t tk_mtx_tag_of_vec (lua_State *L, int i, void **child)
{
  void *p;
  if ((p = tk_dvec_peekopt(L, i)) != NULL) { *child = p; return TK_TAG_F64; }
  if ((p = tk_fvec_peekopt(L, i)) != NULL) { *child = p; return TK_TAG_F32; }
  if ((p = tk_ivec_peekopt(L, i)) != NULL) { *child = p; return TK_TAG_I64; }
  if ((p = tk_svec_peekopt(L, i)) != NULL) { *child = p; return TK_TAG_I32; }
  if ((p = tk_cvec_peekopt(L, i)) != NULL) { *child = p; return TK_TAG_U8; }
  return TK_TAG_NONE;
}

static int tk_mtx_create_lua (lua_State *L)
{
  lua_settop(L, 1);
  luaL_checktype(L, 1, LUA_TTABLE);
  lua_getfield(L, 1, "n_rows");
  uint64_t rows = tk_lua_checkunsigned(L, -1, "n_rows");
  lua_getfield(L, 1, "n_cols");
  uint64_t cols = tk_lua_checkunsigned(L, -1, "n_cols");
  lua_pop(L, 2);
  lua_getfield(L, 1, "data");
  if (!lua_isnil(L, -1)) {
    void *child = NULL;
    tk_tag_t tag = tk_mtx_tag_of_vec(L, -1, &child);
    if (tag == TK_TAG_NONE)
      return tk_lua_verror(L, 3, "mtx", "data", "expected a santoku vec");
    uint64_t n;
    switch (tag) {
      case TK_TAG_I32: n = ((tk_svec_t *) child)->n; break;
      case TK_TAG_I64: n = ((tk_ivec_t *) child)->n; break;
      case TK_TAG_F32: n = ((tk_fvec_t *) child)->n; break;
      case TK_TAG_F64: n = ((tk_dvec_t *) child)->n; break;
      default: n = ((tk_cvec_t *) child)->n; break;
    }
    lua_getfield(L, 1, "bits");
    bool bits = lua_toboolean(L, -1);
    lua_pop(L, 1);
    if (bits) {
      if (tag != TK_TAG_U8)
        return tk_lua_verror(L, 3, "mtx", "bits", "bits layout requires a cvec data vector");
      tag = TK_TAG_BITS;
      if (n < rows * TK_CVEC_BITS_BYTES(cols))
        return tk_lua_verror(L, 3, "mtx", "data", "vector shorter than packed bit matrix");
    } else if (n < rows * cols)
      return tk_lua_verror(L, 3, "mtx", "data", "vector shorter than n_rows * n_cols");
    tk_mtx_push(L, tag, rows, cols, lua_gettop(L), child);
    return 1;
  }
  lua_pop(L, 1);
  lua_getfield(L, 1, "type");
  tk_tag_t tag = lua_isnil(L, -1) ? TK_TAG_F64 : tk_tag_from_string(luaL_checkstring(L, -1));
  lua_pop(L, 1);
  if (tag == TK_TAG_NONE)
    return tk_lua_verror(L, 3, "mtx", "type", "expected one of u8, i32, i64, f32, f64, bits");
  tk_mtx_t *M = tk_mtx_push_new(L, tag, rows, cols);
  memset(tk_mtx_ptr(M), 0, tk_mtx_rowbytes(M) * rows);
  return 1;
}

static int tk_mtx_shape_lua (lua_State *L)
{
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  lua_pushinteger(L, (lua_Integer) M->n_rows);
  lua_pushinteger(L, (lua_Integer) M->n_cols);
  return 2;
}

static int tk_mtx_type_lua (lua_State *L)
{
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  lua_pushstring(L, tk_tag_name(M->tag));
  return 1;
}

static int tk_mtx_data_lua (lua_State *L)
{
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  tk_eph_get(L, 1, M->v);
  return 1;
}

static int tk_mtx_get_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  uint64_t i = tk_mtx_checkidx(L, M, 2, 3);
  if (M->tag == TK_TAG_F32 || M->tag == TK_TAG_F64)
    lua_pushnumber(L, tk_mtx_get1(M, i));
  else
    lua_pushinteger(L, (lua_Integer) tk_mtx_get1(M, i));
  return 1;
}

static int tk_mtx_set_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  uint64_t i = tk_mtx_checkidx(L, M, 2, 3);
  tk_mtx_set1(M, i, luaL_checknumber(L, 4));
  lua_settop(L, 1);
  return 1;
}

static int tk_mtx_fill_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  if (M->tag == TK_TAG_BITS)
    return tk_lua_verror(L, 2, "mtx", "fill not supported for bits layout");
  double x = luaL_checknumber(L, 2);
  uint64_t n = M->n_rows * M->n_cols;
  for (uint64_t i = 0; i < n; i ++)
    tk_mtx_set1(M, i, x);
  lua_settop(L, 1);
  return 1;
}

static int tk_mtx_eq_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_mtx_t *a = tk_mtx_peek(L, 1, "mtx");
  tk_mtx_t *b = tk_mtx_peek(L, 2, "other");
  double eps = t >= 3 ? luaL_checknumber(L, 3) : -1.0;
  if (a->tag != b->tag || a->n_rows != b->n_rows || a->n_cols != b->n_cols) {
    lua_pushboolean(L, false);
    return 1;
  }
  uint64_t n = a->n_rows * a->n_cols;
  bool r = true;
  if (eps >= 0.0 && a->tag != TK_TAG_BITS) {
    for (uint64_t i = 0; i < n; i ++)
      if (fabs(tk_mtx_get1(a, i) - tk_mtx_get1(b, i)) > eps) { r = false; break; }
  } else {
    r = memcmp(tk_mtx_ptr(a), tk_mtx_ptr(b), tk_mtx_rowbytes(a) * a->n_rows) == 0;
  }
  lua_pushboolean(L, r);
  return 1;
}

#define TK_MTX_NUM_DISPATCH(M, EXPR_S, EXPR_I, EXPR_F, EXPR_D) \
  switch ((M)->tag) { \
    case TK_TAG_I32: { tk_svec_t *v = (tk_svec_t *) (M)->v; (void) v; EXPR_S; break; } \
    case TK_TAG_I64: { tk_ivec_t *v = (tk_ivec_t *) (M)->v; (void) v; EXPR_I; break; } \
    case TK_TAG_F32: { tk_fvec_t *v = (tk_fvec_t *) (M)->v; (void) v; EXPR_F; break; } \
    case TK_TAG_F64: { tk_dvec_t *v = (tk_dvec_t *) (M)->v; (void) v; EXPR_D; break; } \
    default: return tk_lua_verror(L, 2, "mtx", "op not supported for this element type"); \
  }

static int tk_mtx_sums_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  bool row = tk_mtx_axis_row(L, 2);
  TK_MTX_NUM_DISPATCH(M,
    row ? (void) tk_svec_rsums(L, v, M->n_cols) : (void) tk_svec_csums(L, v, M->n_cols),
    row ? (void) tk_ivec_rsums(L, v, M->n_cols) : (void) tk_ivec_csums(L, v, M->n_cols),
    row ? (void) tk_fvec_rsums(L, v, M->n_cols) : (void) tk_fvec_csums(L, v, M->n_cols),
    row ? (void) tk_dvec_rsums(L, v, M->n_cols) : (void) tk_dvec_csums(L, v, M->n_cols))
  return 1;
}

static int tk_mtx_maxs_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  bool row = tk_mtx_axis_row(L, 2);
  TK_MTX_NUM_DISPATCH(M,
    row ? (void) tk_svec_rmaxs(L, v, M->n_cols) : (void) tk_svec_cmaxs(L, v, M->n_cols),
    row ? (void) tk_ivec_rmaxs(L, v, M->n_cols) : (void) tk_ivec_cmaxs(L, v, M->n_cols),
    row ? (void) tk_fvec_rmaxs(L, v, M->n_cols) : (void) tk_fvec_cmaxs(L, v, M->n_cols),
    row ? (void) tk_dvec_rmaxs(L, v, M->n_cols) : (void) tk_dvec_cmaxs(L, v, M->n_cols))
  return 1;
}

static int tk_mtx_mins_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  bool row = tk_mtx_axis_row(L, 2);
  TK_MTX_NUM_DISPATCH(M,
    row ? (void) tk_svec_rmins(L, v, M->n_cols) : (void) tk_svec_cmins(L, v, M->n_cols),
    row ? (void) tk_ivec_rmins(L, v, M->n_cols) : (void) tk_ivec_cmins(L, v, M->n_cols),
    row ? (void) tk_fvec_rmins(L, v, M->n_cols) : (void) tk_fvec_cmins(L, v, M->n_cols),
    row ? (void) tk_dvec_rmins(L, v, M->n_cols) : (void) tk_dvec_cmins(L, v, M->n_cols))
  return 1;
}

static int tk_mtx_mags_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  bool row = tk_mtx_axis_row(L, 2);
  TK_MTX_NUM_DISPATCH(M,
    row ? (void) tk_svec_rmagnitudes(L, v, M->n_cols) : (void) tk_svec_cmagnitudes(L, v, M->n_cols),
    row ? (void) tk_ivec_rmagnitudes(L, v, M->n_cols) : (void) tk_ivec_cmagnitudes(L, v, M->n_cols),
    row ? (void) tk_fvec_rmagnitudes(L, v, M->n_cols) : (void) tk_fvec_cmagnitudes(L, v, M->n_cols),
    row ? (void) tk_dvec_rmagnitudes(L, v, M->n_cols) : (void) tk_dvec_cmagnitudes(L, v, M->n_cols))
  return 1;
}

static int tk_mtx_maxargs_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  bool row = tk_mtx_axis_row(L, 2);
  TK_MTX_NUM_DISPATCH(M,
    row ? (void) tk_svec_rmaxargs(L, v, M->n_cols) : (void) tk_svec_cmaxargs(L, v, M->n_cols),
    row ? (void) tk_ivec_rmaxargs(L, v, M->n_cols) : (void) tk_ivec_cmaxargs(L, v, M->n_cols),
    row ? (void) tk_fvec_rmaxargs(L, v, M->n_cols) : (void) tk_fvec_cmaxargs(L, v, M->n_cols),
    row ? (void) tk_dvec_rmaxargs(L, v, M->n_cols) : (void) tk_dvec_cmaxargs(L, v, M->n_cols))
  return 1;
}

static int tk_mtx_minargs_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  bool row = tk_mtx_axis_row(L, 2);
  TK_MTX_NUM_DISPATCH(M,
    row ? (void) tk_svec_rminargs(L, v, M->n_cols) : (void) tk_svec_cminargs(L, v, M->n_cols),
    row ? (void) tk_ivec_rminargs(L, v, M->n_cols) : (void) tk_ivec_cminargs(L, v, M->n_cols),
    row ? (void) tk_fvec_rminargs(L, v, M->n_cols) : (void) tk_fvec_cminargs(L, v, M->n_cols),
    row ? (void) tk_dvec_rminargs(L, v, M->n_cols) : (void) tk_dvec_cminargs(L, v, M->n_cols))
  return 1;
}

static int tk_mtx_argsort_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  bool row = tk_mtx_axis_row(L, 2);
  const char *dir = luaL_optstring(L, 3, "asc");
  bool asc;
  if (strcmp(dir, "asc") == 0) asc = true;
  else if (strcmp(dir, "desc") == 0) asc = false;
  else return tk_lua_verror(L, 2, "mtx", "dir must be \"asc\" or \"desc\"");
  TK_MTX_NUM_DISPATCH(M,
    row ? (asc ? (void) tk_svec_rasc(L, v, M->n_cols) : (void) tk_svec_rdesc(L, v, M->n_cols))
        : (asc ? (void) tk_svec_casc(L, v, M->n_cols) : (void) tk_svec_cdesc(L, v, M->n_cols)),
    row ? (asc ? (void) tk_ivec_rasc(L, v, M->n_cols) : (void) tk_ivec_rdesc(L, v, M->n_cols))
        : (asc ? (void) tk_ivec_casc(L, v, M->n_cols) : (void) tk_ivec_cdesc(L, v, M->n_cols)),
    row ? (asc ? (void) tk_fvec_rasc(L, v, M->n_cols) : (void) tk_fvec_rdesc(L, v, M->n_cols))
        : (asc ? (void) tk_fvec_casc(L, v, M->n_cols) : (void) tk_fvec_cdesc(L, v, M->n_cols)),
    row ? (asc ? (void) tk_dvec_rasc(L, v, M->n_cols) : (void) tk_dvec_rdesc(L, v, M->n_cols))
        : (asc ? (void) tk_dvec_casc(L, v, M->n_cols) : (void) tk_dvec_cdesc(L, v, M->n_cols)))
  return 1;
}

static int tk_mtx_transpose_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  if (M->tag == TK_TAG_BITS) {
    uint64_t in_rb = TK_CVEC_BITS_BYTES(M->n_cols);
    uint64_t out_rb = TK_CVEC_BITS_BYTES(M->n_rows);
    tk_mtx_t *T = tk_mtx_push_new(L, TK_TAG_BITS, M->n_cols, M->n_rows);
    uint8_t *out = (uint8_t *) tk_mtx_ptr(T);
    const uint8_t *in = (const uint8_t *) tk_mtx_ptr(M);
    memset(out, 0, out_rb * M->n_cols);
    for (uint64_t row = 0; row < M->n_rows; row ++)
      for (uint64_t col = 0; col < M->n_cols; col ++)
        if (in[row * in_rb + col / CHAR_BIT] & (1u << (col % CHAR_BIT)))
          out[col * out_rb + row / CHAR_BIT] |= (1u << (row % CHAR_BIT));
    return 1;
  }
  size_t esz = tk_tag_size(M->tag);
  tk_mtx_t *T = tk_mtx_push_new(L, M->tag, M->n_cols, M->n_rows);
  const char *src = (const char *) tk_mtx_ptr(M);
  char *dst = (char *) tk_mtx_ptr(T);
  for (uint64_t r = 0; r < M->n_rows; r ++)
    for (uint64_t c = 0; c < M->n_cols; c ++)
      memcpy(dst + (c * M->n_rows + r) * esz, src + (r * M->n_cols + c) * esz, esz);
  return 1;
}

static int tk_mtx_rows_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  tk_ivec_t *idx = tk_ivec_peek(L, 2, "ids");
  size_t rb = tk_mtx_rowbytes(M);
  tk_mtx_t *O = lua_isnil(L, 3) ? NULL : tk_mtx_peek(L, 3, "out");
  if (O != NULL) {
    if (O->tag != M->tag)
      return tk_lua_verror(L, 2, "mtx", "rows: out type mismatch");
    tk_mtx_reshape(L, O, idx->n, M->n_cols);
    lua_settop(L, 3);
  } else {
    O = tk_mtx_push_new(L, M->tag, idx->n, M->n_cols);
  }
  const char *src = (const char *) tk_mtx_ptr(M);
  char *dst = (char *) tk_mtx_ptr(O);
  for (uint64_t i = 0; i < idx->n; i ++) {
    int64_t r = idx->a[i];
    if (r < 0 || (uint64_t) r >= M->n_rows)
      return tk_lua_verror(L, 2, "mtx", "row id out of range");
    memcpy(dst + i * rb, src + (uint64_t) r * rb, rb);
  }
  return 1;
}

static int tk_mtx_cols_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  if (M->tag == TK_TAG_BITS)
    return tk_lua_verror(L, 2, "mtx", "cols not supported for bits layout");
  tk_ivec_t *idx = tk_ivec_peek(L, 2, "ids");
  size_t esz = tk_tag_size(M->tag);
  tk_mtx_t *O = tk_mtx_push_new(L, M->tag, M->n_rows, idx->n);
  const char *src = (const char *) tk_mtx_ptr(M);
  char *dst = (char *) tk_mtx_ptr(O);
  for (uint64_t r = 0; r < M->n_rows; r ++)
    for (uint64_t i = 0; i < idx->n; i ++) {
      int64_t c = idx->a[i];
      if (c < 0 || (uint64_t) c >= M->n_cols)
        return tk_lua_verror(L, 2, "mtx", "col id out of range");
      memcpy(dst + (r * idx->n + i) * esz, src + (r * M->n_cols + (uint64_t) c) * esz, esz);
    }
  return 1;
}

static int tk_mtx_row_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  uint64_t r = tk_lua_checkunsigned(L, 2, "row");
  if (r >= M->n_rows)
    return tk_lua_verror(L, 2, "mtx", "row out of range");
  size_t rb = tk_mtx_rowbytes(M);
  uint64_t outn = M->tag == TK_TAG_BITS ? rb : M->n_cols;
  void *child = NULL;
  if (!lua_isnil(L, 3)) {
    tk_tag_t otag = M->tag == TK_TAG_BITS ? TK_TAG_U8 : M->tag;
    tk_tag_t got = tk_mtx_tag_of_vec(L, 3, &child);
    if (got != otag)
      return tk_lua_verror(L, 2, "mtx", "row: out vector type mismatch");
    int rc;
    switch (otag) {
      case TK_TAG_I32: rc = tk_svec_ensure((tk_svec_t *) child, outn); if (rc == 0) ((tk_svec_t *) child)->n = outn; break;
      case TK_TAG_I64: rc = tk_ivec_ensure((tk_ivec_t *) child, outn); if (rc == 0) ((tk_ivec_t *) child)->n = outn; break;
      case TK_TAG_F32: rc = tk_fvec_ensure((tk_fvec_t *) child, outn); if (rc == 0) ((tk_fvec_t *) child)->n = outn; break;
      case TK_TAG_F64: rc = tk_dvec_ensure((tk_dvec_t *) child, outn); if (rc == 0) ((tk_dvec_t *) child)->n = outn; break;
      default: rc = tk_cvec_ensure((tk_cvec_t *) child, outn); if (rc == 0) ((tk_cvec_t *) child)->n = outn; break;
    }
    if (rc != 0)
      return tk_lua_verror(L, 2, "mtx", "allocation failed");
    lua_settop(L, 3);
  } else {
    child = tk_mtx_new_child(L, M->tag == TK_TAG_BITS ? TK_TAG_U8 : M->tag, outn);
  }
  char *dst;
  switch (M->tag) {
    case TK_TAG_I32: dst = (char *) ((tk_svec_t *) child)->a; break;
    case TK_TAG_I64: dst = (char *) ((tk_ivec_t *) child)->a; break;
    case TK_TAG_F32: dst = (char *) ((tk_fvec_t *) child)->a; break;
    case TK_TAG_F64: dst = (char *) ((tk_dvec_t *) child)->a; break;
    default: dst = ((tk_cvec_t *) child)->a; break;
  }
  memcpy(dst, (const char *) tk_mtx_ptr(M) + r * rb, rb);
  return 1;
}


static int tk_mtx_hcat_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  if (M->tag == TK_TAG_BITS)
    return tk_lua_verror(L, 2, "mtx", "hcat not supported for bits layout");
  tk_mtx_t *B = tk_mtx_peek(L, 2, "other");
  if (M->tag != B->tag)
    return tk_lua_verror(L, 2, "mtx", "hcat requires matching element types");
  if (M->n_rows != B->n_rows)
    return tk_lua_verror(L, 2, "mtx", "hcat requires matching n_rows");
  size_t esz = tk_tag_size(M->tag);
  uint64_t c1 = M->n_cols, c2 = B->n_cols, nc = c1 + c2;
  tk_mtx_grow(L, M, M->n_rows * nc);
  char *a = (char *) tk_mtx_ptr(M);
  const char *b = (const char *) tk_mtx_ptr(B);
  for (uint64_t r = M->n_rows; r > 0; r --) {
    uint64_t i = r - 1;
    memmove(a + i * nc * esz, a + i * c1 * esz, c1 * esz);
    memcpy(a + (i * nc + c1) * esz, b + i * c2 * esz, c2 * esz);
  }
  M->n_cols = nc;
  lua_settop(L, 1);
  return 1;
}

static int tk_mtx_persist_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  FILE *fh = tk_lua_fopen(L, luaL_checkstring(L, 2), "w");
  char magic[4] = { 'T', 'K', 'm', 'x' };
  uint8_t version = 1;
  uint8_t tag = (uint8_t) M->tag;
  tk_lua_fwrite(L, magic, 4, 1, fh);
  tk_lua_fwrite(L, (char *) &version, 1, 1, fh);
  tk_lua_fwrite(L, (char *) &tag, 1, 1, fh);
  tk_lua_fwrite(L, (char *) &M->n_rows, sizeof(uint64_t), 1, fh);
  tk_lua_fwrite(L, (char *) &M->n_cols, sizeof(uint64_t), 1, fh);
  tk_lua_fwrite(L, (char *) tk_mtx_ptr(M), tk_mtx_rowbytes(M) * M->n_rows, 1, fh);
  tk_lua_fclose(L, fh);
  return 0;
}

static int tk_mtx_load_lua (lua_State *L)
{
  lua_settop(L, 1);
  FILE *fh = tk_lua_fopen(L, luaL_checkstring(L, 1), "r");
  char magic[4];
  uint8_t version, tag8;
  uint64_t rows, cols;
  tk_lua_fread(L, magic, 4, 1, fh);
  if (memcmp(magic, "TKmx", 4) != 0)
    return tk_lua_verror(L, 2, "mtx", "load: bad magic");
  tk_lua_fread(L, (char *) &version, 1, 1, fh);
  if (version != 1)
    return tk_lua_verror(L, 2, "mtx", "load: unsupported version");
  tk_lua_fread(L, (char *) &tag8, 1, 1, fh);
  tk_lua_fread(L, (char *) &rows, sizeof(uint64_t), 1, fh);
  tk_lua_fread(L, (char *) &cols, sizeof(uint64_t), 1, fh);
  tk_mtx_t *M = tk_mtx_push_new(L, (tk_tag_t) tag8, rows, cols);
  tk_lua_fread(L, (char *) tk_mtx_ptr(M), tk_mtx_rowbytes(M) * M->n_rows, 1, fh);
  tk_lua_fclose(L, fh);
  return 1;
}

static int tk_mtx_center_lua (lua_State *L)
{
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  bool apply = lua_gettop(L) >= 2 && !lua_isnil(L, 2);
  switch (M->tag) {
    case TK_TAG_F64:
      if (apply) tk_dvec_mtx_center(L, (tk_dvec_t *) M->v, M->n_cols, tk_dvec_peek(L, 2, "means"), NULL);
      else { tk_dvec_t *mu = NULL; tk_dvec_mtx_center(L, (tk_dvec_t *) M->v, M->n_cols, NULL, &mu); }
      break;
    case TK_TAG_F32:
      if (apply) tk_fvec_mtx_center(L, (tk_fvec_t *) M->v, M->n_cols, tk_fvec_peek(L, 2, "means"), NULL);
      else { tk_fvec_t *mu = NULL; tk_fvec_mtx_center(L, (tk_fvec_t *) M->v, M->n_cols, NULL, &mu); }
      break;
    default:
      return tk_lua_verror(L, 2, "mtx", "center requires f32 or f64");
  }
  if (apply) {
    lua_settop(L, 1);
    return 1;
  }
  return 1;
}

static int tk_mtx_zscore_lua (lua_State *L)
{
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  bool apply = lua_gettop(L) >= 2 && !lua_isnil(L, 2);
  switch (M->tag) {
    case TK_TAG_F64:
      if (apply) tk_dvec_mtx_zscore(L, (tk_dvec_t *) M->v, M->n_cols, tk_dvec_peek(L, 2, "inv_std"), NULL);
      else { tk_dvec_t *is = NULL; tk_dvec_mtx_zscore(L, (tk_dvec_t *) M->v, M->n_cols, NULL, &is); }
      break;
    case TK_TAG_F32:
      if (apply) tk_fvec_mtx_zscore(L, (tk_fvec_t *) M->v, M->n_cols, tk_fvec_peek(L, 2, "inv_std"), NULL);
      else { tk_fvec_t *is = NULL; tk_fvec_mtx_zscore(L, (tk_fvec_t *) M->v, M->n_cols, NULL, &is); }
      break;
    default:
      return tk_lua_verror(L, 2, "mtx", "zscore requires f32 or f64");
  }
  if (apply) {
    lua_settop(L, 1);
    return 1;
  }
  return 1;
}

static int tk_mtx_standardize_lua (lua_State *L)
{
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  bool apply = lua_gettop(L) >= 3 && !lua_isnil(L, 2);
  switch (M->tag) {
    case TK_TAG_F64:
      if (apply)
        tk_dvec_mtx_standardize(L, (tk_dvec_t *) M->v, M->n_cols,
          tk_dvec_peek(L, 2, "means"), tk_dvec_peek(L, 3, "inv_std"), NULL, NULL);
      else {
        tk_dvec_t *mu = NULL, *is = NULL;
        tk_dvec_mtx_standardize(L, (tk_dvec_t *) M->v, M->n_cols, NULL, NULL, &mu, &is);
      }
      break;
    case TK_TAG_F32:
      if (apply)
        tk_fvec_mtx_standardize(L, (tk_fvec_t *) M->v, M->n_cols,
          tk_fvec_peek(L, 2, "means"), tk_fvec_peek(L, 3, "inv_std"), NULL, NULL);
      else {
        tk_fvec_t *mu = NULL, *is = NULL;
        tk_fvec_mtx_standardize(L, (tk_fvec_t *) M->v, M->n_cols, NULL, NULL, &mu, &is);
      }
      break;
    default:
      return tk_lua_verror(L, 2, "mtx", "standardize requires f32 or f64");
  }
  if (apply) {
    lua_settop(L, 1);
    return 1;
  }
  return 2;
}

static int tk_mtx_median_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  switch (M->tag) {
    case TK_TAG_F64: {
      tk_dvec_t *med = NULL;
      tk_dvec_mtx_median(L, (tk_dvec_t *) M->v, M->n_cols, &med);
      break;
    }
    case TK_TAG_F32: {
      tk_fvec_t *med = NULL;
      tk_fvec_mtx_median(L, (tk_fvec_t *) M->v, M->n_cols, &med);
      break;
    }
    default:
      return tk_lua_verror(L, 2, "mtx", "median requires f32 or f64");
  }
  return 2;
}

static int tk_mtx_sign_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  uint64_t n_trunc = lua_type(L, 2) == LUA_TNUMBER
    ? tk_lua_checkunsigned(L, 2, "n_trunc") : M->n_cols;
  if (n_trunc > M->n_cols)
    n_trunc = M->n_cols;
  switch (M->tag) {
    case TK_TAG_F64: tk_dvec_mtx_sign(L, (tk_dvec_t *) M->v, M->n_cols, n_trunc); break;
    case TK_TAG_F32: tk_fvec_mtx_sign(L, (tk_fvec_t *) M->v, M->n_cols, n_trunc); break;
    default: return tk_lua_verror(L, 2, "mtx", "sign requires f32 or f64");
  }
  return 1;
}

static int tk_mtx_normalize_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  if (!tk_mtx_axis_row(L, 2))
    return tk_lua_verror(L, 2, "mtx", "normalize supports axis \"row\" only");
  switch (M->tag) {
    case TK_TAG_F64: tk_dvec_rnorml2(((tk_dvec_t *) M->v)->a, M->n_rows, M->n_cols); break;
    case TK_TAG_F32: tk_fvec_rnorml2(((tk_fvec_t *) M->v)->a, M->n_rows, M->n_cols); break;
    default: return tk_lua_verror(L, 2, "mtx", "normalize requires f32 or f64");
  }
  lua_settop(L, 1);
  return 1;
}

static int tk_mtx_multiply_lua (lua_State *L)
{
  lua_settop(L, 5);
  tk_mtx_t *A = tk_mtx_peek(L, 1, "mtx");
  tk_mtx_t *B = tk_mtx_peek(L, 2, "other");
  bool ta = lua_toboolean(L, 3);
  bool tb = lua_toboolean(L, 4);
  if (A->tag != B->tag)
    return tk_lua_verror(L, 2, "mtx", "multiply requires matching element types");
  if (A->tag != TK_TAG_F32 && A->tag != TK_TAG_F64)
    return tk_lua_verror(L, 2, "mtx", "multiply requires f32 or f64");
  uint64_t m = ta ? A->n_cols : A->n_rows;
  uint64_t k = ta ? A->n_rows : A->n_cols;
  uint64_t kb = tb ? B->n_cols : B->n_rows;
  uint64_t n = tb ? B->n_rows : B->n_cols;
  if (k != kb)
    return tk_lua_verror(L, 2, "mtx", "multiply: inner dimensions do not match");
  tk_mtx_t *C = lua_isnil(L, 5) ? NULL : tk_mtx_peek(L, 5, "out");
  if (C != NULL) {
    if (C->tag != A->tag)
      return tk_lua_verror(L, 2, "mtx", "multiply: out type mismatch");
    tk_mtx_reshape(L, C, m, n);
    lua_settop(L, 5);
  } else {
    C = tk_mtx_push_new(L, A->tag, m, n);
  }
  if (A->tag == TK_TAG_F64)
    tk_dvec_gemm(ta, tb, m, n, k, 1.0,
      ((tk_dvec_t *) A->v)->a, ((tk_dvec_t *) B->v)->a, 0.0, ((tk_dvec_t *) C->v)->a);
  else
    tk_fvec_gemm(ta, tb, m, n, k, 1.0f,
      ((tk_fvec_t *) A->v)->a, ((tk_fvec_t *) B->v)->a, 0.0f, ((tk_fvec_t *) C->v)->a);
  return 1;
}

static int tk_mtx_multiplyv_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_mtx_t *A = tk_mtx_peek(L, 1, "mtx");
  bool tr = lua_toboolean(L, 3);
  uint64_t ylen = tr ? A->n_cols : A->n_rows;
  uint64_t xlen = tr ? A->n_rows : A->n_cols;
  switch (A->tag) {
    case TK_TAG_F64: {
      tk_dvec_t *x = tk_dvec_peek(L, 2, "x");
      if (x->n < xlen)
        return tk_lua_verror(L, 2, "mtx", "multiplyv: vector too short");
      tk_dvec_t *y = lua_isnil(L, 4) ? NULL : tk_dvec_peek(L, 4, "out");
      if (y != NULL) {
        if (tk_dvec_ensure(y, ylen) != 0)
          return tk_lua_verror(L, 2, "mtx", "allocation failed");
        y->n = ylen;
        lua_settop(L, 4);
      } else {
        y = tk_dvec_create(L, ylen);
      }
      tk_dvec_gemv(tr, A->n_rows, A->n_cols, 1.0, ((tk_dvec_t *) A->v)->a, x->a, 0.0, y->a);
      break;
    }
    case TK_TAG_F32: {
      tk_fvec_t *x = tk_fvec_peek(L, 2, "x");
      if (x->n < xlen)
        return tk_lua_verror(L, 2, "mtx", "multiplyv: vector too short");
      tk_fvec_t *y = lua_isnil(L, 4) ? NULL : tk_fvec_peek(L, 4, "out");
      if (y != NULL) {
        if (tk_fvec_ensure(y, ylen) != 0)
          return tk_lua_verror(L, 2, "mtx", "allocation failed");
        y->n = ylen;
        lua_settop(L, 4);
      } else {
        y = tk_fvec_create(L, ylen);
      }
      tk_fvec_gemv(tr, A->n_rows, A->n_cols, 1.0f, ((tk_fvec_t *) A->v)->a, x->a, 0.0f, y->a);
      break;
    }
    default:
      return tk_lua_verror(L, 2, "mtx", "multiplyv requires f32 or f64");
  }
  return 1;
}


static int tk_mtx_from_pairs_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_ivec_t *is = tk_ivec_peek(L, 1, "i");
  tk_ivec_t *js = tk_ivec_peek(L, 2, "j");
  uint64_t ni = tk_lua_checkunsigned(L, 3, "n_i");
  uint64_t nj = tk_lua_checkunsigned(L, 4, "n_j");
  tk_dvec_t *w = t >= 5 && !lua_isnil(L, 5) ? tk_dvec_peek(L, 5, "weights") : NULL;
  if (is->n != js->n || (w != NULL && w->n != is->n))
    return tk_lua_verror(L, 2, "mtx", "from_pairs: input lengths differ");
  tk_mtx_t *M = tk_mtx_push_new(L, TK_TAG_F64, ni, nj);
  double *out = (double *) tk_mtx_ptr(M);
  memset(out, 0, sizeof(double) * ni * nj);
  for (uint64_t x = 0; x < is->n; x ++) {
    int64_t i = is->a[x], j = js->a[x];
    if (i < 0 || (uint64_t) i >= ni || j < 0 || (uint64_t) j >= nj)
      return tk_lua_verror(L, 2, "mtx", "from_pairs: index out of range");
    out[(uint64_t) i * nj + (uint64_t) j] += w ? w->a[x] : 1.0;
  }
  return 1;
}

static inline uint64_t tk_mtx_total_bits (tk_mtx_t *M)
{

  return M->n_rows * tk_mtx_rowbytes(M) * CHAR_BIT;
}

static inline tk_mtx_t *tk_mtx_check_bits_pair (lua_State *L, tk_mtx_t *M)
{
  tk_mtx_t *B = tk_mtx_peek(L, 2, "other");
  if (M->tag != TK_TAG_BITS || B->tag != TK_TAG_BITS)
    tk_lua_verror(L, 2, "mtx", "bit ops require bits layout");
  if (M->n_rows != B->n_rows || M->n_cols != B->n_cols)
    tk_lua_verror(L, 2, "mtx", "bit ops require matching shapes");
  return B;
}

static int tk_mtx_popcount_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  if (M->tag != TK_TAG_BITS)
    return tk_lua_verror(L, 2, "mtx", "popcount requires bits layout");
  lua_pushinteger(L, (lua_Integer) tk_cvec_bits_popcount((const uint8_t *) tk_mtx_ptr(M), tk_mtx_total_bits(M)));
  return 1;
}

static int tk_mtx_hamming_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  tk_mtx_t *B = tk_mtx_check_bits_pair(L, M);
  lua_pushinteger(L, (lua_Integer) tk_cvec_bits_hamming(
    (const uint8_t *) tk_mtx_ptr(M), (const uint8_t *) tk_mtx_ptr(B), tk_mtx_total_bits(M)));
  return 1;
}

static int tk_mtx_band_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  tk_mtx_t *B = tk_mtx_check_bits_pair(L, M);
  tk_cvec_bits_and((uint8_t *) tk_mtx_ptr(M), (const uint8_t *) tk_mtx_ptr(M),
    (const uint8_t *) tk_mtx_ptr(B), tk_mtx_total_bits(M));
  lua_settop(L, 1);
  return 1;
}

static int tk_mtx_bor_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  tk_mtx_t *B = tk_mtx_check_bits_pair(L, M);
  tk_cvec_bits_or((uint8_t *) tk_mtx_ptr(M), (const uint8_t *) tk_mtx_ptr(M),
    (const uint8_t *) tk_mtx_ptr(B), tk_mtx_total_bits(M));
  lua_settop(L, 1);
  return 1;
}

static int tk_mtx_bxor_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  tk_mtx_t *B = tk_mtx_check_bits_pair(L, M);
  tk_cvec_bits_xor((uint8_t *) tk_mtx_ptr(M), (const uint8_t *) tk_mtx_ptr(M),
    (const uint8_t *) tk_mtx_ptr(B), tk_mtx_total_bits(M));
  lua_settop(L, 1);
  return 1;
}

static int tk_mtx_bandnot_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  tk_mtx_t *B = tk_mtx_check_bits_pair(L, M);
  tk_cvec_bits_andnot((uint8_t *) tk_mtx_ptr(M), (const uint8_t *) tk_mtx_ptr(M),
    (const uint8_t *) tk_mtx_ptr(B), tk_mtx_total_bits(M));
  lua_settop(L, 1);
  return 1;
}

static int tk_mtx_flip_interleave_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  if (M->tag != TK_TAG_BITS)
    return tk_lua_verror(L, 2, "mtx", "flip_interleave requires bits layout");
  tk_cvec_bits_flip_interleave((tk_cvec_t *) M->v, M->n_cols);
  M->n_cols *= 2;
  lua_settop(L, 1);
  return 1;
}


static int tk_mtx_topk_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  tk_mtx_t *Q = tk_mtx_peek(L, 2, "queries");
  uint64_t k = tk_lua_checkunsigned(L, 3, "k");
  if (M->tag != Q->tag)
    return tk_lua_verror(L, 2, "mtx", "topk requires matching element types");
  if (M->n_cols != Q->n_cols)
    return tk_lua_verror(L, 2, "mtx", "topk requires matching n_cols");
  tk_lua_require_mod(L, "santoku.csr");
  switch (M->tag) {
    case TK_TAG_F64:
      tk_dvec_mtx_topk(L, (tk_dvec_t *) Q->v, (tk_dvec_t *) M->v, Q->n_rows, M->n_rows, M->n_cols, k);
      break;
    case TK_TAG_F32:
      tk_fvec_mtx_topk(L, (tk_fvec_t *) Q->v, (tk_fvec_t *) M->v, Q->n_rows, M->n_rows, M->n_cols, k);
      break;
    case TK_TAG_BITS:
      tk_cvec_bits_topk(L, (tk_cvec_t *) Q->v, (tk_cvec_t *) M->v, Q->n_rows, M->n_rows, M->n_cols, k);
      break;
    default:
      return tk_lua_verror(L, 2, "mtx", "topk requires f32, f64, or bits");
  }

  int iv = lua_gettop(L), in_ = iv - 1, io = iv - 2;
  tk_ivec_t *off = tk_ivec_peek(L, io, "offsets");
  tk_ivec_t *ids = tk_ivec_peek(L, in_, "ids");
  void *vals = NULL;
  tk_tag_t vtag = TK_TAG_NONE;
  if ((vals = tk_dvec_peekopt(L, iv)) != NULL) vtag = TK_TAG_F64;
  else if ((vals = tk_fvec_peekopt(L, iv)) != NULL) vtag = TK_TAG_F32;
  else if ((vals = tk_ivec_peekopt(L, iv)) != NULL) vtag = TK_TAG_I64;
  else return tk_lua_verror(L, 2, "mtx", "topk: unexpected score vector type");
  tk_csr_push(L, vtag, TK_TAG_I64, M->n_rows, io, off, in_, ids, iv, vals);
  return 1;
}

static int tk_mtx_to_sparse_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_mtx_t *M = tk_mtx_peek(L, 1, "mtx");
  double eps = lua_isnil(L, 2) ? 0.0 : luaL_checknumber(L, 2);
  if (M->tag == TK_TAG_BITS)
    return tk_lua_verror(L, 2, "mtx", "to_sparse not supported for bits layout");
  uint64_t n = M->n_rows * M->n_cols;
  uint64_t total = 0;
  for (uint64_t i = 0; i < n; i ++)
    if (fabs(tk_mtx_get1(M, i)) > eps)
      total ++;
  tk_lua_require_mod(L, "santoku.csr");
  tk_csr_t *outx = lua_isnil(L, 3) ? NULL : tk_csr_peek(L, 3, "out");
  if (outx != NULL) {
    if (outx->tag != M->tag)
      return tk_lua_verror(L, 2, "mtx", "to_sparse: out type mismatch");
    if (outx->ntag != TK_TAG_I64)
      return tk_lua_verror(L, 2, "mtx", "to_sparse: out must have i64 neighbors");
    if (tk_ivec_ensure(outx->offsets, M->n_rows + 1) != 0
      || tk_csr_nbr_ensure(outx, total) != 0)
      return tk_lua_verror(L, 2, "mtx", "allocation failed");
    outx->offsets->n = M->n_rows + 1;
    tk_csr_nbr_setn(outx, total);
    int rc;
    switch (outx->tag) {
      case TK_TAG_I32: rc = tk_svec_ensure((tk_svec_t *) outx->values, total); if (rc == 0) ((tk_svec_t *) outx->values)->n = total; break;
      case TK_TAG_I64: rc = tk_ivec_ensure((tk_ivec_t *) outx->values, total); if (rc == 0) ((tk_ivec_t *) outx->values)->n = total; break;
      case TK_TAG_F32: rc = tk_fvec_ensure((tk_fvec_t *) outx->values, total); if (rc == 0) ((tk_fvec_t *) outx->values)->n = total; break;
      case TK_TAG_F64: rc = tk_dvec_ensure((tk_dvec_t *) outx->values, total); if (rc == 0) ((tk_dvec_t *) outx->values)->n = total; break;
      default: rc = tk_cvec_ensure((tk_cvec_t *) outx->values, total); if (rc == 0) ((tk_cvec_t *) outx->values)->n = total; break;
    }
    if (rc != 0)
      return tk_lua_verror(L, 2, "mtx", "allocation failed");
    outx->n_cols = M->n_cols;
    tk_ivec_t *off = outx->offsets;
    tk_ivec_t *nbr = (tk_ivec_t *) outx->neighbors;
    void *vals = outx->values;
    off->a[0] = 0;
    uint64_t pos = 0;
    for (uint64_t r = 0; r < M->n_rows; r ++) {
      for (uint64_t c = 0; c < M->n_cols; c ++) {
        double v = tk_mtx_get1(M, r * M->n_cols + c);
        if (fabs(v) > eps) {
          nbr->a[pos] = (int64_t) c;
          switch (M->tag) {
            case TK_TAG_I32: ((tk_svec_t *) vals)->a[pos] = (int32_t) v; break;
            case TK_TAG_I64: ((tk_ivec_t *) vals)->a[pos] = (int64_t) v; break;
            case TK_TAG_F32: ((tk_fvec_t *) vals)->a[pos] = (float) v; break;
            case TK_TAG_F64: ((tk_dvec_t *) vals)->a[pos] = v; break;
            default: ((tk_cvec_t *) vals)->a[pos] = (char) (unsigned char) v; break;
          }
          pos ++;
        }
      }
      off->a[r + 1] = (int64_t) pos;
    }
    lua_settop(L, 3);
    return 1;
  }
  tk_ivec_t *off = tk_ivec_create(L, M->n_rows + 1);
  int io = lua_gettop(L);
  tk_ivec_t *nbr = tk_ivec_create(L, total);
  int in_ = lua_gettop(L);
  nbr->n = total;
  void *vals = tk_mtx_new_child(L, M->tag, total);
  int iv = lua_gettop(L);
  off->a[0] = 0;
  uint64_t pos = 0;
  for (uint64_t r = 0; r < M->n_rows; r ++) {
    for (uint64_t c = 0; c < M->n_cols; c ++) {
      double v = tk_mtx_get1(M, r * M->n_cols + c);
      if (fabs(v) > eps) {
        nbr->a[pos] = (int64_t) c;
        switch (M->tag) {
          case TK_TAG_I32: ((tk_svec_t *) vals)->a[pos] = (int32_t) v; break;
          case TK_TAG_I64: ((tk_ivec_t *) vals)->a[pos] = (int64_t) v; break;
          case TK_TAG_F32: ((tk_fvec_t *) vals)->a[pos] = (float) v; break;
          case TK_TAG_F64: ((tk_dvec_t *) vals)->a[pos] = v; break;
          default: ((tk_cvec_t *) vals)->a[pos] = (char) (unsigned char) v; break;
        }
        pos ++;
      }
    }
    off->a[r + 1] = (int64_t) pos;
  }
  tk_csr_push(L, M->tag, TK_TAG_I64, M->n_cols, io, off, in_, nbr, iv, vals);
  return 1;
}

static luaL_Reg tk_mtx_mt_fns[] = {
  { "shape", tk_mtx_shape_lua },
  { "type", tk_mtx_type_lua },
  { "data", tk_mtx_data_lua },
  { "get", tk_mtx_get_lua },
  { "set", tk_mtx_set_lua },
  { "fill", tk_mtx_fill_lua },
  { "eq", tk_mtx_eq_lua },
  { "sums", tk_mtx_sums_lua },
  { "maxs", tk_mtx_maxs_lua },
  { "mins", tk_mtx_mins_lua },
  { "mags", tk_mtx_mags_lua },
  { "maxargs", tk_mtx_maxargs_lua },
  { "minargs", tk_mtx_minargs_lua },
  { "argsort", tk_mtx_argsort_lua },
  { "transpose", tk_mtx_transpose_lua },
  { "rows", tk_mtx_rows_lua },
  { "cols", tk_mtx_cols_lua },
  { "row", tk_mtx_row_lua },
  { "hcat", tk_mtx_hcat_lua },
  { "persist", tk_mtx_persist_lua },
  { "center", tk_mtx_center_lua },
  { "zscore", tk_mtx_zscore_lua },
  { "standardize", tk_mtx_standardize_lua },
  { "median", tk_mtx_median_lua },
  { "sign", tk_mtx_sign_lua },
  { "normalize", tk_mtx_normalize_lua },
  { "multiply", tk_mtx_multiply_lua },
  { "multiplyv", tk_mtx_multiplyv_lua },
  { "to_sparse", tk_mtx_to_sparse_lua },
  { "popcount", tk_mtx_popcount_lua },
  { "hamming", tk_mtx_hamming_lua },
  { "band", tk_mtx_band_lua },
  { "bor", tk_mtx_bor_lua },
  { "bxor", tk_mtx_bxor_lua },
  { "bandnot", tk_mtx_bandnot_lua },
  { "flip_interleave", tk_mtx_flip_interleave_lua },
  { "topk", tk_mtx_topk_lua },
  { NULL, NULL }
};

static luaL_Reg tk_mtx_fns[] = {
  { "create", tk_mtx_create_lua },
  { "load", tk_mtx_load_lua },
  { "from_pairs", tk_mtx_from_pairs_lua },
  { NULL, NULL }
};

int luaopen_santoku_mtx (lua_State *L)
{
  luaL_newmetatable(L, TK_MTX_MT);
  lua_newtable(L);
  luaL_register(L, NULL, tk_mtx_mt_fns);
  lua_setfield(L, -2, "__index");
  lua_pop(L, 1);
  lua_newtable(L);
  luaL_register(L, NULL, tk_mtx_fns);
  return 1;
}
