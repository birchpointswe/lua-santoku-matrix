#include <santoku/cvec.h>
#include <santoku/ivec.h>
#include <santoku/dvec.h>
#include <santoku/mtx.h>
#include <santoku/csr.h>
#include <stdbool.h>
#include <string.h>

static inline uint64_t tk_cvec_bits_rowbytes (lua_State *L, uint64_t n_bits, const char *fn)
{
  if (n_bits == 0)
    tk_lua_verror(L, 2, fn, "n_bits must be positive");
  return TK_CVEC_BITS_BYTES(n_bits);
}

static inline uint64_t tk_cvec_bits_nrows (lua_State *L, tk_cvec_t *B, uint64_t n_bits, const char *fn)
{
  uint64_t rb = tk_cvec_bits_rowbytes(L, n_bits, fn);
  if (B->n % rb != 0)
    tk_lua_verror(L, 2, fn, "size is not a multiple of the row width in bytes");
  return B->n / rb;
}

static inline void tk_cvec_bits_checklen (lua_State *L, tk_cvec_t *B, uint64_t n_rows, uint64_t n_bits, const char *fn)
{
  uint64_t rb = tk_cvec_bits_rowbytes(L, n_bits, fn);
  if (B->n < n_rows * rb)
    tk_lua_verror(L, 2, fn, "buffer is shorter than n_rows rows of n_bits");
}

static inline bool tk_cvec_bits_checkbit (lua_State *L, int i, const char *fn)
{
  if (lua_type(L, i) == LUA_TBOOLEAN)
    return lua_toboolean(L, i);
  lua_Integer v = luaL_checkinteger(L, i);
  if (v != 0 && v != 1)
    tk_lua_verror(L, 2, fn, "bit value must be 0, 1, true or false");
  return v == 1;
}

static inline uint8_t *tk_cvec_bits_at (lua_State *L, tk_cvec_t *B, uint64_t row, uint64_t col, uint64_t n_bits, const char *fn)
{
  uint64_t rb = tk_cvec_bits_rowbytes(L, n_bits, fn);
  if (col >= n_bits)
    tk_lua_verror(L, 2, fn, "col out of range");
  if ((row + 1) * rb > B->n)
    tk_lua_verror(L, 2, fn, "row out of range");
  return (uint8_t *) B->a + row * rb + col / CHAR_BIT;
}

static int tk_cvec_bits_get_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  uint64_t row = tk_lua_checkunsigned(L, 2, "row");
  uint64_t col = tk_lua_checkunsigned(L, 3, "col");
  uint64_t n_bits = tk_lua_checkunsigned(L, 4, "n_bits");
  uint8_t *p = tk_cvec_bits_at(L, B, row, col, n_bits, "bits_get");
  lua_pushinteger(L, (*p >> (col % CHAR_BIT)) & 1u);
  return 1;
}

static int tk_cvec_bits_set_lua (lua_State *L)
{
  lua_settop(L, 5);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  uint64_t row = tk_lua_checkunsigned(L, 2, "row");
  uint64_t col = tk_lua_checkunsigned(L, 3, "col");
  uint64_t n_bits = tk_lua_checkunsigned(L, 4, "n_bits");
  bool v = tk_cvec_bits_checkbit(L, 5, "bits_set");
  uint8_t *p = tk_cvec_bits_at(L, B, row, col, n_bits, "bits_set");
  if (v)
    *p |= (uint8_t) (1u << (col % CHAR_BIT));
  else
    *p &= (uint8_t) ~(1u << (col % CHAR_BIT));
  lua_settop(L, 1);
  return 1;
}

static int tk_cvec_bits_fill_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  bool v = tk_cvec_bits_checkbit(L, 2, "bits_fill");
  memset(B->a, v ? 0xFF : 0x00, B->n);
  lua_settop(L, 1);
  return 1;
}

static int tk_cvec_bits_rows_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  tk_ivec_t *ids = tk_ivec_peek(L, 2, "ids");
  uint64_t n_bits = tk_lua_checkunsigned(L, 3, "n_bits");
  uint64_t n_rows = tk_cvec_bits_nrows(L, B, n_bits, "bits_rows");
  uint64_t rb = TK_CVEC_BITS_BYTES(n_bits);
  for (uint64_t i = 0; i < ids->n; i ++)
    if (ids->a[i] < 0 || (uint64_t) ids->a[i] >= n_rows)
      return tk_lua_verror(L, 2, "bits_rows", "row id out of range");
  tk_cvec_t *out = tk_cvec_create(L, ids->n * rb);
  for (uint64_t i = 0; i < ids->n; i ++)
    memcpy(out->a + i * rb, B->a + (uint64_t) ids->a[i] * rb, rb);
  return 1;
}

static int tk_cvec_bits_cols_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  tk_ivec_t *ids = tk_ivec_peek(L, 2, "ids");
  uint64_t n_bits = tk_lua_checkunsigned(L, 3, "n_bits");
  uint64_t n_rows = tk_cvec_bits_nrows(L, B, n_bits, "bits_cols");
  if (ids->n == 0)
    return tk_lua_verror(L, 2, "bits_cols", "ids must not be empty");
  for (uint64_t i = 0; i < ids->n; i ++)
    if (ids->a[i] < 0 || (uint64_t) ids->a[i] >= n_bits)
      return tk_lua_verror(L, 2, "bits_cols", "col id out of range");
  uint64_t rb = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t ob = TK_CVEC_BITS_BYTES(ids->n);
  tk_cvec_t *out = tk_cvec_create(L, n_rows * ob);
  memset(out->a, 0, n_rows * ob);
  const uint8_t *src = (const uint8_t *) B->a;
  uint8_t *dst = (uint8_t *) out->a;
  for (uint64_t r = 0; r < n_rows; r ++)
    for (uint64_t i = 0; i < ids->n; i ++) {
      uint64_t c = (uint64_t) ids->a[i];
      if (src[r * rb + c / CHAR_BIT] & (1u << (c % CHAR_BIT)))
        dst[r * ob + i / CHAR_BIT] |= (uint8_t) (1u << (i % CHAR_BIT));
    }
  return 1;
}

static int tk_cvec_bits_hcat_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  tk_cvec_t *C = tk_cvec_peek(L, 2, "other");
  uint64_t nb = tk_lua_checkunsigned(L, 3, "nb");
  uint64_t nc = tk_lua_checkunsigned(L, 4, "nc");
  uint64_t n_rows = tk_cvec_bits_nrows(L, B, nb, "bits_hcat");
  if (tk_cvec_bits_nrows(L, C, nc, "bits_hcat") != n_rows)
    return tk_lua_verror(L, 2, "bits_hcat", "row counts differ");
  uint64_t rb = TK_CVEC_BITS_BYTES(nb);
  uint64_t rc = TK_CVEC_BITS_BYTES(nc);
  uint64_t ob = TK_CVEC_BITS_BYTES(nb + nc);
  tk_cvec_t *out = tk_cvec_create(L, n_rows * ob);
  memset(out->a, 0, n_rows * ob);
  uint8_t *dst = (uint8_t *) out->a;
  for (uint64_t r = 0; r < n_rows; r ++) {
    tk_cvec_bits_copy_at(dst + r * ob, 0, (const uint8_t *) B->a + r * rb, nb);
    tk_cvec_bits_copy_at(dst + r * ob, nb, (const uint8_t *) C->a + r * rc, nc);
  }
  return 1;
}

static int tk_cvec_bits_transpose_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  uint64_t n_rows = tk_lua_checkunsigned(L, 2, "n_rows");
  uint64_t n_bits = tk_lua_checkunsigned(L, 3, "n_bits");
  tk_cvec_bits_checklen(L, B, n_rows, n_bits, "bits_transpose");
  uint64_t in_rb = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t out_rb = TK_CVEC_BITS_BYTES(n_rows);
  tk_cvec_t *T = tk_cvec_create(L, n_bits * out_rb);
  uint8_t *out = (uint8_t *) T->a;
  const uint8_t *in = (const uint8_t *) B->a;
  memset(out, 0, n_bits * out_rb);
  for (uint64_t row = 0; row < n_rows; row ++)
    for (uint64_t col = 0; col < n_bits; col ++)
      if (in[row * in_rb + col / CHAR_BIT] & (1u << (col % CHAR_BIT)))
        out[col * out_rb + row / CHAR_BIT] |= (uint8_t) (1u << (row % CHAR_BIT));
  return 1;
}

static int tk_cvec_bits_popcount_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  lua_pushinteger(L, (lua_Integer) tk_cvec_bits_popcount((const uint8_t *) B->a, B->n * CHAR_BIT));
  return 1;
}

static int tk_cvec_bits_hamming_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  tk_cvec_t *Q = tk_cvec_peek(L, 2, "other");
  uint64_t n_bits = tk_lua_checkunsigned(L, 3, "n_bits");
  uint64_t n_rows = tk_cvec_bits_nrows(L, B, n_bits, "bits_hamming");
  uint64_t n_q = tk_cvec_bits_nrows(L, Q, n_bits, "bits_hamming");
  if (n_q != n_rows && n_q != 1)
    return tk_lua_verror(L, 2, "bits_hamming", "other must have one row or the same row count");
  uint64_t rb = TK_CVEC_BITS_BYTES(n_bits);
  uint64_t qstep = n_q == 1 ? 0 : rb;
  tk_dvec_t *out = tk_dvec_create(L, n_rows);
  const uint8_t *a = (const uint8_t *) B->a;
  const uint8_t *b = (const uint8_t *) Q->a;
  for (uint64_t i = 0; i < n_rows; i ++)
    out->a[i] = (double) tk_cvec_bits_hamming_serial(a + i * rb, b + i * qstep, n_bits);
  return 1;
}

static int tk_cvec_bits_topk_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  tk_cvec_t *Q = tk_cvec_peek(L, 2, "queries");
  uint64_t n_bits = tk_lua_checkunsigned(L, 3, "n_bits");
  uint64_t k = tk_lua_checkunsigned(L, 4, "k");
  uint64_t n_rows = tk_cvec_bits_nrows(L, B, n_bits, "bits_topk");
  uint64_t n_q = tk_cvec_bits_nrows(L, Q, n_bits, "bits_topk");
  tk_lua_require_mod(L, "santoku.csr");
  tk_cvec_bits_topk(L, Q, B, n_q, n_rows, n_bits, k);
  int iv = lua_gettop(L), in_ = iv - 1, io = iv - 2;
  tk_ivec_t *off = tk_ivec_peek(L, io, "offsets");
  tk_ivec_t *ids = tk_ivec_peek(L, in_, "ids");
  tk_ivec_t *dists = tk_ivec_peek(L, iv, "distances");
  tk_csr_push(L, TK_TAG_I64, TK_TAG_I64, n_rows, io, off, in_, ids, iv, dists);
  return 1;
}

static inline tk_cvec_t *tk_cvec_bits_peek_pair (lua_State *L, tk_cvec_t *B, const char *fn)
{
  tk_cvec_t *C = tk_cvec_peek(L, 2, "other");
  if (B->n != C->n)
    tk_lua_verror(L, 2, fn, "byte lengths differ");
  return C;
}

static int tk_cvec_bits_and_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  tk_cvec_t *C = tk_cvec_bits_peek_pair(L, B, "bits_and");
  tk_cvec_bits_and((uint8_t *) B->a, (const uint8_t *) B->a, (const uint8_t *) C->a, B->n * CHAR_BIT);
  lua_settop(L, 1);
  return 1;
}

static int tk_cvec_bits_or_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  tk_cvec_t *C = tk_cvec_bits_peek_pair(L, B, "bits_or");
  tk_cvec_bits_or((uint8_t *) B->a, (const uint8_t *) B->a, (const uint8_t *) C->a, B->n * CHAR_BIT);
  lua_settop(L, 1);
  return 1;
}

static int tk_cvec_bits_xor_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  tk_cvec_t *C = tk_cvec_bits_peek_pair(L, B, "bits_xor");
  tk_cvec_bits_xor((uint8_t *) B->a, (const uint8_t *) B->a, (const uint8_t *) C->a, B->n * CHAR_BIT);
  lua_settop(L, 1);
  return 1;
}

static int tk_cvec_bits_andnot_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  tk_cvec_t *C = tk_cvec_bits_peek_pair(L, B, "bits_andnot");
  tk_cvec_bits_andnot((uint8_t *) B->a, (const uint8_t *) B->a, (const uint8_t *) C->a, B->n * CHAR_BIT);
  lua_settop(L, 1);
  return 1;
}

static int tk_cvec_bits_to_dense_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_cvec_t *B = tk_cvec_peek(L, 1, "cvec");
  uint64_t n_rows = tk_lua_checkunsigned(L, 2, "n_rows");
  uint64_t n_bits = tk_lua_checkunsigned(L, 3, "n_bits");
  tk_cvec_bits_checklen(L, B, n_rows, n_bits, "bits_to_dense");
  uint64_t rb = TK_CVEC_BITS_BYTES(n_bits);
  tk_lua_require_mod(L, "santoku.mtx");
  tk_mtx_t *M = tk_mtx_push_new(L, TK_TAG_F32, n_rows, n_bits);
  float *out = (float *) tk_mtx_ptr(M);
  const uint8_t *in = (const uint8_t *) B->a;
  for (uint64_t r = 0; r < n_rows; r ++)
    for (uint64_t c = 0; c < n_bits; c ++)
      out[r * n_bits + c] = (in[r * rb + c / CHAR_BIT] >> (c % CHAR_BIT)) & 1u ? 1.0f : 0.0f;
  return 1;
}

static luaL_Reg tk_cvec_lua_mt_ext2_fns[] =
{
  { "bits_get", tk_cvec_bits_get_lua },
  { "bits_set", tk_cvec_bits_set_lua },
  { "bits_fill", tk_cvec_bits_fill_lua },
  { "bits_rows", tk_cvec_bits_rows_lua },
  { "bits_cols", tk_cvec_bits_cols_lua },
  { "bits_hcat", tk_cvec_bits_hcat_lua },
  { "bits_transpose", tk_cvec_bits_transpose_lua },
  { "bits_popcount", tk_cvec_bits_popcount_lua },
  { "bits_hamming", tk_cvec_bits_hamming_lua },
  { "bits_topk", tk_cvec_bits_topk_lua },
  { "bits_and", tk_cvec_bits_and_lua },
  { "bits_or", tk_cvec_bits_or_lua },
  { "bits_xor", tk_cvec_bits_xor_lua },
  { "bits_andnot", tk_cvec_bits_andnot_lua },
  { "bits_to_dense", tk_cvec_bits_to_dense_lua },
  { NULL, NULL }
};

int luaopen_santoku_cvec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_cvec_lua_fns);
  tk_cvec_create(L, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_cvec_lua_mt_fns);
  luaL_register(L, NULL, tk_cvec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_cvec_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
