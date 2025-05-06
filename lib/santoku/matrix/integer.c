#define tk_base_t int64_t
#define TK_MT "santoku_matrix_integer"
#define TK_OPEN luaopen_santoku_matrix_integer
#include "gen.h"

static inline int tk_matrix_radd (lua_State *L)
{
  lua_settop(L, 5);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  luaL_checktype(L, 3, LUA_TNUMBER);
  luaL_checktype(L, 4, LUA_TNUMBER);
  lua_Integer rowstart = lua_tointeger(L, 2);
  lua_Integer rowend = lua_tointeger(L, 3);
  lua_Number add = lua_tonumber(L, 4);
  if (rowstart > rowend)
    luaL_error(L, "Error in radd: start row is greater than end row");
  size_t idxstart = tk_matrix_index(L, m0, rowstart, 1);
  size_t idxend = tk_matrix_index(L, m0, rowend, m0->columns);
  for (unsigned int i = idxstart; i <= idxend; i++)
    m0->data[i] += add;
  return 0;
}

static inline int tk_matrix_rmult (lua_State *L)
{
  lua_settop(L, 4);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  luaL_checktype(L, 3, LUA_TNUMBER);
  luaL_checktype(L, 4, LUA_TNUMBER);
  lua_Integer rowstart = lua_tointeger(L, 2);
  lua_Integer rowend = lua_tointeger(L, 3);
  lua_Number scal = lua_tonumber(L, 4);
  if (rowstart > rowend)
    luaL_error(L, "Error in rmult: start row is greater than end row");
  size_t idxstart = tk_matrix_index(L, m0, rowstart, 1);
  size_t idxend = tk_matrix_index(L, m0, rowend, m0->columns);
  tk_base_t *x = m0->data;
  for (unsigned int i = idxstart; i <= idxend; i++)
    x[i] *= scal;
  return 0;
}

static inline int tk_matrix_ramax (lua_State *L)
{
  lua_settop(L, 2);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  lua_Integer row = lua_tointeger(L, 2);
  size_t idx = tk_matrix_index(L, m0, row, 1);
  size_t idxval = idx;
  tk_base_t maxval = abs(m0->data[idx]);
  for (size_t i = idx + 1; i < idx + m0->columns; i++) {
    tk_base_t absval = abs(m0->data[i]);
    if (absval > maxval) {
      maxval = absval;
      idxval = i;
    }
  }
  lua_pushnumber(L, m0->data[idx + idxval]);
  lua_pushinteger(L, idxval);
  return 2;
}

static inline int tk_matrix_sums (lua_State *L)
{
  lua_settop(L, 3);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  tk_matrix_t *m1 = tk_matrix_peek(L, 2);
  luaL_checktype(L, 3, LUA_TNUMBER);
  if (m0->columns != m1->columns)
    luaL_error(L, "Error in sums: destination matrix columns don't match source matrix columns");
  lua_Integer rowdest = lua_tointeger(L, 3);
  size_t idxdest = tk_matrix_index(L, m1, rowdest, 1);
  size_t idxsrc = tk_matrix_index(L, m0, 1, 1);
  memcpy(&m1->data[idxdest], &m0->data[idxsrc], sizeof(tk_base_t) * m1->columns);
  for (size_t i = 2; i <= m0->rows; i ++) {
    idxsrc = tk_matrix_index(L, m0, i, 1);
    for (size_t i = 0; i < m0->columns; i++)
      m1->data[idxdest + i] += m0->data[idxsrc + i];
  }
  return 0;
}

static inline int tk_matrix_mmult (lua_State *L)
{
  lua_settop(L, 5);
  tk_matrix_t *a = tk_matrix_peek(L, 1);
  tk_matrix_t *b = tk_matrix_peek(L, 2);
  tk_matrix_t *c = tk_matrix_peek(L, 3);
  bool transpose_a = lua_toboolean(L, 4);
  bool transpose_b = lua_toboolean(L, 5);
  if (!transpose_a && !transpose_b) {
    if (a->columns != b->rows)
      luaL_error(L, "Error in mmult: columns of A don't match rows of B");
    if (a->rows != c->rows)
      luaL_error(L, "Error in mmult: rows of C don't match rows of A");
    if (b->columns != c->columns)
      luaL_error(L, "Error in mmult: columns of C don't match columns of B");
  } else if (transpose_a && !transpose_b) {
    if (a->rows != b->rows)
      luaL_error(L, "Error in mmult: rows of A don't match rows of B");
    if (a->columns != c->rows)
      luaL_error(L, "Error in mmult: rows of C don't match columns of A");
    if (b->columns != c->columns)
      luaL_error(L, "Error in mmult: columns of C don't match columns of B");
  } else if (!transpose_a && transpose_b) {
    if (a->columns != b->columns)
      luaL_error(L, "Error in mmult: columns of A don't match columns of B");
    if (a->rows != c->rows)
      luaL_error(L, "Error in mmult: rows of C don't match rows of A");
    if (b->rows != c->columns)
      luaL_error(L, "Error in mmult: columns of C don't match rows of B");
  } else if (transpose_a && transpose_b) {
    if (a->rows != b->columns)
      luaL_error(L, "Error in mmult: rows of A don't match columns of B");
    if (a->columns != c->columns)
      luaL_error(L, "Error in mmult: columns of C don't match columns of A");
    if (b->rows != c->rows)
      luaL_error(L, "Error in mmult: rows of C don't match rows of B");
  }
  for (size_t i = 0; i < c->rows; i++) {
    for (size_t j = 0; j < c->columns; j++) {
      tk_base_t sum = 0.0;
      for (size_t k = 0; k < a->columns; k++) {
        tk_base_t a_val = transpose_a ? a->data[k * a->columns + i] : a->data[i * a->columns + k];
        tk_base_t b_val = transpose_b ? b->data[j * b->columns + k] : b->data[k * b->columns + j];
        sum += a_val * b_val;
      }
      c->data[i * c->columns + j] = sum;
    }
  }
  return 0;
}

static inline int tk_matrix_magnitude (lua_State *L)
{
  lua_settop(L, 2);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  lua_Integer row = lua_tointeger(L, 2);
  size_t idx = tk_matrix_index(L, m0, row, 1);
  tk_base_t sum = 0.0;
  for (size_t i = 0; i < m0->columns; i++) {
    tk_base_t val = m0->data[idx + i];
    sum += val * val;
  }
  lua_pushnumber(L, sqrt(sum));
  return 1;
}
