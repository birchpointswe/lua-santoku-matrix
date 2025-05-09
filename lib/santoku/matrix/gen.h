#include "lua.h"
#include "lauxlib.h"

#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "ksort.h"

KSORT_INIT_GENERIC(tk_base_t)

typedef struct {
  size_t rows;
  size_t columns;
  size_t values;
  tk_base_t *data;
} tk_matrix_t;

static inline void tk_lua_callmod (lua_State *L, int nargs, int nret, const char *smod, const char *sfn)
{
  lua_getglobal(L, "require");
  lua_pushstring(L, smod);
  lua_call(L, 1, 1);
  lua_pushstring(L, sfn);
  lua_gettable(L, -2);
  lua_remove(L, -2);
  lua_insert(L, - nargs - 1);
  lua_call(L, nargs, nret);
}


static inline void *tk_lua_checkuserdata (lua_State *L, int i, char *mt)
{
  if (mt == NULL && (lua_islightuserdata(L, i) || lua_isuserdata(L, i)))
    return lua_touserdata(L, i);
  void *p = luaL_checkudata(L, -1, mt);
  lua_pop(L, 1);
  return p;
}

static inline uint64_t tk_lua_checkunsigned (lua_State *L, int i)
{
  lua_Integer l = luaL_checkinteger(L, i);
  if (l < 0)
    luaL_error(L, "value can't be negative");
  return (uint64_t) l;
}

static inline unsigned int tk_lua_optunsigned (lua_State *L, int i, unsigned int def)
{
  if (lua_type(L, i) < 1)
    return def;
  return tk_lua_checkunsigned(L, i);
}

static inline tk_matrix_t *_tk_matrix_create (lua_State *L, size_t rows, size_t columns, tk_base_t *data, tk_matrix_t *m0)
{
  size_t values = rows * columns;
  if (m0 == NULL) {
    m0 = malloc(sizeof(tk_matrix_t));
    if (m0 == NULL)
      luaL_error(L, "Error in malloc during matrix create");
    memset(m0, 0, sizeof(tk_matrix_t));
    tk_matrix_t **m0p = (tk_matrix_t **) lua_newuserdata(L, sizeof(tk_matrix_t *));
    luaL_getmetatable(L, TK_MT);
    lua_setmetatable(L, -2);
    *m0p = m0;
  } else if (m0->data && m0->data != data) {
    free(m0->data);
  }
  m0->rows = rows;
  m0->columns = columns;
  m0->values = values;
  m0->data = data != NULL ? data : malloc(sizeof(tk_base_t) * values);
  if (m0->data == NULL)
    luaL_error(L, "Error in malloc during matrix create");
  return m0;
}

static inline tk_matrix_t **tk_matrix_peekp (lua_State *L, int i)
{
  return (tk_matrix_t **) luaL_checkudata(L, i, TK_MT);
}

static inline tk_matrix_t *tk_matrix_peek (lua_State *L, int i)
{
  tk_matrix_t **m0p = tk_matrix_peekp(L, i);
  return *m0p;
}

static inline size_t tk_matrix_index (lua_State *L, tk_matrix_t *m0, size_t row, size_t column)
{
  if (row > m0->rows || row < 1)
    luaL_error(L, "Matrix row index out of bounds");
  if (column > m0->columns || column < 1)
    luaL_error(L, "Matrix column index out of bounds");
  return (row - 1) * m0->columns + column - 1;
}

static inline int tk_matrix_copy (lua_State *L)
{
  lua_settop(L, 5);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  tk_matrix_t *m1 = tk_matrix_peek(L, 2);
  if (m0->columns != m1->columns)
    luaL_error(L, "Error in copy: can't copy between matrices with different column lengths");
  luaL_checktype(L, 3, LUA_TNUMBER);
  luaL_checktype(L, 4, LUA_TNUMBER);
  luaL_checktype(L, 5, LUA_TNUMBER);
  size_t rowstart = lua_tointeger(L, 3);
  size_t rowend = lua_tointeger(L, 4);
  size_t rowdest = lua_tointeger(L, 5);
  if (rowstart > rowend)
    luaL_error(L, "Error in copy: start row is greater than end row");
  if (rowdest + rowend - rowstart > m0->rows)
    luaL_error(L, "Error in copy: copying more rows than space available");
  size_t idxstart = tk_matrix_index(L, m1, rowstart, 1);
  size_t idxend = tk_matrix_index(L, m1, rowend, m1->columns);
  size_t idxdest = tk_matrix_index(L, m0, rowdest, 1);
  memcpy(&m0->data[idxdest], &m1->data[idxstart], sizeof(tk_base_t) * (idxend - idxstart + 1));
  return 0;
}

static inline int tk_matrix_shrink (lua_State *L)
{
  lua_settop(L, 1);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  if (m0->values >= m0->rows * m0->columns) {
    m0->values = m0->rows * m0->columns;
    m0->data = realloc(m0->data, sizeof(tk_base_t) * m0->values);
    if (m0->data == NULL)
      luaL_error(L, "Error in realloc during matrix shrink");
  }
  return 0;
}

struct tk_matrix_rorder_item {
  bool asc;
  size_t column;
  tk_base_t value;
};

static int tk_matrix_rorder_cmp (const void *ap, const void *bp)
{
  struct tk_matrix_rorder_item *a = (struct tk_matrix_rorder_item *) ap;
  struct tk_matrix_rorder_item *b = (struct tk_matrix_rorder_item *) bp;
  tk_base_t r = a->asc
    ? a->value - b->value
    : b->value - a->value;
  return r == 0 ? 0 : r < 0 ? -1 : 1;
}

static inline void tk_matrix_rorder_push_row (
  lua_State *L,
  tk_matrix_t *m0,
  size_t row,
  unsigned int min,
  unsigned int max,
  tk_base_t cutoff,
  bool lt,
  bool asc
) {
  struct tk_matrix_rorder_item items[m0->columns];
  for (size_t column = 1; column <= m0->columns; column ++) {
    items[column - 1].asc = asc;
    items[column - 1].column = column;
    items[column - 1].value = m0->data[tk_matrix_index(L, m0, row, column)];
  }
  qsort(items,
        m0->columns,
        sizeof(struct tk_matrix_rorder_item),
        tk_matrix_rorder_cmp);
  lua_newtable(L);
  for (size_t i = 0; i < m0->columns; i ++) {
    unsigned int n = i + 1;
    if (n > max || (n > min && (lt
        ? items[i].value < cutoff
        : items[i].value > cutoff)))
      break;
    lua_pushinteger(L, n);
    lua_pushinteger(L, items[i].column);
    lua_settable(L, -3);
  }
}

static inline int tk_matrix_rorder (lua_State *L)
{
  lua_settop(L, 5);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  lua_Integer min = tk_lua_optunsigned(L, 2, 0);
  lua_Integer max = tk_lua_optunsigned(L, 3, m0->columns);
  lua_Number cutoff = luaL_optnumber(L, 4, -DBL_MAX);
  bool lt = !lua_toboolean(L, 5);
  bool asc = lua_toboolean(L, 6);
  if (max < min)
    luaL_error(L, "Max can't be less than min");
  lua_newtable(L);
  for (size_t row = 1; row <= m0->rows; row ++) {
    lua_pushinteger(L, row);
    tk_matrix_rorder_push_row(L, m0, row, min, max, cutoff, lt, asc);
    lua_settable(L, -3);
  }
  return 1;
}

static inline int tk_matrix_tabulate (lua_State *L)
{
  lua_settop(L, 1);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  lua_newtable(L);
  for (size_t row = 1; row <= m0->rows; row ++) {
    lua_pushinteger(L, row);
    lua_newtable(L);
    for (size_t column = 1; column <= m0->columns; column ++) {
      size_t idx = tk_matrix_index(L, m0, row, column);
      lua_pushinteger(L, column);
      lua_pushnumber(L, m0->data[idx]);
      lua_settable(L, -3);
    }
    lua_settable(L, -3);
  }
  return 1;
}

static inline int tk_matrix_reshape (lua_State *L)
{
  lua_settop(L, 3);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  luaL_checktype(L, 3, LUA_TNUMBER);
  lua_Integer rows = lua_tointeger(L, 2);
  lua_Integer columns = lua_tointeger(L, 3);
  if (rows < 0)
    luaL_error(L, "Error in reshape: rows less than 0");
  if (columns < 0)
    luaL_error(L, "Error in reshape: columns less than 0");
  m0->rows = rows;
  m0->columns = columns;
  if (rows * columns > m0->values) {
    m0->values = rows * columns;
    m0->data = realloc(m0->data, sizeof(tk_base_t) * m0->values);
    if (m0->data == NULL)
      luaL_error(L, "Error in realloc during matrix reshape");
  }
  return 0;
}

static inline int tk_matrix_create (lua_State *L)
{
  lua_settop(L, 2);
  luaL_checktype(L, 1, LUA_TNUMBER);
  luaL_checktype(L, 2, LUA_TNUMBER);
  size_t rows = lua_tointeger(L, 1);
  size_t columns = lua_tointeger(L, 2);
  _tk_matrix_create(L, rows, columns, NULL, NULL);
  return 1;
}

static inline int tk_matrix_gc (lua_State *L)
{
  lua_settop(L, 1);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  free(m0->data);
  free(m0);
  return 0;
}

static inline int tk_matrix_get (lua_State *L)
{
  lua_settop(L, 3);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  luaL_checktype(L, 3, LUA_TNUMBER);
  lua_Integer row = lua_tointeger(L, 2);
  lua_Integer column = lua_tointeger(L, 3);
  size_t idx = tk_matrix_index(L, m0, row, column);
  lua_pushnumber(L, m0->data[idx]);
  return 1;
}

static inline int tk_matrix_set (lua_State *L)
{
  lua_settop(L, 4);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  luaL_checktype(L, 3, LUA_TNUMBER);
  luaL_checktype(L, 4, LUA_TNUMBER);
  lua_Integer row = lua_tointeger(L, 2);
  lua_Integer column = lua_tointeger(L, 3);
  lua_Number value = lua_tonumber(L, 4);
  size_t idx = tk_matrix_index(L, m0, row, column);
  m0->data[idx] = value;
  return 0;
}

static inline int tk_matrix_exp (lua_State *L)
{
  lua_settop(L, 4);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  luaL_checktype(L, 3, LUA_TNUMBER);
  luaL_checktype(L, 4, LUA_TNUMBER);
  lua_Integer rowstart = lua_tointeger(L, 2);
  lua_Integer rowend = lua_tointeger(L, 3);
  lua_Number exp = lua_tonumber(L, 4);
  if (rowstart > rowend)
    luaL_error(L, "Error in rmult: start row is greater than end row");
  size_t idxstart = tk_matrix_index(L, m0, rowstart, 1);
  size_t idxend = tk_matrix_index(L, m0, rowend, m0->columns);
  for (size_t i = idxstart; i <= idxend; i ++)
    m0->data[i] = pow(m0->data[i], exp);
  return 0;
}

static inline int tk_matrix_rmin (lua_State *L)
{
  lua_settop(L, 2);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  lua_Integer row = lua_tointeger(L, 2);
  size_t idx = tk_matrix_index(L, m0, row, 1);
  size_t mincol = 1;
  tk_base_t minval = m0->data[idx];
  for (size_t i = 2; i <= m0->columns; i ++) {
    if (m0->data[idx + i - 1] < minval) {
      mincol = i;
      minval = m0->data[idx + i - 1];
    }
  }
  lua_pushnumber(L, minval);
  lua_pushinteger(L, mincol);
  return 2;
}

static inline int tk_matrix_rmax (lua_State *L)
{
  lua_settop(L, 2);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  lua_Integer row = lua_tointeger(L, 2);
  size_t idx = tk_matrix_index(L, m0, row, 1);
  size_t maxcol = 1;
  tk_base_t maxval = m0->data[idx];
  for (size_t i = 2; i <= m0->columns; i ++) {
    if (m0->data[idx + i - 1] > maxval) {
      maxcol = i;
      maxval = m0->data[idx + i - 1];
    }
  }
  lua_pushnumber(L, maxval);
  lua_pushinteger(L, maxcol);
  return 2;
}

static inline int tk_matrix_sort (lua_State *L)
{
  lua_settop(L, 2);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  bool unique = lua_toboolean(L, 2);
  if (m0->values < 1)
    return 0;
  size_t values = m0->values;
  tk_base_t *data = m0->data;
  tk_sort(values, data);
  if (!unique)
    return 0;
  size_t write = 1;
  for (size_t read = 1; read < values; read ++)
    if (data[read] != data[write - 1])
      data[write ++] = data[read];
  m0->rows = 1;
  m0->columns = write;
  m0->values = write;
  return 0;
}

static inline int tk_matrix_sum (lua_State *L)
{
  lua_settop(L, 3);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  luaL_checktype(L, 3, LUA_TNUMBER);
  size_t rowstart = lua_tointeger(L, 2);
  size_t rowend = lua_tointeger(L, 3);
  if (rowstart > rowend)
    luaL_error(L, "Error in sum: start row is greater than end row");
  size_t idxstart = tk_matrix_index(L, m0, rowstart, 1);
  size_t idxend = tk_matrix_index(L, m0, rowend, m0->columns);
  tk_base_t sum = 0;
  for (size_t i = idxstart; i <= idxend; i ++)
    sum += m0->data[i];
  lua_pushnumber(L, sum);
  return 1;
}

static inline int tk_matrix_rows (lua_State *L)
{
  lua_settop(L, 1);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  lua_pushinteger(L, m0->rows);
  return 1;
}

static inline int tk_matrix_columns (lua_State *L)
{
  lua_settop(L, 1);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  lua_pushinteger(L, m0->columns);
  return 1;
}

static inline int tk_matrix_shape (lua_State *L)
{
  lua_settop(L, 1);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  lua_pushinteger(L, m0->rows);
  lua_pushinteger(L, m0->columns);
  return 2;
}

static inline int tk_matrix_extend_raw (lua_State *L)
{
  lua_settop(L, 2);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TSTRING);
  size_t size;
  const char *data = lua_tolstring(L, 2, &size);
  if (size % sizeof(tk_base_t) != 0)
    luaL_error(L, "Length of raw string is not a multiple of sizeof(tk_base_t)");
  size_t values = size / sizeof(tk_base_t);
  if (values % m0->columns != 0)
    luaL_error(L, "Length of raw string is not a multiple of matrix columns");
  size_t extend_rows = values / m0->columns;
  size_t extend_rowstart = m0->rows + 1;
  size_t newrows = m0->rows + extend_rows;
  lua_pop(L, 1);
  lua_pushinteger(L, newrows);
  lua_pushinteger(L, m0->columns);
  tk_matrix_reshape(L);
  m0 = tk_matrix_peek(L, 1);
  size_t idxextend = tk_matrix_index(L, m0, extend_rowstart, 1);
  memcpy(&m0->data[idxextend], data, size);
  return 0;
}

static inline int tk_matrix_raw (lua_State *L)
{
  lua_settop(L, 4);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  luaL_checktype(L, 2, LUA_TNUMBER);
  luaL_checktype(L, 3, LUA_TNUMBER);
  size_t rowstart = lua_tointeger(L, 2);
  size_t rowend = lua_tointeger(L, 3);
  if (rowstart > rowend)
    luaL_error(L, "Error in copy: start row is greater than end row");
  size_t idxstart = tk_matrix_index(L, m0, rowstart, 1);
  size_t idxend = tk_matrix_index(L, m0, rowend, m0->columns);
  if (lua_type(L, 4) == LUA_TNIL) {
    lua_pushlstring(L, (char *) &m0->data[idxstart], sizeof(tk_base_t) * (idxend - idxstart + 1));
    return 1;
  }
  const char *fmt = luaL_checkstring(L, 4);
  if (!strcmp(fmt, "u32")) {
    size_t n = idxend - idxstart + 1;
    unsigned int *out = malloc(sizeof(unsigned int) * n);
    if (!out)
      goto err_mem;
    for (size_t i = 0; i < n; i ++) {
      tk_base_t r = m0->data[i];
      if (r < 0) {
        free(out);
        goto err_negative;
      }
      if (r > (tk_base_t) UINT_MAX) {
        free(out);
        goto err_toobig;
      }
      out[i] = r;
    }
    lua_pushlstring(L, (char *) out, sizeof(unsigned int) * n);
    free(out);
    return 1;
  } else {
    goto err_format;
  }
err_format:
  return luaL_error(L, "unexpected format provided to matrix.raw");
err_mem:
  return luaL_error(L, "malloc error in matrix.raw");
err_negative:
  return luaL_error(L, "unexpected negative in matrix.raw");
err_toobig:
  return luaL_error(L, "unexpected negative in matrix.raw");
}

static inline int tk_matrix_view (lua_State *L)
{
  lua_settop(L, 2);
  tk_matrix_t *m0 = tk_matrix_peek(L, 1);
  bool vacate = lua_toboolean(L, 2);
  lua_pushlightuserdata(L, m0->data);
  lua_pushinteger(L, m0->rows);
  lua_pushinteger(L, m0->columns);
  if (vacate) {
    m0->data = NULL;
    m0->rows = 0;
    m0->columns = 0;
  }
  return 3;
}

static inline int tk_matrix_from_raw (lua_State *L)
{
  lua_settop(L, 2);
  luaL_checktype(L, 1, LUA_TSTRING);
  luaL_checktype(L, 2, LUA_TNUMBER);
  size_t columns = lua_tointeger(L, 2);
  size_t size;
  const char *data = lua_tolstring(L, 1, &size);
  if (size % columns != 0)
    luaL_error(L, "Length of raw string is not a multiple of provided column length");
  size_t rows = size / sizeof(tk_base_t) / columns;
  tk_matrix_t *m0 = _tk_matrix_create(L, rows, columns, NULL, NULL);
  memcpy(m0->data, data, size);
  return 1;
}

static inline int tk_matrix_from_view (lua_State *L)
{
  lua_settop(L, 4);
  tk_base_t *data = tk_lua_checkuserdata(L, 1, NULL);
  size_t rows = tk_lua_checkunsigned(L, 2);
  size_t columns = tk_lua_checkunsigned(L, 3);
  tk_matrix_t *m0 = lua_type(L, 4) == LUA_TNIL ? NULL : tk_matrix_peek(L, 4);
  _tk_matrix_create(L, rows, columns, data, m0);
  return 1;
}

static inline int tk_matrix_radd (lua_State *);
static inline int tk_matrix_rmult (lua_State *);
static inline int tk_matrix_ramax (lua_State *);
static inline int tk_matrix_sums (lua_State *);
static inline int tk_matrix_mmult (lua_State *);
static inline int tk_matrix_magnitude (lua_State *);
static luaL_Reg tk_matrix_extra_fns[];

static luaL_Reg tk_matrix_fns[] =
{
  { "create", tk_matrix_create },
  { "from_raw", tk_matrix_from_raw },
  { "raw", tk_matrix_raw },
  { "from_view", tk_matrix_from_view },
  { "view", tk_matrix_view },
  { "extend_raw", tk_matrix_extend_raw },
  { "get", tk_matrix_get },
  { "set", tk_matrix_set },
  { "shape", tk_matrix_shape },
  { "rows", tk_matrix_rows },
  { "columns", tk_matrix_columns },
  { "magnitude", tk_matrix_magnitude },
  { "mmult", tk_matrix_mmult },
  { "rmult", tk_matrix_rmult },
  { "radd", tk_matrix_radd },
  { "copy", tk_matrix_copy },
  { "sums", tk_matrix_sums },
  { "sum", tk_matrix_sum },
  { "exp", tk_matrix_exp },
  { "rmax", tk_matrix_rmax },
  { "rmin", tk_matrix_rmin },
  { "ramax", tk_matrix_ramax },
  { "tabulate", tk_matrix_tabulate },
  { "rorder", tk_matrix_rorder },
  { "reshape", tk_matrix_reshape },
  { "shrink", tk_matrix_shrink },
  { "sort", tk_matrix_sort },
  { NULL, NULL }
};

int TK_OPEN (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_matrix_fns);
  luaL_register(L, NULL, tk_matrix_extra_fns);
  luaL_newmetatable(L, TK_MT);
  lua_pushcfunction(L, tk_matrix_gc);
  lua_setfield(L, -2, "__gc");
  lua_pushvalue(L, -1);
  lua_setfield(L, -3, "mt_matrix");
  lua_pop(L, 1);
  tk_lua_callmod(L, 1, 1, "santoku.matrix.wrap", "wrap");
  return 1;
}
