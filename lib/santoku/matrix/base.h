#ifndef TK_MATRIX_BASE_H
#define TK_MATRIX_BASE_H

typedef struct {
  size_t rows;
  size_t columns;
  size_t values;
  tk_base_t *data;
} tk_matrix_t;

static inline tk_matrix_t **tk_matrix_peekp (lua_State *L, int i)
{
  return (tk_matrix_t **) luaL_checkudata(L, i, TK_MT);
}

static inline tk_matrix_t *tk_matrix_peek (lua_State *L, int i)
{
  tk_matrix_t **m0p = tk_matrix_peekp(L, i);
  return *m0p;
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

#endif
