#ifndef TK_SPANS_H
#define TK_SPANS_H

#include <santoku/shaped.h>
#include <santoku/ivec.h>

#define TK_SPANS_MT "tk_spans_t"



typedef struct {
  uint64_t n_cols;
  tk_ivec_t *offsets;
  tk_ivec_t **cols;
  char **names;
} tk_spans_t;

static inline tk_spans_t *tk_spans_peek (lua_State *L, int i, const char *name)
{
  tk_spans_t *S = (tk_spans_t *) luaL_checkudata(L, i, TK_SPANS_MT);
  if (S == NULL)
    tk_lua_verror(L, 2, name, "expected a spans");
  return S;
}

static inline uint64_t tk_spans_n (tk_spans_t *S)
{
  return S->n_cols > 0 ? S->cols[0]->n : 0;
}

static inline uint64_t tk_spans_docs (tk_spans_t *S)
{
  return S->offsets->n > 0 ? S->offsets->n - 1 : 0;
}

static inline int64_t tk_spans_colidx (tk_spans_t *S, const char *name)
{
  for (uint64_t c = 0; c < S->n_cols; c ++)
    if (strcmp(S->names[c], name) == 0)
      return (int64_t) c;
  return -1;
}

#endif
