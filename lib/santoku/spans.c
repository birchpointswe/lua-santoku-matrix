#include <santoku/iuset.h>
#include <santoku/spans.h>
#include <santoku/pvec/base.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define tk_spans_pair_lt(a, b) ((a).i < (b).i || ((a).i == (b).i && (a).p < (b).p))
KSORT_INIT(tk_spans_pairs, tk_pair_t, tk_spans_pair_lt)

static inline void tk_spans_suppress_unused (void)
{
  (void) ks_mergesort_tk_spans_pairs;
  (void) ks_heapmake_tk_spans_pairs;
  (void) ks_heapsort_tk_spans_pairs;
  (void) ks_ksmall_tk_spans_pairs;
  (void) ks_shuffle_tk_spans_pairs;
}

static inline void tk_spans_free_meta (tk_spans_t *S)
{
  if (S->names != NULL) {
    for (uint64_t c = 0; c < S->n_cols; c ++)
      free(S->names[c]);
    free(S->names);
    S->names = NULL;
  }
  free(S->cols);
  S->cols = NULL;
}

static int tk_spans_gc_lua (lua_State *L)
{
  tk_spans_t *S = (tk_spans_t *) luaL_checkudata(L, 1, TK_SPANS_MT);
  tk_spans_free_meta(S);
  return 0;
}



static tk_spans_t *tk_spans_alloc (lua_State *L, int inames)
{
  uint64_t n_cols = lua_objlen(L, inames);
  if (n_cols == 0)
    tk_lua_verror(L, 2, "spans", "at least one column name required");
  tk_spans_t *S = (tk_spans_t *) lua_newuserdata(L, sizeof(tk_spans_t));
  S->n_cols = n_cols;
  S->offsets = NULL;
  S->cols = (tk_ivec_t **) calloc(n_cols, sizeof(tk_ivec_t *));
  S->names = (char **) calloc(n_cols, sizeof(char *));
  if (S->cols == NULL || S->names == NULL) {
    free(S->cols);
    free(S->names);
    tk_lua_verror(L, 2, "spans", "allocation failed");
  }
  luaL_getmetatable(L, TK_SPANS_MT);
  lua_setmetatable(L, -2);
  tk_eph_init(L, lua_gettop(L));
  for (uint64_t c = 0; c < n_cols; c ++) {
    lua_rawgeti(L, inames, (int) c + 1);
    const char *nm = lua_tostring(L, -1);
    if (nm == NULL)
      tk_lua_verror(L, 2, "spans", "column names must be strings");
    S->names[c] = strdup(nm);
    lua_pop(L, 1);
  }
  return S;
}


static void tk_spans_init_children (lua_State *L, tk_spans_t *S, int is, uint64_t cap)
{
  S->offsets = tk_ivec_create(L, 1);
  S->offsets->a[0] = 0;
  tk_eph_anchor(L, is, lua_gettop(L), S->offsets);
  lua_pop(L, 1);
  for (uint64_t c = 0; c < S->n_cols; c ++) {
    S->cols[c] = tk_ivec_create(L, cap);
    S->cols[c]->n = 0;
    tk_eph_anchor(L, is, lua_gettop(L), S->cols[c]);
    lua_pop(L, 1);
  }
}

static int tk_spans_create_lua (lua_State *L)
{
  lua_settop(L, 1);
  luaL_checktype(L, 1, LUA_TTABLE);
  tk_spans_t *S = tk_spans_alloc(L, 1);
  tk_spans_init_children(L, S, lua_gettop(L), 0);
  return 1;
}

static int tk_spans_push_lua (lua_State *L)
{
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  for (uint64_t c = 0; c < S->n_cols; c ++)
    tk_ivec_push(S->cols[c], luaL_checkinteger(L, (int) c + 2));
  lua_settop(L, 1);
  return 1;
}

static int tk_spans_doc_lua (lua_State *L)
{
  lua_settop(L, 1);
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  tk_ivec_push(S->offsets, (int64_t) tk_spans_n(S));
  lua_settop(L, 1);
  return 1;
}

static int tk_spans_n_lua (lua_State *L)
{
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  lua_pushinteger(L, (lua_Integer) tk_spans_n(S));
  return 1;
}

static int tk_spans_n_docs_lua (lua_State *L)
{
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  lua_pushinteger(L, (lua_Integer) tk_spans_docs(S));
  return 1;
}

static int tk_spans_offsets_lua (lua_State *L)
{
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  tk_eph_get(L, 1, S->offsets);
  return 1;
}

static int tk_spans_col_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  int64_t c = tk_spans_colidx(S, luaL_checkstring(L, 2));
  if (c < 0)
    return tk_lua_verror(L, 2, "spans", "no such column");
  tk_eph_get(L, 1, S->cols[c]);
  return 1;
}

static int tk_spans_names_lua (lua_State *L)
{
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  lua_newtable(L);
  for (uint64_t c = 0; c < S->n_cols; c ++) {
    lua_pushstring(L, S->names[c]);
    lua_rawseti(L, -2, (int) c + 1);
  }
  return 1;
}


static int tk_spans_filter_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  tk_ivec_t *mask = tk_ivec_peek(L, 2, "mask");
  uint64_t n = tk_spans_n(S);
  if (mask->n < n)
    return tk_lua_verror(L, 2, "spans", "mask shorter than span count");
  uint64_t nd = tk_spans_docs(S);
  uint64_t w = 0;
  for (uint64_t d = 0; d < nd; d ++) {
    int64_t lo = S->offsets->a[d], hi = S->offsets->a[d + 1];
    for (int64_t j = lo; j < hi; j ++)
      if (mask->a[j] != 0) {
        for (uint64_t c = 0; c < S->n_cols; c ++)
          S->cols[c]->a[w] = S->cols[c]->a[j];
        w ++;
      }
    S->offsets->a[d + 1] = (int64_t) w;
  }
  for (uint64_t c = 0; c < S->n_cols; c ++)
    S->cols[c]->n = w;
  lua_settop(L, 1);
  return 1;
}


static int tk_spans_docs_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  tk_ivec_t *ids = tk_ivec_peek(L, 2, "ids");
  uint64_t nd = tk_spans_docs(S);
  lua_createtable(L, (int) S->n_cols, 0);
  int inames = lua_gettop(L);
  for (uint64_t c = 0; c < S->n_cols; c ++) {
    lua_pushstring(L, S->names[c]);
    lua_rawseti(L, inames, (int) c + 1);
  }
  tk_spans_t *O = tk_spans_alloc(L, inames);
  int io = lua_gettop(L);
  tk_spans_init_children(L, O, io, 0);
  for (uint64_t i = 0; i < ids->n; i ++) {
    int64_t d = ids->a[i];
    if (d < 0 || (uint64_t) d >= nd)
      return tk_lua_verror(L, 2, "spans", "doc id out of range");
    for (int64_t j = S->offsets->a[d]; j < S->offsets->a[d + 1]; j ++)
      for (uint64_t c = 0; c < S->n_cols; c ++)
        tk_ivec_push(O->cols[c], S->cols[c]->a[j]);
    tk_ivec_push(O->offsets, (int64_t) tk_spans_n(O));
  }
  lua_remove(L, inames);
  return 1;
}


static int tk_spans_append_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  tk_spans_t *B = tk_spans_peek(L, 2, "other");
  if (S->n_cols != B->n_cols)
    return tk_lua_verror(L, 2, "spans", "append requires matching columns");
  for (uint64_t c = 0; c < S->n_cols; c ++)
    if (strcmp(S->names[c], B->names[c]) != 0)
      return tk_lua_verror(L, 2, "spans", "append requires matching column names");
  uint64_t base = tk_spans_n(S);
  uint64_t nb = tk_spans_n(B);
  for (uint64_t c = 0; c < S->n_cols; c ++) {
    if (tk_ivec_ensure(S->cols[c], base + nb) != 0)
      return tk_lua_verror(L, 2, "spans", "allocation failed");
    memcpy(S->cols[c]->a + base, B->cols[c]->a, nb * sizeof(int64_t));
    S->cols[c]->n = base + nb;
  }
  uint64_t bd = tk_spans_docs(B);
  for (uint64_t d = 1; d <= bd; d ++)
    tk_ivec_push(S->offsets, (int64_t) (base + (uint64_t) B->offsets->a[d]));
  lua_settop(L, 1);
  return 1;
}


static int tk_spans_sort_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  int64_t k = tk_spans_colidx(S, luaL_checkstring(L, 2));
  if (k < 0)
    return tk_lua_verror(L, 2, "spans", "no such column");
  uint64_t nd = tk_spans_docs(S);
  uint64_t maxlen = 0;
  for (uint64_t d = 0; d < nd; d ++) {
    uint64_t len = (uint64_t) (S->offsets->a[d + 1] - S->offsets->a[d]);
    if (len > maxlen)
      maxlen = len;
  }
  if (maxlen < 2) {
    lua_settop(L, 1);
    return 1;
  }
  tk_pair_t *pairs = (tk_pair_t *) malloc(maxlen * sizeof(tk_pair_t));
  int64_t *scratch = (int64_t *) malloc(maxlen * sizeof(int64_t));
  if (pairs == NULL || scratch == NULL) {
    free(pairs);
    free(scratch);
    return tk_lua_verror(L, 2, "spans", "allocation failed");
  }
  for (uint64_t d = 0; d < nd; d ++) {
    int64_t lo = S->offsets->a[d];
    uint64_t len = (uint64_t) (S->offsets->a[d + 1] - lo);
    if (len < 2)
      continue;
    for (uint64_t i = 0; i < len; i ++)
      pairs[i] = tk_pair(S->cols[k]->a[lo + (int64_t) i], (int64_t) i);
    ks_introsort(tk_spans_pairs, len, pairs);
    for (uint64_t c = 0; c < S->n_cols; c ++) {
      for (uint64_t i = 0; i < len; i ++)
        scratch[i] = S->cols[c]->a[lo + pairs[i].p];
      memcpy(S->cols[c]->a + lo, scratch, len * sizeof(int64_t));
    }
  }
  free(pairs);
  free(scratch);
  lua_settop(L, 1);
  return 1;
}

static int tk_spans_eq_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_spans_t *a = tk_spans_peek(L, 1, "spans");
  tk_spans_t *b = tk_spans_peek(L, 2, "other");
  bool r = a->n_cols == b->n_cols
    && a->offsets->n == b->offsets->n
    && tk_spans_n(a) == tk_spans_n(b)
    && tk_ivec_eq(a->offsets, b->offsets, 0, a->offsets->n);
  for (uint64_t c = 0; r && c < a->n_cols; c ++)
    r = strcmp(a->names[c], b->names[c]) == 0
      && tk_ivec_eq(a->cols[c], b->cols[c], 0, a->cols[c]->n);
  lua_pushboolean(L, r);
  return 1;
}



static int tk_spans_surfaces_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  luaL_checktype(L, 2, LUA_TTABLE);
  bool lower = lua_toboolean(L, 3);
  int64_t cs = tk_spans_colidx(S, "s");
  int64_t ce = tk_spans_colidx(S, "e");
  if (cs < 0 || ce < 0)
    return tk_lua_verror(L, 2, "spans", "surfaces requires columns \"s\" and \"e\"");
  uint64_t nd = tk_spans_docs(S);
  tk_ivec_t *ids = tk_ivec_create(L, tk_spans_n(S));
  int iids = lua_gettop(L);
  ids->n = 0;
  lua_newtable(L);
  int imap = lua_gettop(L);
  lua_newtable(L);
  int ilist = lua_gettop(L);
  uint64_t n_unique = 0;
  char *buf = NULL;
  size_t bufcap = 0;
  for (uint64_t d = 0; d < nd; d ++) {
    lua_rawgeti(L, 2, (int) d + 1);
    size_t tlen;
    const char *txt = lua_tolstring(L, -1, &tlen);
    if (txt == NULL)
      return tk_lua_verror(L, 2, "spans", "surfaces: texts must be strings");
    for (int64_t j = S->offsets->a[d]; j < S->offsets->a[d + 1]; j ++) {
      int64_t s = S->cols[cs]->a[j], e = S->cols[ce]->a[j];
      if (s < 0) s = 0;
      if (e > (int64_t) tlen) e = (int64_t) tlen;
      size_t len = e > s ? (size_t) (e - s) : 0;
      if (lower) {
        if (len > bufcap) {
          char *nb = (char *) realloc(buf, len);
          if (nb == NULL) {
            free(buf);
            return tk_lua_verror(L, 2, "spans", "allocation failed");
          }
          buf = nb;
          bufcap = len;
        }
        for (size_t i = 0; i < len; i ++)
          buf[i] = (char) tolower((unsigned char) txt[s + (int64_t) i]);
        lua_pushlstring(L, buf, len);
      } else {
        lua_pushlstring(L, txt + s, len);
      }
      lua_pushvalue(L, -1);
      lua_rawget(L, imap);
      if (lua_isnil(L, -1)) {
        lua_pop(L, 1);
        n_unique ++;
        lua_pushvalue(L, -1);
        lua_pushinteger(L, (lua_Integer) (n_unique - 1));
        lua_rawset(L, imap);
        lua_pushvalue(L, -1);
        lua_rawseti(L, ilist, (int) n_unique);
        tk_ivec_push(ids, (int64_t) (n_unique - 1));
        lua_pop(L, 1);
      } else {
        int64_t id = (int64_t) lua_tointeger(L, -1);
        lua_pop(L, 2);
        tk_ivec_push(ids, id);
      }
    }
    lua_pop(L, 1);
  }
  free(buf);
  lua_pushvalue(L, iids);
  lua_pushvalue(L, ilist);
  return 2;
}

static int tk_spans_persist_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_spans_t *S = tk_spans_peek(L, 1, "spans");
  FILE *fh = tk_lua_fopen(L, luaL_checkstring(L, 2), "w");
  char magic[4] = { 'T', 'K', 's', 'p' };
  uint8_t version = 1;
  uint64_t nc = S->n_cols, no = S->offsets->n, n = tk_spans_n(S);
  tk_lua_fwrite(L, magic, 4, 1, fh);
  tk_lua_fwrite(L, (char *) &version, 1, 1, fh);
  tk_lua_fwrite(L, (char *) &nc, sizeof(uint64_t), 1, fh);
  tk_lua_fwrite(L, (char *) &no, sizeof(uint64_t), 1, fh);
  tk_lua_fwrite(L, (char *) &n, sizeof(uint64_t), 1, fh);
  for (uint64_t c = 0; c < nc; c ++) {
    uint64_t l = strlen(S->names[c]);
    tk_lua_fwrite(L, (char *) &l, sizeof(uint64_t), 1, fh);
    tk_lua_fwrite(L, S->names[c], l, 1, fh);
  }
  tk_lua_fwrite(L, (char *) S->offsets->a, sizeof(int64_t) * no, 1, fh);
  for (uint64_t c = 0; c < nc; c ++)
    tk_lua_fwrite(L, (char *) S->cols[c]->a, sizeof(int64_t) * n, 1, fh);
  tk_lua_fclose(L, fh);
  return 0;
}

static int tk_spans_load_lua (lua_State *L)
{
  lua_settop(L, 1);
  FILE *fh = tk_lua_fopen(L, luaL_checkstring(L, 1), "r");
  char magic[4];
  uint8_t version;
  uint64_t nc, no, n;
  tk_lua_fread(L, magic, 4, 1, fh);
  if (memcmp(magic, "TKsp", 4) != 0)
    return tk_lua_verror(L, 2, "spans", "load: bad magic");
  tk_lua_fread(L, (char *) &version, 1, 1, fh);
  if (version != 1)
    return tk_lua_verror(L, 2, "spans", "load: unsupported version");
  tk_lua_fread(L, (char *) &nc, sizeof(uint64_t), 1, fh);
  tk_lua_fread(L, (char *) &no, sizeof(uint64_t), 1, fh);
  tk_lua_fread(L, (char *) &n, sizeof(uint64_t), 1, fh);
  lua_createtable(L, (int) nc, 0);
  int inames = lua_gettop(L);
  for (uint64_t c = 0; c < nc; c ++) {
    uint64_t l;
    tk_lua_fread(L, (char *) &l, sizeof(uint64_t), 1, fh);
    char *nm = (char *) malloc(l);
    if (nm == NULL)
      return tk_lua_verror(L, 2, "spans", "allocation failed");
    tk_lua_fread(L, nm, l, 1, fh);
    lua_pushlstring(L, nm, l);
    free(nm);
    lua_rawseti(L, inames, (int) c + 1);
  }
  tk_spans_t *S = tk_spans_alloc(L, inames);
  int is = lua_gettop(L);
  tk_spans_init_children(L, S, is, n);
  if (tk_ivec_ensure(S->offsets, no) != 0)
    return tk_lua_verror(L, 2, "spans", "allocation failed");
  S->offsets->n = no;
  tk_lua_fread(L, (char *) S->offsets->a, sizeof(int64_t) * no, 1, fh);
  for (uint64_t c = 0; c < nc; c ++) {
    S->cols[c]->n = n;
    tk_lua_fread(L, (char *) S->cols[c]->a, sizeof(int64_t) * n, 1, fh);
  }
  tk_lua_fclose(L, fh);
  lua_remove(L, inames);
  return 1;
}

static luaL_Reg tk_spans_mt_fns[] = {
  { "push", tk_spans_push_lua },
  { "doc", tk_spans_doc_lua },
  { "n", tk_spans_n_lua },
  { "n_docs", tk_spans_n_docs_lua },
  { "offsets", tk_spans_offsets_lua },
  { "col", tk_spans_col_lua },
  { "names", tk_spans_names_lua },
  { "filter", tk_spans_filter_lua },
  { "docs", tk_spans_docs_lua },
  { "append", tk_spans_append_lua },
  { "sort", tk_spans_sort_lua },
  { "eq", tk_spans_eq_lua },
  { "surfaces", tk_spans_surfaces_lua },
  { "persist", tk_spans_persist_lua },
  { NULL, NULL }
};

static luaL_Reg tk_spans_fns[] = {
  { "create", tk_spans_create_lua },
  { "load", tk_spans_load_lua },
  { NULL, NULL }
};

int luaopen_santoku_spans (lua_State *L)
{
  tk_spans_suppress_unused();
  luaL_newmetatable(L, TK_SPANS_MT);
  lua_pushcfunction(L, tk_spans_gc_lua);
  lua_setfield(L, -2, "__gc");
  lua_newtable(L);
  luaL_register(L, NULL, tk_spans_mt_fns);
  lua_setfield(L, -2, "__index");
  lua_pop(L, 1);
  lua_newtable(L);
  luaL_register(L, NULL, tk_spans_fns);
  return 1;
}
