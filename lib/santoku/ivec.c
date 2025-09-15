#include <santoku/ivec.h>
#include <santoku/cvec.h>
#include <santoku/iuset.h>
#include <string.h>

static inline int tk_ivec_bits_rearrange_lua (lua_State *L) {
  lua_settop(L, 3);
  tk_ivec_t *m0 = tk_ivec_peek(L, 1, "set_bits");
  tk_ivec_t *ids = tk_ivec_peek(L, 2, "ids");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "features");
  tk_ivec_bits_rearrange(m0, ids, n_features);
  return 0;
}

static inline int tk_ivec_bits_score_chi2_lua (
  lua_State *L
) {
  lua_settop(L, 6);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  unsigned int n_threads = tk_threads_getn(L, 6, "threads", NULL);
  char *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_isnil(L, 2)) {

  } else if (tk_lua_testuserdata(L, 2, "tk_cvec_t")) {
    tk_cvec_t *cvec = tk_cvec_peek(L, 2, "codes");
    codes = cvec->a;
  } else {
    tk_ivec_t *m1 = tk_ivec_peek(L, 2, "labels");
    labels = m1;
  }
  tk_ivec_bits_score_chi2(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, n_threads);
  return 1;
}

static inline int tk_ivec_bits_score_mi_lua (
  lua_State *L
) {
  lua_settop(L, 6);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  unsigned int n_threads = tk_threads_getn(L, 6, "threads", NULL);
  char *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_isnil(L, 2)) {

  } else if (tk_lua_testuserdata(L, 2, "tk_cvec_t")) {
    tk_cvec_t *cvec = tk_cvec_peek(L, 2, "codes");
    codes = cvec->a;
  } else {
    tk_ivec_t *m1 = tk_ivec_peek(L, 2, "labels");
    labels = m1;
  }
  tk_ivec_bits_score_mi(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, n_threads);
  return 1;
}

static inline int tk_ivec_bits_top_mi_lua (lua_State *L)
{
  lua_settop(L, 7);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  uint64_t top_k = tk_lua_checkunsigned(L, 6, "top_k");
  unsigned int n_threads = tk_threads_getn(L, 7, "threads", NULL);
  char *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_isnil(L, 2)) {

  } else if (tk_lua_testuserdata(L, 2, "tk_cvec_t")) {
    tk_cvec_t *cvec = tk_cvec_peek(L, 2, "codes");
    codes = cvec->a;
  } else {
    tk_ivec_t *m1 = tk_ivec_peek(L, 2, "labels");
    labels = m1;
  }
  tk_ivec_bits_top_mi(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, top_k, n_threads);
  return 2;
}

static inline int tk_ivec_bits_top_chi2_lua (lua_State *L)
{
  lua_settop(L, 7);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 3, "samples");
  uint64_t n_visible = tk_lua_checkunsigned(L, 4, "visible");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 5, "hidden");
  uint64_t top_k = tk_lua_checkunsigned(L, 6, "top_k");
  unsigned int n_threads = tk_threads_getn(L, 7, "threads", NULL);
  char *codes = NULL;
  tk_ivec_t *labels = NULL;
  if (lua_isnil(L, 2)) {

  } else if (tk_lua_testuserdata(L, 2, "tk_cvec_t")) {
    tk_cvec_t *cvec = tk_cvec_peek(L, 2, "codes");
    codes = cvec->a;
  } else {
    tk_ivec_t *m1 = tk_ivec_peek(L, 2, "labels");
    labels = m1;
  }
  tk_ivec_bits_top_chi2(L, set_bits, codes, labels, n_samples, n_visible, n_hidden, top_k, n_threads);
  return 2;
}

static inline int tk_ivec_bits_score_entropy_lua (
  lua_State *L
) {
  lua_settop(L, 4);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 3, "hidden");
  unsigned int n_threads = tk_threads_getn(L, 4, "threads", NULL);
  tk_ivec_bits_score_entropy(L, set_bits, n_samples, n_hidden, n_threads);
  return 1;
}

static inline int tk_ivec_bits_top_entropy_lua (lua_State *L)
{
  lua_settop(L, 5);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_hidden = tk_lua_checkunsigned(L, 3, "hidden");
  uint64_t top_k = tk_lua_checkunsigned(L, 4, "top_k");
  unsigned int n_threads = tk_threads_getn(L, 5, "threads", NULL);
  tk_ivec_bits_top_entropy(L, set_bits, n_samples, n_hidden, top_k, n_threads);
  return 2;
}

static inline int tk_ivec_bits_filter_lua (lua_State *L)
{
  int n_args = lua_gettop(L);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  tk_ivec_t *top_v = lua_isnil(L, 2) ? NULL : tk_ivec_peek(L, 2, "top_v");
  tk_ivec_t *sample_ids = NULL;
  uint64_t n_visible;
  if (n_args == 4) {
    sample_ids = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "sample_ids");
    n_visible = tk_lua_checkunsigned(L, 4, "visible");
  } else {

    n_visible = tk_lua_checkunsigned(L, 3, "visible");
  }
  tk_ivec_bits_filter(set_bits, top_v, sample_ids, n_visible);
  return 0;
}

static inline int tk_ivec_bits_copy_lua (lua_State *L)
{
  int n_args = lua_gettop(L);
  tk_ivec_t *dest = tk_ivec_peek(L, 1, "dest");
  tk_ivec_t *src_bits = tk_ivec_peek(L, 2, "src_bits");
  tk_ivec_t *selected_features = lua_isnil(L, 3) ? NULL : tk_ivec_peek(L, 3, "selected_features");
  tk_ivec_t *sample_ids = NULL;
  uint64_t n_visible;
  uint64_t dest_sample = 0;

  if (n_args == 6) {
    sample_ids = lua_isnil(L, 4) ? NULL : tk_ivec_peek(L, 4, "sample_ids");
    n_visible = tk_lua_checkunsigned(L, 5, "visible");
    dest_sample = tk_lua_checkunsigned(L, 6, "dest_sample");
  } else if (n_args == 5) {
    if (lua_type(L, 4) == LUA_TNUMBER && lua_type(L, 5) == LUA_TNUMBER) {
      n_visible = tk_lua_checkunsigned(L, 4, "visible");
      dest_sample = tk_lua_checkunsigned(L, 5, "dest_sample");
    } else {
      sample_ids = lua_isnil(L, 4) ? NULL : tk_ivec_peek(L, 4, "sample_ids");
      n_visible = tk_lua_checkunsigned(L, 5, "visible");
    }
  } else {
    n_visible = tk_lua_checkunsigned(L, 4, "visible");
  }
  tk_ivec_bits_copy(dest, src_bits, selected_features, sample_ids, n_visible, dest_sample);
  return 0;
}

static inline int tk_ivec_bits_to_cvec_lua (lua_State *L)
{
  lua_settop(L, 4);
  tk_ivec_t *set_bits = tk_ivec_peek(L, 1, "set_bits");
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "features");
  bool flip_interleave = lua_toboolean(L, 4);
  tk_ivec_bits_to_cvec(L, set_bits, n_samples, n_features, flip_interleave);
  return 1;
}

static inline int tk_ivec_bits_from_cvec_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_cvec_t *cvec = tk_cvec_peek(L, 1, "bitmap");
  const char *bm = cvec->a;
  uint64_t n_samples = tk_lua_checkunsigned(L, 2, "samples");
  uint64_t n_features = tk_lua_checkunsigned(L, 3, "features");
  tk_ivec_bits_from_cvec(L, bm, n_samples, n_features);
  return 1;
}


static inline void tk_ivec_bits_extend_cvec_helper (
  lua_State *L,
  tk_ivec_t *base,
  tk_cvec_t *ext_cvec,
  uint64_t n_feat,
  uint64_t n_extfeat
) {
  tk_ivec_t *ext = tk_cvec_bits_to_ivec(L, ext_cvec, n_extfeat);
  tk_ivec_bits_extend(base, ext, n_feat, n_extfeat);
  tk_ivec_destroy(ext);
  lua_remove(L, -1);
}

static inline int tk_ivec_bits_extend_lua (lua_State *L)
{
  int nargs = lua_gettop(L);


  if (nargs == 4) {

    tk_ivec_t *base = tk_ivec_peek(L, 1, "base_bits");
    uint64_t n_feat = tk_lua_checkunsigned(L, 3, "features");
    uint64_t n_extfeat = tk_lua_checkunsigned(L, 4, "extended");
    tk_ivec_t *ext_ivec = tk_ivec_peekopt(L, 2);
    if (ext_ivec) {
      tk_ivec_bits_extend(base, ext_ivec, n_feat, n_extfeat);
    } else {
      tk_cvec_t *ext_cvec = tk_cvec_peek(L, 2, "ext_bits");
      tk_ivec_bits_extend_cvec_helper(L, base, ext_cvec, n_feat, n_extfeat);
    }
  } else if (nargs == 6 || nargs == 7) {
    tk_ivec_t *base = tk_ivec_peek(L, 1, "base_bits");
    tk_ivec_t *aids = tk_ivec_peek(L, 3, "aids");
    tk_ivec_t *bids = tk_ivec_peek(L, 4, "bids");
    uint64_t n_feat = tk_lua_checkunsigned(L, 5, "features");
    uint64_t n_extfeat = tk_lua_checkunsigned(L, 6, "extended");

    bool project = false;
    if (nargs == 7) {
      project = lua_toboolean(L, 7);
    }

    tk_ivec_t *ext_ivec = tk_ivec_peekopt(L, 2);
    if (ext_ivec) {
      tk_ivec_bits_extend_mapped(base, ext_ivec, aids, bids, n_feat, n_extfeat, project);
    } else {

      tk_cvec_t *ext_cvec = tk_cvec_peek(L, 2, "ext_bits");
      tk_ivec_t *ext = tk_cvec_bits_to_ivec(L, ext_cvec, n_extfeat);
      tk_ivec_bits_extend_mapped(base, ext, aids, bids, n_feat, n_extfeat, project);
      tk_ivec_destroy(ext);
      lua_remove(L, -1);
    }
  } else {
    return luaL_error(L, "bits_extend expects 4, 6, or 7 arguments, got %d", nargs);
  }
  return 0;
}

static inline int tk_ivec_copy_pkeys_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_ivec_t *m0 = tk_ivec_peek(L, 1, "dest");
  tk_pvec_t *m1 = tk_pvec_peek(L, 2, "source");
  int64_t start, end, dest;
  if (t == 2) {
    start = 0;
    end = (int64_t) m1->n;
    dest = (int64_t) m0->n;
  } else if (t == 3) {
    start = 0;
    end = (int64_t) m1->n;
    dest = tk_lua_checkinteger(L, 3, "dest");
  } else {
    start = tk_lua_checkinteger(L, 3, "start");
    end = tk_lua_checkinteger(L, 4, "end");
    dest = tk_lua_checkinteger(L, 5, "dest");
  }
  tk_ivec_copy_pkeys(m0, m1, start, end, dest);
  return 0;
}

static inline int tk_ivec_copy_rkeys_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_ivec_t *m0 = tk_ivec_peek(L, 1, "dest");
  tk_rvec_t *m1 = tk_rvec_peek(L, 2, "source");
  int64_t start, end, dest;
  if (t == 2) {
    start = 0;
    end = (int64_t) m1->n;
    dest = (int64_t) m0->n;
  } else if (t == 3) {
    start = 0;
    end = (int64_t) m1->n;
    dest = tk_lua_checkinteger(L, 3, "dest");
  } else {
    start = tk_lua_checkinteger(L, 3, "start");
    end = tk_lua_checkinteger(L, 4, "end");
    dest = tk_lua_checkinteger(L, 5, "dest");
  }
  tk_ivec_copy_rkeys(m0, m1, start, end, dest);
  return 0;
}

static inline int tk_ivec_copy_pvalues_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_ivec_t *m0 = tk_ivec_peek(L, 1, "dest");
  tk_pvec_t *m1 = tk_pvec_peek(L, 2, "source");
  int64_t start, end, dest;
  if (t == 2) {
    start = 0;
    end = (int64_t) m1->n;
    dest = (int64_t) m0->n;
  } else if (t == 3) {
    start = 0;
    end = (int64_t) m1->n;
    dest = tk_lua_checkinteger(L, 3, "dest");
  } else {
    start = tk_lua_checkinteger(L, 3, "start");
    end = tk_lua_checkinteger(L, 4, "end");
    dest = tk_lua_checkinteger(L, 5, "dest");
  }
  tk_ivec_copy_pvalues(m0, m1, start, end, dest);
  return 0;
}

static inline int tk_ivec_copy_rvalues_lua (lua_State *L)
{
  int t = lua_gettop(L);
  tk_ivec_t *m0 = tk_ivec_peek(L, 1, "dest");
  tk_rvec_t *m1 = tk_rvec_peek(L, 2, "source");
  int64_t start, end, dest;
  if (t == 2) {
    start = 0;
    end = (int64_t) m1->n;
    dest = (int64_t) m0->n;
  } else if (t == 3) {
    start = 0;
    end = (int64_t) m1->n;
    dest = tk_lua_checkinteger(L, 3, "dest");
  } else {
    start = tk_lua_checkinteger(L, 3, "start");
    end = tk_lua_checkinteger(L, 4, "end");
    dest = tk_lua_checkinteger(L, 5, "dest");
  }
  tk_ivec_copy_rvalues(m0, m1, start, end, dest);
  return 0;
}

static inline int tk_ivec_set_stats_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_dvec_t *weights = lua_isnil(L, 3) ? NULL : tk_dvec_peek(L, 3, "weights");
  double inter_w, sum_a, sum_b;
  tk_ivec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  lua_pushnumber(L, inter_w);
  lua_pushnumber(L, sum_a);
  lua_pushnumber(L, sum_b);
  return 3;
}

static inline int tk_ivec_set_similarity_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_dvec_t *weights = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_dvec_peek(L, 3, "weights") : NULL;
  tk_ivec_sim_type_t type = TK_IVEC_JACCARD;
  if (nargs >= 4 && lua_isstring(L, 4)) {
    const char *type_str = lua_tostring(L, 4);
    if (strcmp(type_str, "jaccard") == 0) {
      type = TK_IVEC_JACCARD;
    } else if (strcmp(type_str, "overlap") == 0) {
      type = TK_IVEC_OVERLAP;
    } else if (strcmp(type_str, "dice") == 0) {
      type = TK_IVEC_DICE;
    } else if (strcmp(type_str, "tversky") == 0) {
      type = TK_IVEC_TVERSKY;
    }
  }
  double tversky_alpha = (nargs >= 5) ? luaL_checknumber(L, 5) : 0.5;
  double tversky_beta = (nargs >= 6) ? luaL_checknumber(L, 6) : 0.5;
  double sim = tk_ivec_set_similarity(a->a, a->n, b->a, b->n, weights, type, tversky_alpha, tversky_beta);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_ivec_set_jaccard_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_dvec_t *weights = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_dvec_peek(L, 3, "weights") : NULL;
  double inter_w, sum_a, sum_b;
  tk_ivec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  double sim = tk_ivec_set_jaccard(inter_w, sum_a, sum_b);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_ivec_set_overlap_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_dvec_t *weights = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_dvec_peek(L, 3, "weights") : NULL;
  double inter_w, sum_a, sum_b;
  tk_ivec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  double sim = tk_ivec_set_overlap(inter_w, sum_a, sum_b);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_ivec_set_dice_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_dvec_t *weights = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_dvec_peek(L, 3, "weights") : NULL;
  double inter_w, sum_a, sum_b;
  tk_ivec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  double sim = tk_ivec_set_dice(inter_w, sum_a, sum_b);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_ivec_set_tversky_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");







  tk_dvec_t *weights = NULL;
  double alpha = 0.5;
  double beta = 0.5;

  if (nargs >= 3) {
    if (lua_isnumber(L, 3)) {
      alpha = luaL_checknumber(L, 3);
      beta = (nargs >= 4) ? luaL_checknumber(L, 4) : 0.5;
    } else if (!lua_isnil(L, 3)) {
      weights = tk_dvec_peek(L, 3, "weights");
      alpha = (nargs >= 4) ? luaL_checknumber(L, 4) : 0.5;
      beta = (nargs >= 5) ? luaL_checknumber(L, 5) : 0.5;
    }
  }

  double inter_w, sum_a, sum_b;
  tk_ivec_set_stats(a->a, a->n, b->a, b->n, weights, &inter_w, &sum_a, &sum_b);
  double sim = tk_ivec_set_tversky(inter_w, sum_a, sum_b, alpha, beta);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_ivec_set_weights_by_rank_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_ivec_t *features = tk_ivec_peek(L, 1, "features");
  tk_dvec_t *weights = (nargs >= 2 && !lua_isnil(L, 2)) ? tk_dvec_peek(L, 2, "weights") : NULL;
  tk_ivec_t *ranks = (nargs >= 3 && !lua_isnil(L, 3)) ? tk_ivec_peek(L, 3, "ranks") : NULL;
  uint64_t n_ranks = (nargs >= 4) ? tk_lua_checkunsigned(L, 4, "n_ranks") : 1;
  tk_dvec_t *weights_by_rank = tk_dvec_create(L, n_ranks, 0, 0);
  tk_ivec_set_weights_by_rank(features->a, features->n, weights, ranks, n_ranks, weights_by_rank->a);
  weights_by_rank->n = n_ranks;
  return 1;
}

static inline int tk_ivec_set_similarity_by_rank_lua (lua_State *L)
{
  int nargs = lua_gettop(L);
  tk_dvec_t *wacc = tk_dvec_peek(L, 1, "wacc");
  int64_t vsid = tk_lua_checkinteger(L, 2, "vsid");
  tk_dvec_t *q_weights_by_rank = tk_dvec_peek(L, 3, "q_weights_by_rank");
  tk_dvec_t *e_weights_by_rank = tk_dvec_peek(L, 4, "e_weights_by_rank");
  uint64_t n_ranks = tk_lua_checkunsigned(L, 5, "n_ranks");
  int64_t rank_decay_window = (nargs >= 6) ? tk_lua_checkinteger(L, 6, "rank_decay_window") : -1;
  double rank_decay_sigma = (nargs >= 7) ? luaL_checknumber(L, 7) : 10.0;
  double rank_decay_floor = (nargs >= 8) ? luaL_checknumber(L, 8) : 0.01;
  tk_ivec_sim_type_t type = TK_IVEC_JACCARD;
  if (nargs >= 9 && lua_isstring(L, 9)) {
    const char *type_str = lua_tostring(L, 9);
    if (strcmp(type_str, "jaccard") == 0) {
      type = TK_IVEC_JACCARD;
    } else if (strcmp(type_str, "overlap") == 0) {
      type = TK_IVEC_OVERLAP;
    } else if (strcmp(type_str, "dice") == 0) {
      type = TK_IVEC_DICE;
    } else if (strcmp(type_str, "tversky") == 0) {
      type = TK_IVEC_TVERSKY;
    }
  }
  double tversky_alpha = (nargs >= 10) ? luaL_checknumber(L, 10) : 0.5;
  double tversky_beta = (nargs >= 11) ? luaL_checknumber(L, 11) : 0.5;
  double sim = tk_ivec_set_similarity_by_rank(
    wacc, vsid, q_weights_by_rank->a,
    e_weights_by_rank->a, n_ranks,
    rank_decay_window,
    rank_decay_sigma,
    rank_decay_floor, type,
    tversky_alpha,
    tversky_beta);
  lua_pushnumber(L, sim);
  return 1;
}

static inline int tk_ivec_set_find_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_ivec_t *vec = tk_ivec_peek(L, 1, "vector");
  int64_t value = luaL_checkinteger(L, 2);
  int64_t insert_idx = -1;
  int64_t found_idx = tk_ivec_set_find(vec, value, &insert_idx);
  lua_pushinteger(L, found_idx);
  lua_pushinteger(L, insert_idx);
  return 2;
}

static inline int tk_ivec_set_intersect_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_ivec_t *out = NULL;
  if (!lua_isnil(L, 3))
    out = tk_ivec_peek(L, 3, "output");
  tk_ivec_set_intersect(L, a, b, out);
  return out == NULL ? 1 : 0;
}

static inline int tk_ivec_set_union_lua (lua_State *L)
{
  lua_settop(L, 3);
  tk_ivec_t *a = tk_ivec_peek(L, 1, "a");
  tk_ivec_t *b = tk_ivec_peek(L, 2, "b");
  tk_ivec_t *out = NULL;
  if (!lua_isnil(L, 3))
    out = tk_ivec_peek(L, 3, "output");
  tk_ivec_set_union(L, a, b, out);
  return out == NULL ? 1 : 0;
}

static inline int tk_ivec_lookup_lua (lua_State *L)
{
  lua_settop(L, 2);
  tk_ivec_t *indices = tk_ivec_peek(L, 1, "indices");
  tk_ivec_t *source = tk_ivec_peek(L, 2, "source");
  tk_ivec_lookup(indices, source);
  return 0;
}

static luaL_Reg tk_ivec_lua_mt_ext2_fns[] =
{
  { "copy_pkeys", tk_ivec_copy_pkeys_lua },
  { "copy_rkeys", tk_ivec_copy_rkeys_lua },
  { "copy_pvalues", tk_ivec_copy_pvalues_lua },
  { "copy_rvalues", tk_ivec_copy_rvalues_lua },
  { "bits_top_chi2", tk_ivec_bits_top_chi2_lua },
  { "bits_top_mi", tk_ivec_bits_top_mi_lua },
  { "bits_top_entropy", tk_ivec_bits_top_entropy_lua },
  { "bits_score_chi2", tk_ivec_bits_score_chi2_lua },
  { "bits_score_mi", tk_ivec_bits_score_mi_lua },
  { "bits_score_entropy", tk_ivec_bits_score_entropy_lua },
  { "bits_filter", tk_ivec_bits_filter_lua },
  { "bits_copy", tk_ivec_bits_copy_lua },
  { "bits_to_cvec", tk_ivec_bits_to_cvec_lua },
  { "bits_extend", tk_ivec_bits_extend_lua },
  { "bits_rearrange", tk_ivec_bits_rearrange_lua },
  { "set_stats", tk_ivec_set_stats_lua },
  { "set_similarity", tk_ivec_set_similarity_lua },
  { "set_jaccard", tk_ivec_set_jaccard_lua },
  { "set_overlap", tk_ivec_set_overlap_lua },
  { "set_dice", tk_ivec_set_dice_lua },
  { "set_tversky", tk_ivec_set_tversky_lua },
  { "set_weights_by_rank", tk_ivec_set_weights_by_rank_lua },
  { "set_similarity_by_rank", tk_ivec_set_similarity_by_rank_lua },
  { "set_find", tk_ivec_set_find_lua },
  { "set_intersect", tk_ivec_set_intersect_lua },
  { "set_union", tk_ivec_set_union_lua },
  { "lookup", tk_ivec_lookup_lua },
  { NULL, NULL }
};

static luaL_Reg tk_ivec_lua_ext_fns[] =
{
  { "bits_from_cvec", tk_ivec_bits_from_cvec_lua },
  { NULL, NULL }
};

int luaopen_santoku_ivec (lua_State *L)
{
  lua_newtable(L);
  luaL_register(L, NULL, tk_ivec_lua_fns);
  luaL_register(L, NULL, tk_ivec_lua_ext_fns);
  tk_ivec_create(L, 0, 0, 0);
  luaL_getmetafield(L, -1, "__index");
  luaL_register(L, NULL, tk_ivec_lua_mt_fns);
  luaL_register(L, NULL, tk_ivec_lua_mt_ext_fns);
  luaL_register(L, NULL, tk_ivec_lua_mt_ext2_fns);
  lua_pop(L, 2);
  return 1;
}
