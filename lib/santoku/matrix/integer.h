#include <santoku/matrix/integer.conf.h>
#include <santoku/matrix/base.h>

KHASH_INIT(i64, int64_t, char, 0, kh_int64_hash_func, kh_int64_hash_equal)
typedef khash_t(i64) i64_hash_t;

typedef struct { double score; uint64_t v; } tk_matrix_ranked_pair_t;
#define tk_matrix_ranked_pair_gt(a, b) ((a).score > (b).score)
KSORT_INIT(ranked_desc, tk_matrix_ranked_pair_t, tk_matrix_ranked_pair_gt)
KSORT_INIT(u64, uint64_t, ks_lt_generic)

static inline void tk_matrix_push_selected_matrix (
  lua_State *L,
  i64_hash_t *selected
) {

  int64_t *top_v = tk_malloc(L, kh_size(selected) * sizeof(int64_t));
  uint64_t write = 0;
  for (khint_t khi = kh_begin(selected); khi < kh_end(selected); khi ++) {
    if (!kh_exist(selected, khi))
      continue;
    top_v[write ++] = kh_key(selected, khi);
  }
  _tk_matrix_create(L, 1, write, top_v, NULL);
  kh_destroy(i64, selected);
}

static inline tk_matrix_ranked_pair_t *tk_matrix_rankings (
  lua_State *L,
  double *scores,
  uint64_t n_visible,
  uint64_t n_hidden
) {

  tk_matrix_ranked_pair_t *rankings = tk_malloc(L, n_hidden * n_visible * sizeof(tk_matrix_ranked_pair_t));
  for (uint64_t h = 0; h < n_hidden; h ++)
    for (uint64_t v = 0; v < n_visible; v ++)
      rankings[h * n_visible + v] = (tk_matrix_ranked_pair_t) { scores[h * n_visible + v], v };


  for (uint64_t j = 0; j < n_hidden; j ++)
    ks_introsort(ranked_desc, n_visible, rankings + j * n_visible);


  return rankings;
}

static inline void tk_matrix_select_union (
  lua_State *L,
  i64_hash_t *selected,
  tk_matrix_ranked_pair_t *rankings,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k,
  uint64_t trunc
) {

  int kha;
  uint64_t *offsets = tk_malloc(L, n_hidden * sizeof(uint64_t));
  memset(offsets, 0, n_hidden * sizeof(uint64_t));
  while (kh_size(selected) < top_k) {
    bool advanced = false;
    for (uint64_t j = 0; j < n_hidden; j++) {
      tk_matrix_ranked_pair_t *rankings_h = rankings + j * n_visible;
      while (offsets[j] < n_visible) {
        tk_matrix_ranked_pair_t candidate = rankings_h[offsets[j]];
        offsets[j]++;
        if (candidate.v >= (n_visible - trunc))
          continue;
        kh_put(i64, selected, (int64_t) candidate.v, &kha);
        advanced = true;
        break;
      }
      if (kh_size(selected) >= top_k)
        break;
    }
    if (!advanced)
      break;
  }

  free(offsets);
}

static void tk_matrix_top_generic (
  lua_State *L,
  double *scores,
  uint64_t n_visible,
  uint64_t n_hidden,
  uint64_t top_k,
  uint64_t trunc
) {

  i64_hash_t *selected = kh_init(i64);
  tk_matrix_ranked_pair_t *rankings = tk_matrix_rankings(L, scores, n_visible, n_hidden);
  tk_matrix_select_union(L, selected, rankings, n_visible, n_hidden, top_k, trunc);
  free(rankings);

  tk_matrix_push_selected_matrix(L, selected);
}
