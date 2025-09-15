# Now

- Handle sparse labels/codes in bits scoring functions (e.g. accept an ivec of
  labels that are in the list-of-set-bits representation instead of just being a
  single label per sample; this will require extensive changes)
- Heap for top_chi2 and top_mi instead of keeping full list
- Paralellize counting in top_mi and top_chi2

- tk_cvec_t
    - Implement cvec to/from string

- Templatized hash and btree (same peek, pop, create w/out L, register, etc)
- Ensure all APIs exposed (heap, iuset/map extensions)
- Standard peek/test: tk_xvec_check, tk_xvec_test, tk_xvec_fcheck, tk_xvec_ftest

# Consider

- Additional default containers?

    - tk_coset_t: kbtree_t(tk_cvec_t)
    - tk_cumap_t: khash_t(int64_t, tk_cvec_t)
    - tk_comap_t: kbtree_t(int64_t, tk_cvec_t)

    - tk_zomap_t: kbtree_t(tk_cvec_t, int64_t)
    - tk_vumap_t: khash_t(tk_cvec_t, tk_cvec_t)
    - tk_vomap_t: kbtree_t(tk_cvec_t, tk_cvec_t)

    - tk_dumap_t: khash_t(int64_t, double)
    - tk_domap_t: kbtree_t(int64_t, double)
    - tk_rumap_t: khash_t(int64_t, tk_rank_t)
    - tk_romap_t: kbtree_t(int64_t, tk_rank_t)

    - tk_pumap_t: khash_t(int64_t, tk_pair_t)
    - tk_pomap_t: kbtree_t(int64_t, tk_pair_t)
