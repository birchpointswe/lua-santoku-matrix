# Now

- Templatized hash and btree
- Ensure all APIs exposed (heap, iuset/map extensions)
- Standard peek/test: tk_xvec_check, tk_xvec_test, tk_xvec_fcheck, tk_xvec_ftest
- templatize ordered/unordered map/set under a single map.template.h
- Update iuset/ioset/iumap/iomap to have a lua API and to use lua state for cleanup
- Heap for top_chi2 and top_mi instead of keeping full list

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
