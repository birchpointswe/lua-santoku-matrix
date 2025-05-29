# Now

- Integrate into tsetlin, cover all base datastructures used
- Integrate into base toku, ensure functional with emscripten

- cvec: tests, to/from string
- rvec: lua api, especially for heap fns
- ivec: bitmap functions
- dvec: openblas
- dvec: glove clustering test case

- overriable tk_vec_abs/exp so that we can use fabs/labs depending

# Next

- Use santoku-threads with shared/singleton threadpool where openblas is not
  used
    - With HAS_OPENBLAS, openblas is used to accelerate tk_dvec_t operations
    - With HAS_PTHREAD, pthreads are used where openblas is not

- Can we replace dense codes with sparse set-bits format in tk_ivec_t bitmap
  functions?

- Heap for top_chi2 and top_mi instead of keeping full list

# Consider

- Standardized containers?
    - santoku.map/set.u/o
        - tk_cuset_t: khash_t(tk_cvec_t)
        - tk_coset_t: kbtree_t(tk_cvec_t)
        - tk_iuset_t: khash_t(int64_t)
        - tk_ioset_t: kbtree_t(int64_t)
        - tk_cumap_t: khash_t(int64_t, tk_cvec_t)
        - tk_comap_t: kbtree_t(int64_t, tk_cvec_t)
        - tk_iumap_t: khash_t(int64_t, int64_t)
        - tk_iomap_t: kbtree_t(int64_t, int64_t)
        - tk_dumap_t: khash_t(int64_t, double)
        - tk_domap_t: kbtree_t(int64_t, double)
        - tk_rumap_t: khash_t(int64_t, tk_rank_t)
        - tk_romap_t: kbtree_t(int64_t, tk_rank_t)
        - tk_zumap_t: khash_t(tk_cvec_t, int64_t)
        - tk_zomap_t: kbtree_t(tk_cvec_t, int64_t)
        - tk_vumap_t: khash_t(tk_cvec_t, tk_cvec_t)
        - tk_vomap_t: kbtree_t(tk_cvec_t, tk_cvec_t)
        - wrap all klib fns
        - set operations
