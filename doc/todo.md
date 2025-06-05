# Now

- Standard peek/test: tk_xvec_check, tk_xvec_test, tk_xvec_fcheck, tk_xvec_ftest

- usage of lua state in the C api should be limited strictly to create, and
  there should be no destroy function (users must user the lua api)

- templatize ordered/unordered map/set under a single map.template.h

- Update iuset/ioset/iumap/iomap to have a lua API and to use lua state for cleanup

- x_create functions should accept optional int pointer to fill with lua stack
  position (for ref)

- Integrate into toku base or separate library

# Later

- Use these structures across all libraries where possibe

- cvec: tests, to/from string
- rvec: lua api, especially for heap fns
- ivec: bitmap functions
- dvec: glove clustering test case

- Pull ephemeron logic from toku web/python (attach arbitrary lua values weakly
  to other lua values) for linking garbage collection of xvec and xmap/set
  objects to other lua values

# Next

- Intelligent handling of pthread availability. When threads == 0 or 1,
  completely avoid threads (don't just use one spawned thread).

- Can we replace dense codes with sparse set-bits format in tk_ivec_t bitmap
  functions?

- Heap for top_chi2 and top_mi instead of keeping full list

# Consider

- Additional containers?

    - tk_cuset_t: khash_t(tk_cvec_t)
    - tk_coset_t: kbtree_t(tk_cvec_t)

    - tk_cumap_t: khash_t(int64_t, tk_cvec_t)
    - tk_comap_t: kbtree_t(int64_t, tk_cvec_t)
    - tk_zumap_t: khash_t(tk_cvec_t, int64_t)
    - tk_zomap_t: kbtree_t(tk_cvec_t, int64_t)
    - tk_vumap_t: khash_t(tk_cvec_t, tk_cvec_t)
    - tk_vomap_t: kbtree_t(tk_cvec_t, tk_cvec_t)

    - tk_dumap_t: khash_t(int64_t, double)
    - tk_domap_t: kbtree_t(int64_t, double)
    - tk_rumap_t: khash_t(int64_t, tk_rank_t)
    - tk_romap_t: kbtree_t(int64_t, tk_rank_t)
