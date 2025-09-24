# Now

- tk_cvec_t
    - Implement cvec to/from string

- tk_cumap/zumap/etc:
    - Strictly use cvec with proper referencing
        - #define tk_umap_destroy_key/value tk_lua_del_ephemeron
        - #define tk_umap_link_key/value tk_lua_link_ephemeron (like add but with pointer)
    - Support those macros with vector as well

# Next

- Templatized btree

- tk_cvec_t/tk_ivec_t
    - bits_extend: don't create temporary same-format buffer, copy direct

# Consider

- Should cumap/zumap use cvec instead of strings?

- Replace iumap/dumap with puset and ruset, like euset

- Standard peek/test: tk_xvec_check, tk_xvec_test, tk_xvec_fcheck, tk_xvec_ftest

- Replace calls to klib with inlined implementations of vec, sort, hash, and
  tree functions

- Additional default containers?
    - tk_coset_t: kbtree_t(tk_cvec_t)
    - tk_comap_t: kbtree_t(int64_t, tk_cvec_t)
    - tk_zomap_t: kbtree_t(tk_cvec_t, int64_t)
    - tk_vumap_t: khash_t(tk_cvec_t, tk_cvec_t)
    - tk_vomap_t: kbtree_t(tk_cvec_t, tk_cvec_t)
    - tk_dumap_t: khash_t(int64_t, double)
    - tk_domap_t: kbtree_t(int64_t, double)
    - tk_rumap_t: khash_t(int64_t, tk_rank_t)
    - tk_romap_t: kbtree_t(int64_t, tk_rank_t)
    - tk_pomap_t: kbtree_t(int64_t, tk_pair_t)
