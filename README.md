# santoku-matrix

The data layer for santoku: typed numeric containers for ML and retrieval work —
dense matrices (`mtx`), sparse matrices (`csr`), labelled span sets (`spans`), the
typed vectors they're built on, plus disk-backed (mmap) vectors and C-level hash
maps/sets.

This README is a usage guide, not an API reference. **The tests are the spec** —
each section points at the test that exercises the full surface. Read those for the
exhaustive method list; read this for the shape of things and the ergonomics.

## Conventions (shared across the object types)

These hold for `mtx` / `csr` / `spans` unless noted:

- **Zero-copy wrap.** `create({ data = vec, ... })` / `create({ offsets = vec, ... })`
  *adopts* the vecs you pass — storage is shared, not copied. Mutating the vec mutates
  the object (and vice versa).
- **In-place ops return self**, so they chain: `X:normalize()`, `A:hcat(B)`.
- **fit / apply.** Feature transforms return the parameters they learned; pass those back
  to apply the same transform to new data:
  ```lua
  local w = X:standardize()   -- fit: learns + applies, returns weights
  Y:standardize(w)            -- apply: uses the given weights
  ```
  Same pattern for `csr:bns` / `csr:idf` and `mtx:center` / `mtx:standardize`.
- **`out=` reuse.** Gather/transform ops take an optional destination object to avoid
  allocation: `M:rows(ids, out)`, `M:multiplyv(v, false, y)`, `M:to_sparse(nil, X)`.
- **persist / load.** Every object round-trips to a binary file: `X:persist(path)` /
  `csr.load(path)` (likewise `mtx.load`, `spans.load`).
- **Type tags.** `mtx`/`csr` carry an element tag (`i32`/`i64`/`f32`/`f64`/`bits`);
  `csr` additionally tags its column ids as `i32` (svec) or `i64` (ivec). Query with
  `:type()`.

---

## `mtx` — dense matrix  ·  `test/spec/santoku/mtx.lua`

Row-major dense matrix over a typed vector. Allocate-and-zero, or wrap an existing vec.

```lua
local mtx = require("santoku.mtx")
local dvec = require("santoku.dvec")

local M = mtx.create({ data = dvec.create({ 1, 2, 3, 4, 5, 6 }), n_rows = 2, n_cols = 3 })
M:shape()            -- 2, 3
M:get(1, 2)          -- 6
M:sums("row")        -- dvec { 6, 15 }   (axis = "row" | "col")
M:normalize("row")   -- in place, returns self

-- retrieval: brute-force top-k of corpus rows per query row -> P csr
local P = corpus:topk(queries, k)   -- P:neighbors() = corpus ids, P:values() = scores
```

Covers: axis reductions (`sums`/`maxs`/`mins`/`maxargs`/`mags`/`argsort`), `transpose`,
`rows`/`cols`/`row` gather, `hcat`, `center`/`standardize`/`normalize`, `multiply`/
`multiplyv`, `sign`/`median` (→ bit matrices), the `bits` layout (`popcount`/`hamming`/
`band`/...), `topk`, `from_pairs`, and `to_sparse` (↔ `csr`).

## `csr` — sparse matrix  ·  `test/spec/santoku/csr.lua`

Rows of `(column, value)` pairs: `offsets` (length `n_rows+1`), `neighbors` (column ids),
and optional `values` (binary — value 1 — when omitted). This is what tokenizers and
label sets produce.

```lua
local csr = require("santoku.csr")
local ivec, fvec = require("santoku.ivec"), require("santoku.fvec")

-- wrap parts (zero-copy)
local X = csr.create({
  offsets   = ivec.create({ 0, 2, 3, 5 }),
  neighbors = ivec.create({ 0, 2, 1, 0, 3 }),   -- i64 ids; pass an svec for i32 ids
  values    = fvec.create({ 1, 1, 1, 1, 1 }),   -- omit for a binary csr
  n_cols    = 4,
})

-- or build row by row
local B = csr.create({ n_cols = 4, values = "f32" })
B:push(0, 1.5):push(2, 2.5):row()               -- :push(col[, val]) then :row()

X:idf()                  -- feature weighting (fit; also :bns(labels), :standardize)
X:normalize()            -- L2 per row (materializes values if binary)
X:hcat(other)            -- concat feature blocks in place (shifts other's columns)
```

Also: `from_classes` / `from_mask` / `from_bits` constructors; `shape`/`nnz`/`type`/
`offsets`/`neighbors`/`values` accessors; `rows`/`select`/`transpose`/`scale_cols`/
`sumsq_cols`; and the dense/bitmap bridges `to_dense` (→ `mtx`) and `to_bits`.

## `spans` — per-document labelled intervals  ·  `test/spec/santoku/spans.lua`

Per-document sets of integer-keyed records: `offsets` (per-doc, like csr) plus one or more
**named integer columns** (e.g. `s`, `e`, `ty`). Columns are addressed by name via `:col`.

```lua
local spans = require("santoku.spans")

-- build: declare columns, push records, close each doc
local S = spans.create({ "s", "e", "ty" })
S:push(0, 3, 1):push(4, 7, 2):doc()
S:doc()                                  -- empty doc
S:push(2, 5, 1):doc()

-- or wrap existing vecs (zero-copy); any key but "offsets" is a column
local W = spans.create({ offsets = off, s = starts, e = ends, ty = types })

S:n()            -- total spans ;  S:n_docs() ;  S:offsets()
S:col("s")       -- the named column (ivec, by reference)
```

Covers: `filter`/`docs`/`append`/`sort`/`eq`/`surfaces`, and the span algorithms
`enumerate_subspans`, `nms_dp` (weighted interval scheduling), and `union` (with gold
labelling). (santoku-learn extends this metatable with NER helpers like `type_labels`.)

## Vectors — `ivec` `dvec` `fvec` `svec` `cvec` `rvec` `pvec`  ·  `test/spec/santoku/matrix.lua`

Typed, growable arrays — the storage under everything above. Element types: `i64` (ivec),
`f64` (dvec), `f32` (fvec), `i32` (svec), `u8`/bitmap (cvec), and the pair vectors
`rvec` `{i64,f64}` / `pvec` `{i64,i64}`.

```lua
local ivec = require("santoku.ivec")
local v = ivec.create({ 5, 1, 3 })
v:push(9); v:asc()            -- get/set/push/size, sort (asc/desc/uasc/kasc/...)
v:dot(other); v:sum()         -- math on ivec/dvec/fvec
```

Beyond the core: `cvec` is the raw byte vector (blob/string storage, and the backing store for
bit-matrices — the bitmap *operations* live on `mtx`, above); `rvec`/`pvec` provide top-k heaps
(`hmin`/`hmax`); `ivec` carries set-similarity (`set_jaccard`/`set_intersect`/`set_union`/...).

## Disk-backed vectors (mmap)  ·  `test/spec/santoku/mmap.lua`

Any numeric vec can be mmap-backed for bounded-RAM / persistent large arrays:

```lua
local v = fvec.mmap_create(path, n)   -- zeroed, on disk
v:set(0, 1.0); v:mmap_sync()          -- flush
local w = fvec.mmap_open(path)        -- reopen
```

Useful as `out=` storage for big encodes (e.g. writing codes straight to disk).

## Hash maps & sets (C API)

`iuset`/`iumap`/`zumap`/`duset`/`cuset`/... are header-only templated hash containers
(`<santoku/iuset.h>` etc.) used from C extensions — not Lua-facing. See the headers for
the `tk_*_create`/`put`/`get`/`foreach` surface.

## License

MIT License

Copyright 2025 Birch Point SWE

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
