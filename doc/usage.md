# Using santoku-matrix

Worked examples for each data type. See the [README](../README.md) for orientation and the
type map. **The tests are the spec**: every section names the test that covers the rest.

Recurring idioms:
- **wrap = zero-copy.** `create({ data/offsets = vec })` shares storage with `vec`.
- **in-place ops return self.** chain them.
- **fit returns params, apply takes them.** `local w = X:standardize(); Y:standardize(w)`.
- **`out=` reuse.** pass a destination to avoid allocation in hot loops.
- **everything persists.** `obj:persist(path)` / `type.load(path)`.

---

## csr: sparse rows  ·  `test/spec/santoku/csr.lua`

Build by wrapping parallel vecs, or row-by-row:

```lua
local csr, ivec, fvec = require("santoku.csr"), require("santoku.ivec"), require("santoku.fvec")

local X = csr.create({                       -- wrap (zero-copy)
  offsets   = ivec.create({ 0, 2, 3 }),      -- n_rows+1
  neighbors = ivec.create({ 0, 1, 1 }),      -- column ids (svec for i32 ids)
  values    = fvec.create({ 1, 2, 3 }),      -- omit -> binary csr
  n_cols    = 2,
})

local B = csr.create({ n_cols = 4, values = "f32" })   -- builder
B:push(0, 1.5):push(2, 2.5):row()                      -- :push(col[, val]); :row() closes a row
```

Weight / shape features (the fit/apply transforms; binary csr auto-materializes values):

```lua
local w = X:idf()            -- or X:bns(labels_csr), X:standardize()  -- fit -> weights
Y:idf(w)                     -- apply the same weights to held-out data
X:normalize()                -- L2 per row, returns self
X:hcat(other)                -- concat feature blocks in place (shifts other's column ids)
local ss = X:sumsq_cols(ivec.create({ 0, b1, n_cols }))   -- per-block sums (for block scaling)
```

Construct from labels / select / bridge to dense or bits:

```lua
local Y = csr.from_classes(ivec.create({ 2, 0, 1 }))   -- one label per row
local M = X:to_dense()                                  -- -> mtx ;  mtx:to_sparse() is the inverse
local bits = X:to_bits()                                -- packed bitmap ;  csr.from_bits(bits, n, c)
```

## mtx: dense matrix  ·  `test/spec/santoku/mtx.lua`

```lua
local mtx, dvec = require("santoku.mtx"), require("santoku.dvec")

local M = mtx.create({ n_rows = 2, n_cols = 3, type = "f64" })   -- allocate, zeroed
local W = mtx.create({ data = dvec.create({1,2,3,4,5,6}), n_rows = 2, n_cols = 3 })  -- wrap

W:get(1, 2)                       -- 6
W:sums("row")  W:maxargs("col")   -- axis reductions -> a vec   (axis = "row" | "col")
W:transpose()  W:rows(ids)        -- reshapes / gathers (rows/cols/row, optional out=)
W:multiply(other)                 -- matmul ;  W:multiplyv(v[, transpose]) for mat-vec
W:center(); W:standardize(); W:normalize("row")   -- fit/apply, same idiom as csr
```

Retrieval and binary layouts:

```lua
local P = corpus:topk(queries, k)     -- brute-force top-k -> P csr (neighbors=ids, values=scores)
local bits = M:sign()                 -- f-matrix -> sign bitmap (bits-tagged mtx)
-- a bits mtx supports :popcount / :hamming / :band / :transpose
```

## spans: per-document labelled intervals  ·  `test/spec/santoku/spans.lua`

Offsets (per-doc) + named integer columns, addressed by name:

```lua
local spans = require("santoku.spans")

local S = spans.create({ "s", "e", "ty" })     -- declare columns
S:push(0, 3, 1):push(4, 7, 2):doc()            -- :push(values...) per record; :doc() closes a doc
S:doc()                                         -- empty doc
S:push(2, 5, 1):doc()

-- or wrap (zero-copy): any key but "offsets" is a column
local W = spans.create({ offsets = off, s = starts, e = ends, ty = types })

S:n()  S:n_docs()  S:offsets()  S:col("s")     -- col() returns the named ivec by reference
```

Span algorithms:

```lua
local C    = S:enumerate_subspans(max_len, outer_ty)   -- candidate spans within non-outer runs
local keep, cls = S:nms_dp(labels, scores, reject, nl) -- weighted interval scheduling (non-overlap)
local U    = A:union(B, gold)                           -- merge + label against gold
local ids, surfs = S:surfaces(texts, lowercase)        -- intern surface strings
```

## vectors: `ivec dvec fvec svec cvec rvec pvec`  ·  `test/spec/santoku/matrix.lua`

```lua
local ivec = require("santoku.ivec")
local v = ivec.create({ 5, 1, 3 })   -- or create(n) for n zeros
v:push(9); v:set(0, 7); v:get(0)
v:asc()                              -- sort (asc/desc/uasc/kasc/...)
v:dot(other); v:sum()               -- math on ivec/dvec/fvec
```

- `cvec` = raw byte vector (blobs, and storage behind bit-matrices; bitmap *ops* are on `mtx`).
- `rvec` {i64,f64} / `pvec` {i64,i64} = pair vectors with top-k heaps (`hmin`/`hmax`).
- `ivec` carries set similarity (`set_jaccard`/`set_intersect`/`set_union`/...).

## mmap: disk-backed vectors  ·  `test/spec/santoku/mmap.lua`

For arrays too big for RAM, or to persist large results:

```lua
local v = fvec.mmap_create(path, n)   -- zeroed, on disk
v:set(0, 1.0); v:mmap_sync()          -- flush to disk
local w = fvec.mmap_open(path)        -- reopen later
```

An mmap vec can back an `mtx` (`mtx.create({ data = v, ... })`) so a big encode writes straight
to disk via the `out=` parameter.

## Gotchas

- **Aliasing:** because `create` wraps, two objects can share one vec. Mutating one is visible in
  the other (intentional for `out=`/mmap, surprising if unexpected). Copy the vec if you need
  independence.
- **In-place vs new:** `hcat`/`normalize`/`standardize`/`center`/`scale_cols`/`filter`/`sort`/
  `append` mutate and return self; `rows`/`select`/`transpose`/`to_dense`/`topk`/`docs`/`union`/
  `enumerate_subspans` return new objects.
- **Tags must match** for binary ops: `mtx:topk`/`eq` and `csr:hcat` require matching element (and
  neighbor) tags.
