# Santoku Matrix

Templated array-like containers (vectors), hash maps, and sets with type-specific implementations.

## Common Vector Operations

All vector types (ivec, dvec, cvec, rvec, pvec) share these core operations unless otherwise noted in their specific sections.

### Constructor Functions

| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `create` | `[size_or_table]` | `vector` | Creates new vector from size or Lua table |
| `load` | `filename, [is_string]` | `vector` | Loads vector from file or string |
| `from_raw` | `raw_data, size` | `vector` | Creates vector from raw binary data |
| `destroy` | `vector` | `-` | Destroys vector and frees memory |

### Core Methods

| Method | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `copy` | `dest, source, [start], [end], [dest_idx]` | `-` | Copies elements between vectors |
| `persist` | `[filename]` | `[string]` | Saves vector to file or returns binary string |
| `size` | `-` | `integer` | Returns number of elements |
| `capacity` | `-` | `integer` | Returns allocated capacity |
| `resize` | `size` | `-` | Resizes vector to exact size |
| `setn` | `n` | `-` | Sets number of active elements |
| `ensure` | `size` | `-` | Ensures minimum capacity |
| `shrink` | `-` | `-` | Shrinks capacity to match size |
| `clear` | `-` | `-` | Sets size to 0 (keeps capacity) |
| `zero` | `-` | `-` | Sets all elements to 0 |
| `transpose` | `dest, source, cols` | `-` | Transposes matrix representation |

### Element Access

| Method | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `get` | `idx` | `value` | Gets element at index |
| `set` | `idx, value` | `-` | Sets element at index |
| `push` | `value` | `-` | Appends value to end |
| `insert` | `idx, value` | `-` | Inserts value at index |

### Sorting Operations

| Method | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `shuffle` | `-` | `-` | Randomly shuffles elements |
| `asc` | `[start], [end]` | `-` | Sorts in ascending order |
| `desc` | `[start], [end]` | `-` | Sorts in descending order |
| `uasc` | `[start], [end]` | `new_end` | Sorts ascending and removes duplicates |
| `udesc` | `[start], [end]` | `new_end` | Sorts descending and removes duplicates |
| `xasc` | `[start], [end]` | `new_end` | Binary sorts ascending with unique |
| `xdesc` | `[start], [end]` | `new_end` | Binary sorts descending with unique |
| `kasc` | `k, [start], [end]` | `-` | Partial sort to find k smallest |
| `kdesc` | `k, [start], [end]` | `-` | Partial sort to find k largest |

### Mathematical Operations
*Available for ivec and dvec only*

| Method | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `add` | `scalar, [start], [end]` | `-` | Adds scalar to elements |
| `scale` | `scalar, [start], [end]` | `-` | Multiplies elements by scalar |
| `addv` | `vector, [start], [end]` | `-` | Element-wise addition |
| `scalev` | `vector, [start], [end]` | `-` | Element-wise multiplication |
| `abs` | `[start], [end]` | `-` | Takes absolute value of elements |
| `exp` | `[start], [end]` | `-` | Applies exponential function |
| `log` | `[start], [end]` | `-` | Applies natural logarithm |
| `pow` | `exponent, [start], [end]` | `-` | Raises elements to power |
| `dot` | `vector_b` | `number` | Computes dot product |
| `multiply` | `a, b, c, k, [transpose_a], [transpose_b]` | `-` | Matrix multiplication |
| `magnitude` | `-` | `number` | Computes Euclidean magnitude |
| `sum` | `-` | `number` | Sums all elements |
| `csums` | `cols` | `vector` | Column-wise sums |
| `rsums` | `cols` | `vector` | Row-wise sums |
| `min` | `-` | `value, index` | Finds minimum value and index |
| `max` | `-` | `value, index` | Finds maximum value and index |
| `cmins` | `cols` | `vector` | Column-wise minimums |
| `rmins` | `cols` | `vector` | Row-wise minimums |
| `cmaxs` | `cols` | `vector` | Column-wise maximums |
| `rmaxs` | `cols` | `vector` | Row-wise maximums |

### Iteration
*Not available for cvec*

| Method | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `each` | `-` | `iterator` | Returns iterator over values |
| `ieach` | `-` | `iterator` | Returns iterator over index-value pairs |

### Conversion and Utility
*Not available for cvec*

| Method | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `table` | `[start], [end]` | `table` | Converts to Lua table |
| `rtable` | `cols, [start], [end]` | `table` | Converts to row-major table of tables |
| `ctable` | `cols, [start], [end]` | `table` | Converts to column-major table of tables |
| `fill` | `value` | `-` | Fills all elements with value |
| `fill_indices` | `-` | `-` | Fills with sequential indices |
| `raw` | `[format]` | `string` | Returns raw binary representation |

## Vector Type Instances

### `santoku.ivec`
Integer vector module providing dynamic arrays of 64-bit integers.

**Inherits:** All common vector operations including mathematical operations
**Base type:** `int64_t`

#### Module-level Functions (ivec)

| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `bits_from_cvec` | `bitmap, samples, features` | `ivec` | Creates ivec from packed bitmap (cvec) |

#### Additional Operations (ivec-specific)

##### Set Operations

| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `set_jaccard` | `vector_b, [weights]` | `number` | Computes Jaccard similarity |
| `set_overlap` | `vector_b, [weights]` | `number` | Computes overlap coefficient |
| `set_dice` | `vector_b, [weights]` | `number` | Computes Dice coefficient |
| `set_tversky` | `vector_b, alpha, beta, [weights]` | `number` | Computes Tversky index |
| `set_similarity` | `vector_b, type, [alpha], [beta], [weights]` | `number` | Computes specified similarity measure |
| `set_stats` | `vector_b, [weights]` | `inter, sum_a, sum_b` | Computes set statistics |
| `set_intersect` | `vector_b, [output]` | `ivec` | Computes set intersection (assumes sorted) |
| `set_union` | `vector_b, [output]` | `ivec` | Computes set union (assumes sorted) |
| `set_find` | `value` | `index, insert_pos` | Binary search, returns index or -1 and insertion position |

##### Feature Selection

| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `bits_score_chi2` | `set_bits, codes/labels, samples, visible, hidden, [threads]` | `dvec` | Chi-squared feature scores. `codes` can be string, light userdata, or cvec; `labels` can be ivec |
| `bits_score_mi` | `set_bits, codes/labels, samples, visible, hidden, [threads]` | `dvec` | Mutual information scores. `codes` can be string, light userdata, or cvec; `labels` can be ivec |
| `bits_score_entropy` | `codes, samples, hidden, [threads]` | `dvec` | Entropy scores. `codes` can be string, light userdata, or cvec |
| `bits_top_chi2` | `set_bits, codes/labels, samples, visible, hidden, k, [threads]` | `ivec` | Top k features by chi-squared. `codes` can be string, light userdata, or cvec; `labels` can be ivec |
| `bits_top_mi` | `set_bits, codes/labels, samples, visible, hidden, k, [threads]` | `ivec` | Top k features by mutual information. `codes` can be string, light userdata, or cvec; `labels` can be ivec |
| `bits_top_entropy` | `codes, samples, hidden, k, [threads]` | `ivec` | Top k features by entropy. `codes` can be string, light userdata, or cvec |

##### Bit Operations

| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `bits_flip_interleave` | `samples, features` | `-` | Flips and interleaves bit patterns |
| `bits_rearrange` | `ids, features` | `-` | Rearranges bits by feature IDs |
| `bits_extend` | `base_bits, ext_bits, features, extended` | `-` | Extends bit representation |
| `bits_to_cvec` | `set_bits, samples, features, [flip_interleave]` | `cvec` | Converts sparse bit indices to packed bitmap |
| `bits_filter` | `set_bits, top_v, visible` | `-` | Filters elements by selected features |

##### Matrix Sorting Operations

| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `rasc` | `cols` | `ivec` | Row-wise ascending sort indices |
| `rdesc` | `cols` | `ivec` | Row-wise descending sort indices |
| `casc` | `cols` | `ivec` | Column-wise ascending sort indices |
| `cdesc` | `cols` | `ivec` | Column-wise descending sort indices |
| `rmaxargs` | `cols` | `ivec` | Row-wise maximum argument indices |
| `cmaxargs` | `cols` | `ivec` | Column-wise maximum argument indices |
| `rminargs` | `cols` | `ivec` | Row-wise minimum argument indices |
| `cminargs` | `cols` | `ivec` | Column-wise minimum argument indices |
| `rmagnitudes` | `cols` | `dvec` | Row-wise magnitude calculation |
| `cmagnitudes` | `cols` | `dvec` | Column-wise magnitude calculation |

##### Cross-type Operations

| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `copy_pkeys` | `pvec, [start], [end], [dest]` | `-` | Copies keys from pair vector |
| `copy_pvalues` | `pvec, [start], [end], [dest]` | `-` | Copies values from pair vector |
| `copy_rkeys` | `rvec, [start], [end], [dest]` | `-` | Copies keys from rank vector |
| `copy_rvalues` | `rvec, [start], [end], [dest]` | `-` | Copies values from rank vector |

### `santoku.dvec`
Double-precision floating-point vector module.

**Inherits:** All common vector operations including mathematical operations
**Base type:** `double`
**Excluded:** Bit operations, cross-type operations, set operations, feature selection

#### Additional Operations (dvec-specific)

##### Matrix Sorting Operations

| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `rasc` | `cols` | `ivec` | Row-wise ascending sort indices |
| `rdesc` | `cols` | `ivec` | Row-wise descending sort indices |
| `casc` | `cols` | `ivec` | Column-wise ascending sort indices |
| `cdesc` | `cols` | `ivec` | Column-wise descending sort indices |
| `rmaxargs` | `cols` | `ivec` | Row-wise maximum argument indices |
| `cmaxargs` | `cols` | `ivec` | Column-wise maximum argument indices |
| `rminargs` | `cols` | `ivec` | Row-wise minimum argument indices |
| `cminargs` | `cols` | `ivec` | Column-wise minimum argument indices |
| `rmagnitudes` | `cols` | `dvec` | Row-wise magnitude calculation |
| `cmagnitudes` | `cols` | `dvec` | Column-wise magnitude calculation |

### `santoku.cvec`
Character vector module providing byte arrays and bitmap operations.

**Inherits:** Core vector operations and sorting operations only
**Base type:** `unsigned char`
**Excluded:** Mathematical operations, iteration methods, table conversion methods, fill operations

#### Module-level Functions (cvec)

| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `bits_from_ivec` | `set_bits, samples, features` | `cvec` | Creates packed bitmap from sparse bit indices |

#### Additional Operations (cvec-specific)

##### Bitmap Operations

| Method | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `bits_flip_interleave` | `samples, features` | `-` | Flips and interleaves bit patterns in place |
| `bits_to_ivec` | `samples, features` | `ivec` | Converts packed bitmap to sparse bit indices |
| `bits_rearrange` | `ids, features` | `-` | Rearranges bitmap samples by IDs |
| `bits_extend` | `ext, samples, base_features, ext_features` | `-` | Extends bitmap with additional features |
| `bits_filter` | `selected_features, samples, features` | `-` | Filters bitmap to selected features |
| `bits_popcount` | `n_bits` | `integer` | Counts set bits in bitmap |
| `bits_hamming` | `other, n_bits` | `integer` | Hamming distance between bitmaps |
| `bits_hamming_mask` | `other, mask, n_bits` | `integer` | Masked Hamming distance |
| `bits_and` | `out, other, n_bits` | `-` | Bitwise AND operation |
| `bits_or` | `out, other, n_bits` | `-` | Bitwise OR operation |
| `bits_xor` | `out, other, n_bits` | `-` | Bitwise XOR operation |

### `santoku.rvec`
Rank vector module for storing pairs of (integer, double).

**Inherits:** Core vector operations and sorting operations
**Base type:** `struct { int64_t i; double d; }`
**Excluded:** Mathematical operations

#### Additional Operations (rvec-specific)

| Method | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `get` | `idx` | `integer, double` | Gets rank pair at index |
| `set` | `idx, integer, double` | `-` | Sets rank pair at index |
| `push` | `integer, double` | `-` | Appends rank pair |
| `keys` | `[out_ivec]` | `ivec` | Extracts integer keys as ivec. If out_ivec provided, writes to it instead of creating new |
| `values` | `[out_dvec]` | `dvec` | Extracts double values as dvec. If out_dvec provided, writes to it instead of creating new |
| `each` | `-` | `iterator` | Returns iterator over pairs |
| `ieach` | `-` | `iterator` | Returns iterator over indexed pairs |

##### Heap Operations

| Method | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `hmax` | `integer, double, max_size` | `-` | Maintains top-k maximum heap |
| `hmax_init` | `-` | `-` | Initializes as max heap |
| `hmax_push` | `integer, double` | `-` | Pushes to max heap |
| `hmax_pop` | `-` | `integer, double` | Pops from max heap |
| `hmin` | `integer, double, max_size` | `-` | Maintains top-k minimum heap |
| `hmin_init` | `-` | `-` | Initializes as min heap |
| `hmin_push` | `integer, double` | `-` | Pushes to min heap |
| `hmin_pop` | `-` | `integer, double` | Pops from min heap |

### `santoku.pvec`
Pair vector module for storing pairs of (integer, integer).

**Inherits:** Core vector operations and sorting operations
**Base type:** `struct { int64_t i; int64_t p; }`
**Excluded:** Mathematical operations

#### Additional Operations (pvec-specific)

| Method | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `get` | `idx` | `integer, integer` | Gets pair at index |
| `set` | `idx, integer, integer` | `-` | Sets pair at index |
| `push` | `integer, integer` | `-` | Appends pair |
| `keys` | `[out_ivec]` | `ivec` | Extracts first integers as ivec. If out_ivec provided, writes to it instead of creating new |
| `values` | `[out_ivec]` | `ivec` | Extracts second integers as ivec. If out_ivec provided, writes to it instead of creating new |
| `each` | `-` | `iterator` | Returns iterator over pairs |
| `ieach` | `-` | `iterator` | Returns iterator over indexed pairs |

##### Heap Operations

| Method | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `hmax` | `integer, integer, max_size` | `-` | Maintains top-k maximum heap |
| `hmax_init` | `-` | `-` | Initializes as max heap |
| `hmax_push` | `integer, integer` | `-` | Pushes to max heap |
| `hmax_pop` | `-` | `integer, integer` | Pops from max heap |
| `hmin` | `integer, integer, max_size` | `-` | Maintains top-k minimum heap |
| `hmin_init` | `-` | `-` | Initializes as min heap |
| `hmin_push` | `integer, integer` | `-` | Pushes to min heap |
| `hmin_pop` | `-` | `integer, integer` | Pops from min heap |

## Template System

### `santoku/vec/tpl.h`
Template header for generating type-specific vector implementations.

#### Template Configuration Macros
| Macro | Description |
|-------|-------------|
| `tk_vec_name` | Prefix for generated function names |
| `tk_vec_base` | Base type for vector elements |
| `tk_vec_pushbase(L, val)` | Pushes base type to Lua stack |
| `tk_vec_peekbase(L, idx)` | Gets base type from Lua stack |
| `tk_vec_lt(a, b)` | Less-than comparison for sorting |
| `tk_vec_gt(a, b)` | Greater-than comparison for sorting |
| `tk_vec_ltx(a, b)` | Binary less-than comparison |
| `tk_vec_gtx(a, b)` | Binary greater-than comparison |
| `tk_vec_eqx(a, b)` | Binary equality comparison |
| `tk_vec_eq(a, b)` | Standard equality comparison |
| `tk_vec_abs(val)` | Absolute value function |
| `tk_vec_limited` | Define to limit available operations |
| `tk_vec_destroy_item(item)` | Custom destructor for elements |
| `tk_vec_mt` | Metatable name for Lua |

#### Generated Functions

##### Core Operations
| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `tk_vec_name##_create` | `L, size, data*, existing*` | `tk_vec_name##_t*` | Create vector |
| `tk_vec_name##_destroy` | `vec*` | `-` | Destroy vector |
| `tk_vec_name##_peek` | `L, index, name` | `tk_vec_name##_t*` | Get vector from stack |
| `tk_vec_name##_peekopt` | `L, index` | `tk_vec_name##_t*` | Optionally get vector |
| `tk_vec_name##_load` | `L` | `tk_vec_name##_t*` | Load vector from string |
| `tk_vec_name##_persist` | `vec*` | `string` | Serialize vector |
| `tk_vec_name##_raw` | `vec*` | `const char*` | Get raw data pointer |

##### Size Management
| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `tk_vec_name##_size` | `vec*` | `uint64_t` | Get element count |
| `tk_vec_name##_capacity` | `vec*` | `uint64_t` | Get capacity |
| `tk_vec_name##_resize` | `vec*, new_size` | `-` | Resize vector |
| `tk_vec_name##_setn` | `vec*, new_size` | `-` | Set element count |
| `tk_vec_name##_ensure` | `vec*, min_capacity` | `-` | Ensure capacity |
| `tk_vec_name##_shrink` | `vec*` | `-` | Shrink to fit |
| `tk_vec_name##_clear` | `vec*` | `-` | Clear elements |
| `tk_vec_name##_zero` | `vec*` | `-` | Zero all elements |

##### Element Access
| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `tk_vec_name##_get` | `vec*, index` | `tk_vec_base` | Get element |
| `tk_vec_name##_set` | `vec*, index, value` | `-` | Set element |
| `tk_vec_name##_push` | `vec*, value` | `-` | Append element |
| `tk_vec_name##_fill` | `vec*, value` | `-` | Fill with value |
| `tk_vec_name##_fill_indices` | `vec*` | `-` | Fill with indices |

##### Data Manipulation
| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `tk_vec_name##_copy` | `dest*, src*, start, end` | `-` | Copy elements |
| `tk_vec_name##_transpose` | `vec*, cols` | `-` | Transpose as matrix |
| `tk_vec_name##_shuffle` | `vec*` | `-` | Random shuffle |
| `tk_vec_name##_table` | `vec*, L` | `-` | Convert to Lua table |
| `tk_vec_name##_ctable` | `vec*, cols, L` | `-` | To column tables |
| `tk_vec_name##_rtable` | `vec*, cols, L` | `-` | To row tables |

##### Sorting Operations
| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `tk_vec_name##_asc` | `vec*, start, end` | `-` | Sort ascending |
| `tk_vec_name##_desc` | `vec*, start, end` | `-` | Sort descending |
| `tk_vec_name##_uasc` | `vec*, start, end` | `uint64_t` | Sort unique ascending |
| `tk_vec_name##_udesc` | `vec*, start, end` | `uint64_t` | Sort unique descending |
| `tk_vec_name##_xasc` | `vec*, start, end` | `uint64_t` | Binary sort ascending |
| `tk_vec_name##_xdesc` | `vec*, start, end` | `uint64_t` | Binary sort descending |
| `tk_vec_name##_kasc` | `vec*, k, start, end` | `-` | Partial sort ascending |
| `tk_vec_name##_kdesc` | `vec*, k, start, end` | `-` | Partial sort descending |

##### Mathematical Operations (unless tk_vec_limited defined)
| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `tk_vec_name##_magnitude` | `vec*` | `double` | Vector magnitude |
| `tk_vec_name##_sum` | `vec*` | `tk_vec_base` | Sum of elements |
| `tk_vec_name##_csums` | `vec*, cols` | `tk_vec_name##_t*` | Column sums |
| `tk_vec_name##_rsums` | `vec*, cols` | `tk_vec_name##_t*` | Row sums |
| `tk_vec_name##_max` | `vec*` | `tk_vec_base` | Maximum element |
| `tk_vec_name##_cmaxs` | `vec*, cols` | `tk_vec_name##_t*` | Column maxima |
| `tk_vec_name##_rmaxs` | `vec*, cols` | `tk_vec_name##_t*` | Row maxima |
| `tk_vec_name##_min` | `vec*` | `tk_vec_base` | Minimum element |
| `tk_vec_name##_cmins` | `vec*, cols` | `tk_vec_name##_t*` | Column minima |
| `tk_vec_name##_rmins` | `vec*, cols` | `tk_vec_name##_t*` | Row minima |
| `tk_vec_name##_scale` | `vec*, scalar` | `-` | Scale by scalar |
| `tk_vec_name##_add` | `vec*, scalar` | `-` | Add scalar |
| `tk_vec_name##_scalev` | `vec*, vec2*` | `-` | Element-wise scale |
| `tk_vec_name##_addv` | `vec*, vec2*` | `-` | Element-wise add |
| `tk_vec_name##_abs` | `vec*` | `-` | Absolute values |
| `tk_vec_name##_exp` | `vec*` | `-` | Exponential |
| `tk_vec_name##_log` | `vec*` | `-` | Natural logarithm |
| `tk_vec_name##_pow` | `vec*, exponent` | `-` | Power operation |
| `tk_vec_name##_dot` | `vec1*, vec2*` | `double` | Dot product |
| `tk_vec_name##_multiply` | `result*, m1*, m2*, r1, c1, c2` | `-` | Matrix multiply |

### `santoku/vec/ext/tpl.h`

Template header for matrix-specific operations on vectors (treating vector as row-major matrix).

#### Generated Matrix Operations (unless tk_vec_limited defined)

| Function | Arguments | Returns | Description |
|----------|-----------|---------|-------------|
| `tk_vec_name##_rasc` | `L, vec*, cols` | `tk_ivec_t*` | Row-wise ascending sort indices |
| `tk_vec_name##_rdesc` | `L, vec*, cols` | `tk_ivec_t*` | Row-wise descending sort indices |
| `tk_vec_name##_casc` | `L, vec*, cols` | `tk_ivec_t*` | Column-wise ascending sort indices |
| `tk_vec_name##_cdesc` | `L, vec*, cols` | `tk_ivec_t*` | Column-wise descending sort indices |
| `tk_vec_name##_rmagnitudes` | `L, vec*, cols` | `tk_dvec_t*` | Row-wise L2 norms |
| `tk_vec_name##_cmagnitudes` | `L, vec*, cols` | `tk_dvec_t*` | Column-wise L2 norms |
| `tk_vec_name##_rmaxargs` | `L, vec*, cols` | `tk_ivec_t*` | Row-wise argmax indices |
| `tk_vec_name##_cmaxargs` | `L, vec*, cols` | `tk_ivec_t*` | Column-wise argmax indices |
| `tk_vec_name##_rminargs` | `L, vec*, cols` | `tk_ivec_t*` | Row-wise argmin indices |
| `tk_vec_name##_cminargs` | `L, vec*, cols` | `tk_ivec_t*` | Column-wise argmin indices |

#### Lua Bindings
All functions above have corresponding `_lua` suffixed versions for Lua registration:
- `tk_vec_name##_rasc_lua`, `tk_vec_name##_rdesc_lua` - Row sorting
- `tk_vec_name##_casc_lua`, `tk_vec_name##_cdesc_lua` - Column sorting  
- `tk_vec_name##_rmagnitudes_lua`, `tk_vec_name##_cmagnitudes_lua` - Magnitudes
- `tk_vec_name##_rmaxargs_lua`, `tk_vec_name##_cmaxargs_lua` - Argmax
- `tk_vec_name##_rminargs_lua`, `tk_vec_name##_cminargs_lua` - Argmin

### Type Definitions

#### `santoku/vec/base.h`
Base header aggregating all vector type definitions:
- `tk_ivec_t` - Integer vector type
- `tk_dvec_t` - Double vector type
- `tk_cvec_t` - Character vector type
- `tk_rvec_t` - Rank vector type
- `tk_pvec_t` - Pair vector type
- `tk_rank_t` - Rank pair structure `{ int64_t i; double d; }`
- `tk_pair_t` - Integer pair structure `{ int64_t i; int64_t p; }`

## Hash Maps and Sets

### `santoku/iuset.h`
Integer set implementation using hash tables.

| Function | Description |
|----------|-------------|
| `tk_iuset_create()` | Creates new integer set |
| `tk_iuset_destroy(set)` | Destroys set and frees memory |
| `tk_iuset_put(set, key, &ret)` | Adds key to set |
| `tk_iuset_get(set, key)` | Gets iterator for key |
| `tk_iuset_del(set, key)` | Removes key from set |
| `tk_iuset_exist(set, iter)` | Checks if iterator is valid |
| `tk_iuset_key(set, iter)` | Gets key at iterator |
| `tk_iuset_size(set)` | Returns number of elements |
| `tk_iuset_resize(set, size)` | Pre-allocates capacity |
| `tk_iuset_clear(set)` | Removes all elements |
| `tk_iuset_contains(set, key)` | Checks if key exists |
| `tk_iuset_foreach(set, var, code)` | Iterates over all keys |
| `tk_iuset_union(set_a, set_b)` | Unions set_b into set_a |
| `tk_iuset_jaccard(set_a, set_b)` | Computes Jaccard similarity |

### `santoku/duset.h`
Double set implementation using hash tables.

| Macro/Function | Description |
|-------|-------------|
| `tk_duset_create()` | Creates new double set |
| `tk_duset_destroy(set)` | Destroys set and frees memory |
| `tk_duset_put(set, key, &ret)` | Adds key to set |
| `tk_duset_get(set, key)` | Gets iterator for key |
| `tk_duset_del(set, key)` | Removes key from set |
| `tk_duset_exist(set, iter)` | Checks if iterator is valid |
| `tk_duset_key(set, iter)` | Gets key at iterator |
| `tk_duset_size(set)` | Returns number of elements |
| `tk_duset_resize(set, size)` | Pre-allocates capacity |
| `tk_duset_clear(set)` | Removes all elements |
| `tk_duset_contains(set, key)` | Checks if key exists |
| `tk_duset_foreach(set, var, code)` | Iterates over all keys |

### `santoku/cuset.h`
Character set implementation using hash tables.

| Macro/Function | Description |
|-------|-------------|
| `tk_cuset_create()` | Creates new character set |
| `tk_cuset_destroy(set)` | Destroys set and frees memory |
| `tk_cuset_put(set, key, &ret)` | Adds key to set |
| `tk_cuset_get(set, key)` | Gets iterator for key |
| `tk_cuset_del(set, key)` | Removes key from set |
| `tk_cuset_exist(set, iter)` | Checks if iterator is valid |
| `tk_cuset_key(set, iter)` | Gets key at iterator |
| `tk_cuset_size(set)` | Returns number of elements |
| `tk_cuset_resize(set, size)` | Pre-allocates capacity |
| `tk_cuset_clear(set)` | Removes all elements |
| `tk_cuset_contains(set, key)` | Checks if key exists |
| `tk_cuset_foreach(set, var, code)` | Iterates over all keys |

### `santoku/iumap.h`
Integer-to-integer hash map implementation.

| Macro/Function | Description |
|-------|-------------|
| `tk_iumap_create()` | Creates new int-to-int map |
| `tk_iumap_destroy(map)` | Destroys map and frees memory |
| `tk_iumap_put(map, key, &ret)` | Adds/updates key in map |
| `tk_iumap_get(map, key)` | Gets iterator for key |
| `tk_iumap_del(map, key)` | Removes key-value pair |
| `tk_iumap_exist(map, iter)` | Checks if iterator is valid |
| `tk_iumap_key(map, iter)` | Gets key at iterator |
| `tk_iumap_value(map, iter)` | Gets value at iterator |
| `tk_iumap_size(map)` | Returns number of pairs |
| `tk_iumap_resize(map, size)` | Pre-allocates capacity |
| `tk_iumap_clear(map)` | Removes all pairs |
| `tk_iumap_foreach(map, kvar, vvar, code)` | Iterates over all pairs |

### `santoku/dumap.h`
Double-to-double hash map implementation.

| Macro/Function | Description |
|-------|-------------|
| `tk_dumap_create()` | Creates new double-to-double map |
| `tk_dumap_destroy(map)` | Destroys map and frees memory |
| `tk_dumap_put(map, key, &ret)` | Adds/updates key in map |
| `tk_dumap_get(map, key)` | Gets iterator for key |
| `tk_dumap_del(map, key)` | Removes key-value pair |
| `tk_dumap_exist(map, iter)` | Checks if iterator is valid |
| `tk_dumap_key(map, iter)` | Gets key at iterator |
| `tk_dumap_value(map, iter)` | Gets value at iterator |
| `tk_dumap_size(map)` | Returns number of pairs |
| `tk_dumap_resize(map, size)` | Pre-allocates capacity |
| `tk_dumap_clear(map)` | Removes all pairs |
| `tk_dumap_foreach(map, kvar, vvar, code)` | Iterates over all pairs |

### `santoku/pumap.h`
Pair-to-pair hash map implementation for tk_pair_t types.

| Macro/Function | Description |
|-------|-------------|
| `tk_pumap_create()` | Creates new pair-to-pair map |
| `tk_pumap_destroy(map)` | Destroys map and frees memory |
| `tk_pumap_put(map, key, &ret)` | Adds/updates key in map |
| `tk_pumap_get(map, key)` | Gets iterator for key |
| `tk_pumap_del(map, key)` | Removes key-value pair |
| `tk_pumap_exist(map, iter)` | Checks if iterator is valid |
| `tk_pumap_key(map, iter)` | Gets key at iterator |
| `tk_pumap_value(map, iter)` | Gets value at iterator |
| `tk_pumap_size(map)` | Returns number of pairs |
| `tk_pumap_resize(map, size)` | Pre-allocates capacity |
| `tk_pumap_clear(map)` | Removes all pairs |
| `tk_pumap_foreach(map, kvar, vvar, code)` | Iterates over all pairs |

### `santoku/zumap.h`
Size_t-to-size_t hash map implementation.

| Macro/Function | Description |
|-------|-------------|
| `tk_zumap_create()` | Creates new size_t-to-size_t map |
| `tk_zumap_destroy(map)` | Destroys map and frees memory |
| `tk_zumap_put(map, key, &ret)` | Adds/updates key in map |
| `tk_zumap_get(map, key)` | Gets iterator for key |
| `tk_zumap_del(map, key)` | Removes key-value pair |
| `tk_zumap_exist(map, iter)` | Checks if iterator is valid |
| `tk_zumap_key(map, iter)` | Gets key at iterator |
| `tk_zumap_value(map, iter)` | Gets value at iterator |
| `tk_zumap_size(map)` | Returns number of pairs |
| `tk_zumap_resize(map, size)` | Pre-allocates capacity |
| `tk_zumap_clear(map)` | Removes all pairs |
| `tk_zumap_foreach(map, kvar, vvar, code)` | Iterates over all pairs |

