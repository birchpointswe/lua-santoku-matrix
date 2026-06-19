local test = require("santoku.test")
local err = require("santoku.error")
local assert = err.assert
local csr = require("santoku.csr")
local mtx = require("santoku.mtx")
local ivec = require("santoku.ivec")
local dvec = require("santoku.dvec")
local fvec = require("santoku.fvec")
require("santoku.cvec")
local tbl = require("santoku.table")
local teq = tbl.equals

test("csr: wrap parts, accessors", function ()
  local off = ivec.create({ 0, 2, 3, 5 })
  local nbr = ivec.create({ 0, 2, 1, 0, 3 })
  local X = csr.create({ offsets = off, neighbors = nbr, n_cols = 4 })
  local r, c = X:shape()
  assert(r == 3 and c == 4)
  assert(X:nnz() == 5)
  assert(X:type() == "none")
  assert(X:offsets() == off)
  assert(X:neighbors() == nbr)
  assert(X:values() == nil)
end)

test("csr: builder push/row", function ()
  local X = csr.create({ n_cols = 4, values = "f32" })
  X:push(0, 1.5):push(2, 2.5):row()
  X:row()
  X:push(3):row()
  local r, c = X:shape()
  assert(r == 3 and c == 4)
  assert(teq(X:offsets():table(), { 0, 2, 2, 3 }))
  assert(teq(X:neighbors():table(), { 0, 2, 3 }))
  assert(math.abs(X:values():get(1) - 2.5) < 1e-6)
end)

test("csr.from_classes", function ()
  local X = csr.from_classes(ivec.create({ 2, 0, 1, 2 }))
  local r, c = X:shape()
  assert(r == 4 and c == 3)
  assert(teq(X:offsets():table(), { 0, 1, 2, 3, 4 }))
  assert(teq(X:neighbors():table(), { 2, 0, 1, 2 }))
end)

test("csr.from_mask", function ()
  local X = csr.from_mask(ivec.create({ 1, 0, 1, 1, 0 }))
  local r, c = X:shape()
  assert(r == 5 and c == 1)
  assert(teq(X:offsets():table(), { 0, 1, 1, 2, 3, 3 }))
  assert(teq(X:neighbors():table(), { 0, 0, 0 }))
end)

test("csr: to_bits/from_bits roundtrip", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 3, 5 }),
    neighbors = ivec.create({ 0, 2, 1, 0, 3 }),
    n_cols = 4,
  })
  local bits = X:to_bits()
  local Y = csr.from_bits(bits, 3, 4)
  assert(X:eq(Y))
end)

test("csr: to_dense / mtx:to_sparse roundtrip", function ()
  local M = mtx.create({
    data = dvec.create({ 1, 0, 3, 0, 2, 0 }),
    n_rows = 2, n_cols = 3,
  })
  local X = M:to_sparse()
  local r, c = X:shape()
  assert(r == 2 and c == 3)
  assert(teq(X:offsets():table(), { 0, 2, 3 }))
  assert(teq(X:neighbors():table(), { 0, 2, 1 }))
  assert(X:values():get(1) == 3)
  local D = X:to_dense()
  assert(D:eq(M))
end)

test("csr: rows gather (with values)", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 4, 6 }),
    neighbors = ivec.create({ 10, 20, 30, 40, 50, 60 }),
    values = fvec.create({ 1, 2, 3, 4, 5, 6 }),
    n_cols = 100,
  })
  local Y = X:rows(ivec.create({ 0, 2 }))
  assert(teq(Y:offsets():table(), { 0, 2, 4 }))
  assert(teq(Y:neighbors():table(), { 10, 20, 50, 60 }))
  assert(Y:values():get(2) == 5)
end)

test("csr: select columns with remap", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 3, 5 }),
    neighbors = ivec.create({ 0, 1, 2, 1, 3 }),
    values = fvec.create({ 1, 2, 3, 4, 5 }),
    n_cols = 4,
  })
  local Y = X:select(ivec.create({ 1, 3 }))
  local _, c = Y:shape()
  assert(c == 2)
  assert(teq(Y:offsets():table(), { 0, 1, 3 }))
  assert(teq(Y:neighbors():table(), { 0, 0, 1 }))
  assert(Y:values():get(0) == 2)
  assert(Y:values():get(1) == 4)
  assert(Y:values():get(2) == 5)
end)

test("csr: hcat in place with shift", function ()
  local A = csr.create({
    offsets = ivec.create({ 0, 2, 3 }),
    neighbors = ivec.create({ 0, 1, 2 }),
    n_cols = 3,
  })
  local B = csr.create({
    offsets = ivec.create({ 0, 1, 2 }),
    neighbors = ivec.create({ 0, 1 }),
    n_cols = 2,
  })
  assert(A:hcat(B) == A)
  local _, c = A:shape()
  assert(c == 5)
  assert(teq(A:offsets():table(), { 0, 3, 5 }))
  assert(teq(A:neighbors():table(), { 0, 1, 3, 2, 4 }))
  assert(teq(B:neighbors():table(), { 0, 1 }))
end)

test("csr: transpose with values", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 3 }),
    neighbors = ivec.create({ 0, 1, 0 }),
    values = fvec.create({ 1, 2, 3 }),
    n_cols = 2,
  })
  local T = X:transpose()
  local r, c = T:shape()
  assert(r == 2 and c == 2)
  assert(teq(T:offsets():table(), { 0, 2, 3 }))
  assert(teq(T:neighbors():table(), { 0, 1, 0 }))
  assert(T:values():get(0) == 1)
  assert(T:values():get(1) == 3)
  assert(T:values():get(2) == 2)
end)

test("csr: normalize materializes values on binary", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 3 }),
    neighbors = ivec.create({ 0, 1, 1 }),
    n_cols = 2,
  })
  assert(X:normalize() == X)
  assert(X:type() == "f32")
  local v = X:values()
  assert(math.abs(v:get(0) - 1 / math.sqrt(2)) < 1e-6)
  assert(math.abs(v:get(2) - 1) < 1e-6)
end)

test("csr: scale_cols", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 3 }),
    neighbors = ivec.create({ 0, 1, 1 }),
    values = fvec.create({ 2, 3, 4 }),
    n_cols = 2,
  })
  X:scale_cols(fvec.create({ 10, 100 }))
  assert(X:values():get(0) == 20)
  assert(X:values():get(1) == 300)
  assert(X:values():get(2) == 400)
end)

test("csr: sumsq_cols, plain and blocked", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 4 }),
    neighbors = ivec.create({ 0, 2, 1, 2 }),
    values = fvec.create({ 1, 2, 3, 4 }),
    n_cols = 3,
  })
  local ss = X:sumsq_cols()
  assert(math.abs(ss:get(0) - 1) < 1e-6)
  assert(math.abs(ss:get(1) - 9) < 1e-6)
  assert(math.abs(ss:get(2) - 20) < 1e-6)
  local blocks = X:sumsq_cols(ivec.create({ 0, 2, 3 }))
  assert(math.abs(blocks:get(0) - 10) < 1e-6)
  assert(math.abs(blocks:get(1) - 20) < 1e-6)
end)

test("csr: persist/load roundtrip", function ()
  local tmp = ".csr_test.bin"
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 3 }),
    neighbors = ivec.create({ 0, 1, 1 }),
    values = fvec.create({ 1.5, 2.5, 3.5 }),
    n_cols = 2,
  })
  X:persist(tmp)
  local Y = csr.load(tmp)
  os.remove(tmp)
  assert(X:eq(Y))
  local B = csr.from_mask(ivec.create({ 1, 0, 1 }))
  B:persist(tmp)
  local B2 = csr.load(tmp)
  os.remove(tmp)
  assert(B:eq(B2))
end)
