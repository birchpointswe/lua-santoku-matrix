local test = require("santoku.test")
local err = require("santoku.error")
local assert = err.assert
local csr = require("santoku.csr")
local mtx = require("santoku.mtx")
local ivec = require("santoku.ivec")
local dvec = require("santoku.dvec")
local fvec = require("santoku.fvec")
local svec = require("santoku.svec")
require("santoku.cvec")
local tbl = require("santoku.table")
local num = require("santoku.num")
local fs = require("santoku.fs")
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
  assert(num.abs(X:values():get(1) - 2.5) < 1e-6)
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

test("csr: i32 neighbors (svec) wrap + persist roundtrip", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 3, 5 }),
    neighbors = svec.create({ 0, 2, 1, 0, 3 }),
    n_cols = 4,
  })
  local r, c = X:shape()
  assert(r == 3 and c == 4)
  assert(X:nnz() == 5)
  assert(teq(X:neighbors():table(), { 0, 2, 1, 0, 3 }))
  local tmp = ".csr_i32_test.bin"
  X:persist(tmp)
  local Y = csr.load(tmp)
  fs.rm(tmp, true)
  assert(X:eq(Y))
  assert(teq(Y:neighbors():table(), { 0, 2, 1, 0, 3 }))
end)

test("csr: i32 neighbors standardize + hcat", function ()
  local A = csr.create({
    offsets = ivec.create({ 0, 2, 4 }),
    neighbors = svec.create({ 0, 1, 0, 1 }),
    values = fvec.create({ 2, 10, 4, 10 }),
    n_cols = 2,
  })
  local w = A:standardize()
  assert(num.abs(w:get(0) - 1) < 1e-5)
  assert(num.abs(w:get(1) - 0) < 1e-5)
  assert(teq(A:values():table(), { 2, 0, 4, 0 }))
  local B = csr.create({
    offsets = ivec.create({ 0, 1, 2 }),
    neighbors = svec.create({ 0, 1 }),
    values = fvec.create({ 5, 6 }),
    n_cols = 2,
  })
  assert(A:hcat(B) == A)
  local _, c = A:shape()
  assert(c == 4)
  assert(teq(A:neighbors():table(), { 0, 1, 2, 0, 1, 3 }))
end)

test("csr: i32 neighbors builder", function ()
  local X = csr.create({ n_cols = 4, neighbors = "i32", values = "f32" })
  X:push(0, 1.5):push(2, 2.5):row()
  X:push(3):row()
  assert(teq(X:offsets():table(), { 0, 2, 3 }))
  assert(teq(X:neighbors():table(), { 0, 2, 3 }))
  assert(num.abs(X:values():get(1) - 2.5) < 1e-6)
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
  assert(num.abs(v:get(0) - 1 / num.sqrt(2)) < 1e-6)
  assert(num.abs(v:get(2) - 1) < 1e-6)
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
  assert(num.abs(ss:get(0) - 1) < 1e-6)
  assert(num.abs(ss:get(1) - 9) < 1e-6)
  assert(num.abs(ss:get(2) - 20) < 1e-6)
  local blocks = X:sumsq_cols(ivec.create({ 0, 2, 3 }))
  assert(num.abs(blocks:get(0) - 10) < 1e-6)
  assert(num.abs(blocks:get(1) - 20) < 1e-6)
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
  fs.rm(tmp, true)
  assert(X:eq(Y))
  local B = csr.from_mask(ivec.create({ 1, 0, 1 }))
  B:persist(tmp)
  local B2 = csr.load(tmp)
  fs.rm(tmp, true)
  assert(B:eq(B2))
end)

test("csr: standardize fit/apply", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 4 }),
    neighbors = ivec.create({ 0, 1, 0, 1 }),
    values = fvec.create({ 2, 10, 4, 10 }),
    n_cols = 2,
  })
  local w = X:standardize()
  assert(num.abs(w:get(0) - 1) < 1e-5)
  assert(num.abs(w:get(1) - 0) < 1e-5)
  assert(teq(X:values():table(), { 2, 0, 4, 0 }))
  local Y = csr.create({
    offsets = ivec.create({ 0, 2 }),
    neighbors = ivec.create({ 0, 1 }),
    values = fvec.create({ 5, 7 }),
    n_cols = 2,
  })
  assert(Y:standardize(w) == w)
  assert(teq(Y:values():table(), { 5, 0 }))
end)

test("csr: bns fit/apply", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 1, 2, 3, 4 }),
    neighbors = ivec.create({ 0, 0, 1, 1 }),
    n_cols = 2,
  })
  local Y = csr.create({
    offsets = ivec.create({ 0, 1, 2, 2, 2 }),
    neighbors = ivec.create({ 0, 0 }),
    n_cols = 1,
  })
  local w = X:bns(Y)
  assert(#w:table() == 2)
  assert(w:get(1) == 0)
  assert(w:get(0) > 1.5)
  local w0 = w:get(0)
  assert(X:type() == "f32")
  assert(teq(X:values():table(), { w0, w0, 0, 0 }))
  local Z = csr.create({
    offsets = ivec.create({ 0, 2 }),
    neighbors = ivec.create({ 0, 1 }),
    n_cols = 2,
  })
  assert(Z:bns(w) == w)
  assert(teq(Z:values():table(), { w0, 0 }))
end)

test("csr: idf fit", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 1, 2, 4 }),
    neighbors = ivec.create({ 0, 0, 0, 1 }),
    values = fvec.create({ 1, 1, 1, 1 }),
    n_cols = 2,
  })
  local w = X:idf()
  assert(num.abs(w:get(0) - num.log(0.5 / 3.5)) < 1e-5)
  assert(num.abs(w:get(1) - num.log(2.5 / 1.5)) < 1e-5)
end)

test("csr.fuse: union with weighted sum, sorted descending", function ()
  local A = csr.create({
    offsets = ivec.create({ 0, 2 }),
    neighbors = ivec.create({ 10, 20 }),
    values = fvec.create({ 0.9, 0.5 }),
    n_cols = 100,
  })
  local B = csr.create({
    offsets = ivec.create({ 0, 2 }),
    neighbors = ivec.create({ 20, 30 }),
    values = fvec.create({ 0.8, 0.4 }),
    n_cols = 100,
  })
  local Y = csr.fuse(A, B)
  assert(Y:type() == "f64")
  assert(teq(Y:offsets():table(), { 0, 3 }))
  assert(teq(Y:neighbors():table(), { 20, 10, 30 }))
  local v = Y:values()
  assert(num.abs(v:get(0) - 1.3) < 1e-6)
  assert(num.abs(v:get(1) - 0.9) < 1e-6)
  assert(num.abs(v:get(2) - 0.4) < 1e-6)
  for i = 1, Y:nnz() - 1 do
    assert(v:get(i - 1) >= v:get(i))
  end
end)

test("csr.fuse: weights change ordering", function ()
  local A = csr.create({
    offsets = ivec.create({ 0, 2 }),
    neighbors = ivec.create({ 10, 20 }),
    values = fvec.create({ 0.9, 0.5 }),
    n_cols = 100,
  })
  local B = csr.create({
    offsets = ivec.create({ 0, 2 }),
    neighbors = ivec.create({ 20, 30 }),
    values = fvec.create({ 0.8, 0.4 }),
    n_cols = 100,
  })
  local Y = csr.fuse(A, B, { weights = { 1, 10 } })
  assert(teq(Y:neighbors():table(), { 20, 30, 10 }))
  local v = Y:values()
  assert(num.abs(v:get(0) - 8.5) < 1e-5)
  assert(num.abs(v:get(1) - 4.0) < 1e-5)
  assert(num.abs(v:get(2) - 0.9) < 1e-5)
end)

test("csr.fuse: k keeps top k", function ()
  local A = csr.create({
    offsets = ivec.create({ 0, 2 }),
    neighbors = ivec.create({ 10, 20 }),
    values = fvec.create({ 0.9, 0.5 }),
    n_cols = 100,
  })
  local B = csr.create({
    offsets = ivec.create({ 0, 2 }),
    neighbors = ivec.create({ 20, 30 }),
    values = fvec.create({ 0.8, 0.4 }),
    n_cols = 100,
  })
  local Y = csr.fuse(A, B, { k = 2 })
  assert(Y:nnz() == 2)
  assert(teq(Y:offsets():table(), { 0, 2 }))
  assert(teq(Y:neighbors():table(), { 20, 10 }))
end)

test("csr.fuse: rrf ranks are 0-based row positions, rows assumed pre-ranked, missing side adds nothing", function ()
  local A = csr.create({
    offsets = ivec.create({ 0, 2 }),
    neighbors = ivec.create({ 1, 2 }),
    values = fvec.create({ 100, 99 }),
    n_cols = 10,
  })
  local B = csr.create({
    offsets = ivec.create({ 0, 2 }),
    neighbors = ivec.create({ 2, 3 }),
    values = fvec.create({ 0.5, 0.4 }),
    n_cols = 10,
  })
  local S = csr.fuse(A, B)
  assert(S:neighbors():get(0) == 1)
  local Y = csr.fuse(A, B, { mode = "rrf", rrf_k = 1 })
  assert(teq(Y:neighbors():table(), { 2, 1, 3 }))
  local v = Y:values()
  assert(num.abs(v:get(0) - 1.5) < 1e-9)
  assert(num.abs(v:get(1) - 1.0) < 1e-9)
  assert(num.abs(v:get(2) - 0.5) < 1e-9)
end)

test("csr.fuse: multi-row with a row empty on one side", function ()
  local A = csr.create({
    offsets = ivec.create({ 0, 2, 2 }),
    neighbors = ivec.create({ 1, 2 }),
    values = fvec.create({ 0.9, 0.5 }),
    n_cols = 10,
  })
  local B = csr.create({
    offsets = ivec.create({ 0, 1, 3 }),
    neighbors = ivec.create({ 2, 5, 6 }),
    values = fvec.create({ 0.8, 0.3, 0.7 }),
    n_cols = 10,
  })
  local Y = csr.fuse(A, B)
  local r, c = Y:shape()
  assert(r == 2 and c == 10)
  assert(teq(Y:offsets():table(), { 0, 2, 4 }))
  assert(teq(Y:neighbors():table(), { 2, 1, 6, 5 }))
  local v = Y:values()
  assert(num.abs(v:get(0) - 1.3) < 1e-6)
  assert(num.abs(v:get(1) - 0.9) < 1e-6)
  assert(num.abs(v:get(2) - 0.7) < 1e-6)
  assert(num.abs(v:get(3) - 0.3) < 1e-6)
end)

test("csr.fuse: errors on mismatched rows and on valueless input", function ()
  local A = csr.create({
    offsets = ivec.create({ 0, 1 }),
    neighbors = ivec.create({ 0 }),
    values = fvec.create({ 1 }),
    n_cols = 2,
  })
  local B = csr.create({
    offsets = ivec.create({ 0, 1, 2 }),
    neighbors = ivec.create({ 0, 1 }),
    values = fvec.create({ 1, 1 }),
    n_cols = 2,
  })
  assert(not pcall(function () csr.fuse(A, B) end))
  local C = csr.create({
    offsets = ivec.create({ 0, 1 }),
    neighbors = ivec.create({ 1 }),
    n_cols = 2,
  })
  assert(not pcall(function () csr.fuse(A, C) end))
  assert(not pcall(function () csr.fuse(C, A) end))
  assert(not pcall(function () csr.fuse(A, A, { mode = "product" }) end))
end)

