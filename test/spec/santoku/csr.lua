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

test("csr: normalize max divides each row by its max", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 4, 5 }),
    neighbors = ivec.create({ 0, 1, 0, 1, 1 }),
    values = fvec.create({ 2, 4, -3, -1, 0 }),
    n_cols = 2,
  })
  assert(X:normalize("max") == X)
  assert(teq(X:values():table(), { 0.5, 1, -3, -1, 0 }))
  assert(not pcall(function () X:normalize("l1") end))
end)

test("csr: dots writes sampled pair products in place", function ()
  local Q = mtx.create({ data = fvec.create({ 1, 0, 0, 1 }), n_rows = 2, n_cols = 2 })
  local D = mtx.create({ data = fvec.create({ 1, 2, 3, 4, 5, 6 }), n_rows = 3, n_cols = 2 })
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 3 }),
    neighbors = ivec.create({ 2, 0, 1 }),
    values = fvec.create({ 9, 9, 9 }),
    n_cols = 3,
  })
  assert(X:dots(Q, D) == X)
  assert(teq(X:neighbors():table(), { 2, 0, 1 }))
  assert(teq(X:values():table(), { 5, 1, 4 }))
  local B = csr.create({ offsets = ivec.create({ 0, 1 }), neighbors = ivec.create({ 1 }), n_cols = 3 })
  B:dots(Q, D)
  assert(B:type() == "f32" and B:values():get(0) == 3)
  local Q64 = mtx.create({ data = dvec.create({ 1, 1 }), n_rows = 1, n_cols = 2 })
  local D64 = mtx.create({ data = dvec.create({ 2, 3 }), n_rows = 1, n_cols = 2 })
  local C = csr.create({ offsets = ivec.create({ 0, 1 }), neighbors = ivec.create({ 0 }), n_cols = 1 })
  C:dots(Q64, D64)
  assert(C:type() == "f64" and C:values():get(0) == 5)
  assert(not pcall(function () X:dots(Q64, D) end))
  local bad = csr.create({ offsets = ivec.create({ 0, 1 }), neighbors = ivec.create({ 3 }), n_cols = 4 })
  assert(not pcall(function () bad:dots(Q, D) end))
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

test("csr: bm25 fit/apply", function ()
  local function make ()
    return csr.create({
      offsets = ivec.create({ 0, 2, 4, 7, 9, 9 }),
      neighbors = ivec.create({ 0, 3, 1, 3, 1, 2, 3, 2, 3 }),
      values = fvec.create({ 2, 1, 3, 1, 1, 2, 1, 3, 1 }),
      n_cols = 4,
    })
  end
  local X = make()
  local w, avgdl = X:bm25()
  assert(num.abs(avgdl - 3) < 1e-9)
  assert(num.abs(w:get(0) - num.log(3)) < 1e-5)
  assert(num.abs(w:get(1) - num.log(1.4)) < 1e-5)
  assert(num.abs(w:get(2) - num.log(1.4)) < 1e-5)
  assert(w:get(3) > 0 and w:get(3) < 1e-5)
  local v = X:values()
  assert(num.abs(v:get(0) - num.log(3) * 2 * 2.2 / 3.2) < 1e-5)
  assert(num.abs(v:get(1) - 1e-6) < 1e-8)
  assert(num.abs(v:get(2) - num.log(1.4) * 6.6 / 4.5) < 1e-5)
  assert(num.abs(v:get(4) - num.log(1.4) * 2.2 / 2.5) < 1e-5)
  assert(num.abs(v:get(5) - num.log(1.4) * 4.4 / 3.5) < 1e-5)
  assert(num.abs(v:get(7) - num.log(1.4) * 6.6 / 4.5) < 1e-5)
  local Y = csr.create({
    offsets = ivec.create({ 0, 1 }),
    neighbors = ivec.create({ 0 }),
    values = fvec.create({ 1 }),
    n_cols = 4,
  })
  local w2, avgdl2 = Y:bm25(w, avgdl)
  assert(w2 == w)
  assert(num.abs(avgdl2 - 3) < 1e-9)
  assert(num.abs(Y:values():get(0) - num.log(3) * 2.2 / 1.6) < 1e-5)
  local X2 = make()
  local _, avgdl3 = X2:bm25(2.0, 0.0)
  assert(num.abs(avgdl3 - 3) < 1e-9)
  assert(num.abs(X2:values():get(0) - num.log(3) * 2 * 3 / 4) < 1e-5)
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

test("csr: topk scores queries against docs by sparse dot product", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 3, 4 }),
    neighbors = ivec.create({ 0, 1, 1, 0 }),
    values = fvec.create({ 1, 2, 1, 3 }),
    n_cols = 2,
  })
  local Q = csr.create({
    offsets = ivec.create({ 0, 1, 2, 3, 3 }),
    neighbors = ivec.create({ 0, 1, 5 }),
    values = fvec.create({ 1, 2, 1 }),
    n_cols = 6,
  })
  local P = X:topk(Q, 2)
  local r, c = P:shape()
  assert(r == 4 and c == 3)
  assert(teq(P:offsets():table(), { 0, 2, 4, 4, 4 }))
  assert(teq(P:neighbors():table(), { 2, 0, 0, 1 }))
  assert(teq(P:values():table(), { 3, 1, 4, 2 }))
  local P1 = X:topk(Q, 1)
  assert(teq(P1:neighbors():table(), { 2, 0 }))
  local Bad = csr.create({
    offsets = ivec.create({ 0, 1 }),
    neighbors = ivec.create({ -1 }),
    n_cols = 1,
  })
  assert(not pcall(function () X:topk(Bad, 1) end))
end)

test("csr: ndcg, recall and mrr against graded judgments", function ()
  local R = csr.create({
    offsets = ivec.create({ 0, 3, 5, 6 }),
    neighbors = ivec.create({ 5, 3, 9, 1, 2, 7 }),
    values = fvec.create({ 3, 2, 1, 2, 1, 1 }),
    n_cols = 10,
  })
  local J = csr.create({
    offsets = ivec.create({ 0, 3, 4, 4 }),
    neighbors = ivec.create({ 3, 9, 4, 7 }),
    values = fvec.create({ 2, 1, 1, 1 }),
    n_cols = 10,
  })
  local function log2 (x) return num.log(x) / num.log(2) end
  local nd0 = (2 / log2(3) + 1 / log2(4)) / (2 + 1 / log2(3) + 1 / log2(4))
  local nd, ndm = R:ndcg(J, 3)
  assert(nd:size() == 3)
  assert(num.abs(nd:get(0) - nd0) < 1e-9)
  assert(nd:get(1) == 0 and nd:get(2) == 0)
  assert(num.abs(ndm - nd0 / 2) < 1e-9)
  local rc, rcm = R:recall(J, 3)
  assert(num.abs(rc:get(0) - 2 / 3) < 1e-9)
  assert(num.abs(rcm - 1 / 3) < 1e-9)
  local mr, mrm = R:mrr(J, 3)
  assert(mr:get(0) == 0.5 and mr:get(1) == 0)
  assert(mrm == 0.25)
  local _, cut = R:recall(J, 1)
  assert(cut == 0)
  assert(not pcall(function () R:rows(ivec.create({ 0 })):ndcg(J, 3) end))
end)

test("csr: spearman compares two rankings of the same candidates per row", function ()
  local off = ivec.create({ 0, 4, 7, 8, 11 })
  local nbr = ivec.create({ 0, 1, 2, 3, 0, 1, 2, 5, 0, 1, 2 })
  local A = csr.create({ offsets = off, neighbors = nbr, n_cols = 6,
    values = fvec.create({ 1, 2, 3, 4, 1, 2, 3, 7, 1, 1, 2 }) })
  local B = csr.create({ offsets = off:clone(), neighbors = nbr:clone(), n_cols = 6,
    values = fvec.create({ 10, 20, 30, 40, 3, 2, 1, 9, 1, 2, 3 }) })
  local rho = A:spearman(B)
  assert(rho:size() == 4)
  assert(num.abs(rho:get(0) - 1) < 1e-12)
  assert(num.abs(rho:get(1) + 1) < 1e-12)
  assert(rho:get(2) == 0)
  assert(num.abs(rho:get(3) - 1.5 / num.sqrt(3)) < 1e-12)
  local share = A:overlap(B, 2)
  assert(share:get(0) == 1)
  assert(share:get(1) == 0.5)
  assert(share:get(2) == 1)
  local C = csr.create({ offsets = off:clone(), n_cols = 6,
    neighbors = ivec.create({ 0, 1, 2, 4, 0, 1, 2, 5, 0, 1, 2 }) })
  assert(not pcall(function () A:spearman(C) end))
end)

test("csr: unique_cols returns first-seen ids and a remapped copy", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 2, 4 }),
    neighbors = ivec.create({ 7, 3, 7, 9 }),
    values = fvec.create({ 1, 2, 3, 4 }),
    n_cols = 10,
  })
  local ids, Xc = X:unique_cols()
  assert(teq(ids:table(), { 7, 3, 9 }))
  assert(teq(Xc:neighbors():table(), { 0, 1, 0, 2 }))
  assert(teq(Xc:values():table(), { 1, 2, 3, 4 }))
  assert(teq(Xc:offsets():table(), { 0, 2, 4 }))
  local _, c = Xc:shape()
  assert(c == 3)
  for j = 0, 3 do
    assert(ids:get(Xc:neighbors():get(j)) == X:neighbors():get(j))
  end
  assert(teq(X:neighbors():table(), { 7, 3, 7, 9 }))
end)

test("csr: auc against a regression target ranks each column", function ()
  local y = { 3, 1, 4, 1.5, 5 }
  local vals = {}
  for d = 1, #y do
    vals[#vals + 1] = y[d]
    vals[#vals + 1] = -y[d]
    vals[#vals + 1] = 1
  end
  local X = csr.create({
    offsets = ivec.create({ 0, 3, 6, 9, 12, 15 }),
    neighbors = ivec.create({ 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2 }),
    values = fvec.create(vals),
    n_cols = 3,
  })
  local w = X:auc(dvec.create(y))
  assert(w:get(0) > 5)
  assert(num.abs(w:get(0) - w:get(1)) < 1e-5)
  assert(w:get(2) == 0)
end)

