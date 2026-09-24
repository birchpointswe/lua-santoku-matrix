local test = require("santoku.test")
local err = require("santoku.error")
local assert = err.assert
local mtx = require("santoku.mtx")
local ivec = require("santoku.ivec")
local dvec = require("santoku.dvec")
local tbl = require("santoku.table")
local num = require("santoku.num")
local fs = require("santoku.fs")
local teq = tbl.equals

test("mtx: create alloc is zeroed", function ()
  local M = mtx.create({ n_rows = 2, n_cols = 3, type = "f64" })
  local r, c = M:shape()
  assert(r == 2 and c == 3)
  assert(M:type() == "f64")
  for i = 0, 1 do
    for j = 0, 2 do
      assert(M:get(i, j) == 0)
    end
  end
end)

test("mtx: wrap existing vec, data() returns child", function ()
  local v = dvec.create({ 1, 2, 3, 4, 5, 6 })
  local M = mtx.create({ data = v, n_rows = 2, n_cols = 3 })
  assert(M:type() == "f64")
  assert(M:get(1, 2) == 6)
  assert(M:data() == v)
  v:set(5, 60)
  assert(M:get(1, 2) == 60)
end)

test("mtx: wrap infers integer types", function ()
  local v = ivec.create({ 1, 2, 3, 4 })
  local M = mtx.create({ data = v, n_rows = 2, n_cols = 2 })
  assert(M:type() == "i64")
  assert(M:get(0, 1) == 2)
end)

test("mtx: set/fill/eq", function ()
  local A = mtx.create({ n_rows = 2, n_cols = 2 })
  local B = mtx.create({ n_rows = 2, n_cols = 2 })
  A:fill(7)
  B:fill(7)
  assert(A:eq(B))
  B:set(1, 1, 8)
  assert(not A:eq(B))
  assert(A:eq(B, 1.5))
  local C = mtx.create({ n_rows = 2, n_cols = 3 })
  assert(not A:eq(C))
end)

test("mtx: sums/maxs/mins by axis", function ()
  local M = mtx.create({ data = dvec.create({ 1, 2, 3, 4, 5, 6 }), n_rows = 2, n_cols = 3 })
  assert(teq(M:sums("row"):table(), { 6, 15 }))
  assert(teq(M:sums("col"):table(), { 5, 7, 9 }))
  assert(teq(M:maxs("row"):table(), { 3, 6 }))
  assert(teq(M:mins("col"):table(), { 1, 2, 3 }))
end)

test("mtx: sums on i64", function ()
  local M = mtx.create({ data = ivec.create({ 1, 2, 3, 4, 5, 6 }), n_rows = 2, n_cols = 3 })
  assert(teq(M:sums("row"):table(), { 6, 15 }))
end)

test("mtx: maxargs/minargs", function ()
  local M = mtx.create({ data = dvec.create({ 1, 5, 2, 9, 3, 7 }), n_rows = 2, n_cols = 3 })
  assert(teq(M:maxargs("row"):table(), { 1, 0 }))
  assert(teq(M:minargs("row"):table(), { 0, 1 }))
  assert(teq(M:maxargs("col"):table(), { 1, 0, 1 }))
end)

test("mtx: mags", function ()
  local M = mtx.create({ data = dvec.create({ 3, 4, 0, 5, 12, 0 }), n_rows = 2, n_cols = 3 })
  local m = M:mags("row")
  assert(num.abs(m:get(0) - 5) < 1e-10)
  assert(num.abs(m:get(1) - 13) < 1e-10)
end)

test("mtx: argsort", function ()
  local M = mtx.create({ data = dvec.create({ 5, 1, 3, 9, 2, 7 }), n_rows = 2, n_cols = 3 })
  assert(teq(M:argsort("row", "asc"):table(), { 1, 2, 0, 1, 2, 0 }))
  assert(teq(M:argsort("row", "desc"):table(), { 0, 2, 1, 0, 2, 1 }))
end)

test("mtx: transpose", function ()
  local M = mtx.create({ data = dvec.create({ 1, 2, 3, 4, 5, 6 }), n_rows = 2, n_cols = 3 })
  local T = M:transpose()
  local r, c = T:shape()
  assert(r == 3 and c == 2)
  assert(teq(T:data():table(), { 1, 4, 2, 5, 3, 6 }))
end)

test("mtx: rows/cols/row gather", function ()
  local M = mtx.create({ data = dvec.create({ 1, 2, 3, 4, 5, 6, 7, 8, 9 }), n_rows = 3, n_cols = 3 })
  local R = M:rows(ivec.create({ 2, 0 }))
  assert(teq(R:data():table(), { 7, 8, 9, 1, 2, 3 }))
  local C = M:cols(ivec.create({ 2, 1 }))
  assert(teq(C:data():table(), { 3, 2, 6, 5, 9, 8 }))
  local row1 = M:row(1)
  assert(teq(row1:table(), { 4, 5, 6 }))
end)

test("mtx: hcat in place", function ()
  local A = mtx.create({ data = dvec.create({ 1, 2, 5, 6 }), n_rows = 2, n_cols = 2 })
  local B = mtx.create({ data = dvec.create({ 3, 7 }), n_rows = 2, n_cols = 1 })
  local r = A:hcat(B)
  assert(r == A)
  local rr, cc = A:shape()
  assert(rr == 2 and cc == 3)
  assert(teq(A:data():table(), { 1, 2, 3, 5, 6, 7 }))
  assert(teq(B:data():table(), { 3, 7 }))
end)

test("mtx: persist/load roundtrip", function ()
  local tmp = ".mtx_test.bin"
  local M = mtx.create({ data = dvec.create({ 1.5, 2.5, 3.5, 4.5 }), n_rows = 2, n_cols = 2 })
  M:persist(tmp)
  local M2 = mtx.load(tmp)
  fs.rm(tmp, true)
  assert(M:eq(M2))
  assert(M2:type() == "f64")
end)

test("mtx: center fit/apply", function ()
  local M = mtx.create({ data = dvec.create({ 1, 10, 3, 20 }), n_rows = 2, n_cols = 2 })
  local means = M:center()
  assert(num.abs(means:get(0) - 2) < 1e-10)
  assert(num.abs(means:get(1) - 15) < 1e-10)
  assert(teq(M:data():table(), { -1, -5, 1, 5 }))
  local N = mtx.create({ data = dvec.create({ 2, 15 }), n_rows = 1, n_cols = 2 })
  assert(N:center(means) == N)
  assert(teq(N:data():table(), { 0, 0 }))
end)

test("mtx: standardize fit returns means + inv_std", function ()
  local M = mtx.create({ data = dvec.create({ 1, 10, 3, 20 }), n_rows = 2, n_cols = 2 })
  local means, istds = M:standardize()
  assert(means and istds)
  local sums = M:sums("col")
  assert(num.abs(sums:get(0)) < 1e-9)
  assert(num.abs(sums:get(1)) < 1e-9)
end)

test("mtx: normalize row", function ()
  local M = mtx.create({ data = dvec.create({ 3, 4, 0, 0, 5, 12 }), n_rows = 2, n_cols = 3 })
  assert(M:normalize("row") == M)
  local m = M:mags("row")
  assert(num.abs(m:get(0) - 1) < 1e-10)
  assert(num.abs(m:get(1) - 1) < 1e-10)
end)

test("mtx: multiply", function ()
  local A = mtx.create({ data = dvec.create({ 1, 2, 3, 4 }), n_rows = 2, n_cols = 2 })
  local B = mtx.create({ data = dvec.create({ 5, 6, 7, 8 }), n_rows = 2, n_cols = 2 })
  local C = A:multiply(B)
  assert(teq(C:data():table(), { 19, 22, 43, 50 }))
  local Ct = A:multiply(B, false, true)
  assert(teq(Ct:data():table(), { 17, 23, 39, 53 }))
end)

test("mtx: multiplyv", function ()
  local A = mtx.create({ data = dvec.create({ 1, 2, 3, 4, 5, 6 }), n_rows = 2, n_cols = 3 })
  local y = A:multiplyv(dvec.create({ 1, 2, 3 }))
  assert(teq(y:table(), { 14, 32 }))
  local yt = A:multiplyv(dvec.create({ 1, 2 }), true)
  assert(teq(yt:table(), { 9, 12, 15 }))
end)

test("mtx: sign/median produce bitmaps", function ()
  local M = mtx.create({ data = dvec.create({ 1, -1, -2, 2 }), n_rows = 2, n_cols = 2 })
  local bits = M:sign()
  assert(bits ~= nil)
  local bits2, medians = M:median()
  assert(bits2 ~= nil and medians ~= nil)
end)

local function itq_data (n, k)
  local fvec = require("santoku.fvec")
  local x, vals = 12345, {}
  for i = 1, n do
    for j = 1, k do
      x = (x * 1103515245 + 12345) % 2147483648
      vals[(i - 1) * k + j] = (x / 2147483648 - 0.5) / j
    end
  end
  local M = mtx.create({ data = fvec.create(vals), n_rows = n, n_cols = k })
  M:center()
  return M
end

local function orthonormal_cols (W, tol)
  local _, c = W:shape()
  local G = W:multiply(W, true, false)
  for i = 0, c - 1 do
    for j = 0, c - 1 do
      local want = i == j and 1 or 0
      if num.abs(G:get(i, j) - want) > tol then return false end
    end
  end
  return true
end

test("mtx: itq full width rotation is orthonormal and lowers quantization error", function ()
  local M = itq_data(200, 16)
  local W, obj, steps, kept = M:itq({ iterations = 30 })
  local r, c = W:shape()
  assert(r == 16 and c == 16)
  assert(kept == 1)
  assert(orthonormal_cols(W, 1e-2))
  assert(obj:size() == 30 and steps:size() == 30)
  for i = 1, obj:size() - 1 do
    assert(obj:get(i) <= obj:get(i - 1) + 1e-6 * num.abs(obj:get(i - 1)))
  end
  assert(obj:get(obj:size() - 1) < obj:get(0))
end)

test("mtx: itq reduced width keeps top variance and orthonormal columns", function ()
  local M = itq_data(200, 16)
  local W, obj, _, kept = M:itq({ bits = 4, iterations = 20 })
  local r, c = W:shape()
  assert(r == 16 and c == 4)
  assert(kept > 0.5 and kept < 1)
  assert(orthonormal_cols(W, 1e-2))
  assert(obj:get(obj:size() - 1) <= obj:get(0))
  local Wp, objp, _, keptp = M:itq({ bits = 4, rotate = false })
  assert(objp:size() == 0)
  assert(num.abs(keptp - kept) < 1e-6)
  assert(orthonormal_cols(Wp, 1e-2))
end)

test("mtx: itq bits feed exhaustive hamming topk", function ()
  local M = itq_data(200, 16)
  local W = M:itq({ iterations = 20 })
  local F = M:multiply(W)
  local r, c = F:shape()
  local B = mtx.create({ data = F:sign(), n_rows = r, n_cols = c, bits = true })
  local P = B:topk(B, 3)
  local off, dist = P:offsets(), P:values()
  for q = 0, r - 1 do
    assert(dist:get(off:get(q)) == 0)
  end
end)

test("mtx.from_pairs: counts and weights", function ()
  local i = ivec.create({ 0, 0, 1, 2, 2, 2 })
  local j = ivec.create({ 0, 1, 1, 0, 0, 1 })
  local M = mtx.from_pairs(i, j, 3, 2)
  assert(teq(M:data():table(), { 1, 1, 0, 1, 2, 1 }))
  local W = mtx.from_pairs(i, j, 3, 2, require("santoku.dvec").create({ 0.5, 0.5, 1, 2, 2, 3 }))
  assert(teq(W:data():table(), { 0.5, 0.5, 0, 1, 4, 3 }))
end)

test("mtx: bits layout, popcount/hamming/ops", function ()
  local csr = require("santoku.csr")
  local ivec2 = require("santoku.ivec")
  local A = csr.create({
    offsets = ivec2.create({ 0, 2, 3 }),
    neighbors = ivec2.create({ 0, 2, 1 }),
    n_cols = 4,
  })
  local B = csr.create({
    offsets = ivec2.create({ 0, 1, 3 }),
    neighbors = ivec2.create({ 0, 1, 3 }),
    n_cols = 4,
  })
  local MA = mtx.create({ data = A:to_bits(), n_rows = 2, n_cols = 4, bits = true })
  local MB = mtx.create({ data = B:to_bits(), n_rows = 2, n_cols = 4, bits = true })
  assert(MA:type() == "bits")
  assert(MA:popcount() == 3)
  assert(MB:popcount() == 3)
  assert(MA:hamming(MB) == 2)
  local T = MA:transpose()
  local r, c = T:shape()
  assert(r == 4 and c == 2)
  assert(T:popcount() == 3)
  MA:band(MB)
  assert(MA:popcount() == 2)
end)

test("mtx: bits eq and persist roundtrip", function ()
  local csr = require("santoku.csr")
  local ivec2 = require("santoku.ivec")
  local X = csr.create({
    offsets = ivec2.create({ 0, 2, 3 }),
    neighbors = ivec2.create({ 0, 2, 1 }),
    n_cols = 4,
  })
  local M = mtx.create({ data = X:to_bits(), n_rows = 2, n_cols = 4, bits = true })
  local tmp = ".mtx_bits_test.bin"
  M:persist(tmp)
  local M2 = mtx.load(tmp)
  fs.rm(tmp, true)
  assert(M2:type() == "bits")
  assert(M:eq(M2))
end)

test("mtx: topk returns csr", function ()
  local corpus = mtx.create({ data = dvec.create({ 1, 0, 0, 1, 0.5, 0.5 }), n_rows = 3, n_cols = 2 })
  local q = mtx.create({ data = dvec.create({ 1, 0 }), n_rows = 1, n_cols = 2 })
  local P = corpus:topk(q, 2)
  local r, c = P:shape()
  assert(r == 1 and c == 3)
  assert(P:nnz() == 2)
  local ids = P:neighbors()
  assert(ids:get(0) == 0)
end)

test("mtx: out= destinations are reused", function ()
  local M = mtx.create({ data = dvec.create({ 1, 2, 3, 4, 5, 6 }), n_rows = 3, n_cols = 2 })
  local out = mtx.create({ n_rows = 1, n_cols = 1 })
  local R = M:rows(ivec.create({ 2, 0 }), out)
  assert(R == out)
  assert(teq(out:data():table(), { 5, 6, 1, 2 }))
  local buf = dvec.create(0)
  local row = M:row(1, buf)
  assert(row == buf)
  assert(teq(buf:table(), { 3, 4 }))
  local y = dvec.create(0)
  assert(M:multiplyv(dvec.create({ 1, 1 }), false, y) == y)
  assert(teq(y:table(), { 3, 7, 11 }))
  local A = mtx.create({ data = dvec.create({ 1, 2, 3, 4 }), n_rows = 2, n_cols = 2 })
  local B = mtx.create({ data = dvec.create({ 5, 6, 7, 8 }), n_rows = 2, n_cols = 2 })
  local C = mtx.create({ n_rows = 1, n_cols = 1 })
  assert(A:multiply(B, false, false, C) == C)
  assert(teq(C:data():table(), { 19, 22, 43, 50 }))
end)

test("mtx/csr: to_sparse/to_dense out= reuse", function ()
  local csr = require("santoku.csr")
  local M = mtx.create({ data = dvec.create({ 1, 0, 3, 0 }), n_rows = 2, n_cols = 2 })
  local X = csr.create({ n_cols = 2, values = "f64" })
  assert(M:to_sparse(nil, X) == X)
  assert(teq(X:offsets():table(), { 0, 1, 2 }))
  assert(teq(X:neighbors():table(), { 0, 0 }))
  local D = mtx.create({ n_rows = 1, n_cols = 1 })
  assert(X:to_dense(D) == D)
  assert(D:eq(M))
end)

test("mtx: maxs on all-negative rows (rmaxs init regression)", function ()
  local M = mtx.create({ data = dvec.create({ -3, -1, -2, -9, -6, -5 }), n_rows = 2, n_cols = 3 })
  assert(teq(M:maxs("row"):table(), { -1, -5 }))
  local I = mtx.create({ data = ivec.create({ -3, -1, -2, -9, -6, -5 }), n_rows = 2, n_cols = 3 })
  assert(teq(I:maxs("row"):table(), { -1, -5 }))
end)

test("mtx: bits transpose exact positions (port of cvec bits_transpose test)", function ()
  local csr = require("santoku.csr")
  local A = csr.create({
    offsets = ivec.create({ 0, 2, 4, 6 }),
    neighbors = ivec.create({ 0, 2, 1, 3, 0, 1 }),
    n_cols = 4,
  })
  local M = mtx.create({ data = A:to_bits(), n_rows = 3, n_cols = 4, bits = true })
  local T = M:transpose()
  local X = csr.from_bits(T:data(), 4, 3)
  local expected = csr.create({
    offsets = ivec.create({ 0, 2, 4, 5, 6 }),
    neighbors = ivec.create({ 0, 2, 1, 2, 0, 1 }),
    n_cols = 3,
  })
  assert(X:eq(expected))
end)
