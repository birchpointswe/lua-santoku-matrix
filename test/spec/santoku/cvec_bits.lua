local test = require("santoku.test")
local err = require("santoku.error")
local assert = err.assert
local cvec = require("santoku.cvec")
local ivec = require("santoku.ivec")
local csr = require("santoku.csr")
require("santoku.mtx")
local tbl = require("santoku.table")
local teq = tbl.equals

local SETS = { { 0, 3, 9 }, { 1, 3, 8 }, { 2, 9 } }

local function rows_of (sets, n_bits)
  local out = {}
  for r = 1, #sets do
    local row = {}
    for c = 1, n_bits do
      row[c] = 0
    end
    for _, c in ipairs(sets[r]) do
      row[c + 1] = 1
    end
    out[r] = row
  end
  return out
end

local function dump (B, n_rows, n_bits)
  local out = {}
  for r = 0, n_rows - 1 do
    local row = {}
    for c = 0, n_bits - 1 do
      row[c + 1] = B:bits_get(r, c, n_bits)
    end
    out[r + 1] = row
  end
  return out
end

local function make ()
  local B = cvec.create(6)
  B:bits_fill(0)
  for r = 1, #SETS do
    for _, c in ipairs(SETS[r]) do
      B:bits_set(r - 1, c, 10, 1)
    end
  end
  return B
end

test("cvec bits: set and get", function ()
  local B = make()
  assert(B:size() == 6)
  assert(teq(dump(B, 3, 10), rows_of(SETS, 10)))
  assert(B:bits_set(2, 5, 10, true) == B)
  assert(B:bits_get(2, 5, 10) == 1)
  B:bits_set(2, 5, 10, false)
  assert(B:bits_get(2, 5, 10) == 0)
  B:bits_set(0, 3, 10, 0)
  assert(B:bits_get(0, 3, 10) == 0)
  assert(not pcall(function () B:bits_set(0, 0, 10, 2) end))
  assert(not pcall(function () B:bits_get(0, 10, 10) end))
  assert(not pcall(function () B:bits_get(3, 0, 10) end))
  assert(not pcall(function () B:bits_get(0, 0, 0) end))
end)

test("cvec bits: fill writes every byte, padding included", function ()
  local B = cvec.create(6)
  assert(B:bits_fill(1) == B)
  assert(B:bits_popcount() == 48)
  assert(teq(dump(B, 3, 10), rows_of({ { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 },
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } }, 10)))
  B:bits_fill(false)
  assert(B:bits_popcount() == 0)
  B:bits_fill(true)
  assert(B:bits_popcount() == 48)
  assert(not pcall(function () B:bits_fill(2) end))
end)

test("cvec bits: popcount", function ()
  assert(make():bits_popcount() == 8)
end)

test("cvec bits: rows", function ()
  local B = make()
  local R = B:bits_rows(ivec.create({ 2, 0 }), 10)
  assert(R:size() == 4)
  assert(teq(dump(R, 2, 10), rows_of({ SETS[3], SETS[1] }, 10)))
  assert(B:bits_rows(ivec.create(), 10):size() == 0)
  assert(not pcall(function () B:bits_rows(ivec.create({ 3 }), 10) end))
  assert(not pcall(function () cvec.create(5):bits_rows(ivec.create({ 0 }), 10) end))
end)

test("cvec bits: cols", function ()
  local B = make()
  local C = B:bits_cols(ivec.create({ 9, 3, 0 }), 10)
  assert(C:size() == 3)
  assert(teq(dump(C, 3, 3), { { 1, 1, 1 }, { 0, 1, 0 }, { 1, 0, 0 } }))
  assert(not pcall(function () B:bits_cols(ivec.create({ 10 }), 10) end))
  assert(not pcall(function () B:bits_cols(ivec.create(), 10) end))
end)

test("cvec bits: hcat", function ()
  local B = make()
  local C = B:bits_cols(ivec.create({ 9, 3, 0 }), 10)
  local H = B:bits_hcat(C, 10, 3)
  assert(H:size() == 6)
  assert(teq(dump(H, 3, 13), rows_of({ { 0, 3, 9, 10, 11, 12 }, { 1, 3, 8, 11 }, { 2, 9, 10 } }, 13)))
  local F = cvec.create(6):bits_fill(1)
  local G = F:bits_hcat(C, 10, 3)
  assert(teq(dump(G, 3, 13)[2], { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0 }))
  assert(not pcall(function () B:bits_hcat(cvec.create(2), 10, 3) end))
  assert(not pcall(function () B:bits_hcat(cvec.create(5), 10, 8) end))
end)

test("cvec bits: transpose", function ()
  local B = make()
  local T = B:bits_transpose(3, 10)
  assert(T:size() == 10)
  for r = 0, 2 do
    for c = 0, 9 do
      assert(T:bits_get(c, r, 3) == B:bits_get(r, c, 10))
    end
  end
  assert(not pcall(function () cvec.create(5):bits_transpose(3, 10) end))
end)

test("cvec bits: hamming per row, broadcast one row, padding ignored", function ()
  local B = make()
  local Q = B:bits_rows(ivec.create({ 1, 2, 0 }), 10)
  assert(teq(B:bits_hamming(Q, 10):table(), { 4, 5, 3 }))
  local Q1 = B:bits_rows(ivec.create({ 0 }), 10)
  assert(teq(B:bits_hamming(Q1, 10):table(), { 0, 4, 3 }))
  assert(not pcall(function () B:bits_hamming(B:bits_rows(ivec.create({ 0, 1 }), 10), 10) end))
  local X = cvec.create(6):bits_fill(1)
  local Z = cvec.create(6):bits_fill(0)
  for r = 0, 2 do
    for c = 0, 9 do
      Z:bits_set(r, c, 10, 1)
    end
  end
  assert(teq(X:bits_hamming(Z, 10):table(), { 0, 0, 0 }))
end)

test("cvec bits: topk ascending by distance, ties by ascending id", function ()
  local B = make()
  local X = B:bits_topk(B:bits_rows(ivec.create({ 0 }), 10), 10, 2)
  local r, c = X:shape()
  assert(r == 1 and c == 3)
  assert(teq(X:offsets():table(), { 0, 2 }))
  assert(teq(X:neighbors():table(), { 0, 2 }))
  assert(teq(X:values():table(), { 0, 3 }))
  local D = B:bits_rows(ivec.create({ 1, 1, 0, 1 }), 10)
  local Y = D:bits_topk(B:bits_rows(ivec.create({ 1, 0 }), 10), 10, 3)
  assert(teq(Y:offsets():table(), { 0, 3, 6 }))
  assert(teq(Y:neighbors():table(), { 0, 1, 3, 2, 0, 1 }))
  assert(teq(Y:values():table(), { 0, 0, 0, 0, 4, 4 }))
  local W = D:bits_topk(B:bits_rows(ivec.create({ 0 }), 10), 10, 9)
  assert(teq(W:neighbors():table(), { 2, 0, 1, 3 }))
  assert(not pcall(function () B:bits_topk(cvec.create(3), 10, 1) end))
end)

test("cvec bits: and, or, xor, andnot in place", function ()
  local B = make()
  local C = B:bits_rows(ivec.create({ 1, 2, 0 }), 10)
  local A = B:bits_rows(ivec.create({ 0, 1, 2 }), 10)
  assert(A:bits_and(C) == A)
  assert(teq(dump(A, 3, 10), rows_of({ { 3 }, {}, { 9 } }, 10)))
  A = B:bits_rows(ivec.create({ 0, 1, 2 }), 10)
  assert(A:bits_or(C) == A)
  assert(teq(dump(A, 3, 10), rows_of({ { 0, 1, 3, 8, 9 }, { 1, 2, 3, 8, 9 }, { 0, 2, 3, 9 } }, 10)))
  A = B:bits_rows(ivec.create({ 0, 1, 2 }), 10)
  assert(A:bits_xor(C) == A)
  assert(teq(dump(A, 3, 10), rows_of({ { 0, 1, 8, 9 }, { 1, 2, 3, 8, 9 }, { 0, 2, 3 } }, 10)))
  A = B:bits_rows(ivec.create({ 0, 1, 2 }), 10)
  assert(A:bits_andnot(C) == A)
  assert(teq(dump(A, 3, 10), rows_of({ { 0, 9 }, { 1, 3, 8 }, { 2 } }, 10)))
  assert(not pcall(function () A:bits_and(cvec.create(4)) end))
end)

test("cvec bits: to_dense", function ()
  local B = make()
  local M = B:bits_to_dense(3, 10)
  local r, c = M:shape()
  assert(r == 3 and c == 10)
  assert(M:type() == "f32")
  for i = 0, 2 do
    for j = 0, 9 do
      assert(M:get(i, j) == B:bits_get(i, j, 10))
    end
  end
  assert(not pcall(function () cvec.create(5):bits_to_dense(3, 10) end))
end)

test("cvec bits: csr to_bits then bits_to_dense matches csr to_dense", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 3, 6, 8 }),
    neighbors = ivec.create({ 0, 3, 9, 1, 3, 8, 2, 9 }),
    n_cols = 10,
  })
  local B = X:to_bits()
  assert(B:size() == 6)
  assert(B:bits_to_dense(3, 10):eq(X:to_dense("f32")))
  local F = X:to_bits(true)
  assert(F:size() == 9)
  local M = F:bits_to_dense(3, 20)
  local E = make()
  for i = 0, 2 do
    for j = 0, 9 do
      local v = E:bits_get(i, j, 10)
      assert(M:get(i, j) == v)
      assert(M:get(i, j + 10) == 1 - v)
    end
  end
end)

test("cvec bits: csr.from_bits after bits_set", function ()
  local B = make()
  local Y = csr.from_bits(B, 3, 10)
  local r, c = Y:shape()
  assert(r == 3 and c == 10)
  assert(teq(Y:offsets():table(), { 0, 3, 6, 8 }))
  assert(teq(Y:neighbors():table(), { 0, 3, 9, 1, 3, 8, 2, 9 }))
end)

test("cvec bits: csr.from_bits ignores set padding bits", function ()
  local B = cvec.create(6)
  B:bits_fill(1)
  local Y = csr.from_bits(B, 3, 10)
  assert(Y:nnz() == 30)
  assert(Y:neighbors():size() == 30)
  assert(teq(Y:offsets():table(), { 0, 10, 20, 30 }))
end)
