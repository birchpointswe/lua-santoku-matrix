local test = require("santoku.test")

local err = require("santoku.error")
local assert = err.assert

local validate = require("santoku.validate")
local eq = validate.isequal

local tbl = require("santoku.table")
local teq = tbl.equals

local ivec = require("santoku.ivec")
local dvec = require("santoku.dvec")
local csr = require("santoku.csr")

test("typed vectors are zero indexed and mutate in place", function ()
  local v = ivec.create({ 3, 1, 2 })
  v:push(4)
  assert(eq(4, v:size()))
  assert(eq(3, v:get(0)))
  v:asc()
  assert(teq({ 1, 2, 3, 4 }, v:table()))
  assert(teq({ { 1, 2 }, { 3, 4 } }, v:rtable(2)))
end)

test("reductions run in c over the whole buffer", function ()
  local v = dvec.create({ 3, 7, 2, 9, 1 })
  assert(eq(22, v:sum()))
  local val, idx = v:max()
  assert(eq(9, val))
  assert(eq(3, idx))
  assert(eq(32, dvec.create({ 1, 2, 3 }):dot(dvec.create({ 4, 5, 6 }))))
end)

test("sorted integer vectors double as sets", function ()
  local a = ivec.create({ 1, 2, 3, 4 })
  local b = ivec.create({ 3, 4, 5, 6 })
  assert(teq({ 1, 2, 3, 4, 5, 6 }, a:set_union(b):table()))
  assert(teq({ 3, 4 }, a:set_intersect(b):table()))
  assert(eq(0.5, a:set_dice(b)))
end)

test("build a sparse matrix a row at a time", function ()
  local X = csr.create({ n_cols = 4, values = "f32" })
  X:push(0, 1.5):push(2, 2.5):row()
  X:row()
  X:push(3):row()
  local rows, cols = X:shape()
  assert(eq(3, rows))
  assert(eq(4, cols))
  assert(teq({ 0, 2, 2, 3 }, X:offsets():table()))
  assert(teq({ 0, 2, 3 }, X:neighbors():table()))
end)
