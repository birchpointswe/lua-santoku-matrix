local serialize  = require("santoku.serialize") -- luacheck: ignore
local it  = require("santoku.iter")
local imtx = require("santoku.matrix.number")
local nmtx = require("santoku.matrix.integer")
local tbl = require("santoku.table")

for _, mtx in ipairs({ nmtx, imtx }) do

  local m0, m1, m2

  m0 = mtx.create(1, 10)

  for i = 1, 5 do
    mtx.set(m0, 1, i, i)
  end

  for i = 1, 5 do
    assert(mtx.get(m0, 1, i) == i)
  end

  mtx.add(m0, -1)

  for i = 1, 5 do
    assert(mtx.get(m0, 1, i) == i - 1)
  end

  m0 = mtx.create({ 1, 2, 3, 4, 5, 6 }, 2, 4)
  assert(mtx.rows(m0) == 1)
  assert(mtx.columns(m0) == 3)
  assert(mtx.get(m0, 1, 1) == 2)
  assert(mtx.get(m0, 1, 2) == 3)
  assert(mtx.get(m0, 1, 3) == 4)

  m0 = mtx.create({ 1, 2, 3, 4, 5, 6 }, 1, 3)
  m1 = mtx.create({ 1, 2, 3, 4, 5, 6 }, 4, 6)

  assert(mtx.dot(m0, m1) == 32)

  m0 = mtx.create({ 1, 2, 3 })
  m1 = mtx.create({ 4, 5, 6 })

  assert(mtx.rows(m0) == 1)
  assert(mtx.columns(m0) == 3)
  assert(mtx.rows(m1) == 1)
  assert(mtx.columns(m1) == 3)

  m2 = mtx.from_raw(mtx.raw(m0), 3)

  assert(mtx.rows(m2) == 1)
  assert(mtx.columns(m2) == 3)

  m0 = mtx.create({
    { 1, 2, 3 },
    { 4, 5, 6 }
  })

  assert(tbl.equals(mtx.tabulate(m0), {
    { 1, 2, 3 },
    { 4, 5, 6 }
  }))

  m0 = mtx.create({
    { 10, 2, 3 },
    { 4, 5, 6 }
  })

  assert(tbl.equals(mtx.rorder(m0, 1, 3, 0), {
    { 1, 3, 2 },
    { 3, 2, 1 }
  }))

  assert(tbl.equals(mtx.rorder(m0, 1, 2, 4), {
    { 1 },
    { 3, 2 }
  }))

  m0 = mtx.create({ 1, 6, 3, 2, 5, 2, 3, 6, 1, 4 })
  mtx.sort(m0)
  assert(tbl.equals(mtx.tabulate(m0), { { 1, 1, 2, 2, 3, 3, 4, 5, 6, 6 } }))

  m0 = mtx.create({ 1, 6, 3, 2, 5, 2, 3, 6, 1, 4 })
  mtx.sort(m0, false)
  assert(tbl.equals(mtx.tabulate(m0), { { 1, 1, 2, 2, 3, 3, 4, 5, 6, 6 } }))

  m0 = mtx.create({ 1, 6, 3, 2, 5, 2, 3, 6, 1, 4 })
  mtx.sort(m0, true)
  assert(tbl.equals(mtx.tabulate(m0), { { 1, 2, 3, 4, 5, 6 } }))

  if mtx.flip_interleave then
    m0 = mtx.create({ 0, 3 })
    mtx.flip_interleave(m0, 2, 2)
  end

  m0 = mtx.create({ 1, 2, 3, 4 })
  local g = mtx.each(m0)
  assert(g() == 1)
  assert(g() == 2)
  assert(g() == 3)
  assert(g() == 4)
  assert(g() == nil)

  g = it.take(2, mtx.each(mtx.create({ 1, 2, 3, 4 })))
  assert(g() == 1)
  assert(g() == 2)
  assert(g() == nil)

  if mtx.flip_interleave then
    m0 = mtx.create({ 0 })
    m1 = mtx.create({ 1 })
    mtx.add(m1, 2)
    mtx.extend(m0, m1)
    mtx.flip_interleave(m0, 2, 2)
  end

end
