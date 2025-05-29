local serialize  = require("santoku.serialize") -- luacheck: ignore
local it  = require("santoku.iter")
local ivec = require("santoku.ivec")
local dvec = require("santoku.dvec")
local tbl = require("santoku.table")

for _, vec in ipairs({ ivec, dvec }) do

  local m0, m1, m2

  m0 = vec.create(10)

  for i = 0, 4 do
    vec.set(m0, i, i)
  end

  for i = 0, 4 do
    assert(vec.get(m0, i) == i)
  end

  vec.add(m0, -1)

  for i = 0, 4 do
    assert(vec.get(m0, i) == i - 1)
  end

  m0 = vec.create({ 1, 2, 3, 4, 5, 6 })
  assert(vec.size(m0) == 6)
  assert(vec.get(m0, 1) == 2)
  assert(vec.get(m0, 2) == 3)
  assert(vec.get(m0, 3) == 4)

  m0 = vec.create({ 1, 2, 3 })
  m1 = vec.create({ 4, 5, 6 })
  assert(vec.dot(m0, m1) == 32)

  m0 = vec.create({ 1, 2, 3 })
  m1 = vec.create({ 4, 5, 6 })

  assert(vec.size(m0) == 3)
  assert(vec.size(m1) == 3)

  m2 = vec.from_raw(vec.raw(m0), 3)

  assert(vec.size(m2) == 3)

  m0 = vec.create({ 1, 2, 3, 4, 5, 6 })

  assert(tbl.equals(vec.table(m0), {
    1, 2, 3, 4, 5, 6
  }))

  assert(tbl.equals(vec.rtable(m0, 3), {
    { 1, 2, 3 },
    { 4, 5, 6 }
  }))

  m0 = vec.create({
    10, 2, 3,
    4, 5, 6
  })

  assert(tbl.equals(ivec.rtable(vec.rdesc(m0, 3), 3), {
    { 0, 2, 1 },
    { 2, 1, 0 }
  }))

  assert(tbl.equals(ivec.rtable(vec.casc(m0, 3), 2), {
    { 1, 0 },
    { 0, 1 },
    { 0, 1 }
  }))

  m0 = vec.create({ 1, 6, 3, 2, 5, 2, 3, 6, 1, 4 })
  vec.asc(m0)
  assert(tbl.equals(vec.table(m0), { 1, 1, 2, 2, 3, 3, 4, 5, 6, 6 }))

  if vec.flip_interleave then
    m0 = vec.create({ 0, 3 })
    vec.flip_interleave(m0, 2, 2)
  end

  m0 = vec.create({ 1, 2, 3, 4 })
  local g = vec.each(m0)
  assert(g() == 1)
  assert(g() == 2)
  assert(g() == 3)
  assert(g() == 4)
  assert(g() == nil)

  g = it.take(2, vec.each(vec.create({ 1, 2, 3, 4 })))
  assert(g() == 1)
  assert(g() == 2)
  assert(g() == nil)

  if vec.flip_interleave then
    m0 = vec.create({ 0 })
    m1 = vec.create({ 1 })
    vec.add(m1, 2)
    vec.copy(m0, m1)
    vec.flip_interleave(m0, 2, 2)
  end

end
