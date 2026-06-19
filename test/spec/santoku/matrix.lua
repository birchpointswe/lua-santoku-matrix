local test = require("santoku.test")
local err = require("santoku.error")
local assert = err.assert
local dvec = require("santoku.dvec")
local ivec = require("santoku.ivec")
local svec = require("santoku.svec")
local rvec = require("santoku.rvec")
local pvec = require("santoku.pvec")
local tbl = require("santoku.table")
local teq = tbl.equals

test("ivec/dvec/svec: create and basic operations", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create(10)
    assert(v:size() == 10)
    assert(v:capacity() >= 10)

    v:set(0, 100)
    v:set(5, 200)
    assert(v:get(0) == 100)
    assert(v:get(5) == 200)

    v:clear()
    assert(v:size() == 0)
  end
end)

test("ivec/dvec/svec: create from table", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    assert(v:size() == 5)
    assert(v:get(0) == 1)
    assert(v:get(2) == 3)
    assert(v:get(4) == 5)
  end
end)

test("ivec/dvec/svec: resize and setn", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3 })
    v:resize(10)
    assert(v:size() == 10)
    assert(v:get(0) == 1)

    v:setn(5)
    assert(v:size() == 5)
  end
end)

test("ivec/dvec/svec: push", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create()
    v:push(10)
    v:push(20)
    v:push(30)
    assert(v:size() == 3)
    assert(v:get(0) == 10)
    assert(v:get(1) == 20)
    assert(v:get(2) == 30)
  end
end)

test("ivec/dvec/svec: insert", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 3, 4 })
    v:insert(1, 2)
    assert(teq(v:table(), { 1, 2, 3, 4 }))
  end
end)

test("ivec/dvec/svec: copy", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local src = vec.create({ 1, 2, 3, 4, 5 })
    local dst = vec.create()
    dst:copy(src)
    assert(teq(dst:table(), { 1, 2, 3, 4, 5 }))

    dst:clear()
    dst:copy(src, 0)
    assert(teq(dst:table(), { 1, 2, 3, 4, 5 }))

    dst:clear()
    dst:copy(src, 1, 4, 0)
    assert(teq(dst:table(), { 2, 3, 4 }))
  end
end)

test("ivec/dvec/svec: reverse", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    v:reverse()
    assert(teq(v:table(), { 5, 4, 3, 2, 1 }))

    v = vec.create({ 1, 2, 3, 4, 5 })
    v:reverse(1, 4)
    assert(teq(v:table(), { 1, 4, 3, 2, 5 }))
  end
end)

test("ivec/dvec/svec: shuffle (whole)", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 })
    v:shuffle()
    local shuffled = v:table()
    assert(#shuffled == 10)
    local sorted = {}
    for i = 1, 10 do sorted[i] = shuffled[i] end
    table.sort(sorted)
    assert(teq(sorted, { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }))
  end
end)

test("ivec/dvec/svec: shuffle with range", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5, 6 })
    v:shuffle(2, 5)
    assert(v:get(0) == 1)
    assert(v:get(1) == 2)
    assert(v:get(5) == 6)
    local vals = { v:get(2), v:get(3), v:get(4) }
    table.sort(vals)
    assert(teq(vals, { 3, 4, 5 }))
  end
end)

test("ivec/dvec/svec: clear with range", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    v:clear(1, 4)
    assert(v:get(0) == 1)
    assert(v:get(1) == 0)
    assert(v:get(2) == 0)
    assert(v:get(3) == 0)
    assert(v:get(4) == 5)
  end
end)

test("ivec/dvec/svec: zero", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    v:zero()
    for i = 0, 4 do
      assert(v:get(i) == 0)
    end
  end
end)

test("ivec/dvec/svec: fill", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create(5)
    v:fill(42)
    for i = 0, 4 do
      assert(v:get(i) == 42)
    end

    v:fill(99, 1, 4)
    assert(v:get(0) == 42)
    assert(v:get(1) == 99)
    assert(v:get(2) == 99)
    assert(v:get(3) == 99)
    assert(v:get(4) == 42)
  end
end)

test("ivec/dvec/svec: asc sort", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 5, 2, 8, 1, 9 })
    v:asc()
    assert(teq(v:table(), { 1, 2, 5, 8, 9 }))
  end
end)

test("ivec/dvec/svec: desc sort", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 5, 2, 8, 1, 9 })
    v:desc()
    assert(teq(v:table(), { 9, 8, 5, 2, 1 }))
  end
end)

test("ivec/dvec/svec: asc with range", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 5, 2, 8, 9 })
    v:asc(1, 4)
    assert(v:get(0) == 1)
    assert(v:get(1) == 2)
    assert(v:get(2) == 5)
    assert(v:get(3) == 8)
    assert(v:get(4) == 9)
  end
end)

test("ivec/dvec/svec: uasc (unique ascending)", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 3, 1, 2, 1, 3, 2 })
    local new_end = v:uasc()
    v:setn(new_end)
    assert(teq(v:table(), { 1, 2, 3 }))
  end
end)

test("ivec/dvec/svec: kasc (partial sort)", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 9, 5, 2, 8, 1, 7 })
    v:kasc(3)
    assert(v:get(0) <= v:get(1))
    assert(v:get(1) <= v:get(2))
  end
end)

test("ivec/dvec/svec: table conversion", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5, 6 })
    assert(teq(v:table(), { 1, 2, 3, 4, 5, 6 }))
    assert(teq(v:table(0, 3), { 1, 2, 3 }))
  end
end)

test("ivec/dvec/svec: rtable (row-major)", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5, 6 })
    assert(teq(v:rtable(3), {
      { 1, 2, 3 },
      { 4, 5, 6 }
    }))
  end
end)

test("ivec/dvec/svec: ctable (column-major)", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5, 6 })
    assert(teq(v:ctable(3), {
      { 1, 4 },
      { 2, 5 },
      { 3, 6 }
    }))
  end
end)

test("ivec/dvec/svec: add scalar", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    v:add(10)
    assert(teq(v:table(), { 11, 12, 13, 14, 15 }))
  end
end)

test("ivec/dvec/svec: add scalar with range", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    v:add(10, 1, 4)
    assert(teq(v:table(), { 1, 12, 13, 14, 5 }))
  end
end)

test("ivec/dvec/svec: scale", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    v:scale(2)
    assert(teq(v:table(), { 2, 4, 6, 8, 10 }))
  end
end)

test("ivec/dvec/svec: scalev", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v1 = vec.create({ 1, 2, 3 })
    local v2 = vec.create({ 2, 3, 4 })
    v1:scalev(v2)
    assert(teq(v1:table(), { 2, 6, 12 }))
  end
end)

test("ivec/dvec/svec: addv", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v1 = vec.create({ 1, 2, 3 })
    local v2 = vec.create({ 10, 20, 30 })
    v1:addv(v2)
    assert(teq(v1:table(), { 11, 22, 33 }))
  end
end)

test("ivec/dvec/svec: sum", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    assert(v:sum() == 15)
  end
end)

test("ivec/dvec/svec: max", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 3, 7, 2, 9, 1 })
    local val, idx = v:max()
    assert(val == 9)
    assert(idx == 3)
  end
end)

test("ivec/dvec/svec: min", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 3, 7, 2, 9, 1 })
    local val, idx = v:min()
    assert(val == 1)
    assert(idx == 4)
  end
end)

test("ivec/dvec/svec: dot product", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v1 = vec.create({ 1, 2, 3 })
    local v2 = vec.create({ 4, 5, 6 })
    assert(v1:dot(v2) == 32)
  end
end)

test("ivec/dvec/svec: magnitude", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 3, 4 })
    assert(math.abs(v:magnitude() - 5) < 1e-10)
  end
end)

test("ivec/dvec/svec: each iterator", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4 })
    local sum = 0
    for val in v:each() do
      sum = sum + val
    end
    assert(sum == 10)
  end
end)

test("ivec/dvec/svec: ieach iterator", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 10, 20, 30 })
    local indices = {}
    local values = {}
    for i, val in v:ieach() do
      table.insert(indices, i)
      table.insert(values, val)
    end
    assert(teq(indices, { 0, 1, 2 }))
    assert(teq(values, { 10, 20, 30 }))
  end
end)

test("ivec/dvec/svec: find", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 10, 20, 30, 40 })
    assert(v:find(30) == 2)
    assert(v:find(99) == nil)
  end
end)









test("dvec: exp", function ()
  local v = dvec.create({ 0, 1, 2 })
  v:exp()
  assert(math.abs(v:get(0) - 1) < 1e-10)
  assert(math.abs(v:get(1) - math.exp(1)) < 1e-10)
  assert(math.abs(v:get(2) - math.exp(2)) < 1e-10)
end)

test("dvec: log", function ()
  local v = dvec.create({ 1, math.exp(1), math.exp(2) })
  v:log()
  assert(math.abs(v:get(0) - 0) < 1e-10)
  assert(math.abs(v:get(1) - 1) < 1e-10)
  assert(math.abs(v:get(2) - 2) < 1e-10)
end)

test("dvec: pow", function ()
  local v = dvec.create({ 1, 2, 3 })
  v:pow(2)
  assert(teq(v:table(), { 1, 4, 9 }))
end)

test("dvec: abs", function ()
  local v = dvec.create({ -1, 2, -3 })
  v:abs()
  assert(teq(v:table(), { 1, 2, 3 }))
end)





test("dvec: round", function ()
  local v = dvec.create({ 1.4, 1.5, 1.6, -1.4, -1.5, -1.6 })
  v:round()
  assert(v:get(0) == 1)
  assert(v:get(1) == 2)
  assert(v:get(2) == 2)
  assert(v:get(3) == -1)
  assert(v:get(4) == -2)
  assert(v:get(5) == -2)
end)

test("dvec: round with range", function ()
  local v = dvec.create({ 1.4, 2.5, 3.6, 4.7 })
  v:round(1, 3)
  assert(v:get(0) == 1.4)
  assert(v:get(1) == 3)
  assert(v:get(2) == 4)
  assert(v:get(3) == 4.7)
end)

test("dvec: trunc", function ()
  local v = dvec.create({ 1.9, -1.9, 2.1, -2.1 })
  v:trunc()
  assert(v:get(0) == 1)
  assert(v:get(1) == -1)
  assert(v:get(2) == 2)
  assert(v:get(3) == -2)
end)

test("dvec: floor", function ()
  local v = dvec.create({ 1.9, -1.1, 2.0 })
  v:floor()
  assert(v:get(0) == 1)
  assert(v:get(1) == -2)
  assert(v:get(2) == 2)
end)

test("dvec: ceil", function ()
  local v = dvec.create({ 1.1, -1.9, 2.0 })
  v:ceil()
  assert(v:get(0) == 2)
  assert(v:get(1) == -1)
  assert(v:get(2) == 2)
end)

test("dvec: to_ivec", function ()
  local d = dvec.create({ 1.0, 2.5, 3.9, -4.2 })
  local i = d:to_ivec()
  assert(i:get(0) == 1)
  assert(i:get(1) == 2)
  assert(i:get(2) == 3)
  assert(i:get(3) == -4)
end)

test("dvec: round and to_ivec chained", function ()
  local d = dvec.create({ 1.4, 2.5, 3.6 })
  local i = d:round():to_ivec()
  assert(i:get(0) == 1)
  assert(i:get(1) == 3)
  assert(i:get(2) == 4)
end)

test("ivec: to_dvec", function ()
  local i = ivec.create({ 1, 2, 3, -4 })
  local d = i:to_dvec()
  assert(d:get(0) == 1.0)
  assert(d:get(1) == 2.0)
  assert(d:get(2) == 3.0)
  assert(d:get(3) == -4.0)
end)

test("rvec: create and basic operations", function ()
  local v = rvec.create()
  v:push(1, 10.5)
  v:push(2, 20.5)
  v:push(3, 30.5)
  assert(v:size() == 3)
  local i, d = v:get(0)
  assert(i == 1 and d == 10.5)
end)

test("rvec: heap operations (hmax)", function ()
  local heap = rvec.create(5)
  heap:setn(0)
  heap:hmax(1, 100, 5)
  heap:hmax(2, 50, 5)
  heap:hmax(3, 200, 5)
  heap:hmax(4, 25, 5)
  heap:hmax(5, 150, 5)
  heap:hmax(6, 300, 5)
  heap:hmax(7, 75, 5)
  assert(heap:size() == 5)
  heap:asc()
  local _, d1 = heap:get(0)
  local _, d5 = heap:get(4)
  assert(d1 <= d5)
end)

test("rvec: sort ascending", function ()
  local v = rvec.create()
  v:push(1, 50)
  v:push(2, 10)
  v:push(3, 30)
  v:asc()
  local _, d1 = v:get(0)
  local _, d2 = v:get(1)
  local _, d3 = v:get(2)
  assert(d1 <= d2 and d2 <= d3)
end)

test("pvec: create and basic operations", function ()
  local v = pvec.create()
  v:push(1, 100)
  v:push(2, 200)
  assert(v:size() == 2)
  local i, p = v:get(0)
  assert(i == 1 and p == 100)
end)

test("ivec: set operations (jaccard)", function ()
  local v0 = ivec.create({ 1, 2, 3, 4 })
  local v1 = ivec.create({ 3, 4, 5, 6 })
  local j = v0:set_jaccard(v1)
  assert(math.abs(j - 1/3) < 1e-10)
end)

test("ivec: set operations (overlap)", function ()
  local v0 = ivec.create({ 1, 2, 3, 4 })
  local v1 = ivec.create({ 3, 4, 5, 6 })
  assert(v0:set_overlap(v1) == 0.5)
end)

test("ivec: set operations (dice)", function ()
  local v0 = ivec.create({ 1, 2, 3, 4 })
  local v1 = ivec.create({ 3, 4, 5, 6 })
  assert(v0:set_dice(v1) == 0.5)
end)

test("ivec: set operations (tversky)", function ()
  local v0 = ivec.create({ 1, 2, 3, 4 })
  local v1 = ivec.create({ 3, 4, 5, 6 })
  assert(v0:set_tversky(v1, 1, 0) == 0.5)
  assert(v0:set_tversky(v1, 0, 1) == 0.5)
end)

test("ivec: set operations (union)", function ()
  local v0 = ivec.create({ 1, 2, 3, 4 })
  local v1 = ivec.create({ 3, 4, 5, 6 })
  local u = v0:set_union(v1)
  assert(teq({ 1, 2, 3, 4, 5, 6 }, u:table()))
end)

test("ivec: set operations (intersect)", function ()
  local v0 = ivec.create({ 1, 2, 3, 4 })
  local v1 = ivec.create({ 3, 4, 5, 6 })
  local i = v0:set_intersect(v1)
  assert(teq({ 3, 4 }, i:table()))
end)

test("ivec: lookup", function ()
  local indices = ivec.create({ 2, 0, 1, 2 })
  local source = ivec.create({ 100, 200, 300 })
  indices:lookup(source)
  assert(teq(indices:table(), { 300, 100, 200, 300 }))
end)


test("ivec/svec: index (value -> position iumap)", function ()
  for _, vec in ipairs({ ivec, svec }) do
    local v = vec.create({ 10, 20, 30 })
    local m = v:index()
    assert(m:val(m:get(10)) == 0)
    assert(m:val(m:get(20)) == 1)
    assert(m:val(m:get(30)) == 2)
  end
end)

test("ivec/dvec/svec: add_scaled (a[i] += x * i)", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 10, 10, 10 })
    v:add_scaled(1)
    assert(teq(v:table(), { 10, 11, 12 }))
    local w = vec.create({ 5, 5, 5, 5 })
    w:add_scaled(2, 1, 3)
    assert(teq(w:table(), { 5, 7, 9, 5 }))
  end
end)

test("ivec/dvec/svec: eq", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local a = vec.create({ 1, 2, 3 })
    assert(a:eq(vec.create({ 1, 2, 3 })))
    assert(not a:eq(vec.create({ 1, 2, 4 })))
    assert(not a:eq(vec.create({ 1, 2 })))
    assert(a:eq(vec.create({ 9, 2, 3 }), 1, 3))
    assert(not a:eq(vec.create({ 9, 2, 4 }), 1, 3))
  end
end)

test("dvec: eq with eps", function ()
  local a = dvec.create({ 1.0, 2.0 })
  local b = dvec.create({ 1.0005, 1.9995 })
  assert(a:eq(b, 1e-3))
  assert(not a:eq(b, 1e-6))
  assert(a:eq(b, 0, 2, 1e-3))
end)

test("pvec: eq (pair compare)", function ()
  local a = pvec.create()
  a:push(1, 10)
  a:push(2, 20)
  local b = pvec.create()
  b:push(1, 10)
  b:push(2, 20)
  assert(a:eq(b))
  local c = pvec.create()
  c:push(1, 10)
  c:push(2, 21)
  assert(not a:eq(c))
end)

test("ivec/dvec/svec: where", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 5, 0, 7, 0, 3 })
    assert(teq(v:where("gt", 0):table(), { 0, 2, 4 }))
    assert(teq(v:where("eq", 0):table(), { 1, 3 }))
    assert(teq(v:where("ne", 0):table(), { 0, 2, 4 }))
    assert(teq(v:where("le", 3):table(), { 1, 3, 4 }))
    assert(teq(v:where("gt", 0, 1, 4):table(), { 2 }))
  end
end)

test("ivec: where composes with copy gather", function ()
  local v = ivec.create({ 5, 0, 7, 0, 3 })
  local idx = v:where("gt", 0)
  local picked = ivec.create(idx:size())
  picked:copy(v, idx)
  assert(teq(picked:table(), { 5, 7, 3 }))
end)

test("ivec/dvec/svec: fill_segments", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create(5)
    v:zero()
    local off = ivec.create({ 0, 2, 5 })
    local vals = vec.create({ 7, 9 })
    v:fill_segments(off, vals)
    assert(teq(v:table(), { 7, 7, 9, 9, 9 }))
  end
end)

test("ivec/dvec/svec: persist and load", function ()
  local tmp = ".vec_test.bin"
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3, 4 })
    v:persist(tmp)
    local loaded = vec.load(tmp)
    assert(teq(v:table(), loaded:table()))
  end
  os.remove(tmp)
end)

test("ivec/dvec/svec: raw", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create({ 1, 2, 3 })
    local raw = v:raw()
    assert(type(raw) == "string")
    assert(#raw > 0)
  end
end)

test("ivec/dvec/svec: shrink", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create(100)
    v:setn(10)
    v:shrink()
    assert(v:size() == 10)
  end
end)

test("ivec/dvec/svec: ensure", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local v = vec.create()
    v:ensure(100)
    assert(v:capacity() >= 100)
  end
end)

test("iumap: create and operations", function ()
  local iumap = require("santoku.iumap")
  local m = iumap.create(0)
  m:put(10)
  m:setval(m:get(10), 100)
  m:put(20)
  m:setval(m:get(20), 200)
  local count = 0
  for _ in m:each() do
    count = count + 1
  end
  assert(count == 2)
  m:destroy()
end)



test("dvec: fill_indices", function ()
  local v = dvec.create(5)
  v:fill_indices()
  assert(teq(v:table(), { 0, 1, 2, 3, 4 }))
end)

test("ivec: fill_indices", function ()
  local v = ivec.create(5)
  v:fill_indices()
  assert(teq(v:table(), { 0, 1, 2, 3, 4 }))
end)













test("ivec: set_find (binary search)", function ()
  local v = ivec.create({ 10, 20, 30, 40, 50 })
  assert(v:set_find(30) == 2)
  assert(v:set_find(25) == -3)
  assert(v:set_find(10) == 0)
  assert(v:set_find(50) == 4)
  assert(v:set_find(5) == -1)
  assert(v:set_find(55) == -6)
end)

test("ivec: set_insert (sorted insert)", function ()
  local v = ivec.create({ 10, 30, 50 })
  v:set_insert(1, 20)
  assert(teq(v:table(), { 10, 20, 30, 50 }))
  v:set_insert(3, 40)
  assert(teq(v:table(), { 10, 20, 30, 40, 50 }))
end)

test("ivec/dvec/svec: copy_indexed", function ()
  for _, vec in ipairs({ ivec, dvec, svec }) do
    local src = vec.create({ 100, 200, 300, 400, 500 })
    local dst = vec.create(3)
    local indices = ivec.create({ 4, 0, 2 })
    dst:copy(src, indices)
    assert(teq(dst:table(), { 500, 100, 300 }))
  end
end)






test("svec: set operations (jaccard)", function ()
  local v0 = svec.create({ 1, 2, 3, 4 })
  local v1 = svec.create({ 3, 4, 5, 6 })
  local j = v0:set_jaccard(v1)
  assert(math.abs(j - 1/3) < 1e-10)
end)

test("svec: set operations (overlap)", function ()
  local v0 = svec.create({ 1, 2, 3, 4 })
  local v1 = svec.create({ 3, 4, 5, 6 })
  assert(v0:set_overlap(v1) == 0.5)
end)

test("svec: set operations (dice)", function ()
  local v0 = svec.create({ 1, 2, 3, 4 })
  local v1 = svec.create({ 3, 4, 5, 6 })
  assert(v0:set_dice(v1) == 0.5)
end)

test("svec: set operations (tversky)", function ()
  local v0 = svec.create({ 1, 2, 3, 4 })
  local v1 = svec.create({ 3, 4, 5, 6 })
  assert(v0:set_tversky(v1, 1, 0) == 0.5)
  assert(v0:set_tversky(v1, 0, 1) == 0.5)
end)

test("svec: set operations (union)", function ()
  local v0 = svec.create({ 1, 2, 3, 4 })
  local v1 = svec.create({ 3, 4, 5, 6 })
  local u = v0:set_union(v1)
  assert(teq({ 1, 2, 3, 4, 5, 6 }, u:table()))
end)

test("svec: set operations (intersect)", function ()
  local v0 = svec.create({ 1, 2, 3, 4 })
  local v1 = svec.create({ 3, 4, 5, 6 })
  local i = v0:set_intersect(v1)
  assert(teq({ 3, 4 }, i:table()))
end)

test("svec: lookup", function ()
  local indices = svec.create({ 2, 0, 1, 2 })
  local source = svec.create({ 100, 200, 300 })
  indices:lookup(source)
  assert(teq(indices:table(), { 300, 100, 200, 300 }))
end)


test("svec: set_find (binary search)", function ()
  local v = svec.create({ 10, 20, 30, 40, 50 })
  assert(v:set_find(30) == 2)
  assert(v:set_find(25) == -3)
  assert(v:set_find(10) == 0)
  assert(v:set_find(50) == 4)
  assert(v:set_find(5) == -1)
  assert(v:set_find(55) == -6)
end)

test("svec: set_insert (sorted insert)", function ()
  local v = svec.create({ 10, 30, 50 })
  v:set_insert(1, 20)
  assert(teq(v:table(), { 10, 20, 30, 50 }))
  v:set_insert(3, 40)
  assert(teq(v:table(), { 10, 20, 30, 40, 50 }))
end)

test("svec: to_ivec", function ()
  local s = svec.create({ 1, 2, 3, -4 })
  local i = s:to_ivec()
  assert(i:get(0) == 1)
  assert(i:get(1) == 2)
  assert(i:get(2) == 3)
  assert(i:get(3) == -4)
end)

test("svec: to_dvec", function ()
  local s = svec.create({ 1, 2, 3, -4 })
  local d = s:to_dvec()
  assert(d:get(0) == 1.0)
  assert(d:get(1) == 2.0)
  assert(d:get(2) == 3.0)
  assert(d:get(3) == -4.0)
end)

test("ivec: to_svec", function ()
  local i = ivec.create({ 1, 2, 3, -4 })
  local s = i:to_svec()
  assert(s:get(0) == 1)
  assert(s:get(1) == 2)
  assert(s:get(2) == 3)
  assert(s:get(3) == -4)
end)

test("svec: fill_indices", function ()
  local v = svec.create(5)
  v:fill_indices()
  assert(teq(v:table(), { 0, 1, 2, 3, 4 }))
end)




test("ivec: bincount", function ()
  local v = ivec.create({ 0, 2, 2, 1, 2 })
  assert(teq(v:bincount(3):table(), { 1, 1, 3 }))
  local w = dvec.create({ 1, 0.5, 0.5, 2, 1 })
  local b = v:bincount(3, w)
  assert(math.abs(b:get(0) - 1) < 1e-10)
  assert(math.abs(b:get(1) - 2) < 1e-10)
  assert(math.abs(b:get(2) - 2) < 1e-10)
end)
