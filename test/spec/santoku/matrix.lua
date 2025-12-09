local test = require("santoku.test")
local err = require("santoku.error")
local assert = err.assert
local dvec = require("santoku.dvec")
local ivec = require("santoku.ivec")
local rvec = require("santoku.rvec")
local pvec = require("santoku.pvec")
local tbl = require("santoku.table")
local teq = tbl.equals

test("ivec/dvec: create and basic operations", function ()
  for _, vec in ipairs({ ivec, dvec }) do
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

test("ivec/dvec: create from table", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    assert(v:size() == 5)
    assert(v:get(0) == 1)
    assert(v:get(2) == 3)
    assert(v:get(4) == 5)
  end
end)

test("ivec/dvec: resize and setn", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3 })
    v:resize(10)
    assert(v:size() == 10)
    assert(v:get(0) == 1)

    v:setn(5)
    assert(v:size() == 5)
  end
end)

test("ivec/dvec: push", function ()
  for _, vec in ipairs({ ivec, dvec }) do
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

test("ivec/dvec: insert", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 3, 4 })
    v:insert(1, 2)
    assert(teq(v:table(), { 1, 2, 3, 4 }))
  end
end)

test("ivec/dvec: copy", function ()
  for _, vec in ipairs({ ivec, dvec }) do
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

test("ivec/dvec: reverse", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    v:reverse()
    assert(teq(v:table(), { 5, 4, 3, 2, 1 }))

    v = vec.create({ 1, 2, 3, 4, 5 })
    v:reverse(1, 4)
    assert(teq(v:table(), { 1, 4, 3, 2, 5 }))
  end
end)

test("ivec/dvec: shuffle (whole)", function ()
  for _, vec in ipairs({ ivec, dvec }) do
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

test("ivec/dvec: shuffle with range", function ()
  for _, vec in ipairs({ ivec, dvec }) do
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

test("ivec/dvec: clear with range", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    v:clear(1, 4)
    assert(v:get(0) == 1)
    assert(v:get(1) == 0)
    assert(v:get(2) == 0)
    assert(v:get(3) == 0)
    assert(v:get(4) == 5)
  end
end)

test("ivec/dvec: zero", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    v:zero()
    for i = 0, 4 do
      assert(v:get(i) == 0)
    end
  end
end)

test("ivec/dvec: fill", function ()
  for _, vec in ipairs({ ivec, dvec }) do
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

test("ivec/dvec: asc sort", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 5, 2, 8, 1, 9 })
    v:asc()
    assert(teq(v:table(), { 1, 2, 5, 8, 9 }))
  end
end)

test("ivec/dvec: desc sort", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 5, 2, 8, 1, 9 })
    v:desc()
    assert(teq(v:table(), { 9, 8, 5, 2, 1 }))
  end
end)

test("ivec/dvec: asc with range", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 5, 2, 8, 9 })
    v:asc(1, 4)
    assert(v:get(0) == 1)
    assert(v:get(1) == 2)
    assert(v:get(2) == 5)
    assert(v:get(3) == 8)
    assert(v:get(4) == 9)
  end
end)

test("ivec/dvec: uasc (unique ascending)", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 3, 1, 2, 1, 3, 2 })
    local new_end = v:uasc()
    v:setn(new_end)
    assert(teq(v:table(), { 1, 2, 3 }))
  end
end)

test("ivec/dvec: kasc (partial sort)", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 9, 5, 2, 8, 1, 7 })
    v:kasc(3)
    assert(v:get(0) <= v:get(1))
    assert(v:get(1) <= v:get(2))
  end
end)

test("ivec/dvec: table conversion", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5, 6 })
    assert(teq(v:table(), { 1, 2, 3, 4, 5, 6 }))
    assert(teq(v:table(0, 3), { 1, 2, 3 }))
  end
end)

test("ivec/dvec: rtable (row-major)", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5, 6 })
    assert(teq(v:rtable(3), {
      { 1, 2, 3 },
      { 4, 5, 6 }
    }))
  end
end)

test("ivec/dvec: ctable (column-major)", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5, 6 })
    assert(teq(v:ctable(3), {
      { 1, 4 },
      { 2, 5 },
      { 3, 6 }
    }))
  end
end)

test("ivec/dvec: add scalar", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    v:add(10)
    assert(teq(v:table(), { 11, 12, 13, 14, 15 }))
  end
end)

test("ivec/dvec: add scalar with range", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    v:add(10, 1, 4)
    assert(teq(v:table(), { 1, 12, 13, 14, 5 }))
  end
end)

test("ivec/dvec: scale", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    v:scale(2)
    assert(teq(v:table(), { 2, 4, 6, 8, 10 }))
  end
end)

test("ivec/dvec: scalev", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v1 = vec.create({ 1, 2, 3 })
    local v2 = vec.create({ 2, 3, 4 })
    v1:scalev(v2)
    assert(teq(v1:table(), { 2, 6, 12 }))
  end
end)

test("ivec/dvec: addv", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v1 = vec.create({ 1, 2, 3 })
    local v2 = vec.create({ 10, 20, 30 })
    v1:addv(v2)
    assert(teq(v1:table(), { 11, 22, 33 }))
  end
end)

test("ivec/dvec: sum", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5 })
    assert(v:sum() == 15)
  end
end)

test("ivec/dvec: max", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 3, 7, 2, 9, 1 })
    local val, idx = v:max()
    assert(val == 9)
    assert(idx == 3)
  end
end)

test("ivec/dvec: min", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 3, 7, 2, 9, 1 })
    local val, idx = v:min()
    assert(val == 1)
    assert(idx == 4)
  end
end)

test("ivec/dvec: dot product", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v1 = vec.create({ 1, 2, 3 })
    local v2 = vec.create({ 4, 5, 6 })
    assert(v1:dot(v2) == 32)
  end
end)

test("ivec/dvec: magnitude", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 3, 4 })
    assert(math.abs(v:magnitude() - 5) < 1e-10)
  end
end)

test("ivec/dvec: each iterator", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4 })
    local sum = 0
    for val in v:each() do
      sum = sum + val
    end
    assert(sum == 10)
  end
end)

test("ivec/dvec: ieach iterator", function ()
  for _, vec in ipairs({ ivec, dvec }) do
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

test("ivec/dvec: find", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 10, 20, 30, 40 })
    assert(v:find(30) == 2)
    assert(v:find(99) == nil)
  end
end)

test("ivec/dvec: transpose", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local src = vec.create({ 1, 2, 3, 4, 5, 6 })
    local dst = vec.create(6)
    dst:transpose(src, 3)
    assert(teq(dst:table(), { 1, 4, 2, 5, 3, 6 }))
  end
end)

test("ivec/dvec: rsums (row sums)", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5, 6 })
    local sums = v:rsums(3)
    assert(teq(sums:table(), { 6, 15 }))
  end
end)

test("ivec/dvec: csums (column sums)", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4, 5, 6 })
    local sums = v:csums(3)
    assert(teq(sums:table(), { 5, 7, 9 }))
  end
end)

test("ivec/dvec: rmaxs", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 3, 1, 2, 4, 6, 5 })
    local maxs = v:rmaxs(3)
    assert(maxs:size() > 0)
  end
end)

test("ivec/dvec: rmins", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 3, 1, 2, 4, 6, 5 })
    local mins = v:rmins(3)
    assert(mins:size() > 0)
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

test("dvec: center", function ()
  local v = dvec.create({ 1, 2, 3, 4, 5, 6 })
  v:center(2, 3)
  local centered = v:table()
  assert(math.abs(centered[1] - (-1.5)) < 1e-10)
  assert(math.abs(centered[4] - 1.5) < 1e-10)
end)

test("dvec: rnorml2 (L2 row normalization)", function ()
  local v = dvec.create({ 3, 4, 0, 0, 5, 12 })
  v:rnorml2(2, 3)
  local normed = v:table()
  assert(math.abs(normed[1] - 0.6) < 1e-10)
  assert(math.abs(normed[2] - 0.8) < 1e-10)
end)

test("dvec: multiply", function ()
  local a = dvec.create({ 1, 2, 3, 4 })
  local b = dvec.create({ 5, 6, 7, 8 })
  local c = dvec.create(4)
  c:multiply(a, b, 2, false, true)
  assert(c:size() == 4)
end)

test("dvec: scores_elbow methods", function ()
  local v = dvec.create({ 100, 80, 50, 30, 25, 24, 23, 22 })
  local _, idx = v:scores_elbow("max_gap")
  assert(idx >= 1 and idx <= 8)
  _, idx = v:scores_elbow("max_curvature")
  assert(idx >= 1 and idx <= 8)
  _, idx = v:scores_elbow("kneedle")
  assert(idx >= 1 and idx <= 8)
  _, idx = v:scores_elbow("lmethod")
  assert(idx >= 1 and idx <= 8)
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

test("ivec: bits_select", function ()
  local corpus = ivec.create({ 0, 4, 8, 12 })
  local subcorpus = ivec.create()
  corpus:bits_select(nil, ivec.create({ 0, 3 }), 4, subcorpus)
  assert(teq({ 0, 4 }, subcorpus:table()))
end)

test("ivec: bits_to_cvec and bits_top_mi", function ()
  require("santoku.cvec")
  local n_samples = 3
  local n_features = 4
  local n_hidden = 2
  local features = ivec.create({ 0, 1, 4, 6, 9, 11 })
  local labels = ivec.create({ 0, 1, 3, 4 })
  local top_features, weights = features:bits_top_mi(labels, n_samples, n_features, n_hidden, 2)
  assert(top_features:size() == 2)
  assert(weights:size() == 2)
  local w = weights:table()
  assert(w[1] >= w[2])
end)

test("ivec: bits_top_chi2", function ()
  require("santoku.cvec")
  local n_samples = 3
  local n_features = 4
  local n_hidden = 2
  local features = ivec.create({ 0, 1, 4, 6, 9, 11 })
  local labels = ivec.create({ 0, 1, 3, 4 })
  local top_features, weights = features:bits_top_chi2(labels, n_samples, n_features, n_hidden, 2)
  assert(top_features:size() == 2)
  assert(weights:size() == 2)
end)

test("ivec: cvec roundtrip", function ()
  require("santoku.cvec")
  local n_samples = 3
  local n_features = 4
  local features = ivec.create({ 0, 1, 4, 6, 9, 11 })
  local bitmap = features:bits_to_cvec(n_samples, n_features)
  assert(bitmap ~= nil)
end)

test("ivec: scores_elbow", function ()
  local v = ivec.create({ 100, 80, 50, 30, 25, 24, 23, 22 })
  local _, idx = v:scores_elbow("max_gap")
  assert(idx >= 1 and idx <= 8)
  _, idx = v:scores_elbow("max_curvature")
  assert(idx >= 1 and idx <= 8)
end)

test("dvec: mtx_top_variance", function ()
  local matrix = dvec.create({ 1, 2, 10, 2, 3, 20, 3, 4, 30 })
  local top_v, scores = matrix:mtx_top_variance(3, 3, 2)
  assert(top_v:size() == 2)
  assert(scores:size() == 2)
end)

test("ivec/dvec: persist and load", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3, 4 })
    v:persist("/data/data/com.termux/files/home/.vec_test.bin")
    local loaded = vec.load("/data/data/com.termux/files/home/.vec_test.bin")
    assert(teq(v:table(), loaded:table()))
  end
end)

test("ivec/dvec: raw", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create({ 1, 2, 3 })
    local raw = v:raw()
    assert(type(raw) == "string")
    assert(#raw > 0)
  end
end)

test("ivec/dvec: shrink", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local v = vec.create(100)
    v:setn(10)
    v:shrink()
    assert(v:size() == 10)
  end
end)

test("ivec/dvec: ensure", function ()
  for _, vec in ipairs({ ivec, dvec }) do
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

test("ivec: rdesc (row descending argsort)", function ()
  local v = ivec.create({ 10, 2, 3, 4, 5, 6 })
  local result = v:rdesc(3)
  assert(result ~= nil)
end)

test("ivec: casc (column ascending argsort)", function ()
  local v = ivec.create({ 10, 2, 3, 4, 5, 6 })
  local result = v:casc(3)
  assert(result ~= nil)
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

test("dvec: rmagnitudes (row L2 norms)", function ()
  local v = dvec.create({ 3, 4, 0, 5, 12, 0 })
  local mags = v:rmagnitudes(3)
  assert(mags:size() == 2)
  assert(math.abs(mags:get(0) - 5) < 1e-10)
  assert(math.abs(mags:get(1) - 13) < 1e-10)
end)

test("dvec: cmagnitudes (column L2 norms)", function ()
  local v = dvec.create({ 3, 4, 0, 4, 3, 0 })
  local mags = v:cmagnitudes(3)
  assert(mags:size() == 3)
  assert(math.abs(mags:get(0) - 5) < 1e-10)
  assert(math.abs(mags:get(1) - 5) < 1e-10)
  assert(math.abs(mags:get(2) - 0) < 1e-10)
end)

test("dvec: rmags (BLAS row magnitudes)", function ()
  local v = dvec.create({ 3, 4, 0, 5, 12, 0 })
  local mags = v:rmags(3)
  assert(mags:size() == 2)
  assert(math.abs(mags:get(0) - 5) < 1e-10)
  assert(math.abs(mags:get(1) - 13) < 1e-10)
end)

test("dvec: cmags (BLAS column magnitudes)", function ()
  local v = dvec.create({ 3, 4, 0, 4, 3, 0 })
  local mags = v:cmags(3)
  assert(mags:size() == 3)
  assert(math.abs(mags:get(0) - 5) < 1e-10)
  assert(math.abs(mags:get(1) - 5) < 1e-10)
  assert(math.abs(mags:get(2) - 0) < 1e-10)
end)

test("dvec: rmaxargs (row max indices)", function ()
  local v = dvec.create({ 1, 5, 2, 9, 3, 7 })
  local args = v:rmaxargs(3)
  assert(args:size() == 2)
  assert(args:get(0) == 1)
  assert(args:get(1) == 0)
end)

test("dvec: cmaxargs (column max indices)", function ()
  local v = dvec.create({ 1, 5, 2, 9, 3, 7 })
  local args = v:cmaxargs(3)
  assert(args:size() == 3)
  assert(args:get(0) == 1)
  assert(args:get(1) == 0)
  assert(args:get(2) == 1)
end)

test("dvec: rminargs (row min indices)", function ()
  local v = dvec.create({ 1, 5, 2, 9, 3, 7 })
  local args = v:rminargs(3)
  assert(args:size() == 2)
  assert(args:get(0) == 0)
  assert(args:get(1) == 1)
end)

test("dvec: cminargs (column min indices)", function ()
  local v = dvec.create({ 1, 5, 2, 9, 3, 7 })
  local args = v:cminargs(3)
  assert(args:size() == 3)
  assert(args:get(0) == 0)
  assert(args:get(1) == 1)
  assert(args:get(2) == 0)
end)

test("dvec: rasc (row ascending argsort)", function ()
  local v = dvec.create({ 5, 1, 3, 9, 2, 7 })
  local indices = v:rasc(3)
  assert(indices:size() == 6)
  assert(indices:get(0) == 1)
  assert(indices:get(1) == 2)
  assert(indices:get(2) == 0)
  assert(indices:get(3) == 1)
  assert(indices:get(4) == 2)
  assert(indices:get(5) == 0)
end)

test("dvec: rdesc (row descending argsort)", function ()
  local v = dvec.create({ 5, 1, 3, 9, 2, 7 })
  local indices = v:rdesc(3)
  assert(indices:size() == 6)
  assert(indices:get(0) == 0)
  assert(indices:get(1) == 2)
  assert(indices:get(2) == 1)
  assert(indices:get(3) == 0)
  assert(indices:get(4) == 2)
  assert(indices:get(5) == 1)
end)

test("dvec: casc (column ascending argsort, column-major)", function ()
  local v = dvec.create({ 5, 1, 3, 2, 4, 6 })
  local indices = v:casc(3)
  assert(indices:size() == 6)
  assert(indices:get(0) == 1)
  assert(indices:get(1) == 0)
  assert(indices:get(2) == 1)
  assert(indices:get(3) == 0)
  assert(indices:get(4) == 0)
  assert(indices:get(5) == 1)
end)

test("dvec: cdesc (column descending argsort, column-major)", function ()
  local v = dvec.create({ 5, 1, 3, 2, 4, 6 })
  local indices = v:cdesc(3)
  assert(indices:size() == 6)
  assert(indices:get(0) == 0)
  assert(indices:get(1) == 1)
  assert(indices:get(2) == 0)
  assert(indices:get(3) == 1)
  assert(indices:get(4) == 1)
  assert(indices:get(5) == 0)
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

test("ivec/dvec: copy_indexed", function ()
  for _, vec in ipairs({ ivec, dvec }) do
    local src = vec.create({ 100, 200, 300, 400, 500 })
    local dst = vec.create(3)
    local indices = ivec.create({ 4, 0, 2 })
    dst:copy(src, indices)
    assert(teq(dst:table(), { 500, 100, 300 }))
  end
end)

test("dvec: multiplyv (matrix-vector)", function ()
  local A = dvec.create({ 1, 2, 3, 4, 5, 6 })
  local x = dvec.create({ 1, 2, 3 })
  local y = dvec.create(2)
  y:multiplyv(A, x, 3)
  assert(y:size() == 2)
  assert(y:get(0) == 14)
  assert(y:get(1) == 32)
end)

test("dvec: multiplyv with transpose", function ()
  local A = dvec.create({ 1, 2, 3, 4, 5, 6 })
  local x = dvec.create({ 1, 2 })
  local y = dvec.create(3)
  y:multiplyv(A, x, 3, true)
  assert(y:size() == 3)
  assert(y:get(0) == 9)
  assert(y:get(1) == 12)
  assert(y:get(2) == 15)
end)

test("rvec: rankings", function ()
  local scores = dvec.create({ 0.9, 0.1, 0.5, 0.8, 0.2, 0.6 })
  local rankings = rvec.rankings(scores, 3, 2)
  assert(rankings:size() == 6)
end)

test("ivec: from_rvec", function ()
  local r = rvec.create()
  r:push(10, 1.0)
  r:push(20, 2.0)
  r:push(30, 3.0)
  local keys = ivec.from_rvec(r)
  assert(keys:size() == 3)
  assert(keys:get(0) == 10)
  assert(keys:get(1) == 20)
  assert(keys:get(2) == 30)
end)
