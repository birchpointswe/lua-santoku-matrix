local test = require("santoku.test")
local err = require("santoku.error")
local assert = err.assert
local dvec = require("santoku.dvec")
local ivec = require("santoku.ivec")
local fvec = require("santoku.fvec")
local tbl = require("santoku.table")
local teq = tbl.equals

local tmp = os.getenv("PREFIX") and (os.getenv("PREFIX") .. "/tmp") or "/tmp"

test("ivec/dvec/fvec: mmap_create and basic access", function ()
  for i, vec in ipairs({ ivec, dvec, fvec }) do
    local path = tmp .. "/test_mmap_create_" .. i .. ".bin"
    local v = vec.mmap_create(path, 5)
    assert(v:size() == 5)
    for j = 0, 4 do
      assert(v:get(j) == 0)
    end
    v:set(0, 10)
    v:set(2, 30)
    v:set(4, 50)
    assert(v:get(0) == 10)
    assert(v:get(2) == 30)
    assert(v:get(4) == 50)
    os.remove(path)
  end
end)

test("ivec/dvec/fvec: mmap_create, sync, and mmap_open", function ()
  for i, vec in ipairs({ ivec, dvec, fvec }) do
    local path = tmp .. "/test_mmap_cycle_" .. i .. ".bin"
    do
      local v = vec.mmap_create(path, 4)
      v:set(0, 1)
      v:set(1, 2)
      v:set(2, 3)
      v:set(3, 4)
      v:mmap_sync()
    end
    collectgarbage()
    local loaded = vec.mmap_open(path)
    assert(loaded:size() == 4)
    assert(loaded:get(0) == 1)
    assert(loaded:get(1) == 2)
    assert(loaded:get(2) == 3)
    assert(loaded:get(3) == 4)
    os.remove(path)
  end
end)

test("ivec/dvec: mmap table roundtrip", function ()
  for i, vec in ipairs({ ivec, dvec }) do
    local path = tmp .. "/test_mmap_table_" .. i .. ".bin"
    do
      local v = vec.mmap_create(path, 5)
      v:set(0, 10)
      v:set(1, 20)
      v:set(2, 30)
      v:set(3, 40)
      v:set(4, 50)
      v:mmap_sync()
    end
    collectgarbage()
    local loaded = vec.mmap_open(path)
    assert(teq(loaded:table(), { 10, 20, 30, 40, 50 }))
    os.remove(path)
  end
end)

test("ivec/dvec: mmap_open empty file", function ()
  for i, vec in ipairs({ ivec, dvec }) do
    local path = tmp .. "/test_mmap_empty_" .. i .. ".bin"
    do
      local v = vec.mmap_create(path, 0)
      assert(v:size() == 0)
      v:mmap_sync()
    end
    collectgarbage()
    local loaded = vec.mmap_open(path)
    assert(loaded:size() == 0)
    os.remove(path)
  end
end)

test("ivec/dvec: mmap copy to regular vec", function ()
  for i, vec in ipairs({ ivec, dvec }) do
    local path = tmp .. "/test_mmap_copy_" .. i .. ".bin"
    local v = vec.mmap_create(path, 3)
    v:set(0, 100)
    v:set(1, 200)
    v:set(2, 300)
    local regular = vec.create()
    regular:copy(v)
    assert(teq(regular:table(), { 100, 200, 300 }))
    regular:push(400)
    assert(regular:size() == 4)
    os.remove(path)
  end
end)

test("ivec/dvec: mmap fill and scale", function ()
  for i, vec in ipairs({ ivec, dvec }) do
    local path = tmp .. "/test_mmap_ops_" .. i .. ".bin"
    do
      local v = vec.mmap_create(path, 4)
      v:fill(5)
      assert(teq(v:table(), { 5, 5, 5, 5 }))
      v:scale(3)
      assert(teq(v:table(), { 15, 15, 15, 15 }))
      v:mmap_sync()
    end
    collectgarbage()
    local loaded = vec.mmap_open(path)
    assert(teq(loaded:table(), { 15, 15, 15, 15 }))
    os.remove(path)
  end
end)

test("dvec: mmap dot product", function ()
  local p1 = tmp .. "/test_mmap_dot_1.bin"
  local p2 = tmp .. "/test_mmap_dot_2.bin"
  local v1 = dvec.mmap_create(p1, 3)
  local v2 = dvec.mmap_create(p2, 3)
  v1:set(0, 1) v1:set(1, 2) v1:set(2, 3)
  v2:set(0, 4) v2:set(1, 5) v2:set(2, 6)
  assert(v1:dot(v2) == 32)
  os.remove(p1)
  os.remove(p2)
end)

test("dvec: mmap addv between mmap and regular", function ()
  local path = tmp .. "/test_mmap_addv.bin"
  local v1 = dvec.mmap_create(path, 3)
  v1:set(0, 1) v1:set(1, 2) v1:set(2, 3)
  local v2 = dvec.create({ 10, 20, 30 })
  v1:addv(v2)
  assert(teq(v1:table(), { 11, 22, 33 }))
  os.remove(path)
end)
