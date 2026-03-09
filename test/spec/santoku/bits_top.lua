local test = require("santoku.test")
local err = require("santoku.error")
local assert = err.assert
local ivec = require("santoku.ivec")
local csr = require("santoku.csr")
require("santoku.cvec")
local tbl = require("santoku.table")
local teq = tbl.equals

test("cvec:bits_top_idf basic", function ()
  local n_samples = 4
  local n_features = 5
  local off = ivec.create({ 0, 3, 6, 8, 10 })
  local nbr = ivec.create({ 0, 1, 2, 0, 1, 3, 0, 2, 0, 3 })
  local bitmap = csr.to_cvec(off, nbr, n_samples, n_features)
  local top_features, weights = bitmap:bits_top_idf(n_samples, n_features, 3)
  assert(top_features:size() == 3)
  assert(weights:size() == 3)
  local w = weights:table()
  assert(w[1] >= w[2] and w[2] >= w[3], "weights should be descending")
end)

test("cvec:bits_top_idf excludes df=0 and df=n_samples", function ()
  local n_samples = 3
  local n_features = 4
  local off = ivec.create({ 0, 2, 4, 5 })
  local nbr = ivec.create({ 0, 1, 0, 1, 0 })
  local bitmap = csr.to_cvec(off, nbr, n_samples, n_features)
  local top_features, _ = bitmap:bits_top_idf(n_samples, n_features)
  local f = top_features:table()
  for _, feat in ipairs(f) do
    assert(feat ~= 0, "feature 0 (df=3=n_samples) should be excluded")
    assert(feat ~= 2, "feature 2 (df=0) should be excluded")
    assert(feat ~= 3, "feature 3 (df=0) should be excluded")
  end
  assert(teq(f, {1}), "only feature 1 (df=2) should remain")
end)

test("cvec:bits_top_bns with labels basic", function ()
  local n_samples = 6
  local n_features = 4
  local n_classes = 2
  local off = ivec.create({ 0, 2, 4, 6, 8, 10, 12 })
  local nbr = ivec.create({ 0, 1, 0, 1, 0, 1, 2, 3, 2, 3, 2, 3 })
  local bitmap = csr.to_cvec(off, nbr, n_samples, n_features)
  local labels = ivec.create({
    0 * n_classes + 0,
    1 * n_classes + 0,
    2 * n_classes + 0,
    3 * n_classes + 1,
    4 * n_classes + 1,
    5 * n_classes + 1,
  })
  local union_ids, union_weights, class_offsets, class_ids, class_weights = bitmap:bits_top_bns(labels, n_samples, n_features, n_classes, nil, 2, "max")
  assert(union_ids:size() == 2)
  assert(union_weights:size() == 2)
  assert(class_offsets:size() == n_classes + 1)
  assert(class_ids:size() > 0)
  assert(class_ids:size() == class_weights:size())
end)
