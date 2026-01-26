local test = require("santoku.test")
local err = require("santoku.error")
local assert = err.assert
local ivec = require("santoku.ivec")
require("santoku.cvec")
local tbl = require("santoku.table")
local teq = tbl.equals

test("ivec:bits_top_df basic", function ()
  local n_samples = 4
  local n_features = 5
  local bits = ivec.create({
    0 * n_features + 0,
    0 * n_features + 1,
    0 * n_features + 2,
    1 * n_features + 0,
    1 * n_features + 1,
    1 * n_features + 3,
    2 * n_features + 0,
    2 * n_features + 2,
    3 * n_features + 0,
    3 * n_features + 3,
  })
  local top_features, weights = bits:bits_top_df(n_samples, n_features, 3)
  assert(top_features:size() == 3)
  assert(weights:size() == 3)
  local w = weights:table()
  assert(w[1] >= w[2] and w[2] >= w[3], "weights should be descending")
  local f = top_features:table()
  assert(not teq(f, {0}), "feature 0 (df=4) should be excluded")
end)

test("ivec:bits_top_df excludes df=0 and df=n_samples", function ()
  local n_samples = 3
  local n_features = 4
  local bits = ivec.create({
    0 * n_features + 0,
    1 * n_features + 0,
    2 * n_features + 0,
    0 * n_features + 1,
    1 * n_features + 1,
  })
  local top_features = bits:bits_top_df(n_samples, n_features)
  local f = top_features:table()
  for _, feat in ipairs(f) do
    assert(feat ~= 0, "feature 0 (df=3=n_samples) should be excluded")
    assert(feat ~= 2, "feature 2 (df=0) should be excluded")
    assert(feat ~= 3, "feature 3 (df=0) should be excluded")
  end
  assert(teq(f, {1}), "only feature 1 (df=2) should remain")
end)

test("ivec:bits_top_df min_df/max_df filtering", function ()
  local n_samples = 10
  local n_features = 4
  local bits = ivec.create({})
  for s = 0, 0 do bits:push(s * n_features + 0) end
  for s = 0, 2 do bits:push(s * n_features + 1) end
  for s = 0, 4 do bits:push(s * n_features + 2) end
  for s = 0, 7 do bits:push(s * n_features + 3) end
  local top_features, _ = bits:bits_top_df(n_samples, n_features, 4, 0.2, 0.6)
  local f = top_features:table()
  assert(not teq(f, {0}), "feature 0 (df=1, 10%) should be excluded by min_df=0.2")
  assert(not teq(f, {3}), "feature 3 (df=8, 80%) should be excluded by max_df=0.6")
end)

test("cvec:bits_top_df basic", function ()
  local n_samples = 4
  local n_features = 5
  local bits = ivec.create({
    0 * n_features + 0,
    0 * n_features + 1,
    0 * n_features + 2,
    1 * n_features + 0,
    1 * n_features + 1,
    1 * n_features + 3,
    2 * n_features + 0,
    2 * n_features + 2,
    3 * n_features + 0,
    3 * n_features + 3,
  })
  local bitmap = bits:bits_to_cvec(n_samples, n_features)
  local top_features, weights = bitmap:bits_top_df(n_samples, n_features, 3)
  assert(top_features:size() == 3)
  assert(weights:size() == 3)
  local w = weights:table()
  assert(w[1] >= w[2] and w[2] >= w[3], "weights should be descending")
end)

test("cvec:bits_top_df excludes df=0 and df=n_samples", function ()
  local n_samples = 3
  local n_features = 4
  local bits = ivec.create({
    0 * n_features + 0,
    1 * n_features + 0,
    2 * n_features + 0,
    0 * n_features + 1,
    1 * n_features + 1,
  })
  local bitmap = bits:bits_to_cvec(n_samples, n_features)
  local top_features, _ = bitmap:bits_top_df(n_samples, n_features)
  local f = top_features:table()
  for _, feat in ipairs(f) do
    assert(feat ~= 0, "feature 0 (df=3=n_samples) should be excluded")
    assert(feat ~= 2, "feature 2 (df=0) should be excluded")
    assert(feat ~= 3, "feature 3 (df=0) should be excluded")
  end
  assert(teq(f, {1}), "only feature 1 (df=2) should remain")
end)

test("ivec:bits_top_bns with labels basic", function ()
  local n_samples = 6
  local n_features = 4
  local n_classes = 2
  local bits = ivec.create({
    0 * n_features + 0,
    0 * n_features + 1,
    1 * n_features + 0,
    1 * n_features + 1,
    2 * n_features + 0,
    2 * n_features + 1,
    3 * n_features + 2,
    3 * n_features + 3,
    4 * n_features + 2,
    4 * n_features + 3,
    5 * n_features + 2,
    5 * n_features + 3,
  })
  local labels = ivec.create({
    0 * n_classes + 0,
    1 * n_classes + 0,
    2 * n_classes + 0,
    3 * n_classes + 1,
    4 * n_classes + 1,
    5 * n_classes + 1,
  })
  local top_features, weights = bits:bits_top_bns(labels, n_samples, n_features, n_classes, 2)
  assert(top_features:size() == 2)
  assert(weights:size() == 2)
  local w = weights:table()
  assert(w[1] >= w[2], "weights should be descending")
end)

test("ivec:bits_top_bns excludes uninformative features", function ()
  local n_samples = 4
  local n_features = 3
  local n_classes = 2
  local bits = ivec.create({
    0 * n_features + 0,
    1 * n_features + 0,
    2 * n_features + 0,
    3 * n_features + 0,
    0 * n_features + 1,
    2 * n_features + 1,
  })
  local labels = ivec.create({
    0 * n_classes + 0,
    1 * n_classes + 0,
    2 * n_classes + 1,
    3 * n_classes + 1,
  })
  local top_features, _ = bits:bits_top_bns(labels, n_samples, n_features, n_classes)
  local f = top_features:table()
  for _, feat in ipairs(f) do
    assert(feat ~= 0, "feature 0 (appears in all samples) should be excluded")
    assert(feat ~= 2, "feature 2 (appears in no samples) should be excluded")
  end
end)

test("cvec:bits_top_bns with labels basic", function ()
  local n_samples = 6
  local n_features = 4
  local n_classes = 2
  local bits = ivec.create({
    0 * n_features + 0,
    0 * n_features + 1,
    1 * n_features + 0,
    1 * n_features + 1,
    2 * n_features + 0,
    2 * n_features + 1,
    3 * n_features + 2,
    3 * n_features + 3,
    4 * n_features + 2,
    4 * n_features + 3,
    5 * n_features + 2,
    5 * n_features + 3,
  })
  local bitmap = bits:bits_to_cvec(n_samples, n_features)
  local labels = ivec.create({
    0 * n_classes + 0,
    1 * n_classes + 0,
    2 * n_classes + 0,
    3 * n_classes + 1,
    4 * n_classes + 1,
    5 * n_classes + 1,
  })
  local top_features, weights = bitmap:bits_top_bns(labels, n_samples, n_features, n_classes, 2)
  assert(top_features:size() == 2)
  assert(weights:size() == 2)
end)

test("ivec:bits_top_bns_ind with labels", function ()
  local n_samples = 6
  local n_features = 4
  local n_classes = 2
  local bits = ivec.create({
    0 * n_features + 0,
    0 * n_features + 1,
    1 * n_features + 0,
    1 * n_features + 1,
    2 * n_features + 0,
    2 * n_features + 1,
    3 * n_features + 2,
    3 * n_features + 3,
    4 * n_features + 2,
    4 * n_features + 3,
    5 * n_features + 2,
    5 * n_features + 3,
  })
  local labels = ivec.create({
    0 * n_classes + 0,
    1 * n_classes + 0,
    2 * n_classes + 0,
    3 * n_classes + 1,
    4 * n_classes + 1,
    5 * n_classes + 1,
  })
  local ids_union, ids, offsets, weights = bits:bits_top_bns_ind(labels, n_samples, n_features, n_classes, 2)
  assert(ids_union ~= nil, "ids_union should not be nil")
  assert(ids ~= nil, "ids should not be nil")
  assert(offsets ~= nil, "offsets should not be nil")
  assert(offsets:size() == n_classes + 1, "offsets should have n_classes+1 entries")
  assert(weights ~= nil, "weights should not be nil")
  assert(ids:size() == weights:size(), "ids and weights should have same size")
end)

test("cvec:bits_top_bns_ind with labels", function ()
  local n_samples = 6
  local n_features = 4
  local n_classes = 2
  local bits = ivec.create({
    0 * n_features + 0,
    0 * n_features + 1,
    1 * n_features + 0,
    1 * n_features + 1,
    2 * n_features + 0,
    2 * n_features + 1,
    3 * n_features + 2,
    3 * n_features + 3,
    4 * n_features + 2,
    4 * n_features + 3,
    5 * n_features + 2,
    5 * n_features + 3,
  })
  local bitmap = bits:bits_to_cvec(n_samples, n_features)
  local labels = ivec.create({
    0 * n_classes + 0,
    1 * n_classes + 0,
    2 * n_classes + 0,
    3 * n_classes + 1,
    4 * n_classes + 1,
    5 * n_classes + 1,
  })
  local ids_union, ids, offsets, weights = bitmap:bits_top_bns_ind(labels, n_samples, n_features, n_classes, 2)
  assert(ids_union ~= nil, "ids_union should not be nil")
  assert(ids ~= nil, "ids should not be nil")
  assert(offsets ~= nil, "offsets should not be nil")
  assert(offsets:size() == n_classes + 1, "offsets should have n_classes+1 entries")
  assert(weights ~= nil, "weights should not be nil")
  assert(ids:size() == weights:size(), "ids and weights should have same size")
end)

test("ivec:bits_top_df IDF value correctness", function ()
  local n_samples = 100
  local n_features = 3
  local bits = ivec.create({})
  for s = 0, 0 do bits:push(s * n_features + 0) end
  for s = 0, 49 do bits:push(s * n_features + 1) end
  for s = 0, 98 do bits:push(s * n_features + 2) end
  local top_features, weights = bits:bits_top_df(n_samples, n_features)
  local f = top_features:table()
  local w = weights:table()
  assert(f[1] == 0, "feature 0 (df=1) should have highest IDF")
  assert(f[2] == 1, "feature 1 (df=50) should have medium IDF")
  local expected_idf_0 = math.log((n_samples + 1) / (1 + 1))
  local expected_idf_1 = math.log((n_samples + 1) / (50 + 1))
  assert(math.abs(w[1] - expected_idf_0) < 0.01, "IDF for df=1 should be ~" .. expected_idf_0)
  assert(math.abs(w[2] - expected_idf_1) < 0.01, "IDF for df=50 should be ~" .. expected_idf_1)
end)

test("ivec:bits_top_df weights are distinct (not all identical)", function ()
  local n_samples = 50
  local n_features = 10
  local bits = ivec.create({})
  for f = 0, n_features - 1 do
    local df = (f + 1) * 3
    if df < n_samples then
      for s = 0, df - 1 do
        bits:push(s * n_features + f)
      end
    end
  end
  local _, weights = bits:bits_top_df(n_samples, n_features, 5)
  local w = weights:table()
  assert(#w >= 2, "should have at least 2 features")
  assert(w[1] ~= w[2], "IDF weights should not all be identical (regression test)")
end)

test("ivec:bits_top_bns multi-label", function ()
  local n_samples = 4
  local n_features = 4
  local n_labels = 3
  local bits = ivec.create({
    0 * n_features + 0,
    0 * n_features + 1,
    1 * n_features + 0,
    1 * n_features + 2,
    2 * n_features + 1,
    2 * n_features + 3,
    3 * n_features + 2,
    3 * n_features + 3,
  })
  local labels = ivec.create({
    0 * n_labels + 0,
    0 * n_labels + 1,
    1 * n_labels + 0,
    2 * n_labels + 1,
    2 * n_labels + 2,
    3 * n_labels + 2,
  })
  local top_features, weights = bits:bits_top_bns(labels, n_samples, n_features, n_labels, 3)
  assert(top_features ~= nil, "top_features should not be nil")
  assert(weights ~= nil, "weights should not be nil")
  assert(top_features:size() <= 3, "should return at most top_k features")
  assert(top_features:size() == weights:size(), "features and weights should match")
end)

test("ivec:bits_top_bns_ind with codes (cvec)", function ()
  local n_samples = 6
  local n_features = 4
  local n_dims = 2
  local bits = ivec.create({
    0 * n_features + 0,
    0 * n_features + 1,
    1 * n_features + 0,
    1 * n_features + 1,
    2 * n_features + 0,
    2 * n_features + 1,
    3 * n_features + 2,
    3 * n_features + 3,
    4 * n_features + 2,
    4 * n_features + 3,
    5 * n_features + 2,
    5 * n_features + 3,
  })
  local code_bits = ivec.create({
    0 * n_dims + 0,
    1 * n_dims + 0,
    2 * n_dims + 0,
    3 * n_dims + 1,
    4 * n_dims + 1,
    5 * n_dims + 1,
  })
  local codes = code_bits:bits_to_cvec(n_samples, n_dims)
  local ids_union, ids, offsets, weights = bits:bits_top_bns_ind(codes, n_samples, n_features, n_dims, 2)
  assert(ids_union ~= nil, "ids_union should not be nil")
  assert(ids ~= nil, "ids should not be nil")
  assert(offsets ~= nil, "offsets should not be nil")
  assert(offsets:size() == n_dims + 1, "offsets should have n_dims+1 entries")
  assert(weights ~= nil, "weights should not be nil")
  assert(ids:size() == weights:size(), "ids and weights should have same size")
  local o = offsets:table()
  assert(o[1] == 0, "offsets should start at 0")
  for i = 2, #o do
    assert(o[i] >= o[i-1], "offsets should be non-decreasing")
  end
end)

test("cvec:bits_top_bns_ind with codes (cvec)", function ()
  local n_samples = 6
  local n_features = 4
  local n_dims = 2
  local bits = ivec.create({
    0 * n_features + 0,
    0 * n_features + 1,
    1 * n_features + 0,
    1 * n_features + 1,
    2 * n_features + 0,
    2 * n_features + 1,
    3 * n_features + 2,
    3 * n_features + 3,
    4 * n_features + 2,
    4 * n_features + 3,
    5 * n_features + 2,
    5 * n_features + 3,
  })
  local bitmap = bits:bits_to_cvec(n_samples, n_features)
  local code_bits = ivec.create({
    0 * n_dims + 0,
    1 * n_dims + 0,
    2 * n_dims + 0,
    3 * n_dims + 1,
    4 * n_dims + 1,
    5 * n_dims + 1,
  })
  local codes = code_bits:bits_to_cvec(n_samples, n_dims)
  local ids_union, ids, offsets, weights = bitmap:bits_top_bns_ind(codes, n_samples, n_features, n_dims, 2)
  assert(ids_union ~= nil, "ids_union should not be nil")
  assert(ids ~= nil, "ids should not be nil")
  assert(offsets ~= nil, "offsets should not be nil")
  assert(offsets:size() == n_dims + 1, "offsets should have n_dims+1 entries")
  assert(weights ~= nil, "weights should not be nil")
  assert(ids:size() == weights:size(), "ids and weights should have same size")
end)

test("ivec:bits_top_chi2_ind with codes (cvec)", function ()
  local n_samples = 6
  local n_features = 4
  local n_dims = 2
  local bits = ivec.create({
    0 * n_features + 0,
    0 * n_features + 1,
    1 * n_features + 0,
    1 * n_features + 1,
    2 * n_features + 0,
    2 * n_features + 1,
    3 * n_features + 2,
    3 * n_features + 3,
    4 * n_features + 2,
    4 * n_features + 3,
    5 * n_features + 2,
    5 * n_features + 3,
  })
  local code_bits = ivec.create({
    0 * n_dims + 0,
    1 * n_dims + 0,
    2 * n_dims + 0,
    3 * n_dims + 1,
    4 * n_dims + 1,
    5 * n_dims + 1,
  })
  local codes = code_bits:bits_to_cvec(n_samples, n_dims)
  local ids_union, ids, offsets, weights = bits:bits_top_chi2_ind(codes, n_samples, n_features, n_dims, 2)
  assert(ids_union ~= nil, "ids_union should not be nil")
  assert(ids ~= nil, "ids should not be nil")
  assert(offsets ~= nil, "offsets should not be nil")
  assert(offsets:size() == n_dims + 1, "offsets should have n_dims+1 entries")
  assert(weights ~= nil, "weights should not be nil")
  assert(ids:size() == weights:size(), "ids and weights should have same size")
end)
