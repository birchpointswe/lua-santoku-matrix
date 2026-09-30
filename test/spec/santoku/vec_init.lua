local test = require("santoku.test")
local err = require("santoku.error")
local assert = err.assert
local str = require("santoku.string")
local ivec = require("santoku.ivec")
local csr = require("santoku.csr")

test("vec: a vector created from C raises when its module fails to load, then recovers", function ()
  local X = csr.create({
    offsets = ivec.create({ 0, 1 }),
    neighbors = ivec.create({ 0 }),
    n_cols = 8,
  })
  local preload = package.preload["santoku.cvec"]
  package.preload["santoku.cvec"] = function ()
    err.error("cvec unavailable")
  end
  local ok, a, b = err.pcall(function () return X:to_bits() end)
  package.preload["santoku.cvec"] = preload
  package.loaded["santoku.cvec"] = nil
  assert(not ok, "to_bits must raise when santoku.cvec can't load")
  local msg = tostring(a) .. " " .. tostring(b)
  assert(str.find(msg, "cvec unavailable", 1, true), msg)
  local B = X:to_bits()
  assert(B.bits_popcount ~= nil, "a later cvec must get its methods once the module loads")
  assert(B:bits_popcount() == 1)
end)
