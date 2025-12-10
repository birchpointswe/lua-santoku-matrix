local arr = require("santoku.array")
local fs = require("santoku.fs")
arr.flatten = arr.flatten or function(t, d) d = d or 1; local r = {}; local function f(a, l) for i = 1, #a do local v = a[i]; if l > 0 and type(v) == "table" then f(v, l - 1) else r[#r + 1] = v end end end; f(t, d); return r end
local base = fs.runfile("make.common.lua")
base.env.cflags = arr.flatten({ { "-O3", "-march=native", "-fdata-sections", "-ffunction-sections" }, base.env.cflags or {} })
base.env.ldflags = arr.flatten({ { "-O3", "-march=native", "-Wl,--gc-sections" }, base.env.ldflags or {} })
return base
