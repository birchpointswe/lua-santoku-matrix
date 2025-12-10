local arr = require("santoku.array")
local fs = require("santoku.fs")
local sys = require("santoku.system")
arr.flatten = arr.flatten or function(t, d) d = d or 1; local r = {}; local function f(a, l) for i = 1, #a do local v = a[i]; if l > 0 and type(v) == "table" then f(v, l - 1) else r[#r + 1] = v end end end; f(t, d); return r end
local base = fs.runfile("make.common.lua")
base.env.blas_target = "GENERIC"
base.env.primme_cflags = "-O2 -fPIC"
base.env.primme_ldflags = "-O2 -fPIC"
base.env.cflags = arr.flatten({ { "-g3", "-O0", "-fno-inline", "-fno-omit-frame-pointer" }, base.env.cflags or {} })
base.env.ldflags = arr.flatten({ { "-g3", "-O0" }, base.env.ldflags or {} })
return base
