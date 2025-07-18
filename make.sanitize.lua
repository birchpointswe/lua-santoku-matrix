local arr = require("santoku.array")
local fs = require("santoku.fs")
local base = fs.runfile("make.common.lua")
base.env.cflags = arr.extend({ "-fsanitize=address", "-g3", "-O0", "-fno-omit-frame-pointer" }, base.env.cflags)
base.env.ldflags = arr.extend({ "-fsanitize=address", "-g3", "-O0" }, base.env.ldflags)
return base
