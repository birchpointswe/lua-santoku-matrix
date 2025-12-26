local arr = require("santoku.array")
local fs = require("santoku.fs")
local base = fs.runfile("make.common.lua")
base.env.cflags = arr.flatten({ { "-O3", "-fdata-sections", "-ffunction-sections" }, base.env.cflags or {} })
base.env.ldflags = arr.flatten({ { "-O3", "-Wl,--gc-sections" }, base.env.ldflags or {} })
return base
