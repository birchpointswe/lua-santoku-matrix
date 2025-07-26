local env = {

  name = "santoku-matrix",
  version = "0.0.90-1",
  variable_prefix = "TK_MATRIX",
  license = "MIT",
  public = true,

  cflags = {
    "-std=gnu11", "-Wall", "-Wextra",
    "-Wstrict-overflow", "-Wsign-conversion", "-Wsign-compare",
    "-I$(shell luarocks show santoku --rock-dir)/include/",
    "-I$(shell luarocks show santoku-threads --rock-dir)/include/",
  },

  ldflags = {
    "-lm", "-lpthread"
  },

  dependencies = {
    "lua == 5.1",
    "santoku >= 0.0.280-1",
    "santoku-threads >= 0.0.11-1",
  },

  test = {
    dependencies = {
      "luacov == 0.15.0-1",
    }
  },

}

env.homepage = "https://github.com/treadwelllane/lua-" .. env.name
env.tarball = env.name .. "-" .. env.version .. ".tar.gz"
env.download = env.homepage .. "/releases/download/" .. env.version .. "/" .. env.tarball

return {
  type = "lib",
  env = env,
}
